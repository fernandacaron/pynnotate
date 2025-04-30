import logging
import traceback
import pandas as pd
from datetime import datetime
import argparse
import os

from dict import synonym_dict_mtDNA_anim, synonym_dict_mtDNA_plan, synonym_dict_rDNA, synonym_dict_cpDNA
from utils import read_csv_list, read_accession_file, get_official_name, download_and_save_genes, load_synonyms

def get_args():
    parser = argparse.ArgumentParser(
        description="""
            This script downloads and extracts specific genetic features from GenBank based on user-defined accession numbers. 
            It allows the user to specify a list of features to retrieve, handle synonyms for gene names, and manage features that  
            appear multiple times within a genome. The extracted sequences are saved in FASTA format for further analysis.  
            Example usage: python script.py -e user@example.com -a accession_list.txt -f features.csv -r repeated.csv -o output_folder  
        """)

    parser.add_argument("-e", "--email",
                        required=True,
                        help="REQUIRED: User e-mail to access NCBI database.")

    parser.add_argument("-a", "--accession",
                        required=True,
                        help="REQUIRED: The full path to a file containing a list of GenBank accession numbers (.txt or .csv).")

    parser.add_argument("-t", "--type",
                        type=str,
                        help="REQUIRED: Type of DNA (i.e. mtDNA_animals, mtDNA_plants, rDNA, cpDNA, other).",
                        required=True)

    parser.add_argument("-f", "--features",
                        type=str,
                        help="File listing the features to download (.txt or .csv).",
                        required=False)

    parser.add_argument("-r", "--repeated",
                        type=str,
                        help="If any feature appears more than once, this file should list these features with more than one copy (.txt or .csv).",
                        required=False)

    parser.add_argument("-o", "--output",
                        type=str,
                        help="Output folder (will be created if it does not exist)",
                        required=False,
                        default=".")

    parser.add_argument("-s", "--add_synonyms",
                        type=str,
                        help="File containing new synonyms (.txt or .csv).",
                        required=False,
                        default="")

    args = parser.parse_args()
    
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        logging.info(f"Output directory '{args.output}' created successfully.")
    elif not os.path.isdir(args.output):
        logging.error(f"The output path '{args.output}' is not a directory.")
    
    args.features = read_csv_list(args.features) if args.features else ["all"]
    args.repeated = read_csv_list(args.repeated) if args.repeated else []

    return args

def main():
    args = get_args()

    Entrez.email = args.email
    
    ## configuring log file
    date_string = datetime.now().strftime("%d%b%y").lower()
    time_string = datetime.now().strftime("%Hh%M").lower()

    if args.output:
        folder = args.output if isinstance(args.output, str) else args.output[0]
        filename = f"{folder}/pynnotate_{date_string}_{time_string}.log"
    else:
        filename = f"pynnotate_{date_string}_{time_string}.log"

    logging.basicConfig(
        level = logging.INFO,
        format = '%(asctime)s - %(levelname)s - %(message)s',
        handlers = [
            logging.FileHandler(filename),
            logging.StreamHandler()
        ]
    )

    ## synonym dictionary
    if args.add_synonyms:
        if args.type == "mtDNA_animals":
            global synonym_dict_mtDNA_anim
            synonym_dict_mtDNA_anim = load_synonyms(args.add_synonyms, synonym_dict_mtDNA_anim)
            logging.info("Updated synonyms (mtDNA - Animals)!")
        if args.type == "mtDNA_plants":
            global synonym_dict_mtDNA_plan
            synonym_dict_mtDNA_plan = load_synonyms(args.add_synonyms, synonym_dict_mtDNA_plan)
            logging.info("Updated synonyms (mtDNA - Plants)!")
        if args.type == "rDNA":
            global synonym_dict_rDNA
            synonym_dict_rDNA = load_synonyms(args.add_synonyms, synonym_dict_rDNA)
            logging.info("Updated synonyms (rDNA)!")
        if args.type == "cpDNA":
            global synonym_dict_cpDNA
            synonym_dict_cpDNA = load_synonyms(args.add_synonyms, synonym_dict_cpDNA)
            logging.info("Updated synonyms (cpDNA)!")
        if args.type == "other":
            synonym_dict_other = load_synonyms(args.add_synonyms)
            logging.info("New synonyms dictionaty created!")

    ## getting accession numbers to download
    accession_ids = read_accession_file(args.accession)

    data_accession = []

    if "all" in args.features:
        if args.type == "mtDNA_animals":
            feature_names = set(synonym_dict_mtDNA_anim.keys())
        if args.type == "mtDNA_plants":
            feature_names = set(synonym_dict_mtDNA_plan.keys())
        if args.type == "rDNA":
            feature_names = set(synonym_dict_rDNA.keys())
        if args.type == "cpDNA":
            feature_names = set(synonym_dict_cpDNA.keys())
        if args.type == "other":
            feature_names = set(synonym_dict_other.keys())
    else: 
        feature_names = args.features

    for accession in accession_ids:
        species_name = None  
        try:
            if args.type == "mtDNA_animals":
                species_name, features_found = download_and_save_genes(accession, synonym_dict_mtDNA_anim)
            if args.type == "mtDNA_plants":
                species_name, features_found = download_and_save_genes(accession, synonym_dict_mtDNA_plan)
            if args.type == "rDNA":
                species_name, features_found = download_and_save_genes(accession, synonym_dict_rDNA)
            if args.type == "cpDNA":
                species_name, features_found = download_and_save_genes(accession, synonym_dict_cpDNA)
            if args.type == "other":
                species_name, features_found = download_and_save_genes(accession, synonym_dict_other)
            row_data = [species_name, accession]

            for official_name in feature_names:
                if official_name in features_found:
                    row_data.append(features_found[official_name])
                else:
                    row_data.append("NA")  

            data_accession.append(row_data)

        except Exception as e:
            error_info = traceback.format_exc()
            logging.error(f"ERROR processing {accession}: {e}")
            logging.error("Error details:", error_info)

    feature_names_prefix = [f"Acc_{name}" for name in feature_names]
    data_accession_df = pd.DataFrame(data_accession, 
                                     columns = ["Sequence_names",
                                                "General_accession"] + 
                                                feature_names_prefix)

    if args.output:
        folder = args.output if isinstance(args.output, str) else args.output[0]
        filename = f"{folder}/spp_accession_{date_string}.csv"
    else:
        filename = f"spp_accession_{date_string}.csv"

    data_accession_df.to_csv(filename, index = False)

    logging.info(f"Files successfully saved in directory '{folder}'")

if __name__ == "__main__":
    main()