import logging
import traceback
import pandas as pd
from Bio import Entrez, SeqIO
from datetime import datetime
import argparse
import csv
import os

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

    parser.add_argument("-f", "--features",
                        type=str,
                        help="CSV file listing the features to download.",
                        required=False)

    parser.add_argument("-r", "--repeated",
                        type=str,
                        help="CSV file listing features with more than one copy.",
                        required=False)

    parser.add_argument("-o", "--output",
                        type=str,
                        help="Output folder (will be created if it does not exist)",
                        required=False,
                        default=".")

    parser.add_argument("-s", "--add_synonyms",
                        type=str,
                        help="CSV file containing new synonyms.",
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

def read_csv_list(file_path):
   
    if not os.path.exists(file_path):
        logging.warning(f"File not found: {file_path}")
        return []
    
    with open(file_path, "r", newline="", encoding="utf-8") as file:
        reader = csv.reader(file)
        return [row[0].strip() for row in reader if row]

def read_accession_file(file_path):
    
    if not os.path.exists(file_path):
        logging.error(f"Accession file not found: {file_path}")
        return []
    
    accession_list = []
    file_extension = os.path.splitext(file_path)[-1].lower()
    
    if file_extension == ".csv":
        with open(file_path, "r", newline="", encoding="utf-8") as file:
            reader = csv.reader(file)
            accession_list = [row[0].strip() for row in reader if row]
    else:  
        with open(file_path, "r", encoding="utf-8") as file:
            accession_list = [line.strip() for line in file if line.strip()]
    
    return accession_list


def get_official_name(feature_name, synonym_dict):
    for official_name, synonyms in synonym_dict.items():
        if feature_name in synonyms:
            return official_name
    return None

def download_and_save_genes(genbank_id, synonym_dict):

    args = get_args()

    ## download genbank file
    handle = Entrez.efetch(db = "nucleotide", id = genbank_id, rettype = "gb",
                           retmode = "text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    ## getting the species name and replace spaces with "_"
    species_name = record.annotations.get("organism", "unknown_species").replace(" ", "_")

    ## printing which species is running
    logging.info(f"Obtaining sequences for {species_name} ({genbank_id})...")
    
    # creating an object to store the features found
    features_found = {}

    ## creating an object to store the missing features
    missing_features = set(synonym_dict.keys())
    
    ## initial empty dictionary for features with more than one copy
    repeated_locations = {}

    if args.repeated:
        for repeated in args.repeated:
            repeated_locations[repeated] = []

    for feature in record.features:
            
        ## skipping feature if "gene" or "source"
        if feature.type in ["gene", "source"]:  
            continue

        ## if CDS, rRNA, tRNA, then get the feature name in the field 
        ## "gene" ou "product"
        if feature.type == "CDS":
            gene_name = feature.qualifiers.get("gene", feature.qualifiers.get("product", ["unknown_gene"]))[0]
            feature_type = "CDS"
        elif feature.type == "rRNA":
            gene_name = feature.qualifiers.get("product", feature.qualifiers.get("gene", ["unknown_rRNA"]))[0]
            feature_type = "rRNA"
        elif feature.type == "tRNA":
            gene_name = feature.qualifiers.get("product", feature.qualifiers.get("gene", ["unknown_tRNA"]))[0]
            feature_type = "tRNA"
        else:
            continue

        ## checking if the feature name is a known synonym
        official_name = get_official_name(gene_name, synonym_dict)

        if "all" in args.features:
            ## if the name is not a known synonym, skip feature
            if not official_name:
                logging.error(f"SYNONYM unknown for {gene_name}, skipping.")
                continue

        if "all" in args.features or official_name in args.features:
            ## extracting the feature sequence
            feature_seq = feature.extract(record.seq)

            ## creating a name for this sequence
            sequence_id = f">{species_name}_{genbank_id}"
            
            ## storing the sequences
            date_string = datetime.now().strftime("%d%b%y").lower()

            if official_name in repeated_locations:
                start = int(feature.location.start)
                repeated_locations[official_name].append((start, feature_seq))
            else:
                if args.output:
                    folder = args.output if isinstance(args.output, str) else args.output[0]
                    filename = f"{folder}/{official_name}_{date_string}.fasta"
                else:
                    filename = f"{official_name}_{date_string}.fasta"
                
                with open(filename, "a") as output_handle:
                    output_handle.write(f"{sequence_id}\n{feature_seq}\n")

            ## removing the found feature from the missing features list
            missing_features.discard(official_name)

            ## storing the accession number of the feature in the dictionary
            features_found[official_name] = genbank_id

    if args.repeated:
        ## processing the features with more than once copy
        for repeated, locations in repeated_locations.items():
            if len(locations) > 0:
                if len(locations) == 1:
                    start, feature_seq = locations[0]
                    
                    if args.output:
                        folder = args.output if isinstance(args.output, str) else args.output[0]
                        filename = f"{folder}/{repeated}_part1_{date_string}.fasta"
                    else:
                        filename = f"{repeated}_part1_{date_string}.fasta"

                    with open(filename, "a") as output_handle:
                        output_handle.write(f"{sequence_id}\n{feature_seq}\n")
                elif len(locations) == 2:
                    start1, seq1 = locations[0]
                    start2, seq2 = locations[1]

                    if args.output:
                        folder = args.output if isinstance(args.output, str) else args.output[0]
                        filename1 = f"{folder}/{repeated}_part1_{date_string}.fasta"
                        filename2 = f"{folder}/{repeated}_part2_{date_string}.fasta"
                    else:
                        filename1 = f"{repeated}_part1_{date_string}.fasta"
                        filename2 = f"{repeated}_part2_{date_string}.fasta"

                    if start1 < start2:
                        with open(filename1, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq1}\n")
                        
                        with open(filename2, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq2}\n")
                    else:
                        with open(filename1, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq2}\n")
                            
                        with open(filename2, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq1}\n")

    if "all" in args.features:
            ## if any feature is missing, print a warning
            for missing_feature in missing_features:
                logging.warning(f"Warning: '{missing_feature}' not found for {species_name} ({genbank_id})")

    for missing_feature in missing_features:
        if missing_feature in args.features:
            logging.warning(f"Warning: '{missing_feature}' not found for {species_name} ({genbank_id})")

    return species_name, features_found

synonym_dict = {
    "ATP6": ["ATP6", "atp6", "ATPase6", "ATPase 6", "ATP synthase 6", 
             "ATP synthase subunit 6", "ATP synthase F0 subunit 6", 
             "ATPase subunit 6", "MT-ATP6", "mt-Atp6", "mt-atp6",
             "F0-ATP synthase subunit6", "atpase6", "atpase 6", "Atp6",
             "ATP sythase subunit 6", "AT6", "MTATP6"],
    "ATP8": ["ATP8", "atp8", "ATPase8", "ATPase 8", "ATP synthase 8", 
             "ATP synthase subunit 8", "ATP synthase F0 subunit 8", 
             "ATPase subunit 8", "MT-ATP8", "mt-Atp8", "mt-atp8", "Atp8",
             "F0-ATP synthase subunit8", "atpase8", "atpase 8",
             "ATP sythase subunit 8", "AT8", "MTATP8"],
    "COI": ["cytochrome c oxidase subunit 1", "COI", "COX1", "cox1", "CO1", 
            "COXI", "cytochrome c oxidase subunit I", "COX-I", "coi", 
            "MT-CO1", "mt-Co1", "mt-co1", "cytochrome oxidase c subunit 1",
            "cytochrome oxidase subunit 1", "cytochrome oxidase subunit1",
            "Cytochrome c oxidase subunit1", "CO I", "coi", "co1", "coI",
            "coxI", "cytochrome oxidase I", "cytochrome oxidase subunit I",
            "cox I", "Cox1"],
    "COII": ["cytochrome c oxidase subunit 2", "COII", "COX2", "cox2", 
             "COXII", "CO2", "cytochrome c oxidase subunit II", "COX-II", 
             "MT-CO2", "mt-Co2", "mt-co2", "cytochrome oxidase subunit 2", 
             "coxII", "cytochrome oxidase subunit 2", 
             "cytochrome oxidase subunit2", "Cytochrome c oxidase subunit2",
             "CO II", "coii", "co2", "coII", "cytochrome oxidase II", 
             "cytochrome oxidase subunit II", "cox II", "Cox2"],
    "COIII": ["cytochrome c oxidase subunit 3", "COIII", "COX3", "cox3", "CO3", 
              "COXIII", "cytochrome c oxidase subunit III", "COX-III", "MT-CO3",
               "mt-Co3", "mt-co3", "cytochrome oxidase c subunit 3", "coxIII",
               "cytochrome oxidase subunit 3", "cytochrome oxidase subunit3",
               "Cytochrome c oxidase subunit3", "CO III", "coiii",
               "cytochrome c oxidase subunit3", "co3", "coIII",
               "cytochrome oxidase III", "cytochrome oxidase subunit III", 
               "CO3 subunit 3", "cox III", "Cox3"],
    "CYTB": ["CYTB", "cytb", "cob", "Cyt B", "Cyt b", "Cytb", "cyt b", "Cb", 
             "cytochrome b", "CYB", "cytB", "cyb", "MT-CYB", "mt-cyb", "cyt-B",
             "mt-Cytb", "cytochorome b", "Cytochrome b", "ctyb", "COB"],
    "ND1": ["ND1", "NADH1", "NADH dehydrogenase subunit 1", "nad1", "nd1", 
            "MT-ND1", "mt-Nd1", "mt-nd1", "NADH-1", "MTND1", "nadh1",
            "NADH dehydrogenase subunit1", "NAD1", "NADH subunit 1", 
            "NADH dehydrogenase 1", "NADH dehydrogenase subunit I"],
    "ND2": ["ND2", "NADH2", "NADH dehydrogenase subunit 2", "nad2", "nd2", 
            "NADH dehydrogenase subunit II", "NADH subunit 2", "MT-ND2", 
            "mt-Nd2", "mt-nd2", "NAD2", "NADH dehydrogenase subunit2", "nadh2",
            "NAD2", "NADH subunit 2"],
    "ND3": ["ND3", "NADH3", "NADH dehydrogenase subunit 3", "nad3", "nd3", 
            "MT-ND3", "mt-Nd3", "mt-nd3", "NADH dehydrogenase subunit3", 
            "nadh3", "NAD3", "NADH subunit 3", 
            "NADH dehydrogenase subunit III"],
    "ND4": ["ND4", "NADH4", "NADH dehydrogenase subunit 4", "nad4", "nd4", 
            "MT-ND4", "mt-Nd4", "mt-nd4", "NADH dehydrogenase subunit4", 
            "nadh4", "4", "NAD4", "NADH subunit 4",
            "NADH dehydrogenase subunit IV"],
    "ND4L": ["ND4L", "NADH4L", "NADH dehydrogenase subunit 4L", "nad4l", "nd4l",
             "nd4L", "MT-ND4L", "mt-Nd4l", "mt-nd4l", "nad4L", "nadh4L",
             "NADH dehydrogenase subunit 4 L", "NADH dehydrogenase subunit4L",
             "NAD4L", "NADH subunit 4L", "ND4l",
             "NADH dehydrogenase subunit IV L"],
    "ND5": ["ND5", "NADH5", "NADH dehydrogenase subunit 5", "nad5", "nd5", 
            "MT-ND5", "mt-Nd5", "mt-nd5", "nadh5", "nadh5", "NAD5",
            "NADH dehydrogenase subunit5", "NADH subunit 5",
            "NADH dehydrogenase subunit V"],
    "ND6": ["ND6", "NADH6", "NADH dehydrogenase subunit 6", "nad6", "nd6", 
            "MT-ND6", "mt-Nd6", "mt-nd6", "NADH dehydrogenase subunit6",
            "nadh6", "NAD6", "NADH subunit 6",
            "NADH dehydrogenase subunit VI"],
    "12S": ["12S ribosomal RNA", "s-rRNA", "small subunit ribosomal RNA", "12S",
            "12S rrn", "12SrRNA", "MTRNR1", "mt-Rnr1", "mt-rnr1", "mtrnr1", 
            "MT-RNR1", "SSU", "ssu", "rrn12", "ssu rRNA", "12 S ribosomal RNA",
            "small ribosomal RNA subunit RNA", "12S rRNA", "rRNA-12S",
            "12S ribosormal RNA", "12S-rRNA"],
    "16S": ["16S ribosomal RNA", "l-rRNA", "large subunit ribosomal RNA", "16S",
            "16rrn", "16S rrn", "16S rRNA", "16Srrn", "16SrRNA", "lsu", "LSU", 
            "lsu rRNA", "MTRNR2", "mt-Rnr2", "mt-rnr2", "MT-RNR2", "rrn16",
            "large ribosomal RNA subunit RNA", "16 S ribosomal RNA", "l-RNA",
            "16S-rRNA"],
    "tRNA-Ala": ["tRNA-Ala", "tRNA_Ala", "trnA", "trnA-ugc", "trnA TGC"],
    "tRNA-Arg": ["tRNA-Arg", "tRNA_Arg", "trnR", "trnR-ucg", "trnR TCG"],
    "tRNA-Asn": ["tRNA-Asn", "tRNA_Asn", "trnN", "trnN-guu", "trnN GTT"],
    "tRNA-Asp": ["tRNA-Asp", "tRNA_Asp", "trnD", "trnD-guc", "trnD GTC"],
    "tRNA-Cys": ["tRNA-Cys", "tRNA_Cys", "trnC", "trnC-gca", "trnC GCA"],
    "tRNA-Gln": ["tRNA-Gln", "tRNA_Gln", "trnQ", "trnQ-uug", "trnQ TTG"],
    "tRNA-Glu": ["tRNA-Glu", "tRNA_Glu", "trnE", "trnE-uuc", "trnE TTC"],
    "tRNA-Gly": ["tRNA-Gly", "tRNA_Gly", "trnG", "trnG-ucc", "trnG TCC"],
    "tRNA-His": ["tRNA-His", "tRNA_His", "trnH", "trnH-gug", "trnH GTG"],
    "tRNA-Ile": ["tRNA-Ile", "tRNA_Ile", "trnI", "trnI-gau", "trnI GAT"],
    "tRNA-Leu": ["tRNA-Leu", "tRNA_Leu", "trnL", "trnL-uag", "trnL TAG", 
                 "tRNA-Leu (CUN)", "tRNA-Leu (UUR)", "tRNA-Leu(CUN)", 
                 "tRNA-Leu(UUR)"],
    "tRNA-Lys": ["tRNA-Lys", "tRNA_Lys", "trnK", "trnK-uuu", "trnK TTT"],
    "tRNA-Met": ["tRNA-Met", "tRNA_Met", "trnM", "trnM-cau", "trnM CAT"],
    "tRNA-Phe": ["tRNA-Phe", "tRNA_Phe", "trnF-gaa", "trnF GAA"],
    "tRNA-Pro": ["tRNA-Pro", "tRNA_Pro", "trnP-ugg", "trnP TGG", "proline tRNA"],
    "tRNA-Ser": ["tRNA-Ser", "tRNA_Ser", "trnS", "trnS-uga", "trnS GCT", 
                 "tRNA-Ser (UCN)", "tRNA_Ser (AGY)", "tRNA-Ser(UCN)", 
                 "tRNA-Ser(AGY)"],
    "tRNA-Thr": ["tRNA-Thr", "tRNA_Thr", "trnT-ugu", "trnT TGT", 
                 "threonine tRNA"],
    "tRNA-Trp": ["tRNA-Trp", "tRNA_Trp", "trnW", "trnW-uca", "trnW TCA"],
    "tRNA-Tyr": ["tRNA-Tyr", "tRNA_Tyr", "trnY", "trnY-gua", "trnY GTA"],
    "tRNA-Val": ["tRNA-Val", "tRNA_Val", "trnV", "trnV-uac", "trnV TAC"]
}

def load_synonyms(file_path, synonym_dict):
    ## add new synonyms to the dictionary
    
    new_synonyms = {}
    with open(file_path, newline='', encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) > 1:
                key, synonyms = row[0].strip(), [s.strip() for s in row[1:]] 
                new_synonyms[key] = synonyms

    for key, synonyms in new_synonyms.items():
        if key in synonym_dict:
            synonym_dict[key] = list(set(synonym_dict[key] + synonyms)) 
        else:
            synonym_dict[key] = list(set(synonyms))

    return synonym_dict


def main():
    args = get_args()

    Entrez.email = args.email
    
    ## configuring log file
    date_string = datetime.now().strftime("%d%b%y").lower()
    time_string = datetime.now().strftime("%Hh%M").lower()

    if args.output:
        folder = args.output if isinstance(args.output, str) else args.output[0]
        filename = f"{folder}/annopython_{date_string}_{time_string}.log"
    else:
        filename = f"annopython_{date_string}_{time_string}.log"

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
        global synonym_dict 
        synonym_dict = load_synonyms(args.add_synonyms, synonym_dict)
        logging.info("Updated synonyms!")

    ## getting accession numbers to download
    accession_ids = read_accession_file(args.accession)

    data_accession = []

    if "all" in args.features:
        feature_names = set(synonym_dict.keys())
    else: 
        feature_names = args.features

    for accession in accession_ids:
        species_name = None  
        try:
            species_name, features_found = download_and_save_genes(accession, synonym_dict)
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