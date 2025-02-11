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

    parser.add_argument("-t", "--type",
                        type=str,
                        help="REQUIRED: Type of DNA (i.e. mtDNA_animals, mtDNA_plants, rDNA, cpDNA, other).",
                        required=False)

    parser.add_argument("-f", "--features",
                        type=str,
                        help="CSV file listing the features to download.",
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
                prev = official_name_prev
                repeated_locations[official_name].append((prev, feature_seq))
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

        official_name_prev = official_name

    if args.repeated:
    ## processing the features with more than one copy
        for repeated, locations in repeated_locations.items():
            if len(locations) > 0:
                if len(locations) == 1:
                    prev, feature_seq = locations[0]
                    
                    if args.output:
                        folder = args.output if isinstance(args.output, str) else args.output[0]
                        filename = f"{folder}/{repeated}_part1_{date_string}.fasta"
                    else:
                        filename = f"{repeated}_part1_{date_string}.fasta"

                    with open(filename, "a") as output_handle:
                        output_handle.write(f"{sequence_id}\n{feature_seq}\n")
                else:
                    if args.output:
                        folder = args.output if isinstance(args.output, str) else args.output[0]
                    else:
                        folder = ""
                    
                    for i, (prev, seq) in enumerate(locations, start=1):
                        filename = f"{folder}/{repeated}_after_{prev}_{date_string}.fasta" if folder else f"{repeated}_after_{prev}_{date_string}.fasta"
        
                        with open(filename, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq}\n")

    if "all" in args.features:
        ## if any feature is missing, print a warning
        for missing_feature in missing_features:
            logging.warning(f"Warning: '{missing_feature}' not found for {species_name} ({genbank_id})")

    for missing_feature in missing_features:
        if missing_feature in args.features:
            logging.warning(f"Warning: '{missing_feature}' not found for {species_name} ({genbank_id})")

    return species_name, features_found

synonym_dict_mtDNA_anim = {
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

synonym_dict_mtDNA_plan = {
    "26S": ["26S ribosomal RNA", "rrn26", "26S", "26S rrn", "26SrRNA", 
            "26 S ribosomal RNA", "26S rRNA", "rRNA-26S", "26S ribosormal RNA",
            "26S-rRNA"],
    "ATP1": ["atp1", "Atp1", "ATPase subunit A", "ATPase subunit 1", "ATP synthase F1 subunit 1", "ATP synthase F1 subunit alpha"],
    "ATP4": ["atp4", "Atp4", "ATPase subunit 4"],
    "ATP6": ["ATP6", "atp6", "ATPase6", "ATPase 6", "ATP synthase 6", 
             "ATP synthase subunit 6", "ATP synthase F0 subunit 6", 
             "ATPase subunit 6", "MT-ATP6", "mt-Atp6", "mt-atp6",
             "F0-ATP synthase subunit6", "atpase6", "atpase 6", "Atp6",
             "ATP sythase subunit 6", "AT6", "MTATP6"],
    "ATP6-1": ["atp6-1"],
    "ATP6-2": ["atp6-2"],
    "ATP8": ["ATP8", "atp8", "ATPase8", "ATPase 8", "ATP synthase 8", 
             "ATP synthase subunit 8", "ATP synthase F0 subunit 8", 
             "ATPase subunit 8", "MT-ATP8", "mt-Atp8", "mt-atp8", "Atp8",
             "F0-ATP synthase subunit8", "atpase8", "atpase 8",
             "ATP sythase subunit 8", "AT8", "MTATP8"],
    "ATP9": ["atp9", "Atp9", "ATPase subunit 9", "ATP synthase F0 subunit 9"],
    "ccmB": ["ccmB", "ccb2", "ccb206", "cytochrome c biogenesis ccmB", 
             "cytochrome c biogenesis B"],
    "ccmC": ["ccmC", "CcmC", "cytochrome c biogensis protein C", 
             "cytochrome c biogenesis ccmC", "cytochrome c biogenesis C"],
    "ccmFC": ["ccmFC", "CcmFC", "cytochrome c biogenesis ccmF", 
              "cytochrome c biogenesis Fc", "CcmFc"],
    "ccmFn": ["ccmFn", "CcmFn", "cytochrome c biogenesis Fn", "ccmFN"],
    "ccmFN1": ["ccmFN1", "cb6n1", "ccb382", 
               "cytochrome c biogenesis protein FN1", "ccmFN-1"],
    "ccmFN2": ["ccmFN2", "ccb203", "ccb6n2", 
               "cytochrome c biogenesis protein FN2"],
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
    "COII-1": ["cox2-1", "cytochrome c oxidase subunit 2-1"],
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
             "mt-Cytb", "cytochorome b", "Cytochrome b", "ctyb", "COB", 
             "apocytochrome B", "apocytochrome b"],
    "L16": ["rpl16", "ribosomal protein L16"],
    "matR": ["matR", "MatR", "maturase protein", "mat-r", "mat-R"],
    "ND1": ["ND1", "NADH1", "NADH dehydrogenase subunit 1", "nad1", "nd1", 
            "MT-ND1", "mt-Nd1", "mt-nd1", "NADH-1", "MTND1", "nadh1",
            "NADH dehydrogenase subunit1", "NAD1", "NADH subunit 1", 
            "NADH dehydrogenase 1", "NADH dehydrogenase subunit I", "Nad1"],
    "ND2": ["ND2", "NADH2", "NADH dehydrogenase subunit 2", "nad2", "nd2", 
            "NADH dehydrogenase subunit II", "NADH subunit 2", "MT-ND2", 
            "mt-Nd2", "mt-nd2", "NAD2", "NADH dehydrogenase subunit2", "nadh2",
            "NAD2", "NADH subunit 2"],
    "ND3": ["ND3", "NADH3", "NADH dehydrogenase subunit 3", "nad3", "nd3", 
            "MT-ND3", "mt-Nd3", "mt-nd3", "NADH dehydrogenase subunit3", 
            "nadh3", "NAD3", "NADH subunit 3", 
            "NADH dehydrogenase subunit III", "Nad3"],
    "ND4": ["ND4", "NADH4", "NADH dehydrogenase subunit 4", "nad4", "nd4", 
            "MT-ND4", "mt-Nd4", "mt-nd4", "NADH dehydrogenase subunit4", 
            "nadh4", "4", "NAD4", "NADH subunit 4",
            "NADH dehydrogenase subunit IV", "Nad4"],
    "ND4L": ["ND4L", "NADH4L", "NADH dehydrogenase subunit 4L", "nad4l", "nd4l",
             "nd4L", "MT-ND4L", "mt-Nd4l", "mt-nd4l", "nad4L", "nadh4L",
             "NADH dehydrogenase subunit 4 L", "NADH dehydrogenase subunit4L",
             "NAD4L", "NADH subunit 4L", "ND4l",
             "NADH dehydrogenase subunit IV L", "Nad4L"],
    "ND5": ["ND5", "NADH5", "NADH dehydrogenase subunit 5", "nad5", "nd5", 
            "MT-ND5", "mt-Nd5", "mt-nd5", "nadh5", "nadh5", "NAD5",
            "NADH dehydrogenase subunit5", "NADH subunit 5",
            "NADH dehydrogenase subunit V"],
    "ND6": ["ND6", "NADH6", "NADH dehydrogenase subunit 6", "nad6", "nd6", 
            "MT-ND6", "mt-Nd6", "mt-nd6", "NADH dehydrogenase subunit6",
            "nadh6", "NAD6", "NADH subunit 6",
            "NADH dehydrogenase subunit VI", "Nad6"],
    "ND7": ["nad7", "NADH dehydrogenase subunit 7", "Nad7"],
    "ND9": ["ND9", "NADH9", "NADH dehydrogenase subunit 9", "nad9", "nd9", 
            "MT-ND9", "mt-Nd9", "mt-nd9", "nadh9", "nadh9", "NAD9",
            "NADH dehydrogenase subunit9", "NADH subunit 9"],
    "PsaA": ["psaA"],
    "rbcL": ["rbcl", "rbcL"],
    "rpl2": ["rpl2", "ribosomal protein L2"],
    "rpl5": ["rpl5", "ribosomal protein L5"],
    "rpl10": ["rpl10", "RPL10", "RpL10", "Rpl10", "ribosomal protein L10"],
    "rpoB": ["rpoB", "RNA polymerase subunit beta", 
             "RNA polymerase beta subunit"],
    "rps3": ["rps3", "ribosomal protein S3"],
    "rps4": ["rps4", "ribosomal protein S4"],
    "rrn5": ["rrn5", "5S ribosomal RNA", "RRN5"],
    "rps7": ["rps7", "ribosomal protein S7", "RPS7"],
    "rps12": ["rps12", "ribosomal protein S12", "RPS12", "RpS12", "Rps12"],
    "rps14": ["rps14", "ribosomal protein S14", "RPS14", "Rps14"],
    "rrn18": ["rrn18", "RRN18", "18S ribosomal RNA"],
    "S1": ["rps1", "ribosomal protein S1", "RPS1"],
    "S10": ["ribosomal protein S10", "rps10", "RPS10", "Rps10"],
    "S13": ["rps13", "ribosomal protein S13", "RPS13", "Rps13"],
    "S19": ["rps19", "ribosomal protein S19", "RPS19", "Rps19"],
    "SDH3": ["sdh3", "succinate dehydrogenase subunit 3", "SDH3"],
    "tatC": ["tatC", "twin arginine translocation", 
             "twin-arginine translocase subunit TatC", "TATC"],
    "tRNA-Arg": ["tRNA-Arg", "tRNA_Arg", "trnR", "trnR(TCT)"],
    "tRNA-Asn": ["tRNA-Asn", "tRNA_Asn", "trnN", "trnN-guu", "trnN GTT",
                 "trnN-GUU", "trnN(GTT)", "trnN(guu)"],
    "tRNA-Asp": ["tRNA-Asp", "tRNA_Asp", "trnD", "trnD-guc", "trnD GTC",
                 "trnD-GUC", "trnD(GTC)", "trnD(guc)"],
    "tRNA-Cys": ["tRNA-Cys", "tRNA_Cys", "trnC", "trnC-gca", "trnC GCA",
                 "trnC-GCA", "trnC(GCA)", "trnC(gca)"],
    "tRNA-Gln": ["tRNA-Gln", "tRNA_Gln", "trnQ", "trnQ-uug", "trnQ TTG",
                 "trnQ-UUG", "trnQ(TTG)", "trnQ(uug)"],
    "tRNA-Glu": ["tRNA-Glu", "tRNA_Glu", "trnE", "trnE-uuc", "trnE TTC",
                 "trnE-UUC", "trnE(TTC)", "trnE(uuc)"],
    "tRNA-Gly": ["tRNA-Gly", "tRNA_Gly", "trnG", "trnG-ucc", "trnG TCC"],
    "tRNA-His": ["tRNA-His", "tRNA_His", "trnH", "trnH-gug", "trnH GTG",
                 "trnH-GUG", "trnH(GTG)", "trnH(gug)"],
    "tRNA-Ile": ["tRNA-Ile", "tRNA_Ile", "trnI", "trnI-CAU", "trnI(CAT)",
                 "trnI(cau)"],
    "tRNA-Leu": ["tRNA-Leu", "tRNA_Leu", "trnL", "trnL-CAA"],
    "tRNA-Lys": ["tRNA-Lys", "tRNA_Lys", "trnK", "trnK-uuu", "trnK TTT",
                 "trnK-UUU", "trnK(TTT)", "trnK(uuu)"],
    "tRNA-fMet": ["tRNA-fMet", "trnfM-CAU", "trnfM(CAT)", "trnfM(cau)"],
    "tRNA-Met": ["tRNA-Met", "tRNA_Met", "trnM", "trnM-cau", "trnM CAT",
                 "tRNA-met", "trnM-CAU", "trnM(CAT)"],
    "tRNA-Pro": ["tRNA-Pro", "tRNA_Pro", "trnP-ugg", "trnP TGG", "proline tRNA",
                 "trnP-UGG", "trnP", "trnP(TGG)", "trnP(ugg)"],
    "tRNA-Phe": ["tRNA-Phe", "tRNA_Phe", "trnF-gaa", "trnF GAA", "trnF",
                 "trnF(GAA)", "trnF(gaa)"],
    "tRNA-Ser": ["tRNA-Ser", "tRNA_Ser", "trnS", "trnS(GGA)", "trnS(GCT)", 
                 "trnS(uga)", "trnS(gcu)", "trnS(gga)"],
    "tRNA-Trp": ["tRNA-Trp", "tRNA_Trp", "trnW", "trnW-CCA", "trnW(cca)"],
    "tRNA-Tyr": ["tRNA-Tyr", "tRNA_Tyr", "trnY", "trnY-gua", "trnY GTA",
                 "trnY-GUA", "trnY(GTA)", "trnY(gua)"],
    "ycf1": ["ycf1"],
    "ycf2": ["ycf2"]
}

synonym_dict_rDNA = {
    "rRNA_18S": ["18S ribosomal RNA"],
    "ITS1": ["internal transcribed spacer 1", 
             "internal transcribed spacer ITS1"],
    "rRNA_5_8S": ["5.8S ribosomal RNA"],
    "ITS2": ["internal transcribed spacer 2", 
             "internal transcribed spacer ITS2"],
    "rRNA_28S": ["28S ribosomal RNA"]
}

synonym_dict_cpDNA = {
    "rrn4_5": ["rrn4.5", "4.5S ribosomal RNA"],
    "rrn5": ["rrn5", "5S ribosomal RNA"],
    "rrn16": ["rrn16", "16S ribosomal RNA"],
    "rrn23": ["rrn23", "23S ribosomal RNA", "Rrn23"],
    "accD": ["accD", "AccD", "acetyl-CoA carboxylase carboxyltransferase", 
             "acetyl-CoA carboxylase carboxyltransferase beta subunit", 
             "acetyl-CoA carboxylase, carboxyl transferase subunit beta", 
             "beta subunit of acetyl-CoA carboxylase"],
    "atpA": ["atpA", "ATP synthase CF1 alpha subunit"],
    "atpB": ["atpB", "ATP synthase CF1 beta subunit", "ATPase beta chain"],
    "atpE": ["atpE", "ATP synthase CF1 epsilon subunit", 
             "ATP synthase epsilon chain"],
    "atpF": ["atpF", "ATP synthase CF0 subunit I"],
    "atpH": ["atpH", "ATP synthase CF0 subunit III", "ATPase III subunit", 
             "ATP synthase subunit III"],
    "atpI": ["atpI", "ATP synthase CF0 subunit IV"],
    "ccsA": ["ccsA", "ycf5", "YCF5", "CcsA"],
    "cemA": ["cemA", "CemA", "chloroplast envelope membrane protein", 
             "heme-binding protein", "envelope membrane protein", 
             "chloroplast envelope protein", "ycf10"],
    "clpP": ["clpP", "clp protease proteolytic subunit", "ClpP", 
             "ATP-dependent protease proteolytic subunit", 
            "clp protease proteolytic subunit", "clp protease proteolytic"],
    "matK": ["matK", "maturase", "mat", "maturase K", "MATK", "matk"],
    "ndhA": ["ndhA", "NdhA"],
    "ndhB": ["ndhB", "NADH dehydrogenase subunit 2", 
             "NADH-plastoquinone oxidoreductase subunit 2", "NdhB"],
    "ndhC": ["ndhC", "NADH dehydrogenase subunit 3", 
             "NADH-plastoquinone oxidoreductase subunit 3", "NdhC"],
    "ndhD": ["ndhD", "NdhD", "NADH dehydrogenase subunit 4", "ndh4", 
             "NADH-plastoquinone oxidoreductase subunit 4"],
    "ndhE": ["ndhE", "NdhE", "NADH dehydrogenase subunit 4L", "ndh4L", 
             "NADH-plastoquinone oxidoreductase subunit 4L"],
    "ndhF": ["ndhF", "NADH dehydrogenase subunit 5", 
             "NADH-plastoquinone oxidoreductase subunit 5", "NdhF"],
    "ndhG": ["ndhG", "NdhG", "NADH dehydrogenase subunit 6", 
             "NADH-plastoquinone oxidoreductase subunit 6"],
    "ndhH": ["ndhH", "NdhH", "NADH-plastoquinone oxidoreductase subunit 7", 
             "NADH dehydrogenase subunit 7"],
    "ndhI": ["ndhI", "Ndhl"],
    "ndhJ": ["ndhJ", "NdhJ", "NADH dehydrogenase subunit J", 
             "NADH-plastoquinone oxidoreductase subunit J"],
    "ndhK": ["ndhK", "NdhK", "NADH dehydrogenase subunit K", 
             "NADH-plastoquinone oxidoreductase subunit K"],
    "pafI": ["pafI", "ycf3", "photosystem I assembly factor I", 
             "putative chloroplast RF34", "Ycf3"],
    "pafII": ["pafII", "ycf4", "photosystem I assembly protein ycf4", 
              "photosystem I assembly factor II"],
    "petA": ["petA", "PetA", "cytochrome f", "apocytochrome f"],
    "petB": ["petB", "cytochrome b6", "PetB"],
    "petD": ["petD", "cytochrome b6/f complex subunit IV", 
             "cytochrome b6/f complex subunit 4", "PetD", 
             "subunit IV of cytochrome b6/f complex"],
    "petG": ["petG", "PetG", "cytochrome b6/f complex subunit V", 
             "cytochrome B6-F complex subunit 5"],
    "petL": ["petL", "PetL", "cytochrome b6/f complex subunit VI", 
             "subunit VI of cytochrome b6/f complex", "ycf7"],
    "petN": ["petN", "cytochrome b6/f complex subunit VIII", 
             "cytochrome b6f complex subunit VIII", "PetN", "ycf6"],
    "psaA": ["psaA", "photosystem I P700 apoprotein A1", 
             "PSI P700 apoprotein A1", "PsaA"],
    "psaB": ["psaB", "photosystem I P700 apoprotein A2", 
             "PSI P700 apoprotein A2", 
             "photosystem I P700 chlorophyll a apoprotein A2", "PsAB", "PsaB"],
    "psaC": ["psaC", "photosystem I iron-sulfur center", 
             "photosystem I subunit VII", "subunit VII of photosystem I", 
             "PsaC", "PSI C protein", 
             "photosystem I iron-sulfur center subunit VII"],
    "psaI": ["psaI", "PsaI", "photosystem I subunit VIII", 
             "photosystem I reaction center subunit VIII"],
    "psaJ": ["psaJ", "PSI J-protein", "photosystem I subunit IX", 
             "photosystem I reaction center subunit IX", 
             "photosystem I protein J"],
    "psbA": ["psbA", "photosystem II protein D1"],
    "psbB": ["psbB", "photosystem II 47 kDa protein", 
             "photosystem II P680 chlorophyll A apoprotein", "PsbB", 
             "photosystem II CP47 chlorophyll apoprotein"],
    "psb30": ["psb30", "photosystem II protein psb30", "Ycf12", "ycf12", 
              "photosystem II reaction centre protein ycf12"],
    "psbC": ["psbC", "photosystem II CP43 chlorophyll apoprotein"],
    "psbD": ["psbD", "photosystem II protein D2", "PsbD", "PSII D2 protein"],
    "psbE": ["psbE", "cytochrome b559 alpha subunit of photosystem", 
             "cytochrome b559 alpha chain", 
             "photosystem II cytochrome b559 alpha subunit", 
             "photosystem II protein V"],
    "psbF": ["psbF", "PsbF", "photosystem II cytochrome b559 beta subunit",
             "photosystem II protein VI", "PSII subunit VI"],
    "psbI": ["psbI", "photosystem II protein I"],
    "psbJ": ["psbJ", "PsbJ", "photosystem II protein J",
             "J protein of photosystem II"],
    "psbK": ["psbK", "photosystem II protein K", "PSII K protein"],
    "psbL": ["psbL", "psbl", "photosystem II protein L"],
    "psbH": ["psbH", "photosystem II phosphoprotein", 
             "photosystem II reaction center protein H", 
             "photosystem II protein H", "PsbH", "PSII 10 kDa phosphoprotein"],
    "psbM": ["psbM", "photosystem II protein M", "photosystem II M protein",
             "PsbM"],
    "psbN": ["psbN", "photosystem II protein N", 
             "photosystem II reaction center N protein", 
             "N protein of photosystem II", "PSII N protein", "PsbN"],
    "psbT": ["psbT", "photosystem II protein T", "PsbT", 
             "T protein of photosystem II", "ycf8"],
    "psbZ": ["psbZ", "PsbZ", "photosystem II protein Z",
             "photosystem II reaction center Z protein", "ycf9", "Ycf9"],
    "rbcL": ["rbcL", "rbcl", "ribulose 1,5-bisphosphate carboxylase/oxygenase",
             "RuBisCO large subunit", "large subunit of Rubisco"],
    "rpl12": ["rpl12", "ribosomal protein L2", "Ribosomal protein L2"],
    "rpl14": ["rpl14", "ribosomal protein L14"],
    "rpl16": ["rpl16", "ribosomal protein L16", "Rpl16", 
              "Ribosomal protein L16"],
    "rps15": ["rps15", "ribosomal protein S15", "Ribosomal protein S15", 
              "Rps15"],
    "rps19": ["rps19", "ribosomal protein S19", "Ribosomal protein S19", 
              "Rps19"],
    "rpl22": ["rpl22", "ribosomal protein L22", "Rpl22", 
              "Ribosomal protein L22"],
    "rpl12": ["rpl12", "ribosomal protein L2", "Ribosomal protein L2"],
    "rpl36": ["rpl36", "ribosomal protein L16"],
    "rpl20": ["rpl20", "ribosomal protein L20", "Rpl20"],
    "rpl32": ["rpl32", "ribosomal protein L32", "Rpl32", 
              "Ribosomal protein L32"],
    "rpl33": ["rpl33", "ribosomal protein L33"],
    "rps3": ["rps3", "ribosomal protein S3", "Rps3", "Ribosomal protein S3"],
    "rps7": ["rps7", "ribosomal protein S7", "Rps7", "Ribosomal protein S7"],
    "rps8": ["rps8", "ribosomal protein S8", "Rps8", "Ribosomal protein S8"],
    "rps18": ["rps18", "ribosomal protein S18", "Rps18"],
    "rpoA": ["rpoA", "RNA polymerase alpha subunit", "RpoA", 
             "RNA polymerase alpha chain", "RNA polymerase a-subunit"],
    "rpoB": ["rpoB", "RNA polymerase beta subunit"],
    "rpoC1": ["rpoC1", "RNA polymerase beta' subunit"],
    "rpoC2": ["rpoC2", "RNA polymerase beta'' subunit"],
    "rps11": ["rps11", "ribosomal protein S11", "Rps11"],
    "rps12": ["rps12", "ribosomal protein S12", "RPS12"],
    "rps14": ["rps14", "ribosomal protein S14"],
    "rps16": ["rps16", "Rps16", "ribosomal protein S16"],
    "rps2": ["rps2", "ribosomal protein S2"],
    "rps4": ["rps4", "Rps4", "RPS4", "ribosomal protein S4"],
    "trnA": ["trnA-UGC", "tRNA-Ala"],
    "trnC": ["trnC-GCA", "tRNA-Cys"],
    "trnD": ["trnD-GUC", "tRNA-Asp"],
    "trnE": ["trnE-UUC", "tRNA-Glu"],
    "trnF": ["trnF-GAA", "tRNA-Phe"],
    "trnG": ["trnG-UCC", "tRNA-Gly"],
    "trnH": ["trnH-GUG", "tRNA-His"],
    "trnI": ["trnI-CAU", "trnI-GAU", "trnI-GAT", "tRNA-Ile"],
    "trnK": ["trnK-UUU", "tRNA-Lys"],
    "trnL": ["trnL-UAA", "trnL-CAA", "trnL-UAG", "tRNA-Leu"],
    "trnM": ["trnM-CAU", "tRNA-Met"],
    "trnN": ["trnN-GUU", "tRNA-Asn"],
    "trnP": ["trnP-UGG", "tRNA-Pro"],
    "trnQ": ["trnQ-UUG", "tRNA-Gln"],
    "trnR": ["trnR-UCU", "trnR-ACG", "tRNA-Arg"],
    "trnS": ["trnS-GCU", "trnS-UGA", "trnS-GGA", "tRNA-Ser"],
    "trnT": ["trnT-GGU", "trnT-UGU", "tRNA-Thr"],
    "trnV": ["trnV-UAC", "trnV-GAC", "tRNA-Val"],
    "trnW": ["trnW-CCA", "tRNA-Trp"],
    "trnY": ["trnY-GUA", "tRNA-Tyr"],
    "trnfM": ["trnfM-CAU"],
    "ycf1": ["ycf1", "hypothetical chloroplast RF1", "Ycf1"],
    "ycf2": ["ycf2", "Ycf2"],
    "ycf15": ["ycf15", "Ycf15"]
}

def load_synonyms(file_path, syn_dict=None):
    ## load synonyms from a CSV file and update an existing dictionary (if provided). If no dictionary is provided, a new one is created.

    if syn_dict is None:
        syn_dict = {} 
    
    new_synonyms = {}
    with open(file_path, newline='', encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) > 1:
                key, synonyms = row[0].strip(), [s.strip() for s in row[1:]] 
                new_synonyms[key] = synonyms

    for key, synonyms in new_synonyms.items():
        if key in syn_dict:
            syn_dict[key] = list(set(syn_dict[key] + synonyms)) 
        else:
            syn_dict[key] = list(set(synonyms))

    return syn_dict

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