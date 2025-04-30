import logging
import pandas as pd
from Bio import Entrez, SeqIO
from datetime import datetime
import csv
import os

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

    official_name_prev = "NA"

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
                        filename = f"{folder}/{repeated}_{date_string}.fasta"
                    else:
                        filename = f"{repeated}_{date_string}.fasta"

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
