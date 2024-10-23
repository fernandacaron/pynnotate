import logging
import traceback
import pandas as pd
from Bio import Entrez, SeqIO
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description = """Description...""")
    
    #parser.add_argument("-e", "--email",
    #                        required=True,
    #                        help="REQUIRED: User e-mail to acess NCBI database.")
        
    parser.add_argument("-a", "--accession",
                        required = True,
                        help = "REQUIRED: The full path to a file containing the accession numbers.")

    parser.add_argument("--f", "--features",
                        nargs = '+',
                        help = "List the features to download.",
                        required = False,
                        default = "all")

    parser.add_argument("--r", "--repeated",
                        nargs = '+',
                        help = "List of features with more than one copy.",
                        required = False)

    return parser.parse_args()

def read_accession_file(file_path):
    with open(file_path, "r") as file:
        accession_list = [line.strip() for line in file if line.strip()] \
    return accession_list

## finding the general name of a feature based on a synonym
def get_official_name(feature_name):
    for official_name, synonyms in synonym_dict.items():
        if feature_name in synonyms:
            return official_name
    return None

## função para baixar e salvar as sequências das features
def download_and_save_genes(genbank_id):
    
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
        official_name = get_official_name(gene_name)

        if args.features == "all":
            ## if the name is not a known synonym, skip feature
            if not official_name:
                logging.error(f"SYNONYM unknown for {gene_name}, skipping.")
                continue

            if args.features == "all" or gene_name in args.features:
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
                    filename = f"{official_name}_{date_string}.fasta"
                
                    with open(filename, "a") as output_handle:
                        output_handle.write(f"{sequence_id}\n{feature_seq}\n")

                ## removeing the found feature from the missing features list
                missing_features.discard(official_name)

                ## storing the accession number of the feature in the dictionary
                features_found[official_name] = genbank_id

        ## processing the features with more than once copy
        for repeated, locations in repeated_locations.items():
            if len(locations) > 0:
                if len(locations) == 1:
                    start, feature_seq = locations[0]
                    filename = f"{repeated}_part1_{date_string}.fasta"
                    with open(filename, "a") as output_handle:
                        output_handle.write(f"{sequence_id}\n{feature_seq}\n")
                elif len(locations) == 2:
                    start1, seq1 = locations[0]
                    start2, seq2 = locations[1]

                    if start1 < start2:
                        filename1 = f"{repeated}_part1_{date_string}.fasta"
                        with open(filename1, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq1}\n")
                        
                        filename2 = f"{repeated}_part2_{date_string}.fasta"
                        with open(filename2, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq2}\n")
                    else:
                        filename1 = f"{repeated}_part1_{date_string}.fasta"
                        with open(filename1, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq2}\n")
                        
                        filename2 = f"{repeated}_part2_{date_string}.fasta"
                        with open(filename2, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq1}\n")

        ## if any feature is missing, print a warning
        for missing_feature in missing_features:
            logging.warning(f"Warning: '{missing_feature}' not found for {species_name} ({genbank_id})")

        return species_name, features_found

    else: 

        for features in args.features:
            
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
            official_name = get_official_name(gene_name)

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
                filename = f"{official_name}_{date_string}.fasta"
                
                with open(filename, "a") as output_handle:
                    output_handle.write(f"{sequence_id}\n{feature_seq}\n")

            ## removeing the found feature from the missing features list
            missing_features.discard(official_name)

            ## storing the accession number of the feature in the dictionary
            features_found[official_name] = genbank_id

        ## processing the features with more than once copy
        for repeated, locations in repeated_locations.items():
            if len(locations) > 0:
                if len(locations) == 1:
                    start, feature_seq = locations[0]
                    filename = f"{repeated}_part1_{date_string}.fasta"
                    with open(filename, "a") as output_handle:
                        output_handle.write(f"{sequence_id}\n{feature_seq}\n")
                elif len(locations) == 2:
                    start1, seq1 = locations[0]
                    start2, seq2 = locations[1]

                    if start1 < start2:
                        filename1 = f"{repeated}_part1_{date_string}.fasta"
                        with open(filename1, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq1}\n")
                        
                        filename2 = f"{repeated}_part2_{date_string}.fasta"
                        with open(filename2, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq2}\n")
                    else:
                        filename1 = f"{repeated}_part1_{date_string}.fasta"
                        with open(filename1, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq2}\n")
                        
                        filename2 = f"{repeated}_part2_{date_string}.fasta"
                        with open(filename2, "a") as output_handle:
                            output_handle.write(f"{sequence_id}\n{seq1}\n")

        ## if any feature is missing, print a warning
        for missing_feature in missing_features:
            logging.warning(f"Warning: '{missing_feature}' not found for {species_name} ({genbank_id})")

        return species_name, features_found



def main():
    args = get_args()

    ## configuring log file
    date_string = datetime.now().strftime("%d%b%y").lower()
    time_string = datetime.now().strftime("%Hh%M").lower()
    logging.basicConfig(
        level = logging.INFO,
        format = '%(asctime)s - %(levelname)s - %(message)s',
        handlers = [
            logging.FileHandler(f"annotation_{date_string}_{time_string}.log"),
            logging.StreamHandler()
        ]
    )

    ## getting accession numbers to download
    accession_ids = read_accession_file(accession)

    ## synonym dictionary
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
        "tRNA_Ala": ["tRNA-Ala", "trnA", "trnA-ugc", "trnA TGC"],
        "tRNA_Arg": ["tRNA-Arg", "trnR", "trnR-ucg", "trnR TCG"],
        "tRNA_Asn": ["tRNA-Asn", "trnN", "trnN-guu", "trnN GTT"],
        "tRNA_Asp": ["tRNA-Asp", "trnD", "trnD-guc", "trnD GTC"],
        "tRNA_Cys": ["tRNA-Cys", "trnC", "trnC-gca", "trnC GCA"],
        "tRNA_Gln": ["tRNA-Gln", "trnQ", "trnQ-uug", "trnQ TTG"],
        "tRNA_Glu": ["tRNA-Glu", "trnE", "trnE-uuc", "trnE TTC"],
        "tRNA_Gly": ["tRNA-Gly", "trnG", "trnG-ucc", "trnG TCC"],
        "tRNA_His": ["tRNA-His", "trnH", "trnH-gug", "trnH GTG"],
        "tRNA_Ile": ["tRNA-Ile", "trnI", "trnI-gau", "trnI GAT"],
        "tRNA_Leu": ["tRNA-Leu", "trnL", "trnL-uag", "trnL TAG", "tRNA-Leu (CUN)",
                     "tRNA-Leu (UUR)", "tRNA-Leu(CUN)", "tRNA-Leu(UUR)"],
        "tRNA_Lys": ["tRNA-Lys", "trnK", "trnK-uuu", "trnK TTT"],
        "tRNA_Met": ["tRNA-Met", "trnM", "trnM-cau", "trnM CAT"],
        "tRNA_Phe": ["tRNA-Phe", "trnF-gaa", "trnF GAA"],
        "tRNA_Pro": ["tRNA-Pro", "trnP-ugg", "trnP TGG", "proline tRNA"],
        "tRNA_Ser": ["tRNA-Ser", "trnS", "trnS-uga", "trnS GCT", "tRNA-Ser (UCN)",
                     "tRNA-Ser (AGY)", "tRNA-Ser(UCN)", "tRNA-Ser(AGY)"],
        "tRNA_Thr": ["tRNA-Thr", "trnT-ugu", "trnT TGT", "threonine tRNA"],
        "tRNA_Trp": ["tRNA-Trp", "trnW", "trnW-uca", "trnW TCA"],
        "tRNA_Tyr": ["tRNA-Tyr", "trnY", "trnY-gua", "trnY GTA"],
        "tRNA_Val": ["tRNA-Val", "trnV", "trnV-uac", "trnV TAC"]
    }


## criando tabela para armazenar os números de acesso
data_accession = []

## criado um conjunto para rastrear os nomes das features para colocar na tabela
feature_names = set(synonym_dict.keys())

## loop para processar todos os números de acesso
for accession in accession_ids:
    ## usar o nome da espécie como chave no dicionário
    species_name = None  
    try:
        species_name, features_found = download_and_save_genes(accession)
        ## linha na tabela com o nome da espécie e o accession
        row_data = [species_name, accession]

        ## preencher as colunas com os dados das features
        for official_name in feature_names:
            ## se a feature foi encontrada, adicionar o número de acesso
            if official_name in features_found:
                row_data.append(features_found[official_name])
            else:
                row_data.append("NA")  ## senão, adicionar NA

        ## adicionar a linha à lista de dados da tabela
        data_accession.append(row_data)

    except Exception as e:
        error_info = traceback.format_exc()
        logging.error(f"ERRO ao processar {accession}: {e}")
        logging.error("Detalhes do erro:", error_info)

## criando um dataframe usando o dicionário data_accession
feature_names_prefix = [f"Acc_{name}" for name in feature_names]
data_accession_df = pd.DataFrame(data_accession, 
                                 columns = ["Sequence_names", 
                                 "Mitogenome_accession"] + 
                                 feature_names_prefix)

## salvando a tabela em um arquivo csv
data_accession_df.to_csv(f"spp_accession_{date_string}.csv",
                         index = False)


