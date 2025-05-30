import argparse
import json

from terminal.utils import load_config, merge_args_with_config

from main.pynnotate import begin_search

def get_args():
    parser = argparse.ArgumentParser(
        description="""
            Pynnotate is your go-to tool for efficient gene extraction and analysis. 
            It allows the user to specify a list of features to retrieve, handle synonyms for gene names, and manage features that  
            appear multiple times within a genome. 
            Example usage: --config config.yaml
        """)

    parser.add_argument("-c", "--config", type=str, 
                        help="REQUIRED: Path to configuration file",
                        required=True)
    parser.add_argument("-t", "--type",
                        type=str,
                        help="REQUIRED: Type of DNA (i.e. animal_mito, plant_mito, plant_chloro, other)",
                        required=False)
    parser.add_argument("-a", "--accession",
                        type=lambda s: [acc.strip() for acc in a.split(',')],
                        required=False,
                        help="Comma-separated list of GenBank accession numbers to search/download")
    parser.add_argument("-g", "--genes",
                        type=lambda s: [gene.strip() for gene in s.split(',')],
                        help="Comma-separated list of genes to search/download, e.g. COI,CYTB,ATP6",
                        required=False)
    parser.add_argument("--organism",
                        type=str,
                        help="Organism to search/download",
                        required=False)
    parser.add_argument("-p", "--publication",
                        type=str,
                        help="Publication term (title, authors, year)",
                        required=False)
    parser.add_argument("--additional",
                        type=str,
                        help="Any additional search terms, e.g. NOT sp.",
                        required=False)
    parser.add_argument("--mitochondrialgene", 
                        action="store_true", 
                        help="Refine your search terms to Mitochondrial gene")
    parser.add_argument("--mitogenome", 
                        action="store_true", 
                        help="Refine your search terms to Mitogenome")
    parser.add_argument("--chloroplast", 
                        action="store_true", 
                        help="Refine your search terms to Chloroplast")
    parser.add_argument("--annotated", 
                        action="store_true", 
                        help="Delete unannotated records")
    parser.add_argument("--header",
                        type=str,
                        help="Header fields (Genbank fields)",
                        required=False)
    parser.add_argument("--genbankid", 
                        action="store_true", 
                        help="Include GenBank ID in fasta header")
    parser.add_argument("--unique", 
                        action="store_true", 
                        help="Include only 1 individual per species")
    parser.add_argument("--add_synonyms", 
                        type=json.loads, 
                        default={},
                        help="Additional gene name synonyms in JSON format, e.g. '{\"CYTB\": [\"CYTOCHROME B\", \"CYTB\"]}'",
                        required=False)
    parser.add_argument("--minbp",
                        type=str,
                        help="Minimum sequence size (bp)",
                        required=False)
    parser.add_argument("--maxbp",
                        type=str,
                        help="Maximum sequence size (bp)",
                        required=False)
    parser.add_argument("--extraction", 
                        action="store_true", 
                        help="Extract all annotated genes separately")
    parser.add_argument("--overlap", 
                        action="store_true", 
                        help="Fix overlap between extracted genes")
    parser.add_argument("--logmissing", 
                        action="store_true", 
                        help="Generate log of missing genes per sample (useful for mitogenomes)")
    parser.add_argument("-e", "--email",
                        type=str,
                        required=False,
                        help="REQUIRED: User e-mail to access NCBI database")
    parser.add_argument("-o", "--output",
                        type=str,
                        help="Output folder (will be created if it does not exist)",
                        required=False,
                        default=".")

    args = parser.parse_args()

    return args, parser

def main():

    args, parser = get_args()
    config = load_config(args.config)
    final_config = merge_args_with_config(config, args, parser)

    extraction = config.get("extraction", False)

    config = {
        "email": final_config["email"],
        "title": final_config["publication"],
        "organisms": final_config["organism"],
        "genes": final_config["genes"],
        "fields": final_config["header"],
        "output_name": final_config["output"],
        "include_id": final_config["genbankid"],
        "mito": final_config["mitochondrialgene"],
        "mitogenome": final_config["mitogenome"],
        "ids_text": final_config["accession"],
        "chloroplast": final_config["chloroplast"],
        "min_len": final_config["minbp"],
        "max_len": final_config["maxbp"],
        "delete_unverified": final_config["annotated"],
        "extraction": final_config["extraction"],
        "overlap": final_config["overlap"],
        "unique": final_config["unique"],
        "logmissing": final_config["logmissing"],
        "org_type": final_config["type"],
        "add_synonyms": final_config["add_synonyms"],
        "additional": final_config["additional"]
    }
    
    begin_search(config)

if __name__ == "__main__":
    main()