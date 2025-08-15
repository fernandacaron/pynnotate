import argparse
import json

from pynnotate.terminal.utils import load_config, merge_args_with_config

from pynnotate.main.main import begin_search

def get_args():
    parser = argparse.ArgumentParser(
        description="""
            Pynnotate is your go-to tool for efficient gene extraction and analysis. 
            It allows the user to specify a list of features to retrieve, handle synonyms for gene names, and manage features that  
            appear multiple times within a genome. 
            Example usage: --config config.yaml
        """,
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False)

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
    parser.add_argument("--filter_mode", 
                        required=True,
                        choices=["unconstrained", "flexible", "strict"],
                        default="unconstrained",
                        help=(
                            "REQUIRED: Filtering mode: "
                            "'unconstrained' = include all sequences, "
                            "'flexible' = allow multiple sequences per species if they have new genes (e.g. supermatrix), "
                            "'strict' = include only one sequence per species, prioritizing individuals with more genes\n\n"
                            "⚠️ In 'strict' mode, filtering is based on the genes available in the selected gene dictionary and/or any user-provided gene dictionary.\n"
                            "⚠️ When 'unconstrained' mode is used in combination with separate gene extraction, the downloaded sequences will correspond to the gene set in the selected and/or user-provided dictionary."
                            )
                        )
    parser.add_argument("--prioritize", 
                        action="store_true", 
                        help="Prioritize individual with more genes (mitochondrial)")
    parser.add_argument("--add_synonyms", 
                        type=json.loads, 
                        default={},
                        help=("Additional gene name synonyms in JSON format."
                            "Pynnotate already includes an internal dictionary of gene name synonyms to aid extraction. You can provide additional synonyms for genes not automatically recognized. We recommend running the program first to identify any unrecognized gene synonyms. Add any missing synonyms here to improve matching. "
                            "⚠️ ATTENTION: When selecting the genome type and adding synonyms, they will be incorporated into the internal dictionary for that specific genome type. However, if the genome type selected is 'Other', only the synonyms provided by the user will be used."
                        ),
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
                        help="REQUIRED: Output folder",
                        required=True,
                        default=".")
    parser.add_argument("-f", "--folder",
                        type=str,
                        help="Folder name to create inside output folder (will be created if it does not exist)",
                        required=False,
                        default=".")

    parser.add_argument("-h", "--help", 
                        action="store_true",
                        help="Show this help message and exit")

    args = parser.parse_args()

    if args.help:
        help_text = parser.format_help()
        pager = subprocess.Popen(['less', '-R'], stdin=subprocess.PIPE)
        pager.communicate(help_text.encode('utf-8'))
        exit(0)

    return args, parser


def main():

    import queue

    error_queue = queue.Queue()

    args, parser = get_args()
    config = load_config(args.config)
    final_config = merge_args_with_config(config, args, parser)

    selected_mode = final_config["filter_mode"]

    if selected_mode == "strict":
        prioritize_more_genes = True
        unique_species = False
    elif selected_mode == "flexible":
        prioritize_more_genes = False
        unique_species = True
    else: 
        prioritize_more_genes = False
        unique_species = False

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
        "unique": unique_species,
        "prioritize_more_genes": prioritize_more_genes,
        "unconstrained": final_config["filter_mode"] == "unconstrained",
        "logmissing": final_config["logmissing"],
        "org_type": final_config["type"],
        "add_synonyms": final_config["add_synonyms"],
        "additional": final_config["additional"]
    }
    
    begin_search(config)

if __name__ == "__main__":
    main()