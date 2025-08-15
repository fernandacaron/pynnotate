import argparse
import json
import pydoc
import sys
import textwrap

from pynnotate.terminal.utils import load_config, merge_args_with_config
#from terminal.utils import load_config, merge_args_with_config

from pynnotate.main.main import begin_search
#from main.main import begin_search

class SuperCustomHelpFormatter(argparse.HelpFormatter):
    def __init__(self, prog, indent_increment=2, max_help_position=24, width=70):
        super().__init__(prog, indent_increment, max_help_position, width)
        self.width = width
        self.indent = indent_increment

    def _format_info_text(self, text):
        if text:
            indent_prefix = " " * self.indent
            formatted_lines = []
            for line in text.split('\n'):
                wrapped = textwrap.fill(
                    line,
                    width=self.width,
                    initial_indent=indent_prefix,
                    subsequent_indent=indent_prefix
                )
                formatted_lines.append(wrapped)
            return '\n'.join(formatted_lines)
        return text

    def _format_text(self, text):
        return self._format_info_text(text)

    def _format_action_invocation(self, action):
        if not action.option_strings:
            return action.dest
        
        parts = []
        if action.nargs == 0:
            parts.extend(action.option_strings)
        else:
            parts.extend(action.option_strings)
        return ', '.join(parts)

    def _format_action(self, action):
        if action.required:
            required_string = "REQUIRED: "
        else:
            required_string = ""
            
        help_text = action.help if action.help else ""
        help_text = f"{required_string}{help_text}"
        
        help_lines = []
        for line in help_text.split('\n'):
            wrapped_lines = textwrap.wrap(
                line,
                width=self.width,
                subsequent_indent=' ' * (self.indent * 2)
            )
            help_lines.extend(wrapped_lines or [''])

        invocation = self._format_action_invocation(action)
        formatted_help = '\n' + '   ' * (self.indent * 2) + ('\n' + '   ' * (self.indent * 2)).join(help_lines)
        
        return f"{' ' * self.indent}{invocation}{formatted_help}\n"

    def _format_text(self, text):
        if text:
            lines = []
            for line in text.split('\n'):
                wrapped_lines = textwrap.wrap(
                    line,
                    width=self._width,
                    subsequent_indent='  '
                )
                lines.extend(wrapped_lines or [''])
            return '\n'.join(lines) + '\n'
        return text

def get_args():
    parser = argparse.ArgumentParser(
        description="""
            Pynnotate is your go-to tool for efficient gene extraction
            and analysis.\n 
            It allows the user to specify a list of features to
            retrieve, handle synonyms for gene names, and manage
            features that appear multiple times within a genome. \n 
            Example usage: pynnotate --config config.yaml
        """,
        formatter_class=SuperCustomHelpFormatter)

    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    required.add_argument("-c", "--config", type=str, 
                        help="Path to configuration file",
                        required=True)
    optional.add_argument("-t", "--type",
                        type=str,
                        required=False,
                        help="Type of DNA (i.e. animal_mito, plant_mito, plant_chloro, other)")
    optional.add_argument("--filter_mode", 
                        required=False,
                        choices=["unconstrained", "flexible", "strict"],
                        default="unconstrained",
                        help=(
                            "Filtering mode: \n"
                            "  'unconstrained' = include all sequences \n"
                            "  'flexible' = allow multiple sequences per species if they have new genes (e.g. supermatrix)\n"
                            "  'strict' = include only one sequence per species, prioritizing individuals with more genes\n\n"
                            "  ⚠️  In 'strict' mode, filtering is based on the genes available in the selected gene dictionary and/or any user-provided gene dictionary\n"
                            "  ⚠️  When 'unconstrained' mode is used in combination with separate gene extraction, the downloaded sequences will correspond to the gene set in the selected and/or user-provided dictionary"
                            )
                        )
    optional.add_argument("-e", "--email",
                        type=str,
                        required=False,
                        help="User e-mail to access NCBI database")
    optional.add_argument("-o", "--output",
                        type=str,
                        help="Output folder",
                        required=False,
                        default=".")
    optional.add_argument("--accession",
                        type=lambda s: [acc.strip() for acc in a.split(',')],
                        required=False,
                        help="Comma-separated list of GenBank accession numbers to search/download")
    optional.add_argument("--genes",
                        type=lambda s: [gene.strip() for gene in s.split(',')],
                        help="Comma-separated list of genes to search/download (e.g., COI, CYTB, ATP6)",
                        required=False)
    optional.add_argument("--organism",
                        type=str,
                        help="Organism to search/download",
                        required=False)
    optional.add_argument("--publication",
                        type=str,
                        help="Publication term (title, authors, year)",
                        required=False)
    optional.add_argument("--additional",
                        type=str,
                        help="Any additional search terms (e.g., NOT sp.)",
                        required=False)
    optional.add_argument("--mitochondrialgene", 
                        action="store_true",
                        help="Refine your search terms to Mitochondrial gene")
    optional.add_argument("--mitogenome",
                        action="store_true", 
                        help="Refine your search terms to Mitogenome")
    optional.add_argument("--chloroplast",
                        action="store_true", 
                        help="Refine your search terms to Chloroplast")
    optional.add_argument("--annotated",
                        action="store_true", 
                        help="Delete unannotated records")
    optional.add_argument("--header",
                        type=str,
                        help="Header fields (Genbank fields)",
                        required=False)
    optional.add_argument("--genbankid", 
                        action="store_true", 
                        help="Include GenBank ID in fasta header")
    optional.add_argument("--prioritize", 
                        action="store_true",
                        help="Prioritize individual with more genes (mitochondrial)")
    optional.add_argument("--add_synonyms", 
                        type=json.loads,
                        default={},
                        help=("Additional gene name synonyms in JSON format\n\n"
                            "Pynnotate already includes an internal dictionary of gene name synonyms to aid extraction. You can provide additional synonyms for genes not automatically recognized. We recommend running the program first to identify any unrecognized gene synonyms. Add any missing synonyms here to improve matching\n\n"
                            "⚠️  ATTENTION: When selecting the genome type and adding synonyms, they will be incorporated into the internal dictionary for that specific genome type. However, if the genome type selected is 'Other', only the synonyms provided by the user will be used"
                        ),
                        required=False)
    optional.add_argument("--minbp",
                        type=str,
                        help="Minimum sequence size (bp)",
                        required=False)
    optional.add_argument("--maxbp",
                        type=str,
                        help="Maximum sequence size (bp)",
                        required=False)
    optional.add_argument("--extraction", 
                        action="store_true",
                        help="Extract all annotated genes separately")
    optional.add_argument("--overlap", 
                        action="store_true",
                        help="Fix overlap between extracted genes")
    optional.add_argument("--logmissing", 
                        action="store_true",
                        help="Generate log of missing genes per sample (useful for mitogenomes)")
    optional.add_argument("-f", "--folder",
                        type=str,
                        help="Folder name to create inside output folder (will be created if it does not exist)",
                        required=False,
                        default=".")

    # check manual help
    if "--help" in sys.argv or "-h" in sys.argv:
        pydoc.pager(parser.format_help())
        sys.exit(0)

    # parse args normalmente
    args = parser.parse_args()
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
        "folder_name": final_config["folder"],
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