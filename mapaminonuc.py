from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Data import CodonTable
import os

## dictionary for the codon-tables
codon_tables = {
    1: "Standard (Universal)",
    2: "Vertebrate Mitochondrial",
    3: "Yeast Mitochondrial",
    4: "Mold, Protozoan, and Coelenterate Mitochondrial",
    5: "Invertebrate Mitochondrial",
    6: "Ciliate Nuclear",
    9: "Echinoderm Mitochondrial",
    10: "Euplotid Nuclear",
    11: "Bacterial",
    12: "Alternative Yeast Nuclear",
    13: "Ascidian Mitochondrial",
    14: "Flatworm Mitochondrial",
    15: "Blepharisma Nuclear",
    16: "Chlorophyta",
    21: "Fungal Nuclear"
}

## printing the codon-table options
print("Choose a codon-table:")
for id, description in codon_tables.items():
    print(f"{id}: {description}")

## ask for the user input
try:
    user_choice = int(input("Type the codon-table ID that you wish to use: "))
    
    ## checking if the choice is valid
    if user_choice in codon_tables:
        chosen_table = CodonTable.unambiguous_dna_by_id[user_choice]
        print(f"Codon-table chosen: {codon_tables[user_choice]}")
    else:
        print("Invalid ID! Please, choose another listed ID.")
except ValueError:
    print("Please, type a valid number.")

## getting the files
nucleotide_file_path = input("Insert the nucleotide file path: ")
nucleotide_format = os.path.splitext(nucleotide_file_path)[1][1:]

aa_file_path = input("Insert the amino acid file path: ")
aa_format = os.path.splitext(aa_file_path)[1][1:]

nt_seq_records = list(SeqIO.parse(nucleotide_file_path, nucleotide_format))
aa_seq_records = list(SeqIO.parse(aa_file_path, aa_format))

def map_amino_to_nuc_no_alignment(nt_seq, aa_seq):
    codon_map = []
    nt_pos = 0
    
    for aa in aa_seq:
        ## ignoring gaps
        if aa != "-":
            codon = str(nt_seq[nt_pos:nt_pos + 3])
            codon_map.append(codon)            
            nt_pos += 3
    return "".join(codon_map)

def map_amino_to_nuc(nt_seq, aa_seq):
    codon_map = []
    nt_pos = 0
    
    for aa in aa_seq:
        if aa == "-":
            codon_map.append("---")
        else:
            codon = str(nt_seq[nt_pos:nt_pos + 3])
            codon_map.append(codon)
            nt_pos += 3
    return "".join(codon_map) 

## function to reorder sequences
def reorder_nt_sequences(nt_seq_records, aa_seq_records):
    ordered_nt_sequences = []

    # Mapeia os IDs das sequências de aminoácidos para sequências de nucleotídeos
    aa_to_nt_map = {record.id: record.seq for record in nt_seq_records}

    for aa_record in aa_seq_records:
        nt_seq = aa_to_nt_map.get(aa_record.id)

        if nt_seq is None:
            raise ValueError(f"Nucleotide sequence not found for ID: {aa_record.id}")

        ordered_nt_sequences.append((aa_record.id, nt_seq))

    return ordered_nt_sequences

# Reordenando as sequências
ordered_nt_sequences = reorder_nt_sequences(nt_seq_records, aa_seq_records)

# Verificando se os arquivos têm o mesmo número de sequências
if len(ordered_nt_sequences) != len(aa_seq_records):
    raise ValueError("The number of sequences does not match!")

output_filename = "aligned_nucleotide_sequences.fasta"

try:
    user_choice = input("Is the amino acid sequence aligned? (Yes/No) ").strip().lower()
    if user_choice == "yes":
        aligned = True
    elif user_choice == "no":
        aligned = False
    else:
        raise ValueError("Please, type 'Yes' or 'No'.")
except ValueError as e:
    print(e)
    exit()  # Saia do programa se a entrada for inválida

with open(output_filename, "w") as output_handle:
    for (nt_id, ordered_nt_sequence), aa_seq_record in zip(ordered_nt_sequences, aa_seq_records):
        if aligned:
            aligned_nucleotide_seq = map_amino_to_nuc(ordered_nt_sequence, aa_seq_record.seq)
        else:
            aligned_nucleotide_seq = map_amino_to_nuc_no_alignment(ordered_nt_sequence, aa_seq_record.seq)

        output_handle.write(f">{aa_seq_record.id}\n")
        output_handle.write(f"{aligned_nucleotide_seq}\n")

