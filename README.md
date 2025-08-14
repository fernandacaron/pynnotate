# 🧬 Pynnotate

**Pynnotate** is a Python graphical tool (GUI) for searching, downloading, and automatically annotating genetic sequences from GenBank.  
Developed for both advanced researchers and teachers or students new to bioinformatics, phylogeny, and molecular genetics, Pynnotate offers a user-friendly interface that requires no prior programming knowledge.

![Python](https://img.shields.io/badge/Python-3.8+-blue?logo=python)  
[![Licence: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)  

1.  [✨ Features](#-features)  
2.  [🛠️ Requirements](#-requirements)  
3.  [💾 Installation](#-installation)  
4.  [🧪 Usage example](#-usage-example---terminal-version)  
5.  [⚙️ Argument details](#%EF%B8%8F-argument-details)  
6.  [🧾 Generated files](#-generated-files)  
7.  [🤝 Contributing](#-contributing)  
8.  [📣 Citation](#-citation)  

---

## 👥 Authors

**Fernanda de Souza Caron**  
PhD researcher, PPG Ecology and Conservation (UFPR)

**Felipe de Medeiros Magalhães**  
Postdoctoral researcher, Federal University of Paraíba (UFPB)

**Matheus Salles**  
PhD researcher, PPG Zoology (UFPR)

**Fabricius M. C. B. Domingos**  
Researcher and lecturer, PPG Zoology (UFPR)

---

## ✨ Features

- 🔍 Simple search by free terms or specific IDs in GenBank
- 🧠 Automatic gene extraction with synonym grouping
- ✂️ Sequence length filters and options to prioritise samples, ideal for different analysis levels
- Filtering modes:  
  🌐 Unconstrained mode: includes all sequences found  
  🌱 Flexible mode (`unique_species = True`): allows multiple sequences per species if genes differ  
  🔒 Strict mode (`prioritize_more_genes = True`): includes only the best sequence per species, ideal for simple analyses  
- 🧬 Supports mitogenomes, chloroplasts, and nuclear genomes
- 👓 Automatic identification of multiple copies of tRNA-Leu and tRNA-Ser, with grouping by genomic position
- 🖼️ Intuitive graphical interface for configuration, execution, and monitoring without command-line use
- 📂 Complete generation of FASTA files, Excel spreadsheets, and detailed logs, ready for teaching or research

---

## 💾 Installation

### Terminal version

The terminal version of Pynnotate is recommended for users who prefer to run the tool via command line or integrate it into automated pipelines.

1. Clone the GitHub repository:

```bash
git clone https://github.com/fernandacaron/pynnotate.git
cd pynnotate
pip install .
```

> Requirements: Python 3.8+

2. Run the programme in the terminal with:

```bash
python pynnotate.py --help
```

### Graphical version (GUI)

For ease of use, we provide a ready-to-use graphical version packaged as a `.app` file for major operating systems.

1. Go to the Releases page on GitHub
2. Download the installer for your system
3. Install/unpack and run the programme by clicking its icon
4. The graphical interface will open, allowing you to configure and run all functions without the terminal

---

## 🧪 Usage example

### Graphical version

1. Set a gene (e.g., COI) and an organism (e.g., Anura)  
2. Click “💾 Search and download sequences”  
3. The programme will search, download, and extract the data automatically  
4. View the generated files in the chosen location  

### Terminal version

Pynnotate uses a YAML configuration file to simplify option settings. An example file is available in the `examples/` folder of the repository, named `config.yaml`.

Run with the YAML file:

```bash
python pynnotate.py -c examples/config.yaml
```

#### Important notes:

- The YAML file groups all configurations, avoiding the need for multiple command-line arguments.  
- Ensure that all file paths in the YAML are correct.  
- To see all available options and their descriptions, run:

```bash
python pynnotate.py -h
```

---

## ⚙️ Argument Details

Pynnotate is a command-line tool that accepts various arguments to customise the search, download, and extraction of sequences from GenBank. Below is a detailed description of each argument available in the current code.

#### **Mandatory Arguments**

##### `-c` or `--config`

Description: Path to the YAML configuration file containing all options to run Pynnotate.

> Note: The YAML file groups all settings, making it easier to use without multiple command-line arguments. An example is available in the `examples/` folder.

#### **Mandatory Arguments in the YAML File**

To run Pynnotate correctly via terminal, you must provide a YAML configuration file with at least the following required fields:

##### `-e` or `--email`

Description: Your valid email, required by NCBI Entrez for identification and access to GenBank.

##### `-o` or `--output`

Description: Directory where output files will be saved (the folder name can also be provided via the `--folder` argument, but it is not mandatory).

##### `-t` or `--type`

Description: Type of genome/organism to determine the synonym dictionary. Accepted values: *animal\_mito, plant\_mito, plant\_chloro, other*.

##### `--filter-mode`

Description: Defines how sequences will be filtered by species. This parameter is essential to control redundancy and the structure of your dataset.

**Accepted values:**

🌐 Unconstrained: Includes all available sequences regardless of redundancy. Useful when you want to manually explore or curate all records.

🌱 Flexible: Allows multiple sequences per species only if each new sequence adds different genes (e.g., in supermatrix analyses).

🔒 Strict: Includes only one sequence per species, prioritising the one with the highest number of genes present in the main dictionary or in the dictionary provided by the user.

**⚠️ NOTE**: In strict mode, the filter considers the genes listed in the default synonym dictionary and/or the dictionary provided by the user.

**⚠️ NOTE**: When the unconstrained mode is used in combination with separate gene extraction (`--extraction`), all sequences corresponding to the selected genes will be downloaded, even if there are multiple records per species.

**🚨 In addition, you must include either `--accession` or a search term in the query (`--genes`, `--organism`, `--publication` or `--additional`) to indicate the data search:**

##### `-a` or `--accession`

Description: List of GenBank IDs (accessions) to download. Can be null if using a query argument.

> Note: Use only if you want to search for specific IDs instead of using a query.

#### **Optional Arguments (via YAML or Command Line)**

##### `-g` or `--genes`

Description: Comma-separated list of genes to search and download (e.g., COI, CYTB, ATP6).

> Note: Extracts only the listed genes; otherwise, all known genes are extracted.

##### `-organism`

Description: Organisms to search and download (e.g., species, family).

##### `-p` or `--publication`

Description: Publication term (e.g., title, authors, year).

##### `--additional`

Description: Any additional search term (e.g., NOT sp).

##### `--mitochondrialgene`

Description: Refine search terms to "mitochondrial genes".

##### `--mitogenome`

Description: Refine search terms to "mitogenomes".

##### `--chloroplast`

Description: Refine search terms to "chloroplast".

##### `--annotated`

Description: Exclude unannotated records.

##### `--header`

Description: Fields for sequence headers (GenBank fields).

##### `--genbankid`

Description: Include GenBank ID in the fasta headers.

##### `--prioritize`

Description: Prioritise individuals with more genes (valid for mitochondrial genes).

##### `--add_synonyms`

Description: Additional synonyms for gene names in JSON format. Pynnotate already includes an internal dictionary of gene name synonyms to assist in extraction. You can provide additional synonyms for genes not automatically recognised. We recommend running the program first to identify any unrecognised gene synonyms. Add any missing synonyms here to improve matching.

**⚠️ NOTE**: When selecting the genome type and adding synonyms, they will be incorporated into the internal dictionary for that specific genome type. However, if the selected genome type is 'other', only user-provided synonyms will be used.

##### `--min_bp`

Description: Defines the minimum allowed length for a sequence to be retained.

##### `--max_bp`

Description: Defines the maximum allowed length for a sequence to be retained.

##### `--extraction`

Description: Boolean. If True, extracts all genes separately, grouping different individuals/species into the respective files for each gene.

##### `--overlap`

Description: Adjust overlap between extracted genes.

##### `--logmissing`

Description: Generate a log of missing species per sample (useful for mitogenomes).

##### `--folder`

Description: Name of the folder to create inside the output directory (will be automatically created with a predefined name if the argument is not provided).

#### **Other Options**

##### `-h` or `--help`

Description: Shows help with the complete list of arguments and their descriptions.

---

## 🧾 Generated files

After running, Pynnotate creates the following in the specified output directory:

1. `sequences.fasta`: Extracted sequences without gene separation.  
2. `log.txt`: Execution log for debugging and traceability.  
3. `metadata.xlsx`: Metadata from GenBank for each sequence.  
4. `genes_matrix.xlsx`: Presence/absence matrix of genes in downloaded records, with accession numbers.  
5. `genes/`: Folder containing sequences separated by gene.

---

## 🤝 Contributing

Contributions are welcome! This is an open-source project, free for academic purposes.  

To report bugs, request features, or submit improvements, open an issue or pull request.

---

## 📣 Citation

If you use **Pynnotate** in your research, please cite it as:

```
Caron, F. S.*, Magalhães, F. M.*, Salles, M., & Domingos, F. M. B. C. (2025). pynnotate: a flexible tool for retrieving and processing GenBank data in molecular evolution research and education. GitHub: https://github.com/fernandacaron/pynnotate
```
