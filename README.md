# 🧬 Pynnotate

[🇧🇷 Versão em Português](README.pt-BR.md)

**Pynnotate** is a Python graphical tool (GUI) for searching, downloading, and automatically annotating genetic sequences from GenBank. Developed for both advanced researchers and teachers or students new to bioinformatics, phylogeny, and molecular genetics, Pynnotate offers a user-friendly interface that requires no prior programming knowledge.

![Python](https://img.shields.io/badge/Python-3.8+-blue?logo=python)  [![Licence: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)  

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

The terminal version of Pynnotate is recommended for users who prefer to run the tool via command line or integrate it into automated pipelines. This method works on **Windows**, **macOS**, and **Linux**.

> **Requirements:**  
> - Python **3.8 or newer** must be installed.  
> - To check if Python is installed, run:  
>   ```bash
>   python --version
>   ```
>   or  
>   ```bash
>   python3 --version
>   ```
> - If Python is not installed:  
>   - **Windows**: [Download from python.org](https://www.python.org/downloads/windows/) and check “Add Python to PATH” during installation.  
>   - **macOS**: Install via [python.org](https://www.python.org/downloads/macos/) or use Homebrew:  
>     ```bash
>     brew install python
>     ```
>   - **Linux**: Use your package manager, e.g.:  
>     ```bash
>     sudo apt install python3 python3-pip
>     ```

1. Clone the GitHub repository:  

   If you have issues with SSH authentication, use the HTTPS version below (recommended for most users):

   **HTTPS (recommended)**:
   ```bash
   git clone https://github.com/fernandacaron/pynnotate.git
   ```

   **SSH (for advanced users with SSH keys configured)**:
   ```bash
   git clone git@github.com:fernandacaron/pynnotate.git
   ```

2. Navigate to the project folder:

   ```bash
   cd pynnotate
   ```

3. Install Pynnotate:

   ```bash
   pip install -e .
   ```

4. Test if the programme runs:

   ```bash
   pynnotate --help
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
3. Review the complete generated Entrez query and optionally copy it for verification in GenBank
4. Confirm the download; the programme will download and extract the data automatically
5. View the generated files in the chosen location

### Terminal version

Pynnotate uses a YAML configuration file to simplify option settings. An example file is available in the `examples/` folder of the repository, named `config.yaml`.

Run with the YAML file:

```bash
pynnotate --config pynnotate/examples/config.yaml
```

#### Important notes:

- The YAML file groups all configurations, avoiding the need for multiple command-line arguments.  
- Ensure that all file paths in the YAML are correct.  
- To see all available options and their descriptions, run:

```bash
pynnotate --help
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

##### `--api-key` (optional)

Description: NCBI API key recommended for larger searches. The key is used only for Entrez requests and is not written to Pynnotate output logs.

##### `-o` or `--output`

Description: Directory where output files will be saved (the folder name can also be provided via the `--folder` argument, but it is not mandatory).

##### `-t` or `--type`

Description: Type of genome/organism to determine the synonym dictionary. Accepted values: *animal\_mito, plant\_mito, plant\_chloro, other*.

**⚠️ NOTE**: Genome type selection affects gene extraction and filtering. When extraction is disabled, all sequences matching your search will be downloaded regardless of genome type.

##### `-f`, `--filter-mode` or `--filter_mode`

Description: Defines how sequences will be filtered by species. This parameter is essential to control redundancy and the structure of your dataset.

**Accepted values:**

🌐 Unconstrained: Includes all available sequences regardless of redundancy. Useful when you want to manually explore or curate all records.

🌱 Flexible: Groups accessions into individuals internally (`specimen_voucher`, then `isolate`, then GenBank ID), ranks individuals by recognized gene coverage, and adds further individuals only when they contribute genes still missing for the species. Separate gene extraction is required. Because this mode builds species-level supermatrices, every per-gene FASTA header contains only the standardized organism name and `genes_matrix.xlsx` always contains one row per organism. The fields selected for the complete FASTA and the GenBank-ID checkbox affect only `sequences.fasta`. A second `gene_provenance` worksheet records the GenBank accession, isolate, and specimen voucher used for each species/gene entry. If an organism name is unavailable, the accession is used as a logged fallback.

🔒 Strict: Includes only one sequence per species, prioritising the one with the highest number of genes present in the main dictionary or in the dictionary provided by the user.

**⚠️ NOTE**: In strict mode, the filter considers the genes listed in the default synonym dictionary and/or the dictionary provided by the user.

**⚠️ NOTE**: When the unconstrained mode is used in combination with separate gene extraction (`--extraction`), all sequences corresponding to the selected genes will be downloaded, even if there are multiple records per species.

**🚨 In addition, you must include either `--accession` or a search term in the query (`--genes`, `--organism`, `--publication` or `--additional`) to indicate the data search:**

##### `--accession`

Description: List of GenBank IDs (accessions) to download. Can be null if using a query argument.

> Note: Use only if you want to search for specific IDs instead of using a query.

#### **Optional Arguments (via YAML or Command Line)**

##### `--genes`

Description: Comma-separated list of genes to search and download (e.g., COI, CYTB, ATP6).

> Note: Extracts only the listed genes; otherwise, all known genes are extracted.

##### `-organism`

Description: Organisms to search and download (e.g., species, family).

##### `--publication`

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

**⚠️ NOTE**: Gene extraction will be limited to the selected synonym dictionary. For example, when selecting 'plant_chloro', only chloroplast genes will be extracted.

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
6. `synonyms_used.json`: Exact effective synonym dictionary and database version used in the run.

## 🔄 Updating the official synonym database

In Step 1, click **Check for updates** only when you want to contact the official
Pynnotate repository. A single small update covers the Animal mitochondrial,
Plant mitochondrial, and Chloroplast dictionaries. Pynnotate never checks on
startup, continues to work offline, validates the downloaded file before using
it, and can return to the database included with the app through **Restore
included**. User-added synonyms remain separate and are applied after the
official dictionary.

Maintainer and local-test instructions are in
[`SYNONYM_DATABASE.md`](SYNONYM_DATABASE.md).

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
