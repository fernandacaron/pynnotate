# 🧬 Pynnotate

**Pynnotate** is a Python tool designed to streamline the annotation and extraction of specific gene features from GenBank files. It supports user-defined synonym dictionaries for gene names and automates the generation of per-feature FASTA files and summary tables.

![Python](https://img.shields.io/badge/Python-3.8+-blue?logo=python)
![License](https://img.shields.io/github/license/fernandacaron/pynnotate)
![Status](https://img.shields.io/badge/status-in%20development-orange)
![Version](https://img.shields.io/github/v/release/fernandacaron/pynnotate?logo=github)

## 📋 Table of Contents

1.  [🚀 Installation](#installation)
2.  [💡 Features](#features)
3.  [🧪 Example usage](#example-usage)
4.  [⚙️ Argument Details](#argument-details)
5.  [📘 Documentation](#documentation)
6.  [🤝 Contributing](#contributing)
7.  [📣 Citation](#citation)

---

## 🚀 Installation

You can install the latest version directly from GitHub:

```bash
git clone https://github.com/fernandacaron/pynnotate.git
cd pynnotate
pip install .
```

> Requirements: Python 3.8+, Biopython

---

## 💡 Features

- 🔍 Automated detection and extraction of specific features (e.g., COI, ATP6, 16S)
- 🧠 Flexible support for feature name synonyms
- 📂 Generation of per-feature `.fasta` files for nucleotide and amino acid sequences
- 📊 Summary table reporting the presence/absence of features across species
- 🔄 Compatible with multi-record GenBank files

---

## 🧪 Example usage

To run Pynnotate, you need to provide several arguments via the command line. Here's a comprehensive example using the provided example files:

```bash
python main.py \
    -e seu_email@exemplo.com \
    -a examples/accession.csv \
    -t mtDNA_animals \
    -f examples/features.csv \
    -r examples/repeated.csv \
    -s examples/synonyms.csv \
    -o output_folder
```

### Important Notes:

Make sure the paths to your input files are correct relative to where you're running the python main.py command.
The -f, -r, and -s arguments are optional; include them only if you have files for specific features, repeated features and additional synonyms, respectively.
Pynnotate will create the output_folder if it doesn't already exist.

For more details, run:

```bash
python main.py -h
```

This will display the help message with a full list of available options.

---

## ⚙️ Argument Details

Pynnotate is a command-line tool that accepts several arguments to customize the GenBank data retrieval and feature extraction process. Below is a detailed explanation of each argument:

### **Required Arguments**

These arguments must be provided for Pynnotate to run successfully.

#### `-e` or `--email`

* **Description:** Your valid email address. This is required by the NCBI Entrez API to identify users and is essential for accessing the GenBank database.
* **Usage:**
    ```bash
    -e your_email@example.com
    ```
* **Example:**
    ```bash
    -e john.doe@university.edu
    ```
* **Importance:** Providing a valid email is crucial; otherwise, you may encounter errors or be blocked from accessing the NCBI Entrez API.

#### `-a` or `--accession`

* **Description:** The path to a file containing a list of GenBank accession numbers. Each accession number should be on a new line in the file. Pynnotate uses these accession numbers to retrieve the corresponding GenBank records. The file can be either a plain text file (`.txt`) or a CSV file (`.csv`) with one accession number per row.
* **Usage:**
    ```bash
    -a path/to/accession_list.txt
    ```
    or
    ```bash
    -a path/to/accession_list.csv
    ```
* **Examples:**
    * `accession_list.txt`:
        ```text
        AC12345.1
        AC12346.2
        AC12347.3
        ```
    * `accession_list.csv`:
        ```csv
        AC12345.1
        AC12346.2
        AC12347.3
        ```
* **Note:** Ensure that the file path is correct and the file exists.

#### `-t` or `--type`

* **Description:** Specifies the type of DNA being processed. This argument is crucial because it determines which predefined synonym dictionary Pynnotate will use to identify gene features.
* **Supported Types:**
    * `mtDNA_animals`: For animal mitochondrial DNA.
    * `mtDNA_plants`: For plant mitochondrial DNA.
    * `rDNA`: For ribosomal DNA.
    * `cpDNA`: For chloroplast DNA.
    * `other`:  If none of the above types apply, you can use 'other' and provide your own synonym dictionary.
* **Usage:**
    ```bash
    -t mtDNA_animals
    ```
* **Caution:** Using an incorrect type may result in features not being identified correctly.

### **Optional Arguments**

These arguments are not required but can be used to customize the script's behavior.

#### `-f` or `--features`

* **Description:** The path to a file listing the specific features you want to extract from the GenBank records. If this argument is not provided, Pynnotate will extract all known features (using the specified synonym dictionary). The file should contain one feature name per line. It can be a `.txt` or `.csv` file.
* **Usage:**
    ```bash
    -f path/to/features.txt
    ```
    or
    ```bash
    -f path/to/features.csv
    ```
* **Examples:**
    * `features.txt`:
        ```text
        COI
        ATP6
        16S
        ```
    * `features.csv`:
        ```csv
        COI
        ATP6
        16S
        ```
* **Flexibility:** This argument allows you to extract only the features relevant to your analysis, reducing processing time and output file sizes.

#### `-r` or `--repeated`

* **Description:** The path to a file listing features that may appear multiple times within a genome.  This is important for handling cases where you need to extract each instance of a feature (e.g., multiple copies of a tRNA). The file should contain one feature name per line (.txt or .csv).
* **Usage:**
    ```bash
    -r path/to/repeated.txt
    ```
    or
    ```bash
    -r path/to/repeated.csv
    ```
* **Example:**
    * `repeated.txt`:
        ```text
        tRNA-Leu
        ```
    * `repeated.csv`:
        ```csv
        tRNA-Leu
        ```
* **Advanced Usage:** Pynnotate will create separate output files for each instance of these repeated features, naming them in a way that indicates their order or context.

#### `-o` or `--output`

* **Description:** The path to the directory where the output files (FASTA files and summary CSV) will be saved. If this argument is not provided, the output files will be saved in the current working directory. If the specified directory does not exist, Pynnotate will create it.
* **Usage:**
    ```bash
    -o path/to/output_directory
    ```
* **Example:**
    ```bash
    -o results/
    ```
* **Organization:** Using this argument helps keep your output organized, especially when processing multiple datasets.

#### `-s` or `--add_synonyms`

* **Description:** The path to a file containing additional synonyms for gene features. This allows you to extend or override the default synonym dictionaries used by Pynnotate. The file should be a CSV with the official gene name in the first column and its synonyms in subsequent columns.
* **Usage:**
    ```bash
    -s path/to/synonyms.csv
    ```
* **Example:**
    * `synonyms.csv`:
        ```csv
        COI,COX1,CO1,cytochrome oxidase I
        ND4,NADH4,NADH dehydrogenase subunit 4
        ```
* **Customization:** This argument is useful when working with non-standard gene names or when the default dictionaries do not cover all the synonyms you need.

By understanding these arguments, you can effectively use Pynnotate to extract and analyze genetic features from GenBank data.

---

## 📘 Documentation

Detailed documentation will be available soon. Stay tuned!

---

## 🤝 Contributing

Contributions are welcome! To report bugs, request features, or submit improvements, please open an issue or a pull request. 

---

## 📣 Citation

If you use **Pynnotate** in your research, please cite it as follows:

```
Caron, F. S., Salles, M., & Domingos, F. M. B. C. (2025). pynnotate: a flexible tool for retrieving and processing genetic feature annotations from GenBank. GitHub: https://github.com/fernandacaron/pynnotate
```
