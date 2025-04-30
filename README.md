# 🧬 Pynnotate

**Pynnotate** is a Python tool designed to streamline the annotation and extraction of specific gene features from GenBank files. It supports user-defined synonym dictionaries for gene names and automates the generation of per-feature FASTA files and summary tables.

![Python](https://img.shields.io/badge/Python-3.8+-blue?logo=python)
![License](https://img.shields.io/github/license/fernandacaron/pynnotate)
![Status](https://img.shields.io/badge/status-in%20development-orange)
![Version](https://img.shields.io/github/v/release/fernandacaron/pynnotate?logo=github)

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

#### Explanation of Arguments:

-e or --email: Your email address (required to access GenBank).
-a or --accession: Path to the file containing GenBank accession numbers. In this example, we're using examples/accession.csv.
-t or --type: The type of DNA being processed (e.g., mtDNA_animals, mtDNA_plants, rDNA, cpDNA, or other). Here, we're using mtDNA_animals.
-f or --features: (Optional) Path to the file listing specific features to extract. Here, we're using examples/features.csv. If not provided, all features are extracted.
-r or --repeated: (Optional) Path to the file listing features that appear multiple times. Here, we're using examples/repeated.csv.
-s or --add_synonyms: (Optional) Path to the file containing additional synonyms. Here, we're using examples/synonyms.csv.
-o or --output: (Optional) Directory where output files will be saved. If not provided, files are saved in the current directory. In this case, a folder named output_folder will be created (if it doesn't exist).

#### Important Notes:

Make sure the paths to your input files are correct relative to where you're running the python main.py command.
The -f, -r, and -s arguments are optional; include them only if you have files for specific features, repeated features and additional synonyms, respectively.
Pynnotate will create the output_folder if it doesn't already exist.

For more details, run:

```bash
python main.py -h
```

This will display the help message with a full list of available options.

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
