# Pynnotate

A Python tool for annotation processing.

## Installation

Clone the repository and install with:

```bash
git clone https://github.com/fernandacaron/pynnotate.git
cd pynnotate
pip install .

# 🧬 Pynnotate

**Pynnotate** is a Python tool designed to streamline the annotation and extraction of specific gene features from GenBank files. It supports user-defined synonym dictionaries for gene names and automates the generation of per-feature FASTA files and summary tables.

![Python](https://img.shields.io/badge/Python-3.8+-blue?logo=python)
![License](https://img.shields.io/github/license/SEU-USUARIO/pynnotate)
![Status](https://img.shields.io/badge/status-in%20development-orange)

---

## 🚀 Installation

You can install the latest version directly from GitHub:

```bash
pip install git+https://github.com/SEU-USUARIO/pynnotate.git
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

```python
from pynnotate import process_genbank

# Input GenBank file
gb_path = "your_file.gb"

# Dictionary of feature name synonyms
synonyms_dict = {
    "COI": ["COI", "COX1", "CO1"],
    "ATP6": ["ATP6"],
    "ATP8": ["ATP8"]
}

# Run annotation and export results
process_genbank(
    gb_path,
    feature_synonyms=synonyms_dict,
    nucleotide_output="nucleotides.fasta",
    aa_output="proteins.fasta",
    summary_table="summary.csv"
)
```

---

## 📘 Documentation

Detailed documentation will be available soon. Stay tuned!

---

## 🤝 Contributing

Contributions are welcome! To report bugs, request features, or submit improvements, please open an issue or a pull request. See the `CONTRIBUTING.md` file for guidelines.

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).

---

## 📣 Citation

If you use **Pynnotate** in your research, please cite it as follows:

```
CARON, F. M. et al. Pynnotate: A feature annotation and extraction tool for GenBank files. GitHub: https://github.com/SEU-USUARIO/pynnotate
```
