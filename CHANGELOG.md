# Changelog

## 1.0.2rc1

- Added a pre-download query review window showing the complete generated Entrez query, with scrolling and a one-click clipboard copy option for independent verification in GenBank.
- Added an elapsed-time counter to the graphical download window and recorded the final duration in both the completion message and `log.txt`.
- Strict mode now ranks individuals by recognized target-gene count, then cumulative extracted target-gene length, and finally a stable GenBank-accession tie-breaker.
- Replaced the ambiguous `Unwanted gene according to the parameters specified` log entry with distinct reasons for missing annotations, unrecognized aliases, non-target genes, and extraction failures.
- Renamed the ambiguous `DUPLICATED ALIASES` log section to `REPEATED GENE ANNOTATIONS` and clarified that repeated features are not yet classified as redundant annotations or distinct biological copies.
- Added a manual **Check for updates** control below the Step 1 synonym editor.
- One small versioned download updates the Animal mitochondrial, Plant mitochondrial, and Chloroplast dictionaries together.
- Kept the bundled database as an offline fallback and added **Restore included**.
- Added schema, revision, HTTPS, size, and SHA-256 validation plus atomic installation.
- Added **Suggest a synonym**, which opens a pre-filled GitHub issue.
- Added synonym database version and hashes to `log.txt`.
- Added `synonyms_used.json` to every result set for reproducibility.
- Added offline tests for full updates, checksum rejection, and restoration.

## 1.0.1

Detailed release notes: [RELEASE_NOTES_1.0.1.md](RELEASE_NOTES_1.0.1.md)

- Fixed the Windows gene-synonym editor modal freeze.
- Fixed missing-gene logging so metadata workbooks are still generated.
- Delayed output-directory creation until the download is confirmed.
- Added optional NCBI API-key support to the GUI, CLI, and YAML configuration.
- Added the GenBank `geo_loc_name` source qualifier to FASTA-header selection and metadata workbooks.
- Fixed `Workbook is not defined` when generating the multi-gene sample matrix.
- Collapsed repeated `specimen_voucher` or `isolate` values into one individual row in the gene matrix, regardless of filter mode.
- Preserved and logged multiple accessions found for the same individual and gene instead of silently overwriting them.
- Decoupled Flexible-mode individual identification from FASTA header choices.
- Flexible mode now ranks individuals by canonical recognized-gene coverage before filling missing genes for each species.
- Flexible-mode per-gene FASTA files now always use organism-only species labels and omit GenBank IDs automatically, while the complete FASTA still respects the selected fields and GenBank-ID checkbox.
- The Flexible-mode gene matrix now always contains one row per organism, matching the species labels used in per-gene FASTA files.
- Added a `gene_provenance` worksheet with the GenBank accession, isolate, and specimen voucher used for every species/gene entry.
- Added a logged GenBank-ID fallback for the rare case in which an organism name is unavailable.
- Added explicit individual identity provenance to non-Flexible gene matrices when only non-individual header fields are selected.
- Made separate gene extraction mandatory when Flexible mode is selected.
- Fixed animal mitochondrial gene-overlap detection and logging.
- Prevented provisional taxa such as `sp. 1` and `sp. 2` from being collapsed.
- Fixed conflicting gene aliases for 12S/18S and chloroplast `rpl36`.
- Added retries that stop incomplete downloads instead of reporting false success.
- Included graphical resources in Python packages and PyInstaller bundles.
- Added separate Windows, macOS, and Linux build specifications.
- Added both Windows folder and single-executable build specifications.
- Reconciled the legacy compile notes with the current resource paths and outputs.
- Added offline regression tests and reproducible build instructions.
