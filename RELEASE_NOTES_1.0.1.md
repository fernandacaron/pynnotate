# Pynnotate 1.0.1

Release date: 2026-08-22

Pynnotate 1.0.1 is a maintenance and reliability release focused on fixing
cross-platform GUI issues, improving GenBank downloads and metadata output,
making gene matrices auditable, and clarifying the behavior of Flexible mode
for species-level supermatrices.

The existing Strict and Unconstrained filtering strategies remain available.
Flexible mode received an intentional behavior update, described below.

## Highlights

- Fixed the gene-synonym editor freeze observed in Windows builds.
- Fixed missing-gene logging and Excel workbook generation.
- Added optional NCBI API-key support.
- Added the GenBank `geo_loc_name` qualifier.
- Improved individual and species grouping in `genes_matrix.xlsx`.
- Added explicit reporting of multiple accessions assigned to the same matrix
  taxon and gene.
- Defined Flexible mode as an organism-level supermatrix workflow.
- Added reproducible PyInstaller configurations for Windows, macOS, and Linux.
- Added an offline regression test suite.

## Bug fixes

### Windows synonym editor

- Removed the modal-window behavior that could freeze the application when the
  user opened or saved the additional gene-synonym editor on Windows.
- Synonym editing now uses non-blocking callbacks and returns the saved JSON to
  the main interface safely.

### Missing-gene log and Excel files

- Fixed `name 'Workbook' is not defined` when **Generate log of missing genes
  per sample** was enabled.
- Enabling the missing-gene log no longer prevents `metadata.xlsx` or
  `genes_matrix.xlsx` from being generated.

### Output-directory creation

- Output folders are now created only after the user confirms the download.
- Cancelling a search after the result count is displayed no longer leaves an
  empty results directory behind.

### GenBank download reliability

- Added retry handling for incomplete or interrupted GenBank responses.
- A batch that cannot be retrieved after the configured retries now stops with
  an explicit error instead of reporting a false successful download.
- API-key values are forwarded consistently to both Entrez searches and
  downloads.

### Gene overlap correction

- Fixed animal mitochondrial overlap detection and trimming.
- Corrected overlap events and the number of removed base pairs are now
  recorded in `log.txt`.

### Taxon and synonym handling

- Prevented provisional taxa such as `sp. 1` and `sp. 2` from being collapsed
  into the same standardized organism name.
- Corrected conflicting aliases involving 12S/18S and chloroplast `rpl36`.

## New and improved metadata support

### NCBI API key

- Added an optional NCBI API-key field to the graphical interface.
- Added `--api-key` / `--api_key` support to the command-line interface.
- Added the optional `api_key` setting to the YAML configuration example.
- Searches continue to work without an API key; the field is optional.

### `geo_loc_name`

- Added the GenBank `geo_loc_name` source qualifier to the selectable FASTA
  header fields.
- Added `geo_loc_name` to metadata workbook output.

## Gene matrix improvements

### Individual-aware grouping

- Repeated `specimen_voucher` or `isolate` values are collapsed into the same
  individual row in non-Flexible gene matrices.
- Individual identity follows this priority when required:
  `specimen_voucher` > `isolate` > GenBank accession.
- When the selected output fields do not identify an individual, the matrix
  can include `individual_id` and `individual_id_source` for provenance.

### Multiple-accession conflicts

- Multiple GenBank accessions assigned to the same matrix taxon and gene are
  preserved instead of being silently overwritten.
- Conflicting accessions are stored together in the matrix cell and reported
  under `[MATRIX CONFLICTS]` in `log.txt`.
- This is particularly useful in Unconstrained mode, where redundant GenBank
  and RefSeq records are intentionally retained for later user curation.
- Pynnotate reports these cases but does not automatically decide which record
  should be removed.

### Workbook usability

- Added clearer header formatting, filters, frozen headers, and practical
  column widths to gene-matrix workbooks.
- These presentation changes do not alter the filtering logic of Strict or
  Unconstrained mode.

## Flexible mode: organism-level supermatrices

Flexible mode is now explicitly defined as a species-level workflow. It can
combine complementary genes from different identified individuals belonging to
the same organism.

### Selection strategy

- Records are identified internally using `specimen_voucher`, then `isolate`,
  then GenBank accession.
- Individuals are ranked by the number of recognized selected genes and the
  amount of usable sequence data.
- Additional individuals from the same organism are retained only when they
  contribute genes that are still missing for that organism.
- **Separate gene extraction** is mandatory in Flexible mode because gene
  recognition is required to apply this strategy.

### Stable FASTA headers

- FASTA files inside the `genes` directory always use the standardized
  organism name as the sequence header.
- GenBank accessions are omitted from Flexible-mode per-gene FASTA headers,
  even when **Include GenBank ID** is enabled.
- The header fields selected by the user and the GenBank-ID checkbox continue
  to control the complete `sequences.fasta` file only.
- If an organism name is unavailable, the GenBank accession is used as a
  logged fallback.

This prevents unique accession numbers from turning every gene record into a
different taxon during downstream concatenation.

### Species matrix and provenance

- The first sheet of `genes_matrix.xlsx` contains one row per organism and one
  column per extracted gene, following the compact matrix layout used in
  previous versions.
- A second sheet named `gene_provenance` records the origin of every selected
  species/gene entry using these columns:

  - `organism`
  - `gene`
  - `GenBank_ID`
  - `isolate`
  - `specimen_voucher`

This keeps the supermatrix species-based while allowing every gene to be
traced back to its original accession and individual metadata.

## Logging and auditing

`log.txt` now provides clearer reporting of:

- sequences found before filtering;
- sequences and species retained;
- individuals represented;
- selected FASTA header fields;
- whether GenBank IDs were included in the complete FASTA;
- Flexible-mode organism-only header behavior;
- missing organism-name fallbacks;
- missing genes per sample;
- unknown and duplicated gene aliases;
- corrected mitochondrial overlaps;
- records removed by filtering;
- multiple-accession matrix conflicts.

## Graphical interface

- Added the NCBI API-key input to Step 4.
- Updated warnings and option states for gene extraction and overlap
  correction.
- Flexible mode automatically keeps separate gene extraction enabled.
- Clarified that the selected header fields apply to the complete
  `sequences.fasta`, while Flexible per-gene FASTAs use organism-only labels.
- Improved packaged icon and image resource discovery.

## Command-line and configuration updates

- Added optional API-key support.
- Kept GUI, CLI, and YAML configuration keys synchronized.
- Added validation for incompatible option combinations, including Flexible
  mode without separate gene extraction.
- Preserved compatibility with existing 1.0.0 configurations; the new API-key
  field is optional.

## Packaging and build system

- Added dedicated PyInstaller specifications for:

  - Windows folder build: `Pynnotate_win.spec`
  - Windows single executable: `Pynnotate_win_onefile.spec`
  - macOS application bundle: `Pynnotate_mac.spec`
  - Linux folder build: `Pynnotate_linux.spec`

- Included graphical resources in both Python packages and PyInstaller builds.
- Added reproducible build dependencies in `requirements-build.txt`.
- Added English and Portuguese build instructions in `BUILD.md` and
  `BUILD.pt-BR.md`.
- Reconciled the legacy compilation notes with the current package structure
  and resource paths.

PyInstaller applications must still be compiled on their target operating
system: Windows builds on Windows, macOS builds on macOS, and Linux builds on
Linux.

## Code quality and validation

- Removed unnecessary runtime dependencies from the installation list.
- Added offline regression tests covering the principal fixed behaviors.
- Added compilation and indentation checks for the Python modules and
  PyInstaller specifications.
- Validated gene-matrix grouping against the supplied 905-row example matrix.

The current regression suite contains 12 tests covering cancellation behavior,
API-key forwarding, Flexible-mode extraction and headers, matrix grouping and
provenance, missing-gene logging, `geo_loc_name`, overlap correction, and
provisional taxon names.

## Upgrade notes

- Users of Strict and Unconstrained mode can continue using their existing
  filtering workflows.
- Unconstrained mode intentionally retains all compatible records; consult
  `[MATRIX CONFLICTS]` before downstream concatenation or manual curation.
- Flexible-mode users should expect per-gene FASTA headers and the first gene
  matrix sheet to be fixed at the organism level.
- Existing YAML files remain valid. Add `api_key` only if desired.
- Rebuild the executable for each target operating system after updating the
  source.

## Recommended post-upgrade checks

1. Open and save a small custom synonym dictionary.
2. Cancel a confirmed search and verify that no empty output folder is created.
3. Run a small search with the missing-gene log enabled.
4. Confirm that `sequences.fasta`, `metadata.xlsx`, `genes_matrix.xlsx`,
   `log.txt`, and the `genes` directory are generated.
5. In Flexible mode, verify organism-only headers in the per-gene FASTAs and
   inspect the `gene_provenance` sheet.
6. In Unconstrained mode, inspect any `[MATRIX CONFLICTS]` entries before
   downstream analysis.
7. Repeat the synonym-editor test on Windows before publishing the Windows
   executable.
