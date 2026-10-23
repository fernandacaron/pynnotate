# Building Pynnotate 1.0.2rc1

Portuguese instructions: [BUILD.pt-BR.md](BUILD.pt-BR.md).

PyInstaller applications must be built on the operating system they target.
Use a clean virtual environment for each build.

## Prepare the environment

```bash
python -m venv .venv
```

Activate it on macOS/Linux:

```bash
source .venv/bin/activate
```

Activate it on Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install the runtime and build dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt -r requirements-build.txt
```

On Linux, Tk may also need to be installed by the operating system (for
example, `python3-tk` on Debian/Ubuntu).

## Build

Windows, folder build (recommended for testing):

```powershell
python -m PyInstaller --clean --noconfirm Pynnotate_win.spec
```

This creates `dist\Pynnotate\Pynnotate.exe`. Keep the complete
`dist\Pynnotate` folder together.

Windows, single executable (release option):

```powershell
python -m PyInstaller --clean --noconfirm Pynnotate_win_onefile.spec
```

This creates `dist\Pynnotate.exe`. One-file mode starts more slowly because it
extracts the application to a temporary directory every time it is opened.

macOS:

```bash
python -m PyInstaller --clean --noconfirm Pynnotate_mac.spec
```

The macOS build uses the architecture of the active Python installation. Build
with an ARM Python for Apple Silicon and an x86_64 Python for Intel, or use a
universal2 Python and set `target_arch="universal2"` in the macOS spec.

Linux:

```bash
python -m PyInstaller --clean --noconfirm Pynnotate_linux.spec
```

The generated application is placed in `dist/`. For folder builds, test and
distribute the complete `Pynnotate` folder (or `Pynnotate.app` on macOS).

## macOS local testing and DMG

Open the application:

```bash
open dist/Pynnotate.app
```

If macOS quarantines an application that you built locally, clear the
quarantine attribute only from that local test build:

```bash
xattr -cr dist/Pynnotate.app
```

After validating the application, an unsigned test DMG can be created with:

```bash
hdiutil create -volname "Pynnotate" \
  -srcfolder dist/Pynnotate.app \
  -ov -format UDZO Pynnotate-mac.dmg
```

For public distribution, sign and notarize the `.app` instead of asking users
to bypass Gatekeeper.

## Before releasing

1. Open the gene-synonym editor and save a small JSON dictionary.
2. Search by two manual GenBank accessions, first cancelling and then
   confirming the download.
3. Repeat with the missing-gene log enabled.
4. Confirm that `sequences.fasta`, `metadata.xlsx`, `log.txt`, and the `genes`
   directory are generated as expected.
5. Test once without an API key and once with an NCBI API key.
6. In Step 1, click **Check for updates** and confirm that a network failure preserves the current database.
7. Confirm that `synonyms_used.json` is generated with the other result files.
