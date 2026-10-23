# Synonym database maintenance and testing

Beginner instructions in Portuguese: [TUTORIAL_ATUALIZAR_SINONIMOS.md](TUTORIAL_ATUALIZAR_SINONIMOS.md).

Pynnotate ships with `pynnotate/data/synonyms_database.json`. The file contains
the four official dictionaries (`animal_mito`, `animal_nuclear`, `plant_mito`,
and `plant_chloro`) and remains available when the computer is offline. The
`other` genome type continues to use only synonyms entered by the user.

The update button reads `pynnotate/data/synonyms_manifest.json` from the main
GitHub branch. It installs a candidate only when its schema and version are
valid and its raw-file SHA-256 matches the manifest. Installation uses an
atomic file replacement in the user's application-data directory; it never
modifies the executable.

## Publishing an official update

1. Edit `pynnotate/data/synonyms_database.json`.
2. Increase `database_revision` by exactly one and change
   `dictionary_version` to a human-readable version/date.
3. Run `python scripts/rebuild_synonym_manifest.py`.
4. Run `python -m unittest discover -s tests -v`.
5. Commit and publish both JSON files together on the `main` branch.

The manifest URL compiled into the app is:

`https://raw.githubusercontent.com/fernandacaron/pynnotate/main/pynnotate/data/synonyms_manifest.json`

Users with an older database revision will then be offered the whole updated
base. Updates are never downloaded without the user's click and confirmation.

## Local macOS test before publishing

Create a harmless revision containing `PYNNOTATE-TEST-ALIAS`:

```bash
python scripts/create_test_synonym_update.py /tmp/pynnotate-synonym-test
cd /tmp/pynnotate-synonym-test
python -m http.server 8765
```

Leave that terminal open. In a second terminal, from the project root, launch
the source GUI against the localhost manifest:

```bash
PYNNOTATE_SYNONYM_MANIFEST_URL=http://127.0.0.1:8765/synonyms_manifest.json \
python -m pynnotate.pynnotate_graphic
```

For a compiled macOS application:

```bash
PYNNOTATE_SYNONYM_MANIFEST_URL=http://127.0.0.1:8765/synonyms_manifest.json \
dist/Pynnotate.app/Contents/MacOS/Pynnotate
```

In Step 1, click **Check for updates**, confirm installation, and verify that
the status changes to `local-test-update`. Then click **Restore included**.

Local HTTP is accepted only for `localhost`, `127.0.0.1`, and `::1`; normal
update servers must use HTTPS.
