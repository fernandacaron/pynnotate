#!/usr/bin/env python3
"""Validate synonyms_database.json and refresh its published checksum manifest."""

import hashlib
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from pynnotate.main.synonym_updates import validate_synonym_database


DATA_DIR = PROJECT_ROOT / "pynnotate" / "data"
DATABASE_PATH = DATA_DIR / "synonyms_database.json"
MANIFEST_PATH = DATA_DIR / "synonyms_manifest.json"
DATABASE_URL = (
    "https://raw.githubusercontent.com/fernandacaron/pynnotate/"
    "main/pynnotate/data/synonyms_database.json"
)


def main():
    database_bytes = DATABASE_PATH.read_bytes()
    database = validate_synonym_database(json.loads(database_bytes.decode("utf-8")))
    manifest = {
        "schema_version": 1,
        "database_revision": database["database_revision"],
        "dictionary_version": database["dictionary_version"],
        "database_url": DATABASE_URL,
        "sha256": hashlib.sha256(database_bytes).hexdigest(),
    }
    MANIFEST_PATH.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Validated revision {database['database_revision']}")
    print(f"Updated {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
