#!/usr/bin/env python3
"""Build the versioned synonym database and its update manifest."""

import hashlib
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from pynnotate.main.alias_maps import (
    alias_map_animal,
    alias_map_chloroplast,
    alias_map_mito_plant,
    alias_map_nuclear_animal,
)


DATA_DIR = PROJECT_ROOT / "pynnotate" / "data"
DATABASE_PATH = DATA_DIR / "synonyms_database.json"
MANIFEST_PATH = DATA_DIR / "synonyms_manifest.json"
DATABASE_URL = (
    "https://raw.githubusercontent.com/fernandacaron/pynnotate/"
    "main/pynnotate/data/synonyms_database.json"
)


def build_database():
    return {
        "schema_version": 1,
        "database_revision": 1,
        "dictionary_version": "2026.08.23",
        "dictionaries": {
            "animal_mito": alias_map_animal,
            "animal_nuclear": alias_map_nuclear_animal,
            "plant_mito": alias_map_mito_plant,
            "plant_chloro": alias_map_chloroplast,
        },
    }


def main():
    database = build_database()
    database_bytes = (
        json.dumps(database, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    DATABASE_PATH.write_bytes(database_bytes)

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
    print(f"Wrote {DATABASE_PATH}")
    print(f"Wrote {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
