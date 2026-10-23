#!/usr/bin/env python3
"""Create a localhost-only synonym update for manual GUI testing."""

import argparse
import copy
import hashlib
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from pynnotate.main.synonym_updates import get_bundled_database


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", nargs="?", default="/tmp/pynnotate-synonym-test")
    args = parser.parse_args()
    output_dir = Path(args.output).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    database = copy.deepcopy(get_bundled_database())
    database["database_revision"] += 1
    database["dictionary_version"] = "local-test-update"
    database["dictionaries"]["animal_mito"]["COI"].append("PYNNOTATE-TEST-ALIAS")
    database_bytes = (
        json.dumps(database, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    (output_dir / "synonyms_database.json").write_bytes(database_bytes)
    manifest = {
        "schema_version": 1,
        "database_revision": database["database_revision"],
        "dictionary_version": database["dictionary_version"],
        "database_url": "synonyms_database.json",
        "sha256": hashlib.sha256(database_bytes).hexdigest(),
    }
    (output_dir / "synonyms_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(output_dir)


if __name__ == "__main__":
    main()
