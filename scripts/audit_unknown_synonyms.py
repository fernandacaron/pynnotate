#!/usr/bin/env python3
"""Audit Pynnotate [UNKNOWN ALIASES] logs without changing the database.

The script uses the same case-insensitive NFKC normalization as Pynnotate,
checks every bundled genome dictionary, and writes review-friendly CSV, JSON,
and Markdown files.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import unicodedata
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from difflib import SequenceMatcher
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = PROJECT_ROOT / "pynnotate" / "data" / "synonyms_database.json"
GENOME_TYPES = ("animal_mito", "animal_nuclear", "plant_mito", "plant_chloro")

TRNA_ONE_LETTER = {
    "A": "tRNA-Ala",
    "R": "tRNA-Arg",
    "N": "tRNA-Asn",
    "D": "tRNA-Asp",
    "C": "tRNA-Cys",
    "Q": "tRNA-Gln",
    "E": "tRNA-Glu",
    "G": "tRNA-Gly",
    "H": "tRNA-His",
    "I": "tRNA-Ile",
    "L": "tRNA-Leu",
    "K": "tRNA-Lys",
    "M": "tRNA-Met",
    "F": "tRNA-Phe",
    "P": "tRNA-Pro",
    "S": "tRNA-Ser",
    "T": "tRNA-Thr",
    "W": "tRNA-Trp",
    "Y": "tRNA-Tyr",
    "V": "tRNA-Val",
}


@dataclass
class UnknownAlias:
    alias: str
    records: tuple[str, ...]


def normalize(text: str) -> str:
    """Match Pynnotate's comparison: Unicode NFKC, trim, uppercase."""

    return unicodedata.normalize("NFKC", text).strip().upper()


def relaxed(text: str) -> str:
    """Comparison key used only to suggest punctuation/spacing variants."""

    return re.sub(r"[^A-Z0-9]+", "", normalize(text))


def parse_log(path: Path) -> tuple[dict[str, str], list[UnknownAlias]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    metadata = {}
    metadata_patterns = {
        "genome_type": r"^Genome type chosen:\s*(.+?)\s*$",
        "database_version": r"^Synonym database version:\s*(.+?)\s*$",
        "database_revision": r"^Synonym database revision:\s*(.+?)\s*$",
    }
    for key, pattern in metadata_patterns.items():
        match = re.search(pattern, text, flags=re.MULTILINE)
        if match:
            metadata[key] = match.group(1)

    marker = "[UNKNOWN ALIASES]"
    if marker not in text:
        return metadata, []

    section = text.split(marker, 1)[1]
    aliases = []
    line_pattern = re.compile(r"^(.*?) \(from records: (.*)\)$")
    for line in section.splitlines()[1:]:
        if line.startswith("["):
            break
        match = line_pattern.match(line)
        if not match:
            continue
        records = tuple(item.strip() for item in match.group(2).split(",") if item.strip())
        aliases.append(UnknownAlias(match.group(1).strip(), records))
    return metadata, aliases


def load_database(path: Path) -> dict:
    database = json.loads(path.read_text(encoding="utf-8"))
    dictionaries = database.get("dictionaries")
    if not isinstance(dictionaries, dict) or not dictionaries:
        raise ValueError("The synonym database has no 'dictionaries' object.")
    return database


def build_indexes(dictionaries: dict) -> dict[str, dict[str, set[str]]]:
    indexes = {}
    for genome_type, genes in dictionaries.items():
        index = defaultdict(set)
        for gene, aliases in genes.items():
            for alias in [gene, *aliases]:
                index[normalize(alias)].add(gene)
        indexes[genome_type] = dict(index)
    return indexes


def build_relaxed_index(genes: dict) -> dict[str, set[tuple[str, str]]]:
    index = defaultdict(set)
    for gene, aliases in genes.items():
        for alias in [gene, *aliases]:
            index[relaxed(alias)].add((gene, alias))
    return dict(index)


def infer_gene(alias: str, target_genes: dict) -> tuple[str, float, str]:
    """Return a conservative candidate; never applies it automatically."""

    value = normalize(alias)
    if "PSEUDO" in value:
        return "", 0.0, "Pseudogene annotation; review instead of adding automatically."

    rrna_rules = (
        ("12S", (r"12\s*S", r"RRN\s*12", r"^RRNS$", r"^RNS$", r"^S[- _]?RRNA$", r"SMALL .*RIBOSOMAL")),
        ("16S", (r"16\s*S", r"RRN\s*16", r"^RRNL$", r"^RNL$", r"^L[- _]?RRNA$", r"LARGE .*RIBOSOMAL")),
    )
    for gene, patterns in rrna_rules:
        if gene in target_genes and any(re.search(pattern, value) for pattern in patterns):
            return gene, 0.99, "Recognized rRNA nomenclature pattern."

    # Standard mitochondrial tRNA symbols use TRN + one-letter amino-acid code.
    # Do not interpret the generic word TRNA as tRNA-Ala.
    if value.startswith("TRNSER") and "tRNA-Ser" in target_genes:
        return "tRNA-Ser", 0.97, "Recognized spelled-out TRNSER gene symbol."
    match = re.match(r"^TRN([A-Z])(?:$|[0-9 (\-_])", value)
    if match and not value.startswith("TRNA"):
        gene = TRNA_ONE_LETTER.get(match.group(1), "")
        if gene in target_genes:
            return gene, 0.97, "Recognized TRN one-letter amino-acid code."

    return "", 0.0, "No conservative biological rule matched."


def fuzzy_candidate(alias: str, genes: dict) -> tuple[str, float, str]:
    query = relaxed(alias)
    scored = []
    for gene, aliases in genes.items():
        best_alias = ""
        best_score = 0.0
        for existing in [gene, *aliases]:
            score = SequenceMatcher(None, query, relaxed(existing)).ratio()
            if score > best_score:
                best_alias, best_score = existing, score
        scored.append((best_score, gene, best_alias))
    scored.sort(reverse=True)
    if not scored:
        return "", 0.0, ""
    best_score, best_gene, best_alias = scored[0]
    second_score = scored[1][0] if len(scored) > 1 else 0.0
    if best_score >= 0.88 and best_score - second_score >= 0.06:
        return best_gene, round(best_score, 3), f"Similar to existing alias '{best_alias}'."
    return "", round(best_score, 3), "No unique high-confidence textual match."


def classify_alias(alias: str, target: str, dictionaries: dict, indexes: dict) -> dict:
    key = normalize(alias)
    target_matches = indexes[target].get(key, set())
    if target_matches:
        genes = sorted(target_matches)
        return {
            "status": "already_in_target",
            "suggested_gene": "; ".join(genes),
            "confidence": 1.0,
            "found_in": target,
            "reason": "Exact match after Pynnotate case-insensitive normalization.",
        }

    relaxed_matches = build_relaxed_index(dictionaries[target]).get(relaxed(alias), set())
    relaxed_genes = sorted({gene for gene, _ in relaxed_matches})
    if len(relaxed_genes) == 1:
        existing = sorted({existing for _, existing in relaxed_matches})[0]
        return {
            "status": "suggested_high_confidence",
            "suggested_gene": relaxed_genes[0],
            "confidence": 0.99,
            "found_in": "",
            "reason": f"Differs from existing alias '{existing}' only by spacing/punctuation.",
        }

    gene, confidence, reason = infer_gene(alias, dictionaries[target])
    if gene:
        return {
            "status": "suggested_high_confidence",
            "suggested_gene": gene,
            "confidence": confidence,
            "found_in": "",
            "reason": reason,
        }
    if "PSEUDO" in normalize(alias):
        return {
            "status": "review_pseudogene",
            "suggested_gene": "",
            "confidence": 0.0,
            "found_in": "",
            "reason": reason,
        }

    # A valid symbol can use a different canonical gene name in another genome
    # dictionary (for example TRNT vs tRNA-Thr). Prefer a strong target-specific
    # biological rule, but retain exact cross-dictionary matches as evidence.
    other_matches = []
    for genome_type, index in indexes.items():
        if genome_type == target:
            continue
        for other_gene in sorted(index.get(key, set())):
            other_matches.append((genome_type, other_gene))
    if other_matches:
        unique_genes = sorted({other_gene for _, other_gene in other_matches})
        return {
            "status": "exists_in_other_dictionary",
            "suggested_gene": "; ".join(unique_genes),
            "confidence": 1.0 if len(unique_genes) == 1 else 0.7,
            "found_in": "; ".join(f"{kind}:{other_gene}" for kind, other_gene in other_matches),
            "reason": "Exact alias exists, but not in the target genome dictionary.",
        }

    gene, confidence, reason = fuzzy_candidate(alias, dictionaries[target])
    if gene:
        return {
            "status": "suggested_review",
            "suggested_gene": gene,
            "confidence": confidence,
            "found_in": "",
            "reason": reason,
        }
    return {
        "status": "unresolved",
        "suggested_gene": "",
        "confidence": confidence,
        "found_in": "",
        "reason": reason,
    }


def coverage_by_dictionary(aliases: list[UnknownAlias], indexes: dict) -> dict[str, dict[str, int]]:
    coverage = {}
    for genome_type, index in indexes.items():
        matching = [row for row in aliases if normalize(row.alias) in index]
        coverage[genome_type] = {
            "aliases": len(matching),
            "records": sum(len(row.records) for row in matching),
        }
    return coverage


def write_outputs(
    output_dir: Path,
    logs: list[Path],
    target: str,
    database: dict,
    rows: list[dict],
    coverage: dict,
    log_types: list[str],
    target_was_explicit: bool,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    total_records = sum(row["record_count"] for row in rows)
    best_type = max(coverage, key=lambda key: coverage[key]["records"])
    best = coverage[best_type]
    target_coverage = coverage[target]["records"]
    likely_wrong_target = bool(
        total_records
        and best_type != target
        and best["records"] / total_records >= 0.75
        and target_coverage / total_records <= 0.25
        and not target_was_explicit
    )

    csv_path = output_dir / "unknown_aliases_audit.csv"
    fields = [
        "alias", "record_count", "status", "suggested_gene", "confidence",
        "found_in", "reason", "sample_records", "source_logs",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)

    aliases_to_add = defaultdict(list)
    if not likely_wrong_target:
        for row in rows:
            if row["status"] == "suggested_high_confidence" and row["suggested_gene"]:
                aliases_to_add[row["suggested_gene"]].append(row["alias"])

    proposal = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_logs": [str(path) for path in logs],
        "log_genome_types": sorted(set(filter(None, log_types))),
        "target_genome": target,
        "database_revision": database.get("database_revision"),
        "dictionary_version": database.get("dictionary_version"),
        "likely_wrong_target": likely_wrong_target,
        "note": (
            "No aliases were proposed because another bundled dictionary matches most occurrences. "
            "Check the genome type and rerun with an explicit --genome-type."
            if likely_wrong_target else
            "Review every proposal before editing synonyms_database.json. This file does not modify the database."
        ),
        "aliases_to_add": {gene: sorted(set(values), key=normalize) for gene, values in sorted(aliases_to_add.items())},
        "review_items": [
            {
                "alias": row["alias"],
                "record_count": row["record_count"],
                "status": row["status"],
                "suggested_gene": row["suggested_gene"],
                "confidence": row["confidence"],
                "reason": row["reason"],
            }
            for row in rows
            if row["status"] not in {"already_in_target", "suggested_high_confidence"}
        ],
    }
    (output_dir / "synonym_proposals.json").write_text(
        json.dumps(proposal, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    counts = defaultdict(int)
    record_counts = defaultdict(int)
    for row in rows:
        counts[row["status"]] += 1
        record_counts[row["status"]] += row["record_count"]
    mismatch = bool(log_types) and any(kind != target for kind in log_types)

    lines = [
        "# Synonym audit report",
        "",
        f"- Target dictionary: `{target}`",
        f"- Aliases audited: {len(rows)}",
        f"- Record occurrences represented: {total_records}",
        f"- Database revision: {database.get('database_revision', 'unknown')}",
        f"- Dictionary version: `{database.get('dictionary_version', 'unknown')}`",
    ]
    if log_types:
        lines.append(f"- Genome type recorded in log(s): `{', '.join(sorted(set(log_types)))}`")
    if mismatch:
        lines.extend((
            "",
            "> **Attention:** the target dictionary differs from the genome type recorded in at least one log.",
        ))
    if likely_wrong_target:
        lines.extend((
            "",
            f"> **Safety stop:** `{best_type}` recognizes most occurrences while `{target}` recognizes very few. `aliases_to_add` was intentionally left empty. Check the search genome type or rerun with an explicit `--genome-type`.",
        ))
    lines.extend(("", "## Classification", ""))
    for status in sorted(counts):
        lines.append(f"- `{status}`: {counts[status]} aliases / {record_counts[status]} record occurrences")

    lines.extend(("", "## Exact coverage of the unknown list by bundled dictionaries", ""))
    for genome_type in GENOME_TYPES:
        if genome_type not in coverage:
            continue
        item = coverage[genome_type]
        percentage = (100 * item["records"] / total_records) if total_records else 0
        lines.append(f"- `{genome_type}`: {item['aliases']} aliases / {item['records']} occurrences ({percentage:.1f}%)")
    if total_records and best["records"] / total_records >= 0.75:
        lines.extend((
            "",
            f"> The dictionary `{best_type}` already recognizes {best['records']} of {total_records} occurrences ({100 * best['records'] / total_records:.1f}%). Check whether the correct genome type was selected before adding synonyms.",
        ))

    lines.extend((
        "",
        "## How to use the files",
        "",
        "1. Filter `unknown_aliases_audit.csv` by `status` and inspect the record examples.",
        "2. Review `aliases_to_add` in `synonym_proposals.json`; nothing is applied automatically.",
        "3. Add only approved aliases to `pynnotate/data/synonyms_database.json`.",
        "4. Increase `database_revision` and `dictionary_version`.",
        "5. Run `python scripts/rebuild_synonym_manifest.py` and commit both JSON files.",
        "",
        "Case alone never requires a new alias: Pynnotate matching is case-insensitive.",
    ))
    (output_dir / "audit_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit unknown aliases from one or more Pynnotate log.txt files."
    )
    parser.add_argument("logs", nargs="+", type=Path, help="Pynnotate log.txt file(s)")
    parser.add_argument(
        "--genome-type", choices=GENOME_TYPES,
        help="Dictionary to audit. Defaults to the type recorded in the first log.",
    )
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--output-dir", type=Path, default=Path("synonym_audit"))
    args = parser.parse_args()

    database = load_database(args.database)
    dictionaries = database["dictionaries"]
    indexes = build_indexes(dictionaries)

    combined = {}
    log_types = []
    for log_path in args.logs:
        metadata, aliases = parse_log(log_path)
        if metadata.get("genome_type"):
            log_types.append(metadata["genome_type"])
        for item in aliases:
            key = normalize(item.alias)
            entry = combined.setdefault(key, {"alias": item.alias, "records": set(), "logs": set()})
            entry["records"].update(item.records)
            entry["logs"].add(log_path.name)

    if not combined:
        parser.error("No [UNKNOWN ALIASES] entries were found in the supplied log(s).")
    target = args.genome_type or (log_types[0] if log_types else "animal_mito")
    if target not in dictionaries:
        parser.error(f"Dictionary '{target}' is not present in the database.")

    unknowns = [UnknownAlias(value["alias"], tuple(sorted(value["records"]))) for value in combined.values()]
    rows = []
    for item in sorted(unknowns, key=lambda row: normalize(row.alias)):
        result = classify_alias(item.alias, target, dictionaries, indexes)
        rows.append({
            "alias": item.alias,
            "record_count": len(item.records),
            **result,
            "sample_records": "; ".join(item.records[:5]),
            "source_logs": "; ".join(sorted(combined[normalize(item.alias)]["logs"])),
        })

    coverage = coverage_by_dictionary(unknowns, indexes)
    write_outputs(
        args.output_dir,
        args.logs,
        target,
        database,
        rows,
        coverage,
        log_types,
        target_was_explicit=args.genome_type is not None,
    )
    print(f"Audited {len(rows)} aliases for {target}")
    print(f"Outputs written to {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
