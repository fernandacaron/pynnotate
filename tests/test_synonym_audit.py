import json
import tempfile
import unittest
from pathlib import Path

from scripts.audit_unknown_synonyms import (
    build_indexes,
    classify_alias,
    normalize,
    parse_log,
    write_outputs,
)


class SynonymAuditTests(unittest.TestCase):
    def setUp(self):
        self.dictionaries = {
            "animal_mito": {
                "12S": ["12S RIBOSOMAL RNA"],
                "16S": ["16S RIBOSOMAL RNA"],
                "tRNA-Cys": ["TRNC", "TRNC(GCA)"],
            },
            "plant_mito": {"rrn18": ["18S RIBOSOMAL RNA"]},
            "plant_chloro": {"rbcL": ["RBCL"]},
        }
        self.indexes = build_indexes(self.dictionaries)

    def test_normalization_is_case_insensitive_and_nfkc(self):
        self.assertEqual(normalize(" 12s ribosomal rna "), "12S RIBOSOMAL RNA")

    def test_alias_in_another_dictionary_is_not_proposed_as_new(self):
        result = classify_alias(
            "12s ribosomal rna", "plant_mito", self.dictionaries, self.indexes
        )
        self.assertEqual(result["status"], "exists_in_other_dictionary")
        self.assertEqual(result["suggested_gene"], "12S")

    def test_trna_symbol_gets_conservative_suggestion(self):
        result = classify_alias("TRNC(TGC)", "animal_mito", self.dictionaries, self.indexes)
        self.assertEqual(result["status"], "suggested_high_confidence")
        self.assertEqual(result["suggested_gene"], "tRNA-Cys")

    def test_target_specific_trna_rule_wins_over_other_dictionary_name(self):
        dictionaries = {
            **self.dictionaries,
            "plant_chloro": {"trnT": ["TRNT"]},
        }
        dictionaries["animal_mito"] = {
            **dictionaries["animal_mito"],
            "tRNA-Thr": ["TRNA-THR"],
        }
        result = classify_alias(
            "TRNT", "animal_mito", dictionaries, build_indexes(dictionaries)
        )
        self.assertEqual(result["status"], "suggested_high_confidence")
        self.assertEqual(result["suggested_gene"], "tRNA-Thr")

    def test_pseudogene_is_review_only(self):
        result = classify_alias(
            "PSEUDO TRNS(AGY)", "animal_mito", self.dictionaries, self.indexes
        )
        self.assertEqual(result["status"], "review_pseudogene")

    def test_log_parser_reads_metadata_and_records(self):
        content = """Genome type chosen: plant_mito
Synonym database revision: 1

[UNKNOWN ALIASES] Some aliases were not detected:
12S RIBOSOMAL RNA (from records: ABC123, DEF456)

[REPEATED GENE ANNOTATIONS]
"""
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "log.txt"
            path.write_text(content, encoding="utf-8")
            metadata, aliases = parse_log(path)
        self.assertEqual(metadata["genome_type"], "plant_mito")
        self.assertEqual(aliases[0].records, ("ABC123", "DEF456"))

    def test_likely_wrong_implicit_target_suppresses_proposals(self):
        database = {
            "database_revision": 1,
            "dictionary_version": "test",
            "dictionaries": self.dictionaries,
        }
        rows = [{
            "alias": "TRNC(TGC)",
            "record_count": 10,
            "status": "suggested_high_confidence",
            "suggested_gene": "rrn18",
            "confidence": 0.99,
            "found_in": "",
            "reason": "test",
            "sample_records": "ABC123",
            "source_logs": "log.txt",
        }]
        coverage = {
            "animal_mito": {"aliases": 1, "records": 10},
            "plant_mito": {"aliases": 0, "records": 0},
            "plant_chloro": {"aliases": 0, "records": 0},
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            write_outputs(
                Path(temp_dir), [Path("log.txt")], "plant_mito", database,
                rows, coverage, ["plant_mito"], target_was_explicit=False,
            )
            proposal = json.loads(
                (Path(temp_dir) / "synonym_proposals.json").read_text(encoding="utf-8")
            )
        self.assertTrue(proposal["likely_wrong_target"])
        self.assertEqual(proposal["aliases_to_add"], {})


if __name__ == "__main__":
    unittest.main()
