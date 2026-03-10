import unittest
import pandas as pd
import asyncio

from BioTools.protein_annotation import (
    _parse_pdb_data,
    _parse_uniprot_data,
    protein_results_to_dataframe,
    get_proteins_info,
)


class TestProteinAnnotationParsers(unittest.TestCase):
    def test_parse_uniprot_data_extracts_expected_fields(self):
        data = {
            "gene": [{"name": {"value": "TP53"}}],
            "organism": {"taxonomy": 9606},
            "comments": [
                {"type": "FUNCTION", "text": [{"value": "Tumor suppressor activity."}]}
            ],
            "dbReferences": [
                {"type": "GO", "id": "GO:0003677", "properties": {"term": "DNA binding"}},
                {"type": "GO", "id": "GO:0006355", "properties": {"term": "Regulation of transcription"}},
            ],
            "sequence": {"sequence": "MEEPQSDPSV"},
        }

        result = _parse_uniprot_data(data)

        self.assertEqual(result["Gene"], "TP53")
        self.assertEqual(result["TaxID"], 9606)
        self.assertEqual(result["Annotation"], "Tumor suppressor activity.")
        self.assertEqual(result["Sequence"], "MEEPQSDPSV")
        self.assertIn("GO:0003677", result["GO_terms"])

    def test_parse_pdb_data_deduplicates_ids(self):
        data = {
            "P04637": [
                {"pdb_id": "1TUP"},
                {"pdb_id": "1TUP"},
                {"pdb_id": "2OCJ"},
            ]
        }
        pdb_ids = _parse_pdb_data(data, "P04637")
        self.assertEqual(set(pdb_ids), {"1TUP", "2OCJ"})

    def test_protein_results_to_dataframe_flattens_nested_fields(self):
        results = [
            {
                "UniProtID": "P04637",
                "Gene": "TP53",
                "TaxID": 9606,
                "Annotation": "Tumor suppressor activity.",
                "GO_terms": {"GO:0003677": "DNA binding"},
                "Sequence": "MEEPQSDPSV",
                "PDB": ["1TUP", "2OCJ"],
            }
        ]

        df = protein_results_to_dataframe(results, set_index=True, flatten_nested=True)

        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn("GO_ids", df.columns)
        self.assertIn("GO_labels", df.columns)
        self.assertIn("PDB_ids", df.columns)
        self.assertEqual(df.loc["P04637", "GO_ids"], "GO:0003677")
        self.assertEqual(df.loc["P04637", "GO_labels"], "DNA binding")

    def test_protein_results_to_dataframe_handles_empty_results(self):
        df = protein_results_to_dataframe([])
        self.assertTrue(df.empty)

    def test_get_proteins_info_validates_max_concurrent(self):
        with self.assertRaises(ValueError):
            asyncio.run(get_proteins_info([], max_concurrent=0))

    def test_get_proteins_info_validates_request_timeout(self):
        with self.assertRaises(ValueError):
            asyncio.run(get_proteins_info([], request_timeout=0))


if __name__ == "__main__":
    unittest.main()
