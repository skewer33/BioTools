import unittest

from BioTools.protein_annotation import _parse_pdb_data, _parse_uniprot_data


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


if __name__ == "__main__":
    unittest.main()
