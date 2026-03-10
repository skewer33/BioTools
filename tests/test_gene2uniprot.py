import unittest

from BioTools.gene2uniprot import _find_UID


class TestFindUID(unittest.TestCase):
    def test_returns_swiss_prot_id(self):
        payload = {"hits": [{"uniprot": {"Swiss-Prot": "P04637"}}]}
        self.assertEqual(_find_UID(payload), "P04637")

    def test_falls_back_to_trembl_and_takes_first_item(self):
        payload = {"hits": [{"uniprot": {"Swiss-Prot": None, "TrEMBL": ["A0A024RBG1"]}}]}
        self.assertEqual(_find_UID(payload), "A0A024RBG1")

    def test_returns_none_for_missing_hits(self):
        payload = {"hits": []}
        self.assertIsNone(_find_UID(payload))


if __name__ == "__main__":
    unittest.main()
