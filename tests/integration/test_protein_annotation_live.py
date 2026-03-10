import asyncio
import os
import unittest

import aiohttp

from BioTools.protein_annotation import get_proteins_info


# 1) Add IDs here for probing live API.
PROBE_UNIPROT_IDS = [
    "P04637",  # TP53
    "P00533",  # EGFR
    "P38398"
]


EXPECTED_GENE_BY_UID = {
    "P04637": "TP53",
    "P00533": "EGFR",
    "P38398": "BRCA1",

}



def _run_live_protein_info(uniprot_ids):
    try:
        return asyncio.run(
            get_proteins_info(
                uniprot_ids,
                max_concurrent=5,
                request_timeout=15.0,
                per_request_retries=2,
                per_request_retry_delay=0.5,
            )
        )
    except (aiohttp.ClientError, OSError) as exc:
        raise unittest.SkipTest(f"Live API is unavailable from this environment: {exc}")


@unittest.skipUnless(
    os.getenv("RUN_LIVE_TESTS") == "1",
    "Live tests are disabled. Set RUN_LIVE_TESTS=1 to enable.",
)
class TestProteinAnnotationLive(unittest.TestCase):
    def test_probe_uid_to_gene_mapping(self):
        if not PROBE_UNIPROT_IDS:
            self.skipTest("PROBE_UNIPROT_IDS is empty.")

        results, _ = _run_live_protein_info(PROBE_UNIPROT_IDS)

        actual_map = {item.get("UniProtID"): item.get("Gene") for item in results}
        print("\nResolved mapping from UniProt API:")
        for uid in PROBE_UNIPROT_IDS:
            print(f"  {uid}: {actual_map.get(uid)}")

        self.assertIsInstance(results, list)

    def test_expected_uid_and_gene_regression(self):
        if not EXPECTED_GENE_BY_UID:
            self.skipTest("EXPECTED_GENE_BY_UID is empty. Fill it after probing.")

        ids = list(EXPECTED_GENE_BY_UID.keys())
        results, _ = _run_live_protein_info(ids)
        actual_map = {item.get("UniProtID"): item.get("Gene") for item in results}

        missing_ids = [uid for uid in ids if uid not in actual_map]
        self.assertEqual(missing_ids, [], f"Missing UniProtIDs in results: {missing_ids}")

        mismatches = {}
        for uid, expected_gene in EXPECTED_GENE_BY_UID.items():
            actual_gene = actual_map.get(uid)
            if actual_gene != expected_gene:
                mismatches[uid] = {"expected": expected_gene, "actual": actual_gene}

        self.assertEqual(mismatches, {}, f"Gene mismatches: {mismatches}")


if __name__ == "__main__":
    unittest.main()
