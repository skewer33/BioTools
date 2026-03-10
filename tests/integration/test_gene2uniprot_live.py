import asyncio
import os
import unittest

import aiohttp

from BioTools.gene2uniprot import gene2uniprotid


PROBE_GENES = [
    "TP53",
    "EGFR",
    "BRCA1",
    'LOH2',
    'KEK'
]

EXPECTED_UIDS = {
    "TP53": "P04637",
    "EGFR": "P00533",
    "BRCA1": "P38398",
    "LOH2": None,
    "KEK": None,
}

EXPECTED_ERRORS = ['LOH2', 'KEK']


def _run_live_gene2uniprot(genes, taxid=9606):
    try:
        return asyncio.run(
            gene2uniprotid(
                genes,
                taxid=taxid,
                max_cycle=3,
                retry_delay=0.5,
            )
        )
    except (aiohttp.ClientError, OSError) as exc:
        raise unittest.SkipTest(f"Live API is unavailable from this environment: {exc}")


@unittest.skipUnless(
    os.getenv("RUN_LIVE_TESTS") == "1",
    "Live tests are disabled. Set RUN_LIVE_TESTS=1 to enable.",
)
class TestGene2UniProtLive(unittest.TestCase):
    def test_probe_genes_print_mapping(self):
        if not PROBE_GENES:
            self.skipTest("PROBE_GENES is empty.")

        found, errors = _run_live_gene2uniprot(PROBE_GENES)

        print("\nResolved mapping from mygene.info:")
        for gene in PROBE_GENES:
            print(f"  {gene}: {found.get(gene)}")
        if errors:
            print(f"Unresolved genes: {errors}")

        self.assertIsInstance(found, dict)
        self.assertIsInstance(errors, list)

    def test_expected_mapping_regression(self):
        if not EXPECTED_UIDS:
            self.skipTest("EXPECTED_UIDS is empty. Fill it after probing.")

        genes = list(EXPECTED_UIDS.keys())
        found, errors = _run_live_gene2uniprot(genes)

        self.assertEqual(errors, EXPECTED_ERRORS, f"Some genes were not resolved: {errors}")

        mismatches = {}
        for gene, expected_uid in EXPECTED_UIDS.items():
            actual_uid = found.get(gene)
            if actual_uid != expected_uid:
                mismatches[gene] = {"expected": expected_uid, "actual": actual_uid}

        self.assertEqual(mismatches, {}, f"UID mismatches: {mismatches}")



if __name__ == "__main__":
    unittest.main()
