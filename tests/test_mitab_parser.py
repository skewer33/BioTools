import unittest

import pandas as pd

from BioTools.MITAB_parser import Check_Value, MITAB_parser


def make_test_df():
    """Build a small self-contained MITAB-like DataFrame (PSI-MI TAB 2.7 format)."""
    return pd.DataFrame({
        '# ID(s) interactor A': ['uniprotkb:P68133', 'uniprotkb:Q6ZQX7-4'],
        'ID(s) interactor B': ['uniprotkb:Q9UQM7', 'uniprotkb:P68133'],
        'Alt. ID(s) interactor A': ['uniprotkb:P02568', 'intact:EBI-25830459'],
        'Alt. ID(s) interactor B': ['uniprotkb:Q9UL21', 'uniprotkb:P02568'],
        'Alias(es) interactor A': [
            'psi-mi:acts_human(display_short)|uniprotkb:ACTA1(gene name)|uniprotkb:ACTA(gene name synonym)',
            'psi-mi:q6zqx74(display_short)|uniprotkb:C17orf97(gene name synonym)|uniprotkb:LIAT1(gene name)',
        ],
        'Alias(es) interactor B': [
            'psi-mi:kcc2a_human(display_short)|uniprotkb:CAMK2A(gene name)',
            'psi-mi:acts_human(display_short)|uniprotkb:ACTA1(gene name)',
        ],
        'Interaction detection method(s)': [
            'psi-mi:"MI:0398"(two hybrid pooling approach)',
            'psi-mi:"MI:0006"(anti bait coimmunoprecipitation)|psi-mi:"MI:0019"(coimmunoprecipitation)',
        ],
        'Publication Identifier(s)': [
            'intact:EBI-25827495|doi:10.1016/j.celrep.2020.108050|pubmed:32814053',
            '-',
        ],
        'Taxid interactor A': ['taxid:9606(human)|taxid:9606(Homo sapiens)', None],
        'Taxid interactor B': ['taxid:9606(human)', 'taxid:9606(human)'],
    })


class TestCheckValue(unittest.TestCase):
    def test_accepts_valid_value(self):
        Check_Value("human", {"human", "mouse"}, "species")

    def test_raises_for_invalid_value(self):
        with self.assertRaises(Exception):
            Check_Value("yeast", {"human", "mouse"}, "species")


class TestMITABParser(unittest.TestCase):
    def setUp(self):
        self.df = make_test_df()
        self.parser = MITAB_parser(self.df, parsing_data=['protein_id', 'taxid', 'publications', 'detection_method'])

    def test_constructor_normalizes_hash_column(self):
        # '# ID(s) interactor A' must become '#ID(s) interactor A' internally
        self.assertIn('#ID(s) interactor A', self.parser.df.columns)
        self.assertNotIn('# ID(s) interactor A', self.parser.df.columns)

    def test_constructor_rejects_unknown_data(self):
        with self.assertRaises(Exception):
            MITAB_parser(self.df, parsing_data=['bogus'])

    def test_constructor_rejects_missing_column(self):
        bad = self.df.drop(columns=['Publication Identifier(s)'])
        with self.assertRaises(Exception):
            MITAB_parser(bad, parsing_data=['publications'])

    def test_default_parsing_data(self):
        p = MITAB_parser(self.df)
        self.assertEqual(p.required_data, ['protein_id'])

    def test_uid_gene_extraction(self):
        res = self.parser.get_UID_Gene_from_mitab()
        self.assertEqual(res.loc[0, 'UniProtID_A'], 'P68133')
        self.assertEqual(res.loc[0, 'UniProtID_B'], 'Q9UQM7')
        self.assertEqual(res.loc[0, 'Gene_A'], 'ACTA1')   # from '(gene name)' annotation
        self.assertEqual(res.loc[0, 'Gene_B'], 'CAMK2A')
        # isoform accession must be recognized
        self.assertEqual(res.loc[1, 'UniProtID_A'], 'Q6ZQX7-4')
        self.assertEqual(res.loc[1, 'Gene_A'], 'LIAT1')

    def test_taxid_extraction(self):
        res = self.parser.get_taxid_from_mitab()
        self.assertEqual(res.loc[0, 'taxid_A'], '9606')
        self.assertEqual(res.loc[0, 'taxid_B'], '9606')
        self.assertIsNone(res.loc[1, 'taxid_A'])  # NaN cell -> None

    def test_publications_extraction(self):
        res = self.parser.get_publications_from_mitab()
        pubs = res['Publications'].iloc[0]
        self.assertEqual(pubs['pubmed'], ['32814053'])
        self.assertEqual(pubs['doi'], ['10.1016/j.celrep.2020.108050'])
        # '-' placeholder must not raise and must produce an empty dict
        self.assertEqual(res['Publications'].iloc[1], {})

    def test_detection_method_extraction(self):
        res = self.parser.get_detection_method_from_mitab()
        self.assertEqual(res.loc[0, 'Detection_method_id'], 'MI:0398')
        self.assertEqual(res.loc[0, 'Detection_method_name'], 'two hybrid pooling approach')
        # multiple methods are joined with '|'
        self.assertEqual(res.loc[1, 'Detection_method_id'], 'MI:0006|MI:0019')
        self.assertEqual(res.loc[1, 'Detection_method_name'], 'anti bait coimmunoprecipitation|coimmunoprecipitation')

    def test_parse_combines_all_columns(self):
        res = self.parser.parse()
        expected = ['UniProtID_A', 'UniProtID_B', 'Gene_A', 'Gene_B',
                    'taxid_A', 'taxid_B', 'Publications',
                    'Detection_method_id', 'Detection_method_name']
        self.assertEqual(list(res.columns), expected)
        self.assertEqual(len(res), len(self.df))

    def test_df_not_mutated_by_extraction(self):
        p = MITAB_parser(self.df, parsing_data=['protein_id', 'taxid', 'publications', 'detection_method'])
        columns_before = list(p.df.columns)  # columns after constructor normalization
        p.get_UID_Gene_from_mitab()
        p.get_taxid_from_mitab()
        p.get_publications_from_mitab()
        p.get_detection_method_from_mitab()
        self.assertEqual(list(p.df.columns), columns_before)


if __name__ == "__main__":
    unittest.main()
