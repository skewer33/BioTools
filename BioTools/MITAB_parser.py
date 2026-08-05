import re
import pandas as pd


def Check_Value(val: str | float | int, valid_values: set, valname: str, message='Wrong value123'):
    """
    Function for check correctness of input value

    Parameters
    ----------
    val : str, float, int
        input value
    valid_values : set
        set of valid values
    valname : str
        group name of valid_values set. Or name of val variable
    message : str
        Error message

    Returns
    -------
    None
    """
    if val not in valid_values:
        if message == 'Wrong value123':
            message = f'Wrong value of "{valname}" variable! Choose one of {valid_values}'
        raise Exception(message)


class MITAB_parser():

    MITAB_columns = {'#ID(s) interactor A', 'ID(s) interactor B', 'Alt. ID(s) interactor A',
       'Alt. ID(s) interactor B', 'Alias(es) interactor A',
       'Alias(es) interactor B', 'Interaction detection method(s)',
       'Publication 1st author(s)', 'Publication Identifier(s)',
       'Taxid interactor A', 'Taxid interactor B', 'Interaction type(s)',
       'Source database(s)', 'Interaction identifier(s)',
       'Confidence value(s)', 'Expansion method(s)',
       'Biological role(s) interactor A', 'Biological role(s) interactor B',
       'Experimental role(s) interactor A',
       'Experimental role(s) interactor B', 'Type(s) interactor A',
       'Type(s) interactor B', 'Xref(s) interactor A', 'Xref(s) interactor B',
       'Interaction Xref(s)', 'Annotation(s) interactor A',
       'Annotation(s) interactor B', 'Interaction annotation(s)',
       'Host organism(s)', 'Interaction parameter(s)', 'Creation date',
       'Update date', 'Checksum(s) interactor A', 'Checksum(s) interactor B',
       'Interaction Checksum(s)', 'Negative', 'Feature(s) interactor A',
       'Feature(s) interactor B', 'Stoichiometry(s) interactor A',
       'Stoichiometry(s) interactor B', 'Identification method participant A',
       'Identification method participant B',
       'Biological effect(s) interactor A',
       'Biological effect(s) interactor B', 'Causal regulatory mechanism',
       'Causal statement'}

    # All data types that can be extracted from a MITAB table.
    valid_parsing_data = {'protein_id', 'taxid', 'publications', 'detection_method'}

    # Columns required to extract each data type for interactors A and B.
    # Note: 'publications' and 'detection_method' are interaction-level columns,
    # so they are listed once (in data_path_A) and reused for validation.
    data_path_A = {'protein_id': ['#ID(s) interactor A', 'Alt. ID(s) interactor A', 'Alias(es) interactor A'],
                   'taxid': ['Taxid interactor A'],
                   'publications': ['Publication Identifier(s)'],
                   'detection_method': ['Interaction detection method(s)']}

    data_path_B = {'protein_id': ['ID(s) interactor B', 'Alt. ID(s) interactor B', 'Alias(es) interactor B'],
                   'taxid': ['Taxid interactor B'],
                   'publications': ['Publication Identifier(s)'],
                   'detection_method': ['Interaction detection method(s)']}

    def __init__(self, df, parsing_data = ['protein_id']):
        # Normalize a common MITAB header variant: '# ID(s) interactor A' -> '#ID(s) interactor A'
        self.df = df.rename(columns={'# ID(s) interactor A': '#ID(s) interactor A'})
        # instructions for parsing
        self.get_data = {'protein_id': self.get_UID_Gene_from_mitab,
                                'taxid': self.get_taxid_from_mitab,
                                'publications': self.get_publications_from_mitab,
                                'detection_method': self.get_detection_method_from_mitab}

        self._validate_required_data(parsing_data)
        self.required_data = list(parsing_data)
        self._check_columns()

    def _validate_required_data(self, required_data):
        for datatype in required_data:
            valid_values = self.valid_parsing_data
            Check_Value(datatype, valid_values, valname='',
                        message=f"Valid members of 'required_data' list is {valid_values}.\n\t   Example: required_data=['protein_id', 'detection_method']")

    def _check_columns(self):
        valid_columns = set(self.df.columns)
        necessary_cols = {item for k in self.required_data
                          for item in list(self.data_path_A[k]) + list(self.data_path_B[k])}  # get set of necessary columns for required data
        for column in necessary_cols:
            Check_Value(column, valid_columns, valname='',
                    message=f"Your MITAB Table doesn't contain '{column}'.\nFor current required_data it must contain at least:\n{necessary_cols}.")

    def parse(self):
        """
        Run all requested extraction methods and return their results combined
        horizontally (one DataFrame with all requested columns).
        """
        result_parts = [self.get_data[key]() for key in self.required_data]
        return pd.concat(result_parts, axis=1)

    # ------------------------------------------------------------------
    #  Regex helpers
    # ------------------------------------------------------------------
    _db_pattern = re.compile(r'(uniprotkb|psi-mi|entrez gene/locuslink):([^|]+)')
    # UniProt accessions: 6 or 10 characters, first letter rules + optional isoform suffix '-N'
    _uniprot_pattern = re.compile(r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?$')
    _gene_pattern = re.compile(r'^[A-Za-z][A-Za-z0-9\-]*$')  # Starts with a letter

    @staticmethod
    def _split_entry(value):
        """
        Split a single MITAB db:value(annotation) entry into (raw_value, annotation).
        Example: 'ACTA1(gene name)' -> ('ACTA1', 'gene name')
        """
        m = re.match(r'^(.*?)\(([^()]*)\)$', value)
        if m:
            return m.group(1).strip().strip('"'), m.group(2)
        return value.strip().strip('"'), ''

    @classmethod
    def _parse_entries(cls, s):
        """
        Return a list of (db, raw_value, annotation) tuples parsed from a
        '|'-separated MITAB cell, keeping only recognized databases.
        """
        result = []
        for db, value in cls._db_pattern.findall(str(s)):
            raw, ann = cls._split_entry(value)
            result.append((db, raw, ann))
        return result

    @classmethod
    def _find_uniprot(cls, s):
        for db, raw, _ann in cls._parse_entries(s):
            if db == 'uniprotkb' and cls._uniprot_pattern.match(raw):
                return raw
        return None

    @classmethod
    def _find_gene(cls, s):
        gene_db_priority = {'uniprotkb': 0, 'entrez gene/locuslink': 1, 'psi-mi': 2}
        gene_annotation_priority = {'gene name': 0, 'gene name synonym': 1, 'orf name': 2}
        best, best_rank = None, None
        for db, raw, ann in cls._parse_entries(s):
            if not cls._gene_pattern.match(raw):
                continue
            if cls._uniprot_pattern.match(raw):
                continue  # UniProt accession, not a gene name
            rank = (gene_db_priority.get(db, 99), gene_annotation_priority.get(ann, 99))
            if best_rank is None or rank < best_rank:
                best_rank, best = rank, raw
        return best

    # ------------------------------------------------------------------
    #  Extraction methods
    # ------------------------------------------------------------------
    def get_UID_Gene_from_mitab(self):

        '''
        This function takes a MITAB file and returns a DataFrame with UniProt IDs and Gene names for each interactor.
        The function looks for the following columns in the MITAB file:
        - '#ID(s) interactor A'
        - 'ID(s) interactor B'
        - 'Alt. ID(s) interactor A'
        - 'Alt. ID(s) interactor B'
        - 'Alias(es) interactor A'
        - 'Alias(es) interactor B'

        Parameters:
        df : pd.DataFrame
            The input DataFrame containing the MITAB data.
        Returns:
        result : pd.DataFrame
            A DataFrame containing the UniProt IDs and Gene names for each interactor.
            The columns are:
                - 'UniProtID_A'
                - 'UniProtID_B'
                - 'Gene_A'
                - 'Gene_B'
        '''

        target_columns = {'#ID(s) interactor A': 'ID_A',
                          'ID(s) interactor B': 'ID_B',
                          'Alt. ID(s) interactor A': 'Alt_A',
                          'Alt. ID(s) interactor B': 'Alt_B',
                          'Alias(es) interactor A': 'Alias_A',
                          'Alias(es) interactor B': 'Alias_B'}

        # local copy -> do not destroy self.df (other extraction methods may need it)
        df = self.df.loc[:, list(target_columns.keys())].rename(columns=target_columns)

        # data processing (search priority: ID -> Alt. ID -> Alias)
        result = pd.DataFrame({
            'UniProtID_A': (df['ID_A'].apply(self._find_uniprot)
                            .combine_first(df['Alt_A'].apply(self._find_uniprot))
                            .combine_first(df['Alias_A'].apply(self._find_uniprot))),
            'UniProtID_B': (df['ID_B'].apply(self._find_uniprot)
                            .combine_first(df['Alt_B'].apply(self._find_uniprot))
                            .combine_first(df['Alias_B'].apply(self._find_uniprot))),
            'Gene_A': (df['ID_A'].apply(self._find_gene)
                       .combine_first(df['Alt_A'].apply(self._find_gene))
                       .combine_first(df['Alias_A'].apply(self._find_gene))),
            'Gene_B': (df['ID_B'].apply(self._find_gene)
                       .combine_first(df['Alt_B'].apply(self._find_gene))
                       .combine_first(df['Alias_B'].apply(self._find_gene))),
        })
        return result


    def get_publications_from_mitab(self):

        target_col = 'Publication Identifier(s)'

        def parse_identifiers(input_str):
            # Создаем пустой словарь для результатов
            result = {}

            # Пропускаем пустые/NaN ячейки и заглушки '-'
            if pd.isna(input_str):
                return result
            s = str(input_str).strip()
            if not s or s == '-':
                return result

            # Разбиваем строку на элементы по разделителю "|"
            for item in s.split("|"):
                # Разделяем каждый элемент на ключ и значение по первому вхождению ":"
                parts = item.split(":", 1)
                if len(parts) != 2:
                    continue

                key = parts[0].strip()
                value = parts[1].strip()

                # Добавляем значение в список соответствующего ключа
                if key:
                    result.setdefault(key, []).append(value)

            return result

        result = pd.DataFrame({
            'Publications': self.df[target_col].apply(parse_identifiers),
        })

        return result

    def get_taxid_from_mitab(self, taxid_type='digits'):
        '''
        This function takes a MITAB file and returns a DataFrame with taxid for each interactor.
        The function looks for the following columns in the MITAB file:
        - 'Taxid interactor A'
        - 'Taxid interactor B'

        Parameters:
        df : pd.DataFrame
            The input DataFrame containing the MITAB data.
        taxid_type : str
            The type of taxid to extract. It can be one of the following:
                - 'digits': only digits
                - 'text': text in brackets
                - 'full': everything after taxid:
                Example: 'taxid:333284("Hepatitis C virus genotype 1b (isolate Con1)")'
            If 'digits': 333284
            If 'text': "Hepatitis C virus genotype 1b (isolate Con1"
            If 'full': taxid:333284("Hepatitis C virus genotype 1b (isolate Con1)"

        Returns:
        result : pd.DataFrame
            A DataFrame containing the taxid for each interactor.
            The columns are:
                - 'taxid_A'
                - 'taxid_B'
        '''
        if taxid_type not in ['digits', 'text', 'full']:
            raise ValueError("taxid_type must be one of ['digits', 'text', 'full']")

        target_columns = ['Taxid interactor A', 'Taxid interactor B']
        taxid_patterns = {
            'digits': r'taxid:(\d+)',  # only digits
            'text': r'taxid:\d+\s*\(([^)]+)',  # text in brackets
            'full': r'taxid:(.*)'  # everything after taxid:
        }
        pattern = re.compile(taxid_patterns[taxid_type])

        def find_taxid(row, pattern):
            if pd.isna(row):
                return None
            match = re.search(pattern, str(row))
            if match:
                return match.group(1)
            else:
                return None

        # data processing (self.df is left untouched)
        result = pd.DataFrame({
            'taxid_A': self.df[target_columns[0]].apply(lambda x: find_taxid(x, pattern)),
            'taxid_B': self.df[target_columns[1]].apply(lambda x: find_taxid(x, pattern))
        })

        return result

    def get_detection_method_from_mitab(self):
        '''
        This function extracts the interaction detection method(s) from the
        'Interaction detection method(s)' column.

        The cell contains one or more '|'-separated PSI-MI ontology terms, e.g.:
            psi-mi:"MI:0398"(two hybrid pooling approach)

        Parameters:
        df : pd.DataFrame
            The input DataFrame containing the MITAB data.

        Returns:
        result : pd.DataFrame
            A DataFrame with the extracted methods. The columns are:
                - 'Detection_method_id': PSI-MI identifier(s), e.g. 'MI:0398'
                - 'Detection_method_name': human-readable name(s), e.g. 'two hybrid pooling approach'
            If a row contains several methods, they are joined with '|'.
        '''
        target_col = 'Interaction detection method(s)'
        method_pattern = re.compile(r'psi-mi:"?MI:(\d+)"?\s*\(([^)]*)\)')

        def parse_methods(input_str):
            if pd.isna(input_str):
                return None, None
            matches = method_pattern.findall(str(input_str))
            if not matches:
                return None, None
            ids = '|'.join('MI:%s' % num for num, _name in matches)
            names = '|'.join(name for _num, name in matches)
            return ids, names

        result = pd.DataFrame({
            'Detection_method_id': self.df[target_col].apply(lambda x: parse_methods(x)[0]),
            'Detection_method_name': self.df[target_col].apply(lambda x: parse_methods(x)[1]),
        })

        return result
