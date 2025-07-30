import re
import pandas as pd


def Check_Value(val:[str, float, int], valid_values:set, valname:str, message='Wrong value123'):
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
    
    valid_parsing_data = {'protein_id', 'taxid', 'publications'}
    
    def __init__(self, df, parsing_data = ['UniProtID']):
        self.df = df
        # instructions for parsing
        self.get_data = {'protein_id': self.get_UID_Gene_from_mitab, 
                                'taxid': self.get_taxid_from_mitab,
                                'publications': self.get_publications_from_mitab}
        
        self._validate_required_data(parsing_data)
        self.required_data = parsing_data
        self._check_columns()
        
    def _validate_required_data(self, required_data):
        for datatype in required_data:
            valid_values = set(self.data_path_A.keys())
            Check_Value(datatype, valid_values, valname='', 
                        message=f"Valid members of 'required_data' list is {valid_values}.\n\t   Example: required_data=['UniProtID', 'GO_terms']")
    
    def _check_columns(self):
        valid_columns = set(self.data.columns)
        necessary_cols = {item for k in self.required_data for item in list(self.data_path_A[k].keys()) + list(self.data_path_B[k].keys())} # get set of necessary columns for required data
        for column in necessary_cols:
            Check_Value(column, valid_columns, valname='',
                    message=f"Your MITAB Table isn`t contain '{column}'.\nFor current required_data it must contain at least:\n{necessary_cols}.")
            
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
        
        target_columns = ['#ID(s) interactor A', 'ID(s) interactor B', 'Alias(es) interactor A', 'Alias(es) interactor B']
        
        
        # Regular expression for DB pattern
        db_pattern = re.compile(r'(uniprotkb|psi-mi|entrez gene/locuslink):([^|]+)')
        
        # UniProtID and Gene Patterns
        uniprot_pattern = re.compile(r'^[A-Z][A-Z0-9]{5,9}$')  # 6-10 символов, первая буква
        gene_pattern = re.compile(r'^[A-Za-z][A-Za-z0-9\-]*$')  # Начинается с буквы

        # DB priority for Gene
        gene_db_priority = ['uniprotkb', 'psi-mi', 'entrez gene/locuslink']

        # function for finding UniProtID
        def find_uniprot(s):
            for db, value in db_pattern.findall(str(s)):
                if db == 'uniprotkb' and uniprot_pattern.match(value):
                    return value
            return None

        # function for finding Gene
        def find_gene(s):
            entries = db_pattern.findall(str(s))
            # search by DB priority
            for db in gene_db_priority:
                for entry_db, value in entries:
                    if entry_db == db and gene_pattern.match(value):
                        return value
            return None

        # transform data
        renamed_columns = ['ID_A', 'ID_B', 'Alias_A', 'Alias_B']
        self.df = self.df.loc[:, target_columns].rename(columns={c_old:c_new for c_old, c_new in zip(target_columns, renamed_columns)})

        # data processing
        result = pd.DataFrame({
            # UniProtID from ID and Alias
            'UniProtID_A': self.df['ID_A'].apply(find_uniprot).combine_first(self.df['Alias_A'].apply(find_uniprot)),
            'UniProtID_B': self.df['ID_B'].apply(find_uniprot).combine_first(self.df['Alias_B'].apply(find_uniprot)),
            
            # Gene from ID and Alias (ID priority)
            'Gene_A': self.df['ID_A'].apply(find_gene).combine_first(self.df['Alias_A'].apply(find_gene)),
            'Gene_B': self.df['ID_B'].apply(find_gene).combine_first(self.df['Alias_B'].apply(find_gene))
        })
        return result
    
    
    def get_publication_from_mitab(self):

        target_col = 'Publication Identifier(s)'

        def parse_identifiers(input_str):
            # Создаем пустой словарь для результатов
            result = {}
            
            # Разбиваем строку на элементы по разделителю "|"
            for item in input_str.split("|"):
                # Разделяем каждый элемент на ключ и значение по первому вхождению ":"
                key, value = re.split(r":", item, 1)
                
                # Удаляем возможные пробелы вокруг ключа и значения
                key = key.strip()
                value = value.strip()
                
                # Добавляем значение в список соответствующего ключа
                if key in result:
                    result[key].append(value)
                else:
                    result[key] = [value]
            
            return result
        
        result = pd.DataFrame({
            'Publications': self.df[target_col].apply(lambda x: parse_identifiers(x)),
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
            match = re.search(pattern, row)
            if match:
                return match.group(1)
            else:
                return None
        
        
        # transform data
        self.df = self.df.loc[:, target_columns]
        
        # data processing
        result = pd.DataFrame({
            'taxid_A': self.df[target_columns[0]].apply(lambda x: find_taxid(x, pattern)),
            'taxid_B': self.df[target_columns[1]].apply(lambda x: find_taxid(x, pattern))
        })
        
        return result