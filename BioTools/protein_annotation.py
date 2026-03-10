import aiohttp
import asyncio
from tqdm.asyncio import tqdm
import pandas as pd

UNIPROT_API_URL = "https://www.ebi.ac.uk/proteins/api/proteins/"
PDB_API_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/"

async def _get_uniprot_data(uniprot_id, session):
    url = f"{UNIPROT_API_URL}{uniprot_id}"
    try:
        async with session.get(url, headers={"Accept": "application/json"}) as response:
            if response.status == 200:
                return await response.json()
            else:
                #print(f"Ошибка при запросе UniProt ID {uniprot_id}: {response.status}")
                return None
    except Exception as e:
        #print(f"Ошибка соединения для UniProt ID {uniprot_id}: {str(e)}")
        return None

def _parse_uniprot_data(data):
    
    result = {
        "Gene": "N/A",
        "TaxID": "N/A",
        "Annotation": "N/A",
        "GO_terms": {},
        "Sequence": "N/A"
    }
    
    #get gene name
    
    result['Gene'] = data.get("gene", [{}])[0].get("name", {}).get("value", "N/A")
    #get taxid
    result['TaxID'] = data.get("organism", {}).get("taxonomy", "N/A")
    
    # get annotation
    for comment in data.get("comments", []):
        if comment.get("type") == "FUNCTION":
            try:
                texts = comment.get("text", [])
                if texts:
                    result["Annotation"] = texts[0]["value"]
                    break
            except (KeyError, IndexError):
                continue
    
    go_terms = []
    go_descripion = []
    for db_reference in data.get("dbReferences", []):
        if db_reference.get("type") == "GO":
            go_terms.append(db_reference.get("id"))
            go_descripion.append(db_reference.get("properties", {}).get("term"))
    result['GO_terms'] = dict(zip(go_terms, go_descripion))
    
    # Извлечение информации о последовательности
    sequence_info = data.get("sequence", {})
    if sequence_info:
        result["Sequence"] = sequence_info.get("sequence", "N/A")
    return result

async def _get_pdb_structures(uniprot_id, session, error_ids):
    url = f"{PDB_API_URL}{uniprot_id}"
    try:
        async with session.get(url) as response:
            if response.status == 200:
                return await response.json()
            else:
                #print(f"Ошибка при запросе PDB для UniProt ID {uniprot_id}: {response.status}")
                error_ids['PDB'].append(uniprot_id)
                return None
    except Exception as e:
        #print(f"Ошибка соединения для PDB {uniprot_id}: {str(e)}")
        error_ids['PDB'].append(uniprot_id)
        return None

def _parse_pdb_data(data, uniprot_id):
    pdb_structures = []
    if data and uniprot_id in data:
        for entry in data[uniprot_id]:
            pdb_structures.append(entry.get("pdb_id"))
    return list(set(pdb_structures))


def protein_results_to_dataframe(results, set_index=True, flatten_nested=True, list_sep=";"):
    """
    Convert results from get_proteins_info to a pandas DataFrame.

    Parameters
    ----------
    results : list[dict]
        List returned by get_proteins_info.
    set_index : bool
        If True and UniProtID exists, set it as DataFrame index.
    flatten_nested : bool
        If True, add convenient flat string columns for GO and PDB data.
    list_sep : str
        Separator for flattened list-like values.

    Returns
    -------
    pandas.DataFrame
        Tabular representation of protein results.
    """
    df = pd.DataFrame(results)

    if flatten_nested and not df.empty:
        if "GO_terms" in df.columns:
            df["GO_ids"] = df["GO_terms"].apply(
                lambda x: list_sep.join(x.keys()) if isinstance(x, dict) and x else ""
            )
            df["GO_labels"] = df["GO_terms"].apply(
                lambda x: list_sep.join(str(v) for v in x.values()) if isinstance(x, dict) and x else ""
            )
        if "PDB" in df.columns:
            df["PDB_ids"] = df["PDB"].apply(
                lambda x: list_sep.join(str(v) for v in x) if isinstance(x, list) and x else ""
            )

    if set_index and "UniProtID" in df.columns:
        df = df.set_index("UniProtID", drop=False)

    return df

async def _get_protein_info(uniprot_id, session, error_ids):
    uniprot_data = await _get_uniprot_data(uniprot_id, session)
    if not uniprot_data:
        error_ids['UniProtID'].append(uniprot_id)
        return None
    
    try:
        result = _parse_uniprot_data(uniprot_data)
    except Exception as e:
        #print(f"Ошибка парсинга данных UniProt для {uniprot_id}: {str(e)}")
        error_ids['ParseError'].append(uniprot_id)
        return None
    
    # Проверка и сбор ошибок (кроме PDB, которые обрабатываются отдельно)
    for k, v in result.items():
        if v == "N/A" or v == [] or v == {}:
            error_ids[k].append(uniprot_id)
    
    pdb_data = await _get_pdb_structures(uniprot_id, session, error_ids)
    result['PDB'] = _parse_pdb_data(pdb_data, uniprot_id) if pdb_data else []
    result['UniProtID'] = uniprot_id
    
    return result

async def get_proteins_info(uniprot_ids, max_concurrent=10, return_dataframe=False, flatten_nested=True):
    """
    Asynchronously retrieves protein information for a list of UniProt IDs.
    This function queries the UniProt and PDB APIs for each UniProt ID provided, retrieves the corresponding protein information,
    and returns a list of dictionaries containing the protein information, along with a dictionary of error IDs.
    
    Parameters
    ----------
    uniprot_ids : list
        A list of UniProt IDs to query.
    max_concurrent : int
        The maximum number of concurrent requests to the APIs. Default is 10.
    return_dataframe : bool
        If True, return pandas DataFrame instead of list of dictionaries.
    flatten_nested : bool
        Used only when return_dataframe=True. Adds flat GO/PDB columns.
        
    Returns
    -------
    tuple
        A tuple containing:
        - results: list or pandas.DataFrame
            Protein information per UniProt ID.
        - error_ids: dict
            A dictionary where each key is a category of error and each value is a list of UniProt IDs that encountered that error.
    """
    
    error_ids = {
        'UniProtID': [],
        'PDB': [],
        'Gene': [],
        'GO_terms': [],
        'TaxID': [],
        'Annotation': [],
        'ParseError': []
    }
    
    connector = aiohttp.TCPConnector(limit=max_concurrent)
    async with aiohttp.ClientSession(connector=connector) as session:
        tasks = [_get_protein_info(uid, session, error_ids) for uid in tqdm(uniprot_ids)]
        results = await tqdm.gather(
            *tasks, 
            desc="Fetching protein data",
            unit="protein")
    
    valid_results = [res for res in results if res is not None]
    
    print(f"{len(valid_results)} UniProtID were successfully processed")
    print(f"{len(error_ids['UniProtID'])} UniProtID not found")

    if return_dataframe:
        return protein_results_to_dataframe(valid_results, flatten_nested=flatten_nested), error_ids

    return valid_results, error_ids
