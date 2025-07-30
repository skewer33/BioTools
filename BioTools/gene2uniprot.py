import aiohttp
import asyncio
import time


async def _async_query(genename, taxid):
    #url = f"http://mygene.info/v3/query?q=symbol:{genename}&species_facet_filter={taxid}&fields=uniprot&species={taxid}"
    url = f"http://mygene.info/v3/query?q={genename}&species_facet_filter={taxid}&fields=uniprot&species={taxid}"
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            json_res = await response.json()
            uniprot_id = _find_UID(json_res)
            if uniprot_id is None:
                return [genename, None]
            else:
                return [genename, uniprot_id]


def _find_UID(json_res):
    """
    This function takes the result of mygene.info/v3/query in json format and returns the Uniprot ID (Swiss-Prot or TrEMBL) of the first hit, if available.

    Parameters
    ----------
    json_res : json
        The result of mygene.info/v3/query in json format

    Returns
    -------
    str
        The Uniprot ID (Swiss-Prot or TrEMBL) of the first hit, if available
    """
    try:
        UID = json_res['hits'][0]['uniprot'].get('Swiss-Prot')
        if UID is None:
            UID = json_res['hits'][0]['uniprot'].get('TrEMBL')
    except:
        UID = None
    if isinstance(UID, list):
        UID = UID[0]
    return UID

async def _async_request(genenames:list, taxid):
    
    """
    Asynchronously retrieves UniProt IDs for a list of gene names.

    This function queries the mygene.info API for each gene name provided, retrieves the corresponding UniProt ID,
    and returns a dictionary mapping each gene name to its UniProt ID, along with a list of gene names that encountered errors.

    Parameters
    ----------
    genenames : list
        A list of gene names to query.

    Returns
    -------
    tuple
        A tuple containing:
        - results_dict: dict
            A dictionary where each key is a gene name and each value is a dictionary containing the UniProt ID and an error flag.
        - error_list: list
            A list of gene names for which the query did not successfully retrieve a UniProt ID.
    """
    tasks = [_async_query(genename, taxid) for genename in genenames]
    results = await asyncio.gather(*tasks)
    UniprotIDs = {res[0]: res[1] for res in results if res[1] is not None}
    error_list = [res[0] for res in results if res[1] is None]
    
    return UniprotIDs, error_list

async def gene2uniprotid(genenames:list, taxid = 9606, cycles = 2):
    """
    Retrieves UniProt IDs for a list of gene names.

    This function uses asynchronous requests to query the mygene.info API for each gene name provided,
    retrieves the corresponding UniProt ID, and returns a dictionary mapping each gene name to its UniProt ID.
    It also returns a list of gene names for which the queries did not successfully retrieve a UniProt ID.

    Parameters
    ----------
    genenames : list
        A list of gene names to query.
    taxid : int, optional
        The taxonomic ID of the species. Default is 9606 (Homo sapiens).
    cycles : int, optional
        The number of times to retry the queryv if not all UniProt IDs were successfully retrieved. Default is 2.

    Returns
    -------
    tuple
        A tuple containing:
        - UniProtId_dict: dict
            A dictionary where each key is a gene name and the value is its corresponding UniProt ID.
        - error_genes: list
            A list of gene names for which the query did not successfully retrieve a UniProt ID.
    """
    if cycles < 1:
        raise ValueError("Number of cycles must be at least 1")
    UniProtId_dict, error_genes = await _async_request(genenames, taxid=taxid)
    for i in range(cycles-1):
        time.sleep(5)
        if len(error_genes) > 0:
            UniProtId_dict_2, error_genes = await _async_request(error_genes, taxid=taxid)
            UniProtId_dict = {**UniProtId_dict, **UniProtId_dict_2}
        else:
            break
    
    print(f'{len(UniProtId_dict)} genes succesfully converted to UniprotIDs')
    print(f'{len(error_genes)} genes not converted')
    
    return UniProtId_dict, error_genes
