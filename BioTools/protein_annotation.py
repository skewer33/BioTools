import asyncio

import aiohttp
import pandas as pd
from tqdm.asyncio import tqdm

UNIPROT_API_URL = "https://www.ebi.ac.uk/proteins/api/proteins/"
PDB_API_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/"


def _init_error_ids():
    return {
        "UniProtID": [],
        "PDB": [],
        "Gene": [],
        "GO_terms": [],
        "TaxID": [],
        "Annotation": [],
        "Sequence": [],
        "ParseError": [],
        "UnhandledError": [],
    }


async def _get_uniprot_data(uniprot_id, session, per_request_retries=2, per_request_retry_delay=0.5):
    url = f"{UNIPROT_API_URL}{uniprot_id}"
    for attempt in range(per_request_retries + 1):
        try:
            async with session.get(url, headers={"Accept": "application/json"}) as response:
                if response.status == 200:
                    return await response.json()
                if attempt < per_request_retries:
                    await asyncio.sleep(per_request_retry_delay)
                    continue
                return None
        except (aiohttp.ClientError, asyncio.TimeoutError):
            if attempt < per_request_retries:
                await asyncio.sleep(per_request_retry_delay)
                continue
            return None
    return None


def _parse_uniprot_data(data):
    result = {
        "Gene": "N/A",
        "TaxID": "N/A",
        "Annotation": "N/A",
        "GO_terms": {},
        "Sequence": "N/A",
    }

    result["Gene"] = data.get("gene", [{}])[0].get("name", {}).get("value", "N/A")
    result["TaxID"] = data.get("organism", {}).get("taxonomy", "N/A")

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
    go_description = []
    for db_reference in data.get("dbReferences", []):
        if db_reference.get("type") == "GO":
            go_terms.append(db_reference.get("id"))
            go_description.append(db_reference.get("properties", {}).get("term"))
    result["GO_terms"] = dict(zip(go_terms, go_description))

    sequence_info = data.get("sequence", {})
    if sequence_info:
        result["Sequence"] = sequence_info.get("sequence", "N/A")
    return result


async def _get_pdb_structures(uniprot_id, session, error_ids, per_request_retries=2, per_request_retry_delay=0.5):
    url = f"{PDB_API_URL}{uniprot_id}"
    for attempt in range(per_request_retries + 1):
        try:
            async with session.get(url) as response:
                if response.status == 200:
                    return await response.json()
                if attempt < per_request_retries:
                    await asyncio.sleep(per_request_retry_delay)
                    continue
                error_ids["PDB"].append(uniprot_id)
                return None
        except (aiohttp.ClientError, asyncio.TimeoutError):
            if attempt < per_request_retries:
                await asyncio.sleep(per_request_retry_delay)
                continue
            error_ids["PDB"].append(uniprot_id)
            return None
    error_ids["PDB"].append(uniprot_id)
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


async def _get_protein_info(
    uniprot_id,
    session,
    error_ids,
    per_request_retries=2,
    per_request_retry_delay=0.5,
):
    uniprot_data = await _get_uniprot_data(
        uniprot_id,
        session,
        per_request_retries=per_request_retries,
        per_request_retry_delay=per_request_retry_delay,
    )
    if not uniprot_data:
        error_ids["UniProtID"].append(uniprot_id)
        return None

    try:
        result = _parse_uniprot_data(uniprot_data)
    except Exception:
        error_ids["ParseError"].append(uniprot_id)
        return None

    for k, v in result.items():
        if v == "N/A" or v == [] or v == {}:
            error_ids.setdefault(k, []).append(uniprot_id)

    pdb_data = await _get_pdb_structures(
        uniprot_id,
        session,
        error_ids,
        per_request_retries=per_request_retries,
        per_request_retry_delay=per_request_retry_delay,
    )
    result["PDB"] = _parse_pdb_data(pdb_data, uniprot_id) if pdb_data else []
    result["UniProtID"] = uniprot_id

    return result


async def get_proteins_info(
    uniprot_ids,
    max_concurrent=10,
    return_dataframe=False,
    flatten_nested=True,
    request_timeout=20.0,
    per_request_retries=2,
    per_request_retry_delay=0.5,
):
    """
    Asynchronously retrieves protein information for a list of UniProt IDs.

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
    request_timeout : float
        Request timeout in seconds. Default is 20.0.
    per_request_retries : int
        Number of retries for each request. Default is 2.
    per_request_retry_delay : float
        Delay between per-request retries in seconds. Default is 0.5.

    Returns
    -------
    tuple
        - results: list or pandas.DataFrame
            Protein information per UniProt ID.
        - error_ids: dict
            A dictionary where each key is a category of error and each value is a list of UniProt IDs that encountered that error.
    """
    if max_concurrent < 1:
        raise ValueError("max_concurrent must be at least 1")
    if request_timeout <= 0:
        raise ValueError("request_timeout must be > 0")
    if per_request_retries < 0:
        raise ValueError("per_request_retries must be >= 0")

    error_ids = _init_error_ids()

    timeout = aiohttp.ClientTimeout(total=request_timeout)
    connector = aiohttp.TCPConnector(limit=max_concurrent)
    semaphore = asyncio.Semaphore(max_concurrent)

    async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:

        async def _bounded(uid):
            async with semaphore:
                return await _get_protein_info(
                    uid,
                    session,
                    error_ids,
                    per_request_retries=per_request_retries,
                    per_request_retry_delay=per_request_retry_delay,
                )

        tasks = [_bounded(uid) for uid in tqdm(uniprot_ids)]
        results = await tqdm.gather(
            *tasks,
            desc="Fetching protein data",
            unit="protein",
            return_exceptions=True,
        )

    valid_results = []
    for uid, res in zip(uniprot_ids, results):
        if isinstance(res, Exception):
            error_ids["UnhandledError"].append(uid)
            continue
        if res is not None:
            valid_results.append(res)

    print(f"{len(valid_results)} UniProtID were successfully processed")
    print(f"{len(error_ids['UniProtID'])} UniProtID not found")

    if return_dataframe:
        return protein_results_to_dataframe(valid_results, flatten_nested=flatten_nested), error_ids

    return valid_results, error_ids
