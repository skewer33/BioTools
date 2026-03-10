import asyncio

import aiohttp

MYGENE_API_URL = "https://mygene.info/v3/query"


async def _async_query(genename, taxid, session, per_request_retries=2, per_request_retry_delay=0.5):
    url = f"{MYGENE_API_URL}?q={genename}&species_facet_filter={taxid}&fields=uniprot&species={taxid}"

    for attempt in range(per_request_retries + 1):
        try:
            async with session.get(url) as response:
                if response.status != 200:
                    if attempt < per_request_retries:
                        await asyncio.sleep(per_request_retry_delay)
                        continue
                    return [genename, None]

                json_res = await response.json()
                uniprot_id = _find_UID(json_res)
                return [genename, uniprot_id]
        except (aiohttp.ClientError, asyncio.TimeoutError):
            if attempt < per_request_retries:
                await asyncio.sleep(per_request_retry_delay)
                continue
            return [genename, None]

    return [genename, None]


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
        uid = json_res["hits"][0]["uniprot"].get("Swiss-Prot")
        if uid is None:
            uid = json_res["hits"][0]["uniprot"].get("TrEMBL")
    except Exception:
        uid = None
    if isinstance(uid, list):
        uid = uid[0]
    return uid


async def _async_request(
    genenames: list,
    taxid,
    max_concurrent: int = 20,
    request_timeout: float = 20.0,
    per_request_retries: int = 2,
    per_request_retry_delay: float = 0.5,
):
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
    timeout = aiohttp.ClientTimeout(total=request_timeout)
    connector = aiohttp.TCPConnector(limit=max_concurrent)
    semaphore = asyncio.Semaphore(max_concurrent)

    async with aiohttp.ClientSession(timeout=timeout, connector=connector) as session:

        async def _bounded_query(gene):
            async with semaphore:
                return await _async_query(
                    gene,
                    taxid,
                    session,
                    per_request_retries=per_request_retries,
                    per_request_retry_delay=per_request_retry_delay,
                )

        tasks = [_bounded_query(genename) for genename in genenames]
        results = await asyncio.gather(*tasks)

    uniprot_ids = {res[0]: res[1] for res in results if res[1] is not None}
    error_list = [res[0] for res in results if res[1] is None]

    return uniprot_ids, error_list


async def gene2uniprotid(
    genenames: list,
    taxid: int = 9606,
    max_cycle: int = 50,
    retry_delay: float = 5.0,
    max_concurrent: int = 20,
    request_timeout: float = 20.0,
    per_request_retries: int = 2,
    per_request_retry_delay: float = 0.5,
):
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
        Taxonomy ID (default 9606, human).
    max_cycle : int, optional
        Maximum number of cycles (including the first), default 50.
    retry_delay : float, optional
        Pause between cycles in seconds, default 5.0.
    max_concurrent : int, optional
        Maximum number of parallel HTTP requests, default 20.
    request_timeout : float, optional
        Request timeout in seconds, default 20.0.
    per_request_retries : int, optional
        Retries for each gene request, default 2.
    per_request_retry_delay : float, optional
        Delay between per-request retries in seconds, default 0.5.

    Returns
    -------
    tuple[dict, list]
        - UniProtId_dict: {gene: uniprot_id}
        - error_genes: list of genes for which ID was not found.
    """
    if max_cycle < 1:
        raise ValueError("max_cycle must be at least 1")
    if max_concurrent < 1:
        raise ValueError("max_concurrent must be at least 1")
    if request_timeout <= 0:
        raise ValueError("request_timeout must be > 0")
    if per_request_retries < 0:
        raise ValueError("per_request_retries must be >= 0")

    cycle = 1
    uniprot_id_dict, error_genes = await _async_request(
        genenames,
        taxid=taxid,
        max_concurrent=max_concurrent,
        request_timeout=request_timeout,
        per_request_retries=per_request_retries,
        per_request_retry_delay=per_request_retry_delay,
    )

    found_previous_cycle = 0
    found_current_cycle = len(uniprot_id_dict)

    while cycle < max_cycle and error_genes:
        if found_previous_cycle == 0:
            ratio = float("inf") if found_current_cycle > 0 else 1.0
        else:
            ratio = found_current_cycle / found_previous_cycle

        if ratio == 1.0:
            break

        await asyncio.sleep(retry_delay)
        cycle += 1

        found_previous_cycle = found_current_cycle
        new_dict, error_genes = await _async_request(
            error_genes,
            taxid=taxid,
            max_concurrent=max_concurrent,
            request_timeout=request_timeout,
            per_request_retries=per_request_retries,
            per_request_retry_delay=per_request_retry_delay,
        )
        uniprot_id_dict.update(new_dict)
        found_current_cycle = len(uniprot_id_dict)

    print(f"{len(uniprot_id_dict)} genes successfully converted to UniProtIDs")
    print(f"{len(error_genes)} genes not converted")

    return uniprot_id_dict, error_genes
