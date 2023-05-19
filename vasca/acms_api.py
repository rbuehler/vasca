"""
API tools for the Ampel Catalog Matching Service
"""

import json

import backoff
import requests
from loguru import logger

API_BASEURL = "https://ampel.zeuthen.desy.de"
API_CATALOGMATCH_URL = API_BASEURL + "/api/catalogmatch"


@backoff.on_exception(
    backoff.expo,
    (requests.exceptions.RequestException, json.JSONDecodeError),
    max_time=60,
)
def get_acms_cat_info(queryurl=None):
    """
    Retrieve catalog information from the Ampel Catalog Matching Service API.

    Parameters
    ----------
    queryurl : str, optional
        The URL to retrieve catalog information from the API. Defaults to the
        catalog endpoint URL, which is
        https://ampel.zeuthen.desy.de/api/catalogmatch/catalogs/.

    Returns
    -------
    dict or None
        A dictionary with catalog information, where the keys are the catalog names
        and the values are the corresponding catalog information. If the API response
        does not contain a list of catalogs, None is returned.

    Raises
    ------
    ValueError
        If the API response contains an inconsistent list of catalog names (i.e., some
        names are duplicated).

    Notes
    -----
    This function uses an exponential backoff strategy for retrying the API request
    in case of `requests.exceptions.RequestException` or `json.JSONDecodeError`.
    The maximum time for retrying is set to 60 seconds.
    """

    if queryurl is None:
        queryurl = f"{API_CATALOGMATCH_URL}/catalogs/"

    response = requests.get(queryurl)

    product = response.json()

    if isinstance(product, list):
        if len(product) > 0:
            # Catalog name list
            cat_names = [item["name"] for item in product]

            # Check if names are unique
            n_names = len(cat_names)
            if len(set(cat_names)) != n_names:
                raise ValueError(f"Inconsistent catalog name list: \n{cat_names}")

            # Dictionary with catalog info
            cat_info = {name: info for name, info in zip(cat_names, product)}
    else:
        return None

    logger.debug("Fetching available catalogs from Ampel Catalog Matching Service.")
    return cat_info


@backoff.on_exception(backoff.expo, requests.exceptions.RequestException, max_time=600)
def acms_xmatch_query(
    catalog: str,
    catalog_type: str,
    ra_deg: float,
    dec_deg: float,
    search_radius_arcsec: float = 10,
    search_type: str = "all",
    cache_dir: str = None,
    crash_hard: bool = False,
):
    """
    Method for querying catalogs via the Ampel API.
    'catalog' must be the name of a supported catalog, e.g., SDSS_spec, PS1,
    NEDz_extcats... For a full list of catalogs, refer to
    https://ampel.zeuthen.desy.de/api/catalogmatch/docs.

    Parameters
    ----------
    catalog : str
        The name of the catalog to query.
    catalog_type : str
        The type of catalog to use. Supported types are "extcats" and "catsHTM".
    ra_deg : float
        The right ascension in degrees.
    dec_deg : float
        The declination in degrees.
    search_radius_arcsec : float, optional
        The search radius in arcseconds. Default is 10 arcseconds.
    search_type : str, optional
        The search type. Supported types are "all" (default) and "nearest".
    cache_dir : str, optional
        The directory to cache the query results. Default is None.
    crash_hard : bool, optional
        Whether to raise a `requests.exceptions.RequestException` when encountering
        an error. If False, None is returned instead. Default is False.

    Returns
    -------
    dict or None
        The query result as a dictionary. None is returned if an error occurs and
        ``crash_hard`` is False.

    Raises
    ------
    ValueError
        If ``catalog_type`` or ``search_type`` are not among the supported types
    requests.exceptions.RequestException
        If an error occurs and ``crash_hard`` is True.

    Notes
    -----
    This method uses the exponential backoff strategy to retry the API request in case
    of `requests.exceptions.RequestException`. The maximum time for retrying is set to
    600 seconds.

    """

    known_catalog_types = ["extcats", "catsHTM"]
    known_search_types = ["all", "nearest"]

    if catalog_type not in known_catalog_types:
        raise ValueError(
            f"Expected parameter 'catalog_type' in {known_catalog_types}, "
            f"got '{catalog_type}'."
        )
    if search_type not in known_search_types:
        raise ValueError(
            f"Expected parameter 'search_type' in {known_search_types}, "
            f"got '{search_type}'."
        )

    queryurl_catalogmatch = API_CATALOGMATCH_URL + "/cone_search/" + search_type

    # Creates a json body to post
    headers = {"accept": "application/json", "Content-Type": "application/json"}
    # Creates query
    query = {
        "ra_deg": ra_deg,
        "dec_deg": dec_deg,
        "catalogs": [
            {"name": catalog, "rs_arcsec": search_radius_arcsec, "use": catalog_type}
        ],
    }
    # Posts query
    response = requests.post(url=queryurl_catalogmatch, json=query, headers=headers)

    # Catch HTML status codes
    html_err = 503  # 503 Service Unavailable
    if response.status_code == html_err:
        logger.warning(f"RequestException: HTML Status {html_err}. ")
        if response.headers:
            logger.warning(f"Response header: {response.headers}. ")
        if not crash_hard:
            logger.warning("Returning None")
            return None
        else:
            raise requests.exceptions.RequestException

    # Fetch json payload
    try:
        res = response.json()[0]
    except Exception as e:
        logger.warning(f"Could not extract json. Exception {e}. ")
        if response.headers:
            logger.warning(f"Response header: {response.headers}. ")
        logger.warning("Returning None")
        if not crash_hard:
            logger.warning("Returning None")
            return None
        else:
            raise requests.exceptions.RequestException

    return res
