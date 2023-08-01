"""
VASCA x-matching functions.
"""
from time import sleep

import numpy as np
from astropy import units as uu
from astropy.coordinates import SkyCoord
from astropy.table import Table
from loguru import logger
from tqdm import tqdm

import vasca.xmatch.acms_api as acmsapi


def xmatch(api, catalog, radius, coords, ids=None, dropna=True, hide_progress=True):
    """
    Cross-match function using different APIs.

    Parameters
    ----------
    api : str
        The name of the API to use for cross-matching.
    catalog : str
        The name of the catalog to cross-match with.
    radius : float
        The search radius in arcseconds.
    coords : astropy.table.Table
        Table containing the coordinates of the sources to cross-match.
        It must have "ra" and "dec" columns specifying the right ascension
        and declination in degrees, respectively.
    ids : list or None, optional
        List of identifiers for the sources. If provided, it must have the same length
        as `coords`. Default is None.
    dropna : bool, optional
        Whether to drop rows where no match was found within the search radius.
        Default is True.
    hide_progress : bool, optional
        Whether to hide the progress bar. Default is True.

    Returns
    -------
    astropy.table.Table
        The table containing the cross-matched results.

    Raises
    ------
    ValueError
        If the provided API name is not supported.

    Notes
    -----
    This function provides a unified interface for performing cross-matching using
    different APIs. Currently, only the 'ampel' API is supported. More APIs may be
    implemented in the future.
    """

    known_apis = [
        "ampel",
        # "astroquery_mast",
        # "astroquery_something_else",
    ]  # more APIs to be implemented in the future

    # Checks if api is supported
    if api not in known_apis:
        raise ValueError(f"Expected parameter 'api' in {known_apis}, got '{api}'")

    # Selects x-matching method
    if api == "ampel":
        xmatch_method = xmatch_ampel

    # Run x-matching
    logger.debug(f"Running x-matching via {api} API")
    tt_matched = xmatch_method(
        catalog=catalog,
        search_radius_arcsec=radius,
        coords=coords,
        ids=ids,
        hide_progress=hide_progress,
    )

    # Removes rows where no match was found within radius
    if dropna:
        tt_matched = tt_matched[~np.isnan(tt_matched["vasca_dist"])]

    return tt_matched


def xmatch_ampel(
    catalog,
    search_radius_arcsec,
    coords,
    ids=None,
    hide_progress=True,
):
    """
    This function performs a nearest neighbor search using the Ampel Catalog Matching
    Service. It queries the Ampel API to find the nearest matching objects in the
    specified catalog within the given search radius around the input coordinates.
    The cross-matched results are returned as an Astropy Table.

    The Ampel Catalog Matching Service provides a collection of catalogs that can be
    used for cross-matching. You can refer to the documentation for the list of
    available catalogs.
    (https://ampelproject.github.io/astronomy/ztf/index#catalog-matching)

    Parameters
    ----------
    catalog : str
        The name of the catalog to perform the cross-match with.
    search_radius_arcsec : float
        The search radius in arcseconds.
    coords : astropy.table.Table
        Table containing the coordinates of the sources to cross-match.
        It must have "ra" and "dec" columns specifying the right ascension
        and declination in degrees, respectively.
    ids : list or None, optional
        List of identifiers for the sources. If provided, it must have the same length
        as `coords`. Default is None.
    hide_progress : bool, optional
        Whether to hide the progress bar. Default is True.

    Returns
    -------
    astropy.table.Table
        The table containing the cross-matched results. The columns include "vasca_id",
        "vasca_ra", "vasca_dec", "vasca_dist" (distance to matched object), and
        additional columns specific to the catalog.

    Raises
    ------
    ValueError
        If the provided catalog name is unknown.
    """
    # Get list of available catalogs from Ampel Catalog Matching Service
    acms_cat_info = acmsapi.get_acms_cat_info()

    # Checks if catalog name exits
    if catalog not in acms_cat_info:
        raise ValueError(
            f"Unknown catalog name '{catalog}'. "
            f"Available (ACMS) catalogs:\n{sorted(list(acms_cat_info.keys()))}"
        )

    # Ampel-hosted catalogs
    cat_info = acms_cat_info[catalog]

    # Gets catalog columns
    cat_cols = [item["name"] for item in cat_info["columns"]]

    # Creates the columns for the output table
    # Always starts with VASCA ID,  VASCA coordinates and distance to matched object.
    # Catalog-specific are appended columns after
    cols = ["vasca_id", "vasca_ra", "vasca_dec", "vasca_dist", *cat_cols]

    # Get x-matching distances
    data = list()
    n_srcs = len(coords)
    # Loops over input sources
    for i, coord in tqdm(
        enumerate(coords),
        total=n_srcs,
        desc=f"X-match {catalog}",
        disable=hide_progress,
    ):
        # Coordinates to SkyyCoord object
        coord = SkyCoord(ra=coord["ra"], dec=coord["dec"], unit=uu.deg)

        # Store data row-wise
        row_data = dict.fromkeys(cols, np.nan)

        # Fill default data
        row_data.update(
            {
                "vasca_id": ids[i] if id is not None else i,
                "vasca_ra": coord.ra,
                "vasca_dec": coord.dec,
            }
        )

        # x-match coordinate
        res = acmsapi.acms_xmatch_query(
            catalog=cat_info["name"],
            catalog_type=cat_info["use"],
            ra_deg=coord.ra.value,
            dec_deg=coord.dec.value,
            search_radius_arcsec=search_radius_arcsec,
            search_type="nearest",
            cache_dir=None,
            crash_hard=True,
        )

        # Collects data if response is valid
        # If not, the distance and catalog columns will only hold NaN values
        if isinstance(res, dict):
            row_data.update({"vasca_dist": res["dist_arcsec"], **res["body"]})

        data.append(row_data)

        # Mitigates server-side rate limit
        sleep(0.06)

    # Creates astropy table
    tt_matched = Table(rows=data)

    return tt_matched
