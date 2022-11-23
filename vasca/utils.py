#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for VASCA
"""

from itertools import zip_longest

import numpy as np
from astropy import units as uu
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.nddata import Cutout2D
from astropy.time import Time

# Global variable liking observator+obsfilter to a field ID addon
# The number of Id adon letter has to be three
# See get_region_field_id funtion below.
dd_obs_id_add = {"GALEXNUV": "GNU", "GALEXFUV": "GFU"}


def get_region_field_id(obs_field_id, observaory, obs_filter):
    """
    Return region field id, which also includes observatory and filter identifier.

    Parameters
    ----------
    obs_field_id : int
        Field ID for the given observatory.
    observaory : str
        Observatory name.
    obs_filter : TYPE
        Observation filter.

    Returns
    -------
    str
        Region field identifier.

    """
    return dd_obs_id_add[str(observaory) + str(obs_filter)] + str(obs_field_id)


# def get_field_file_name(field_id, observatory, obs_filter):
#     """
#     Helper function to create and load fields with uniform naming.

#     Parameters
#     ----------
#     field_id : int
#         Field ID.
#     observatory : str
#         Observatory of the field.
#     obs_filter : str
#         Observation filter of the field.

#     Returns
#     -------
#     str
#         Field default file name

#     """

#     return (
#         "field_"
#         + str(field_id)
#         + "_"
#         + str(observatory)
#         + "_"
#         + str(obs_filter)
#         + ".fits"
#     )


def extr_value(inputlist, upper=False):
    """
    Computes the extremum value in a list of numeric lists
    upper = False: minimum (default)
    upper = True: maximum
    """
    # if nested lists have different length
    # pads with np.nan to bring to uniform length
    padded_array = np.array(list(zip_longest(*inputlist, fillvalue=np.nan))).T
    return np.nanmax(padded_array) if upper else np.nanmin(padded_array)


def get_hist_bins(data, bin_size, is_list=False):
    """
    Generates list of bin edges according to data (min/max) and bin_size.
    """
    # get minimum and maximum rounded to integer
    if is_list:
        vmin = np.floor(extr_value(data))
        vmax = np.ceil(extr_value(data, upper=True))
    else:
        vmin = np.floor(np.min(data))
        vmax = np.ceil(np.max(data))
    # generate bin array
    bins = np.arange(vmin, vmax + bin_size, bin_size)
    return bins


def sky_sep2d(coords, seperation_limit=180):
    """
    computes 2D distances between all possible pairs of coordinates
    that are closer together than the separation limit
    (default: 180 -> all-sky)

    if the all-sky case is set, the number of unique distance values
    is determined as follows:
    for the number of coordinates n = len(coords) a number of
    pairwise distances values N = 1/2 * (n-1)*n is returned
    """
    n_coords = len(coords)  # number of coordinates
    seperation_limit = seperation_limit * uu.degree  # maximum separation

    # computes the 2D distences between a pair of coordinates A and B
    # (returned as list, together with the indeces of the pairs)
    # both permutations of coordinates A and B are returned
    # this yields redundent results for the distance C = dist(A,B) = dist(B,A)
    idx1, idx2, sep2d, _ = search_around_sky(coords, coords, seperation_limit)

    # removing redundent results:

    # reshapes distances to square matrix
    n_coords = len(coords)
    sep2d_square = np.zeros((n_coords, n_coords))
    sep2d_square[(idx1, idx2)] = sep2d

    # takes the lower triangular matrix
    # without the main diagonal elements
    sep2d_tril = np.tril(sep2d_square, k=-1)

    # returns the unique distance values (flattened triangular matrix)
    nonzero_idx = np.flatnonzero(sep2d_tril)
    sep2d_flat = sep2d_tril.ravel()[nonzero_idx]
    return sep2d_flat


def get_time_delta(tt_visit_dates, unit=None):
    # set astropy unit, defaults to seconds
    if unit is None:
        unit = uu.second

    # return zero if input has only one element
    if len(tt_visit_dates) < 2:
        return 0
    else:
        # convert to array of astropy Time objects
        times = Time(tt_visit_dates, format="mjd")
        # calculate time deltas
        dtimes_sec = np.array([dt.to_value(unit) for dt in np.diff(times)])
        return dtimes_sec


def get_time_delta_mean(tt_visit_dates, unit=None, deviation=True):
    # set astropy unit, defaults to seconds
    if unit is None:
        unit = uu.second

    # return zero if input has only one element
    if len(tt_visit_dates) < 2:
        return 0 if not deviation else (0, 0)
    else:
        # convert to array of astropy Time objects
        times = Time(tt_visit_dates, format="mjd")
        # calculate mean value and standard deviation of time deltas
        dtimes_mean = np.mean(np.diff(times)).to_value(unit)
        dtimes_std = np.std([dt.to_value(unit) for dt in np.diff(times)])

        return dtimes_mean if not deviation else (dtimes_mean, dtimes_std)


def table_to_array(table):
    """
    Converts an astropy Table object into a regular numpy ndarray.
    Caveat: All table columns must have the same type.
    """
    # Table to structured ndarray
    x = np.array(table)
    dtype = x.dtype[0]  # Type of first column

    # Checks consistency of data types
    assert all(
        [x.dtype[i] == dtype for i in range(len(x.dtype.names))]
    ), f"Expected same dtype '{dtype}' for all columns. {x.dtype}"

    # Creates view (not a copy) to return a regular numpy array
    return x.view((dtype, len(x.dtype.names)))


def get_cutout(field, position, size):
    """
    Wrapper function to create a :py:class:`~astropy.nddata.Cutout2D` from a
    :class:`~vasca.field.BaseField` input.

    Parameters
    ----------
    field : :class:`~vasca.field.BaseField`
        Field data object from which the reference image data and
        corresponding WCS are used to create the cutout.
    position : :py:class:`tuple`, :py:class:`~astropy.coordinates.SkyCoord`
        The position of the cutout array’s center with respect to the data array.
        The position can be specified either as a (x, y) tuple of pixel coordinates
        or a :py:class:`~astropy.coordinates.SkyCoord`,
        in which case wcs is a required input.
    size : :py:class:`astropy.units.Quantity`
        The size of the cutout array along each axis. For more inforamtion
        see documentaiont of :py:class:`~astropy.nddata.Cutout2D`.

    Returns
    -------
    :py:class:`~astropy.nddata.Cutout2D`
    """
    return Cutout2D(field.ref_img.data, position, size, field.ref_wcs)


def get_cutout_bounds(cutout, out_frame="icrs"):
    """
    Returns :py:class:`~astropy.coordinates.SkyCoord` object
    defining a rectangular cutout with upperight and lower left pixels
    in coordinate system of ``out_frame``.

    Parameters
    ----------
    cutout : :py:class:`~astropy.nddata.Cutout2D`
        Cutout object from which to derive the bounds.
    out_frame : :py:class:`str`, :py:class:`~astropy.coordinates.BaseCoordinateFrame` \
    class or instance, or `:py:class:`~astropy.coordinates.SkyCoord` instance, optional
        Coordinate frame of returned :py:class:`~astropy.coordinates.SkyCoord` object.

    Returns
    -------
    :py:class:`~astropy.coordinates.SkyCoord`
    """
    # Corner pixels

    # lower left
    x_ll = 0
    y_ll = 0
    # upper right
    x_ur = cutout.data.shape[1]
    y_ur = cutout.data.shape[0]

    # Convert to world coordinates
    ll = cutout.wcs.pixel_to_world(x_ll, y_ll)
    ur = cutout.wcs.pixel_to_world(x_ur, y_ur)

    # Separation
    # print(ll.separation(ur).arcsec)

    return SkyCoord([ll.ra, ur.ra], [ll.dec, ur.dec], frame=ll.frame).transform_to(
        out_frame
    )


def get_cutout_mask(tt_mcat, cutout_bounds, frame="icrs"):
    """
    Returns a bool array to mask a coordinate dataset for given cutout bounds.

    Parameters
    ----------
    tt_mcat : :py:class:`~astropy.table.Table`
        Table to be masked
    cutout_bounds : :py:class:`~astropy.coordinates.SkyCoord`

    Returns
    -------
    :py:class:`~numpy.ndarray`
    """
    # Handle column naming according to specified frame
    x, y = None, None
    if frame == "icrs" or frame == "fk5":
        x = "ra"
        y = "dec"
        mask_cutout_ra = (tt_mcat[x] <= cutout_bounds[0].ra) * (
            tt_mcat[x] >= cutout_bounds[1].ra
        )
        mask_cutout_dec = (tt_mcat[y] >= cutout_bounds[0].dec) * (
            tt_mcat[y] <= cutout_bounds[1].dec
        )
    elif frame == "galctic":
        x = "lon"
        y = "lat"
        mask_cutout_ra = (tt_mcat[x] <= cutout_bounds[0].l) * (
            tt_mcat[x] >= cutout_bounds[1].l
        )
        mask_cutout_dec = (tt_mcat[y] >= cutout_bounds[0].b) * (
            tt_mcat[y] <= cutout_bounds[1].b
        )

    return mask_cutout_ra * mask_cutout_dec
