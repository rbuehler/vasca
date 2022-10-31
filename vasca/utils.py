#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for VASCA
"""
from itertools import zip_longest

import numpy as np
from astropy import units as uu
from astropy.coordinates import search_around_sky
from astropy.time import Time


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
