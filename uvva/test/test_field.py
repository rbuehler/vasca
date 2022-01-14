#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 07:32:11 2022

@author: buehler
"""

import sys
from pprint import pprint

import astropy.units as uu
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from loguru import logger

from uvva.field import BaseField, Field, GALEXField


@pytest.mark.parametrize(
    "input, expected",
    [
        (
            {
                "field_id": None,
                "field_name": None,
                "ra": None,
                "dec": None,
                "observatory": None,
                "obsfilter": None,
            },
            {
                "field_id": None,
                "field_name": None,
                "ra": None,
                "dec": None,
                "observatory": None,
                "obsfilter": None,
            },
        ),
        (
            {
                "field_id": 42,
                "field_name": "foo",
                "ra": 23,
                "dec": 23,
                "observatory": None,
                "obsfilter": None,
            },
            {
                "field_id": 1 * uu.dimensionless_unscaled,
                "field_name": "foo",
                "ra": 1 * uu.deg,
                "dec": 1 * uu.deg,
                "observatory": None,
                "obsfilter": None,
            },
        ),
    ],
)
def test_base_field_input_type(input, expected):
    par_keys = ["field_id", "field_name", "ra", "dec", "observatory", "obsfilter"]
    bf = BaseField(**input)
    assert all([isinstance(bf.__dict__[key], type(expected[key])) for key in par_keys])


@pytest.mark.parametrize(
    "input, expected",
    [
        ({"ra": 23, "dec": 23}, SkyCoord(1, 1, unit="deg")),
        ({"ra": None, "dec": 23 * uu.deg}, None),
        ({"ra": None, "dec": None}, None),
    ],
)
def test_base_field_center_type(input, expected):
    bf = BaseField(**input)
    assert isinstance(bf.center, type(expected))


def test_base_field_print_info():
    ff = BaseField()
    ff.info()
    print(ff)


def test_base_field_io():
    ff = BaseField()
    ff.write_to_fits()
    ff_olfstr = ff.__str__
    ff.load_from_fits()
    assert ff_olfstr == ff.__str__


def main():
    # logging
    logger.enable("uvva")

    # Tests to excecute
    # test_base_field_print_info()
    test_base_field_io()


if __name__ == "__main__":
    main()
