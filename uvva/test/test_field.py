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


def main():
    # logging
    # replacing the standard handler
    # logger.remove(0)
    # setup custom handler to make use of
    # additional fields defined by uvva
    fmt = (
        "<g>{time:YYYY-MM-DD HH:mm:ss.SSS}</> | "
        "<lvl>{level: ^8}</> | "
        "<c>{name}</>:<c>{line}</> in <c>{extra[classname]}</>:<c>{function}</> | "
        "<lvl>{message}</>"
    )
    logger.add(sys.stdout, format=fmt)
    # activate logging (is deactivated by import of uvva)
    logger.enable("uvva")

    # test BaseField
    ff = BaseField()
    print("\n******** INFO *********")
    ff.info()

    print("\n******** PRINT ********")
    print(ff)

    # test GALEX use case
    gf = GALEXField(parobs_id=42, field_name="gTest", ra=np.float64(5))
    print(gf.tt_field.info)
    pprint(gf.tt_field.meta)

    # # tt_field is the old tt_coadd
    gf = Field(6388191295067652096)
    print("GALEX field info:")


if __name__ == "__main__":
    main()
