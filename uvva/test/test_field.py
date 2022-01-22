#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pprint import pprint

import astropy.units as uu
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from loguru import logger

from uvva.field import BaseField, Field, GALEXField


@pytest.fixture
def new_field():
    bf = BaseField()
    bf.tt_field = bf.get_table(
        np.asarray([42, "name", 2, 3, "obs", "filter"]), "base_field:tt_field"
    )

    return bf


def test_set_field_attr_type(new_field):
    new_field.set_field_attr()
    expected = {
        "id": int(1),
        "name": "name",
        "ra": 1 * uu.deg,
        "dec": 1 * uu.deg,
        "observatory": "obs",
        "obsfilter": "filter",
        "center": SkyCoord(1, 1, unit="deg"),
    }
    assert all(
        [isinstance(new_field.__dict__[key], type(expected[key])) for key in expected]
    )


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

    # Tests to execute
    test_base_field_print_info()
    test_base_field_io()


if __name__ == "__main__":
    main()
