#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import sys
from pprint import pprint

import astropy.units as uu
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.time import Time
from loguru import logger

from uvva.field import BaseField, Field, GALEXField
from uvva.uvva_table import UVVATable


@pytest.fixture
def new_field():
    field_data = [
        2605053246158405632,
        "PS_COSMOS_MOS23",
        151.00379402274802,
        2.20171000810559,
        "GALEX",
        "NUV",
    ]
    visits_data = [
        [2605053108543291392, 54510.66210648148, 54510.68184027778, 1705.0],
        [2605053108576845824, 54520.248125, 54520.2578125, 837.05],
    ]
    bf = BaseField()
    bf.tt_field = UVVATable.from_template(np.asarray(field_data), "base_field:tt_field")
    bf.tt_visits = UVVATable.from_template(
        np.asarray(visits_data), "base_field:tt_visits"
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
        "n_visits": 1,
        "t_exp_sum": 42.0 * uu.s,
        "t_start": Time(54520.248125, format="mjd"),
        "t_stop": Time(54520.248125, format="mjd"),
    }
    assert all(
        [isinstance(new_field.__dict__[key], type(expected[key])) for key in expected]
    )


def test_base_field_print_info(new_field):
    new_field = BaseField()
    new_field.info()
    print(new_field)


def test_base_field_io(new_field):
    new_field.write_to_fits()
    olfstr = new_field.__str__
    new_field.load_from_fits()
    assert olfstr == new_field.__str__


def test_base_field_io_alt(tmp_path, new_field):
    # create temporary directory to store fits data in
    d = tmp_path / "fits_out_dir"
    d.mkdir()
    # get path pointing to the fits file as string
    file_path = (d / "uvva_tables_output.fits").resolve()
    new_field.info()
    new_field.write_to_fits(file_path)
    # test if a fits file exists
    assert not glob.glob(f"{d.resolve}/*.fits")


def main():
    # logging
    logger.enable("uvva")

    bf = BaseField()
    test_base_field_io(bf)
    test_base_field_print_info(bf)
    test_base_field_io(bf)


if __name__ == "__main__":
    main()
