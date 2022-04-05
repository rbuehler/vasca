#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob

import astropy.units as uu
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.time import Time
from loguru import logger

from uvva.field import BaseField, GALEXField
from uvva.resource_manager import ResourceManager


@pytest.fixture
def galex_test_field_from_archive():
    field_id = 6388191295067652096
    filter = "NUV"
    with ResourceManager() as rm:
        test_resource_path = rm.get_path("test_resources", "uvva")
        data_path = f"{test_resource_path}/{field_id}"
        visits_data_path = f"{test_resource_path}/GALEX_visits_list.fits"
    gf = GALEXField.from_archive(
        obs_id=field_id,
        filter=filter,
        data_path=data_path,
        visits_data_path=visits_data_path,
    )
    return gf


def test_galex_field_from_archive(galex_test_field_from_archive):
    gf = galex_test_field_from_archive
    assert gf.field_id == 6388191295067652096


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
        [2605053108543291392, 54510.66210648148, 1705.0],
        [2605053108576845824, 54520.248125, 837.05],
    ]
    bf = BaseField()
    bf.add_table(np.asarray(field_data), "base_field:tt_field")
    bf.add_table(np.asarray(visits_data), "base_field:tt_visits")
    return bf


def test_set_field_attr_type(new_field):
    new_field.set_field_attr()
    expected = {
        "field_id": int(1),
        "name": "name",
        "ra": 1 * uu.deg,
        "dec": 1 * uu.deg,
        "observatory": "obs",
        "obsfilter": "filter",
        "center": SkyCoord(1, 1, unit="deg"),
        "n_visits": 1,
        "time_bin_size_sum": 42.0 * uu.s,
        "time_start": Time(54520.248125, format="mjd"),
        "time_stop": Time(54520.248125, format="mjd"),
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
    new_field.write_to_hdf5()
    save_str = new_field.__str__
    new_field.load_from_fits()
    load_str_fits = new_field.__str__
    new_field.load_from_hdf5()
    load_str_hdf5 = new_field.__str__
    assert save_str == load_str_fits


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
    logger.level("DEBUG")

    bf = BaseField()
    test_base_field_io(bf)
    test_base_field_print_info(bf)
    test_base_field_io(bf)


if __name__ == "__main__":
    main()
