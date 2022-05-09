#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import os

import astropy.units as uu
import numpy as np
import pytest
import yaml
from astropy.coordinates import SkyCoord
from astropy.time import Time
from loguru import logger

from uvva import uvva_pipe
from uvva.field import BaseField, GALEXField
from uvva.resource_manager import ResourceManager


@pytest.fixture
def test_paths(tmp_path):
    """
    Returns dictionary with paths to static & cached test resources
    as well as a path to a common temporary directory.
    """
    paths = dict()

    # cached test resources
    with ResourceManager() as rm:
        paths["resource_root"] = rm.get_path("test_resources", "uvva")

    paths["galex_visits"] = f"{paths['resource_root']}/GALEX_visits_list.fits"
    paths["pipeline_cfg"] = f"{paths['resource_root']}/uvva_test_cfg.yaml"

    # temporary directory
    d = tmp_path
    paths["temp_path"] = f"/{d.resolve()}"

    return paths


@pytest.fixture
def galex_test_field_from_MAST_online(test_paths):
    field_id = 6381787756527353856  # AIS_309_1_28 2 visits (Crab pulsar)
    obs_filter = "NUV"

    data_path = test_paths["temp_path"]
    visits_data_path = test_paths["galex_visits"]

    gf = GALEXField.from_MAST(
        obs_id=field_id,
        obs_filter=obs_filter,
        data_path=data_path,
        visits_data_path=visits_data_path,
    )
    return gf


@pytest.fixture
def galex_test_field_from_MAST_offline(test_paths):
    field_id = 6388191295067652096  # NGC4993-GW170817 2 visits
    obs_filter = "NUV"

    test_resource_path = test_paths["resource_root"]
    data_path = f"{test_resource_path}/{field_id}"
    visits_data_path = f"{test_resource_path}/GALEX_visits_list.fits"
    gf = GALEXField.from_MAST(
        obs_id=field_id,
        obs_filter=obs_filter,
        data_path=data_path,
        visits_data_path=visits_data_path,
    )
    return gf


@pytest.fixture(
    params=[
        "galex_test_field_from_MAST_offline",
        "galex_test_field_from_MAST_online",
    ],
    ids=["NGC4993-GW170817 (cached)", "AIS_309_1_28 (downloaded)"],
)
def galex_test_field_from_MAST(request):
    return request.getfixturevalue(request.param)


def test_galex_field_from_MAST(galex_test_field_from_MAST):
    gf = galex_test_field_from_MAST
    assert gf.field_id in [6388191295067652096, 6381787756527353856]


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
        "obs_filter": "obs_filter",
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
    assert save_str == load_str_fits == load_str_hdf5


def test_base_field_io_alt(test_paths, new_field):
    # store fits data in temporary directory
    d = f"{test_paths['temp_path']}/fits_out_dir"
    os.mkdir(d)
    # get path to the fits file
    file_path = f"{d}/uvva_tables_output.fits"
    # writ out uvva field file
    new_field.info()  # debugging
    new_field.write_to_fits(file_path)
    # test if a fits file exists
    assert os.path.isfile(file_path)


def test_pipeline(test_paths):
    # get pipeline config
    with open(test_paths["pipeline_cfg"]) as f:
        cfg = yaml.safe_load(f)

    # temporary output directory
    pipeline_out = f"{test_paths['temp_path']}/pipe_out"
    os.mkdir(pipeline_out)

    # edit paths
    cfg["general"]["out_dir"] = pipeline_out
    cfg["ressources"]["load_kwargs"]["data_path"] = test_paths["resource_root"]
    cfg["ressources"]["load_kwargs"]["visits_data_path"] = test_paths["galex_visits"]

    # Write modified config file
    cfg_mod = f"{pipeline_out}/uvva_test_cfg_modified.yaml"
    with open(cfg_mod, "w") as f:
        yaml.dump(cfg, f)

    # run pipeline
    cfg = uvva_pipe.set_config(cfg_mod)
    uvva_pipe.run(cfg)


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
