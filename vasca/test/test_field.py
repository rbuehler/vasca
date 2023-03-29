#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil

import astropy.units as uu
import pytest
from astropy.coordinates import SkyCoord
from astropy.time import Time
from loguru import logger

from vasca import vasca_pipe
from vasca.region import Region
from vasca.field import BaseField, GALEXField
from vasca.resource_manager import ResourceManager
import vasca.visualization as vvis


@pytest.fixture
def test_paths(tmp_path):
    """
    Returns dictionary with paths to static & cached test resources
    as well as a path to a common temporary directory.
    """
    paths = dict()

    # cached test resources
    with ResourceManager() as rm:
        paths["resource_root"] = rm.get_path("test_resources", "vasca")

    paths["galex_visits"] = f"{paths['resource_root']}/GALEX_visits_list.fits"
    paths["pipeline_cfg"] = f"{paths['resource_root']}/vasca_test_cfg.yaml"

    # temporary directory
    d = tmp_path
    paths["temp_path"] = f"/{d.resolve()}"

    return paths


@pytest.fixture
def galex_test_field_from_MAST_online(test_paths):
    gfield_id = 6381787756527353856  # AIS_309_1_28 2 visits (Crab pulsar)
    obs_filter = "NUV"

    data_path = test_paths["temp_path"]
    visits_data_path = test_paths["galex_visits"]

    gf = GALEXField.from_MAST(
        obs_id=gfield_id,
        obs_filter=obs_filter,
        data_path=data_path,
        visits_data_path=visits_data_path,
        write=False,
    )
    return gf


@pytest.fixture
def galex_test_field_from_MAST_offline(test_paths):
    gfield_id = 6388191295067652096  # NGC4993-GW170817 2 visits
    obs_filter = "NUV"

    test_resource_path = test_paths["resource_root"]
    data_path = f"{test_resource_path}/{gfield_id}"
    visits_data_path = f"{test_resource_path}/GALEX_visits_list.fits"
    gf = GALEXField.from_MAST(
        obs_id=gfield_id,
        obs_filter=obs_filter,
        data_path=data_path,
        visits_data_path=visits_data_path,
        write=False,
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
    assert gf.field_id in ["GNU6388191295067652096", "GNU6381787756527353856"]


@pytest.fixture
def new_field():

    field_data = {
        "field_id": [2605053246158405632],
        "name": ["PS_COSMOS_MOS23"],
        "ra": [151.00379402274802],
        "dec": [2.20171000810559],
        "observatory": ["GALEX"],
        "obs_filter": ["NUV"],
    }

    visits_data = {
        "vis_id": [2605053108543291392, 2605053108576845824],
        "time_bin_start": [54510.66210648148, 54520.248125],
        "time_bin_size": [1705.0, 837.05],
    }

    bf = BaseField()
    bf.add_table(field_data, "base_field:tt_field")
    bf.add_table(visits_data, "base_field:tt_visits")
    return bf


def test_set_field_attr_type(new_field):
    new_field.set_field_attr()
    expected = {
        "field_id": "str",
        "name": "name",
        "ra": 1 * uu.deg,
        "dec": 1 * uu.deg,
        "observatory": "obs",
        "obs_filter": "obs_filter",
        "center": SkyCoord(1, 1, unit="deg"),
        "nr_vis": 1,
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


def test_base_field_io(test_paths, new_field):
    # store fits data in temporary directory
    d = f"{test_paths['temp_path']}/fits_out_dir"
    os.mkdir(d)
    # get path to the fits file
    file_path = f"{d}/vasca_tables_output"  # note: no file extension

    new_field.write_to_fits(f"{file_path}.fits")
    new_field.write_to_hdf5(f"{file_path}.hdf5")
    save_str = new_field.__str__
    new_field.load_from_fits(f"{file_path}.fits")
    load_str_fits = new_field.__str__
    new_field.load_from_hdf5(f"{file_path}.hdf5")
    load_str_hdf5 = new_field.__str__
    assert save_str == load_str_fits == load_str_hdf5


def test_base_field_io_alt(test_paths, new_field):
    # store fits data in temporary directory
    d = f"{test_paths['temp_path']}/fits_out_dir"
    os.mkdir(d)
    # get path to the fits file
    file_path = f"{d}/vasca_tables_output.fits"
    # writ out vasca field file
    new_field.info()  # debugging
    new_field.write_to_fits(file_path)
    # test if a fits file exists
    assert os.path.isfile(file_path)


def test_pipeline(test_paths):
    # get pipeline config
    vasca_cfg = vasca_pipe.set_config(test_paths["pipeline_cfg"])

    # temporary output directory
    pipeline_out = f"{test_paths['temp_path']}/pipe_out"
    os.mkdir(pipeline_out)

    # edit paths
    vasca_cfg["general"]["out_dir_base"] = pipeline_out
    vasca_cfg["ressources"]["field_kwargs"]["data_path"] = test_paths["resource_root"]
    vasca_cfg["ressources"]["field_kwargs"]["visits_data_path"] = test_paths[
        "galex_visits"
    ]

    # run pipeline
    vasca_pipe.run(vasca_cfg)

    # Test visualizations
    # Load region
    region_name = "CrabGW"
    region_fname = pipeline_out + "/" + region_name + "/region_" + region_name + ".fits"
    rg = Region()
    rg.load_from_fits(region_fname)

    # Plot skypmap
    fd = rg.fields[rg.tt_fields[0]["field_id"]]
    fig, ax = vvis.plot_field_sky_map(fd)
    ax = vvis.plot_sky_sources(rg.tt_sources, tt_det=rg.tt_detections)

    # Plot light curve
    sel = rg.tt_sources["nr_det"] > 1
    fig_lc, ax_lc = vvis.plot_light_curve(
        rg, rg_src_ids=rg.tt_sources[sel]["rg_src_id"][0]
    )
    # delete field data from test resource directory
    # this forces to download the data from mast,
    # i.e., tests also the fallback from load method "MAST_LOCAL" to "MAST_REMOTE"
    # -> Todo: write dedicated test for the various load methods
    field_data_path = f"{test_paths['resource_root']}/6381787756527353856"
    if os.path.isdir(field_data_path):
        shutil.rmtree(field_data_path)


def main():
    # logging
    logger.enable("vasca")
    logger.level("DEBUG")

    bf = BaseField()
    test_base_field_io(bf)
    test_base_field_print_info(bf)
    test_base_field_io(bf)


if __name__ == "__main__":
    main()
