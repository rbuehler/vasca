#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil

import pytest

from vasca import vasca_pipe
from vasca.region import Region
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


def test_pipeline_vis(test_paths):
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
    fd = rg.get_field(field_id=rg.tt_fields[0]["field_id"])
    fig, ax = vvis.plot_field_sky_map(fd)
    vvis.plot_sky_sources(rg.tt_sources, tt_det=rg.tt_detections)

    # Plot light curve
    sel = rg.tt_sources["nr_det"] > 1
    fig_lc, ax_lc = vvis.plot_light_curve(
        rg, rg_src_ids=rg.tt_sources[sel]["rg_src_id"][0]
    )

    # Plot all fields
    rg.add_coverage_hp(nside=4096, coord_sys="icrs")
    vvis.plot_region_sky_mollview(rg, var="nr_vis")
    vvis.plot_region_sky_gnomeview(rg, rg.tt_fields[0]["ra"], rg.tt_fields[0]["dec"])

    # Plot pipeline diagnostics
    vvis.plot_pipe_diagnostic(rg, "tt_sources", "scatter")
    vvis.plot_pipe_diagnostic(rg, "tt_sources", "hist")
    vvis.plot_pipe_diagnostic(fd, "tt_detections", "hist", fig_size=(8, 10))

    # delete field data from test resource directory
    # this forces to download the data from mast,
    # i.e., tests also the fallback from load method "MAST_LOCAL" to "MAST_REMOTE"
    # -> Todo: write dedicated test for the various load methods
    field_data_path = f"{test_paths['resource_root']}/6381787756527353856"
    if os.path.isdir(field_data_path):
        shutil.rmtree(field_data_path)
