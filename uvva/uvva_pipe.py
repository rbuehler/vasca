#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that runs the UVVA pipeline.
"""


import argparse
import os
import sys
from itertools import zip_longest

import numpy as np
import healpy as hpy
import matplotlib.pyplot as plt
import yaml
from loguru import logger

from uvva.region import Region

# from multiprocessing import Pool
from multiprocessing import Pool

from collections import OrderedDict


def set_config(cfg_file):
    """
    Setup pipeline configuration file from yaml

    Parameters
    ----------
    cfg_file : str
        yaml configuration file name

    Returns
    -------
    uvva_cfg : dict
        UVVA pipeline configuration dictionary
    """
    with open(cfg_file) as file:
        global uvva_cfg
        uvva_cfg = yaml.safe_load(file)  # yaml.

    # Set output directory
    if uvva_cfg["general"]["out_dir_base"] == "CWD":
        uvva_cfg["general"]["out_dir_base"] = os.getcwd()

    # Store uvva_cfg file name in uvva_cfg dictionary
    uvva_cfg["cfg_file"] = cfg_file

    return uvva_cfg


# TODO: there are some GALEX specific variables, will need to be adapted to other missions
def diagnostic(tc, table_name, plot_type):

    det_vars = OrderedDict()
    if plot_type == "hist":
        # Detections diagnostic
        if table_name == "tt_detections":
            det_vars["s2n"] = {"logx": True}
            det_vars["mag"] = {"range": [14.5, 26.5]}
            det_vars["mag_err"] = {"range": [0.0, 3.0]}
            det_vars["r_fov"] = {"range": [0.0, 0.7]}
            det_vars["point_src_prob"] = {}
            det_vars["artifacts"] = {"histtype": "step"}
            fig, axs = plt.subplots(3, 2, figsize=(18, 12), squeeze=False)
        elif table_name == "tt_sources":
            det_vars["nr_det"] = {}
            det_vars["nr_uls"] = {}
            det_vars["mag_mean"] = {}
            det_vars["mag_rchiq"] = {"logx": True, "range": [-3, 3]}
            det_vars["mag_dmax"] = {}
            det_vars["mag_dmax_sig"] = {"logx": True, "range": [-3, 2]}
            det_vars["mag_var"] = {"logx": True, "range": [-3, 0]}
            det_vars["ul_weight"] = {"bins": 100}
            fig, axs = plt.subplots(2, 4, figsize=(22, 12), squeeze=False)
        else:
            logger.warning("Diegnostic for table '{table_name}' not defined")

    elif plot_type == "scatter":
        # Detections diagnostic
        if table_name == "tt_detections":
            det_vars[("s2n", "mag")] = {
                "invert_yaxis": True,
                "xlim": [1, 100],
                "ylim": [14.5, 26.5],
            }
            det_vars[("r_fov", "mag")] = {"invert_yaxis": True, "ylim": [15.5, 26.5]}
            det_vars[("point_src_prob", "mag")] = {
                "invert_yaxis": True,
                "ylim": [14.5, 26.5],
            }
            det_vars[("artifacts", "mag")] = {
                "invert_yaxis": True,
                "ylim": [14.5, 26.5],
            }
            det_vars[("r_fov", "artifacts")] = {}
            det_vars[("mag_err", "mag")] = {
                "xlim": [0.01, 3],
                "invert_yaxis": True,
                "ylim": [14.5, 26.5],
            }
            fig, axs = plt.subplots(3, 2, figsize=(14, 12), squeeze=False)
        elif table_name == "tt_sources":
            det_vars[("mag_rchiq", "mag_dmax_sig")] = {"xscale": "log"}
            det_vars[("mag_rchiq", "ul_weight")] = {"xscale": "log"}
            det_vars[("mag_dmax_sig", "ul_weight")] = {}
            det_vars[("mag_dmax_sig", "mag_mean")] = {
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            det_vars[("ul_weight", "mag_mean")] = {
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            det_vars[("mag_rchiq", "mag_mean")] = {
                "xscale": "log",
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            det_vars[("mag_var", "mag_mean")] = {
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            det_vars[("nr_uls", "mag_mean")] = {
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            fig, axs = plt.subplots(2, 4, figsize=(22, 12), squeeze=False)
        else:
            logger.warning("Diegnostic for table '{table_name}' not defined")
    else:
        logger.warning("Plot type '{plot_type}' unknown")

    axs = axs.flatten()
    ax_ctr = 0
    for var, plot_arg in det_vars.items():
        if plot_type == "hist":
            tc.plot_hist(table_name, var, axs[ax_ctr], **plot_arg)
        elif plot_type == "scatter":
            tc.plot_scatter(table_name, var[0], var[1], axs[ax_ctr], **plot_arg)
        ax_ctr += 1

    plt.tight_layout()
    plt.legend()
    return fig


def set_logger():
    """
    Setup logger. Gets configuration from global config variable set with set_config

    Returns
    -------
    None.

    """
    log_dir = (
        uvva_cfg["general"]["out_dir_base"] + "/" + uvva_cfg["general"]["name"] + "/"
    )
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    log_cfg = {
        "handlers": [
            {
                "sink": sys.stdout,
                "format": "<green>{time:YYYY-MM-DD HH:mm:ss.SSSS}</green> "
                "<cyan>{name}</cyan>:<cyan>{line}</cyan> |"
                "<level>{level}:</level> {message}",
                "level": uvva_cfg["general"]["log_level"],
                "colorize": True,
                "backtrace": True,
                "diagnose": True,
            },
            {
                "sink": log_dir + uvva_cfg["general"]["log_file"],
                "serialize": True,
                "backtrace": True,
                "diagnose": True,
            },
        ],
    }
    logger.configure(**log_cfg)
    logger.enable("uvva")

    logger.info("Runing '" + __file__ + "'")
    logger.debug("Config. file: '" + uvva_cfg["cfg_file"] + "'")
    logger.debug("Output log. file: '" + log_cfg["handlers"][1]["sink"] + "'")


def run_field(field):
    """
    Run analysis on a single field

    Parameters
    ----------
    field : uvva.field
        Field to run analysis on.

    Returns
    -------
    field : uvva.field
        Modified field with results

    """

    logger.info("Analysing field:" + str(field.field_id))

    # Create directory structure for fields
    field_dir = (
        uvva_cfg["general"]["out_dir_base"]
        + "/"
        + uvva_cfg["general"]["name"]
        + "/fields/"
        + str(field.field_id)
        + "/"
    )

    # Create folder
    if not os.path.exists(field_dir):
        os.makedirs(field_dir)
    if not os.path.exists(field_dir + "lcs/"):
        os.makedirs(field_dir + "lcs/")

    # Apply selections
    field.select_rows(uvva_cfg["selection"]["det_quality"])

    # Run clustering
    field.cluster_meanshift(
        uvva_cfg["cluster"]["add_upper_limits"], **uvva_cfg["cluster"]["meanshift"]
    )

    # Source selection
    field.select_rows(uvva_cfg["selection"]["src_quality"])
    field.select_rows(uvva_cfg["selection"]["src_variability"])

    # Plot results
    fig_sky = field.plot_sky(plot_detections=True)
    field.write_to_fits(field_dir + "field_" + str(field.field_id) + ".fits")
    if uvva_cfg["general"]["hd_img_out"]:
        fig_sky.savefig(
            field_dir + "sky_map_hr_" + str(field.field_id) + ".pdf", dpi=3000
        )
    else:
        fig_sky.savefig(
            field_dir + "sky_map_hr_" + str(field.field_id) + ".pdf", dpi=150
        )

    # Make field diagnostocs
    diags = [
        ("detections", "hist"),
        ("sources", "hist"),
        ("detections", "scatter"),
        ("sources", "scatter"),
    ]
    for diag in diags:
        fig_diag = diagnostic(field, "tt_" + diag[0], diag[1])
        fig_diag.savefig(
            field_dir + diag[1] + "_" + diag[0] + "_" + str(field.field_id) + ".png",
            dpi=150,
        )
        plt.close(fig_diag)

    # Draw selected  lightcurve
    sel_srcs = field.tt_sources["sel"]
    if sel_srcs.sum() > 0:
        all_srcs_ids = field.tt_sources[sel_srcs]["src_id"].data

        # Loop over list in chuks of 14
        for src_ids_chunk in zip_longest(
            *([np.nditer(all_srcs_ids)]) * 14, fillvalue=-1
        ):
            fig_lc = plt.figure(figsize=(10, 10))
            src_ids_chunk = np.array(src_ids_chunk, dtype=np.int64).flatten()
            src_ids_chunk = np.delete(src_ids_chunk, np.where(src_ids_chunk == -1))
            field.plot_light_curve(src_ids_chunk, ylim=[24.5, 15.5])
            plt.tight_layout()
            srcs_name = "_".join([str(elem) for elem in src_ids_chunk])

            fig_lc.savefig(
                field_dir + "lcs/lc_" + str(field.field_id) + "_" + srcs_name + ".png",
                dpi=150,
            )
            plt.close(fig_lc)

    return field


def run(uvva_cfg):
    """
    Runs the UVVA pipeline

    Parameters
    ----------
    uvva_cfg : dict
        UVVA pipeline configuration dictionary

    Returns
    -------
    None.

    """

    # Setup logger
    set_logger()

    # Load region fields
    rg = Region.load_from_config(uvva_cfg)
    rg.add_table_from_fields("tt_visits")

    # Setup output directors
    region_dir = (
        uvva_cfg["general"]["out_dir_base"] + "/" + uvva_cfg["general"]["name"] + "/"
    )

    # Write our helpix coverage maps
    hp_vis, hp_exp = rg.add_coverage_hp(nside=4096)

    hpy.fitsfunc.write_map(
        region_dir + "/region_" + uvva_cfg["general"]["name"] + "_coverage_hp.fits",
        [hp_vis, hp_exp],
        coord="C",
        column_names=["nr_vis", "exposure"],
        overwrite=True,
        partial=True,
    )

    # Run each field in a separate process in parallel
    with Pool(uvva_cfg["general"]["nr_cpus"]) as pool:
        pool_return = pool.map(run_field, list(rg.fields.values()))

    # update region fields
    for field in pool_return:
        rg.fields[field.field_id] = field

    rg.add_table_from_fields("tt_ref_sources")
    rg.add_table_from_fields("tt_sources")
    rg.add_table_from_fields("ta_sources_lc")
    # rg.add_table_from_fields("tt_detections", only_selected=True)

    # Write out regions
    rg.write_to_fits(
        file_name=region_dir + "/region_" + uvva_cfg["general"]["name"] + ".fits"
    )

    # Make diagnostics
    fig_diag_rgsrcs_hist = diagnostic(rg, "tt_sources", "hist")
    fig_diag_rgsrcs_scat = diagnostic(rg, "tt_sources", "scatter")
    fig_diag_rgsrcs_hist.savefig(
        region_dir
        + "/region_"
        + uvva_cfg["general"]["name"]
        + "_diagnostic_srcs_hist.png",
        dpi=150,
    )
    fig_diag_rgsrcs_scat.savefig(
        region_dir
        + "/region_"
        + uvva_cfg["general"]["name"]
        + "_diagnostic_srcs_scat.png",
        dpi=150,
    )

    # Write used config file
    yaml_out_name = region_dir + "/uvva_ran_cfg.yaml"
    with open(yaml_out_name, "w") as yaml_file:
        yaml.dump(uvva_cfg, yaml_file)


def run_from_file(file_name="./uvva_cfg.yaml"):
    set_config(file_name)
    run(uvva_cfg)


if __name__ == "__main__":

    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--cfg", type=str, default=os.getcwd() + "/uvva_cfg.yaml")
    args = parser.parse_args()
    run_from_file(args.cfg)
