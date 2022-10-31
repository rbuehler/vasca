#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that runs the VASCA pipeline.
"""

import argparse
import os
import sys
from itertools import zip_longest

# from multiprocessing import Pool
from multiprocessing import Pool

import healpy as hpy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import yaml
from loguru import logger

from vasca.region import Region
import vasca.visualization as vvis


def set_config(cfg_file):
    """
    Setup pipeline configuration file from yaml

    Parameters
    ----------
    cfg_file : str
        yaml configuration file name

    Returns
    -------
    vasca_cfg : dict
        VASCA pipeline configuration dictionary
    """
    with open(cfg_file) as file:
        vasca_cfg = yaml.safe_load(file)  # yaml.

    # Set output directory
    if vasca_cfg["general"]["out_dir_base"] == "CWD":
        vasca_cfg["general"]["out_dir_base"] = os.getcwd()

    # Store vasca_cfg file name in vasca_cfg dictionary
    vasca_cfg["cfg_file"] = cfg_file

    return vasca_cfg


def set_logger(vasca_cfg):
    """
    Setup logger. Gets configuration from global config variable set with set_config

    Returns
    -------
    None.

    """
    log_dir = (
        vasca_cfg["general"]["out_dir_base"] + "/" + vasca_cfg["general"]["name"] + "/"
    )
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    # create log file name
    log_file_name = "log_" + vasca_cfg["general"]["name"] + ".txt"
    if vasca_cfg["general"]["log_file"] != "default":
        log_file_name = vasca_cfg["general"]["log_file"]

    log_cfg = {
        "handlers": [
            {
                "sink": sys.stdout,
                "format": "<green>{time:YYYY-MM-DD HH:mm:ss.SSSS}</green> "
                "<cyan>{name}</cyan>:<cyan>{line}</cyan> |"
                "<level>{level}:</level> {message}",
                "level": vasca_cfg["general"]["log_level"],
                "colorize": True,
                "backtrace": True,
                "diagnose": True,
            },
            {
                "sink": log_dir + log_file_name,
                "serialize": True,
                "backtrace": True,
                "diagnose": True,
            },
        ],
    }
    logger.configure(**log_cfg)
    logger.enable("vasca")

    logger.info("Runing '" + __file__ + "'")
    logger.debug("Config. file: '" + vasca_cfg["cfg_file"] + "'")
    logger.debug("Output log. file: '" + log_cfg["handlers"][1]["sink"] + "'")


def run_field(field, vasca_cfg):
    """
    Run analysis on a single field

    Parameters
    ----------
    field : vasca.field
        Field to run analysis on.

    Returns
    -------
    field : vasca.field
        Modified field with results

    """
    logger.info("Analysing field:" + str(field.field_id))

    # Create directory structure for fields
    field_dir = (
        vasca_cfg["general"]["out_dir_base"]
        + "/"
        + vasca_cfg["general"]["name"]
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
    field.select_rows(vasca_cfg["selection"]["det_quality"])

    # Run clustering
    field.cluster_meanshift(
        vasca_cfg["cluster"]["add_upper_limits"], **vasca_cfg["cluster"]["meanshift"]
    )

    # Source selection
    field.select_rows(vasca_cfg["selection"]["src_quality"])
    field.select_rows(vasca_cfg["selection"]["src_variability"])

    # Write out field
    field.write_to_fits(field_dir + "field_" + str(field.field_id) + ".fits")

    # Plot sky maps
    fig_sky = vvis.plot_field_sky(field, plot_detections=True)
    if vasca_cfg["general"]["hd_img_out"]:
        fig_sky.savefig(
            field_dir + "sky_map_hr_" + str(field.field_id) + ".png",
            dpi=2500,
        )
    else:
        fig_sky.savefig(
            field_dir + "sky_map_lr_" + str(field.field_id) + ".png",
            dpi=150,
        )
    plt.close(fig_sky)

    # Make field diagnostocs
    diags = [
        ("detections", "hist"),
        ("sources", "hist"),
        ("detections", "scatter"),
        ("sources", "scatter"),
    ]
    for diag in diags:
        fig_diag = vvis.plot_pipe_diagnostic(field, "tt_" + diag[0], diag[1])
        fig_diag.savefig(
            field_dir + diag[1] + "_" + diag[0] + "_" + str(field.field_id) + ".png",
            dpi=150,
        )
        plt.close(fig_diag)

    # Draw selected  light curves
    sel_srcs = field.tt_sources["sel"]
    if sel_srcs.sum() > 0:
        all_srcs_ids = field.tt_sources[sel_srcs]["fd_src_id"].data

        # Loop over list in chuks of 14
        for fd_src_ids_chunk in zip_longest(
            *([np.nditer(all_srcs_ids)]) * 14, fillvalue=-1
        ):
            fig_lc = plt.figure(figsize=(10, 10))
            fd_src_ids_chunk = np.array(fd_src_ids_chunk, dtype=np.int64).flatten()
            fd_src_ids_chunk = np.delete(
                fd_src_ids_chunk, np.where(fd_src_ids_chunk == -1)
            )
            field.plot_light_curve(fd_src_ids_chunk, ylim=[24.5, 15.5])
            plt.tight_layout()
            srcs_name = "_".join([str(elem) for elem in fd_src_ids_chunk])

            fig_lc.savefig(
                field_dir + "lcs/lc_" + str(field.field_id) + "_" + srcs_name + ".png",
                dpi=150,
            )
            plt.close(fig_lc)

    # Remove some items which are not further needed to free memory
    # del field.__dict__["tt_detections"]
    # field._table_names.remove("tt_detections")
    # field.ref_img = None
    # field.ref_wcs = None

    return field


def run(vasca_cfg):
    """
    Runs the VASCA pipeline

    Parameters
    ----------
    vasca_cfg : dict
        VASCA pipeline configuration dictionary

    Returns
    -------
    None.

    """

    # Setup logger
    set_logger(vasca_cfg)

    # Set matplotlib backend
    if vasca_cfg["general"]["mpl_backend"] != "SYSTEM":
        logger.info(
            f'setting matplotlib backend to {vasca_cfg["general"]["mpl_backend"]}'
        )
        matplotlib.use(vasca_cfg["general"]["mpl_backend"])

    # Load region fields
    rg = Region.load_from_config(vasca_cfg)
    rg.add_table_from_fields("tt_visits")

    # Setup output directors
    region_dir = (
        vasca_cfg["general"]["out_dir_base"] + "/" + vasca_cfg["general"]["name"] + "/"
    )

    # Write our healpix coverage maps
    hp_vis, hp_exp = rg.add_coverage_hp(nside=4096)

    if vasca_cfg["general"]["hp_coverage_out"]:
        hpy.fitsfunc.write_map(
            region_dir
            + "/region_"
            + vasca_cfg["general"]["name"]
            + "_coverage_hp.fits",
            [hp_vis, hp_exp],
            coord="C",
            column_names=["nr_vis", "exposure"],
            dtype=[np.float32, np.float32],
            overwrite=True,
            partial=True,
        )
    # Run each field in a separate process in parallel
    with Pool(vasca_cfg["general"]["nr_cpus"]) as pool:
        pool_return = pool.starmap(
            run_field,
            [(field, vasca_cfg) for field in rg.fields.values()],
        )
    pool.join()

    # update region fields
    for field in pool_return:
        rg.fields[field.field_id] = field

    rg.add_table_from_fields("tt_ref_sources")
    rg.add_table_from_fields("tt_sources")
    rg.add_table_from_fields("ta_sources_lc")

    if vasca_cfg["general"]["save_dets"] == "selected":
        rg.add_table_from_fields("tt_detections", only_selected=True)
    elif vasca_cfg["general"]["save_dets"] == "all":
        rg.add_table_from_fields("tt_detections", only_selected=False)
    else:
        logger.info("Not saving detection table into region")

    # Write out regions
    rg.write_to_fits(
        file_name=region_dir + "/region_" + vasca_cfg["general"]["name"] + ".fits"
    )

    # Make diagnostics
    fig_diag_rgsrcs_hist = vvis.plot_pipe_diagnostic(rg, "tt_sources", "hist")
    fig_diag_rgsrcs_hist.savefig(
        region_dir
        + "/region_diagnostic_srcs_hist_"
        + vasca_cfg["general"]["name"]
        + ".png",
        dpi=150,
    )
    plt.close(fig_diag_rgsrcs_hist)

    fig_diag_rgsrcs_scat = vvis.plot_pipe_diagnostic(rg, "tt_sources", "scatter")
    fig_diag_rgsrcs_scat.savefig(
        region_dir
        + "/region_diagnostic_srcs_scat_"
        + vasca_cfg["general"]["name"]
        + ".png",
        dpi=150,
    )
    plt.close(fig_diag_rgsrcs_scat)

    # Write used config file
    yaml_out_name = (
        region_dir + "/vasca_ran_cfg_" + vasca_cfg["general"]["name"] + ".yaml"
    )
    with open(yaml_out_name, "w") as yaml_file:
        yaml.dump(vasca_cfg, yaml_file)


def run_from_file():
    """
    Run VASCA pipeline from a configuration file passed with the ArgumentParser

    Returns
    -------
    None.

    """

    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--cfg", type=str, default=os.getcwd() + "/vasca_cfg.yaml")
    args = parser.parse_args()

    print(
        50 * "-"
        + f"\n Running vasca_pipe.py with configuration file:\n {args.cfg}\n"
        + 50 * "-"
    )
    vasca_cfg = set_config(args.cfg)
    run(vasca_cfg)


if __name__ == "__main__":
    run_from_file()
