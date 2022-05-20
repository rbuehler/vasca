#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that runs the UVVA pipeline.
"""


import argparse
import os
import sys

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
    cfg : dict
        UVVA pipeline configuration dictionary
    """
    with open(cfg_file) as file:
        global cfg
        cfg = yaml.safe_load(file)  # yaml.

    # Set output directory
    if cfg["general"]["out_dir_base"] == "CWD":
        cfg["general"]["out_dir_base"] = os.getcwd()

    # Store cfg file name in cfg dictionary
    cfg["cfg_file"] = cfg_file

    return cfg


# TODO: there are some GALEX specific variables, will need to be adapted to other missions
def diagnostic(tc, table_name):

    # Detections diagnostic
    det_vars = OrderedDict()
    if table_name == "tt_detections":
        det_vars["s2n"] = {"logx": True}
        det_vars["mag"] = {"range": [13.0, 25.0]}
        det_vars["mag_err"] = {"range": [0.0, 3.0]}
        det_vars["r_fov"] = {"range": [0.0, 0.7]}
        det_vars["point_src_prob"] = {}
        det_vars["artifacts"] = {"histtype": "step"}
        fig, axs = plt.subplots(3, 2, figsize=(18, 12), squeeze=False)
    elif table_name == "tt_sources":
        det_vars["nr_det"] = {}
        det_vars["nr_det_meas"] = {}
        det_vars["mag_mean"] = {}
        det_vars["mag_rchiq"] = {"logx": True, "range": [-3, 3]}
        det_vars["mag_dmax"] = {}
        det_vars["mag_dmax_sig"] = {"logx": True, "range": [-3, 2]}
        det_vars["nr_ul_mean"] = {}
        det_vars["mag_var"] = {"logx": True, "range": [-3, 0]}
        fig, axs = plt.subplots(2, 4, figsize=(22, 12), squeeze=False)
    else:
        logger.warning("Table '{table_name}' does not exist")

    axs = axs.flatten()
    ax_ctr = 0
    for var, hist_arg in det_vars.items():
        tc.plot_hist(table_name, var, axs[ax_ctr], **hist_arg)
        ax_ctr += 1

    plt.tight_layout()
    plt.legend()
    return fig


def set_logger():
    """
    Setup logger. Gets configuration from global config cariable cfg set with set_config

    Returns
    -------
    None.

    """
    log_dir = cfg["general"]["out_dir_base"] + "/" + cfg["general"]["name"] + "/"
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    log_cfg = {
        "handlers": [
            {
                "sink": sys.stdout,
                "format": "<green>{time:YYYY-MM-DD HH:mm:ss.SSSS}</green> "
                "<cyan>{name}</cyan>:<cyan>{line}</cyan> |"
                "<level>{level}:</level> {message}",
                "level": cfg["general"]["log_level"],
                "colorize": True,
                "backtrace": True,
                "diagnose": True,
            },
            {
                "sink": log_dir + cfg["general"]["log_file"],
                "serialize": True,
                "backtrace": True,
                "diagnose": True,
            },
        ],
    }
    logger.configure(**log_cfg)
    logger.enable("uvva")

    logger.info("Runing '" + __file__ + "'")
    logger.debug("Config. file: '" + cfg["cfg_file"] + "'")
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
        cfg["general"]["out_dir_base"]
        + "/"
        + cfg["general"]["name"]
        + "/"
        + str(field.field_id)
        + "/"
    )

    # Create folder
    if not os.path.exists(field_dir):
        os.makedirs(field_dir)

    # Apply selections
    field.select_rows(cfg["selection"]["det_quality"])

    # Run clustering
    field.cluster_meanshift(
        cfg["cluster"]["add_upper_limits"], **cfg["cluster"]["meanshift"]
    )

    # Source selection
    field.select_rows(cfg["selection"]["src_quality"])
    field.select_rows(cfg["selection"]["src_variability"])

    # Plot results
    fig_sky = field.plot_sky(plot_detections=True)

    fig_diag_sel = diagnostic(field, "tt_detections")
    fig_diag_srcs = diagnostic(field, "tt_sources")

    fig_lc = plt.figure(figsize=(10, 4))
    field.plot_light_curve(range(0, 10), ylim=[25.5, 13.5])
    plt.tight_layout()

    # Write field out
    field.write_to_fits(field_dir + "field_" + str(field.field_id) + ".fits")
    if cfg["general"]["hd_img_out"]:
        fig_sky.savefig(
            field_dir + "sky_map_hr_" + str(field.field_id) + ".pdf", dpi=3000
        )
    else:
        fig_sky.savefig(
            field_dir + "sky_map_hr_" + str(field.field_id) + ".pdf", dpi=150
        )
    fig_diag_sel.savefig(
        field_dir + str(field.field_id) + "_diagnostic_det.pdf", dpi=150
    )
    fig_diag_srcs.savefig(
        field_dir + str(field.field_id) + "_diagnostic_srcs.pdf", dpi=150
    )
    fig_lc.savefig(field_dir + str(field.field_id) + "_lc.pdf", dpi=150)

    return field


def run(cfg):
    """
    Runs the UVVA pipeline

    Parameters
    ----------
    cfg : dict
        UVVA pipeline configuration dictionary

    Returns
    -------
    None.

    """

    # Setup logger
    set_logger()

    # Load region fields
    rg = Region.load_from_config(cfg)

    # Run each field in a separate process in parallel
    with Pool(cfg["general"]["nr_cpus"]) as pool:
        pool_return = pool.map(run_field, list(rg.fields.values()))

    # update region fields
    for field in pool_return:
        rg.fields[field.field_id] = field

    rg.add_table_from_fields("tt_ref_sources")
    rg.add_table_from_fields("tt_sources")
    rg.add_table_from_fields("tt_detections", only_selected=True)

    fig_diag_rgsrcs = diagnostic(rg, "tt_sources")

    # Write out regions
    region_dir = cfg["general"]["out_dir_base"] + "/" + cfg["general"]["name"] + "/"
    rg.write_to_fits(
        file_name=region_dir + "/region_" + cfg["general"]["name"] + ".fits"
    )
    fig_diag_rgsrcs.savefig(
        region_dir + "/region_" + cfg["general"]["name"] + "_diagnostic_srcs.pdf",
        dpi=150,
    )

    # Write used config file
    yaml_out_name = region_dir + "/uvva_ran_cfg.yaml"
    with open(yaml_out_name, "w") as yaml_file:
        yaml.dump(cfg, yaml_file)


if __name__ == "__main__":

    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--cfg", type=str, default=os.getcwd() + "/uvva_cfg.yaml")
    args = parser.parse_args()

    # Get UVVA configuration file
    set_config(args.cfg)
    run(cfg)
