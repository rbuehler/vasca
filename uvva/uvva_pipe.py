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

from uvva.field import GALEXField
from uvva.region import Region

# from multiprocessing import Pool
from multiprocessing import Pool


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

    # Write out regions
    region_dir = cfg["general"]["out_dir_base"] + "/" + cfg["general"]["name"] + "/"
    rg.write_to_fits(
        file_name=region_dir + "/region_" + cfg["general"]["name"] + ".fits"
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
