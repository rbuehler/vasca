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
        cfg = yaml.safe_load(file)  # yaml.

    # Set output directory
    if cfg["general"]["out_dir_base"] == "CWD":
        cfg["general"]["out_dir_base"] = os.getcwd()

    # Store cfg file name in cfg dictionary
    cfg["cfg_file"] = cfg_file

    return cfg


def set_logger(cfg):
    """
    Setup logger

    Parameters
    ----------
    cfg : dict
        UVVA pipeline configuration dictionary

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
    set_logger(cfg)

    # Load region fields
    rg = Region.load_from_config(cfg)

    for field_id, gf in rg.fields.items():
        logger.info("Analysing field:" + str(field_id))

        # Create directory structure for fields
        field_dir = (
            cfg["general"]["out_dir_base"]
            + "/"
            + cfg["general"]["name"]
            + "/"
            + str(field_id)
            + "/"
        )

        # Create folder
        if not os.path.exists(field_dir):
            os.makedirs(field_dir)

        # Apply selections
        gf.select_rows("tt_detections", cfg["selection"]["detections"])

        # Run clustering
        gf.cluster_meanshift(
            cfg["cluster"]["add_upper_limits"], **cfg["cluster"]["meanshift"]
        )

        # Plot results
        fig_sky = gf.plot_sky(plot_detections=True)

        fig_lc = plt.figure(figsize=(10, 4))
        gf.plot_light_curve(range(0, 10), ylim=[25.5, 13.5])
        plt.tight_layout()

        # Write field out
        gf.write_to_fits(field_dir + "field_" + str(field_id) + ".fits")
        if cfg["general"]["hd_img_out"]:
            fig_sky.savefig(
                field_dir + "sky_map_hr_" + str(field_id) + ".pdf", dpi=3000
            )
        else:
            fig_sky.savefig(field_dir + "sky_map_hr_" + str(field_id) + ".pdf", dpi=150)

        fig_lc.savefig(field_dir + str(field_id) + "_lc.pdf", dpi=150)

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
    cfg = set_config(args.cfg)

    run(cfg)
