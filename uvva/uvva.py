#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that runs the UVVA pipeline.
"""


import sys
import os
import argparse
import yaml
from loguru import logger
from uvva.region import Region
from uvva.field import GALEXField
import matplotlib.pyplot as plt


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
        cfg = yaml.safe_load(file)     # yaml.

    # Set output directory
    if cfg["general"]["out_dir"] == "CWD":
        cfg["general"]["out_dir"] = os.getcwd()

    # Store cfg file name in cfg dictionary
    cfg["cfg_file"] = cfg_file

    print(cfg)

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
    log_cfg = {
        "handlers": [
            {"sink": sys.stdout,
             "format": "<green>{time:YYYY-MM-DD HH:mm:ss.SSSS}</green> "
             "<cyan>{name}</cyan>:<cyan>{line}</cyan> |"
             "<level>{level}:</level> {message}",
             "level": cfg["general"]["log_level"], "colorize":True,
             "backtrace":True, "diagnose":True},
            {"sink": cfg["general"]["out_dir"]+"/"+cfg["general"]["log_file"],
             "serialize": True, "backtrace":True, "diagnose":True},
        ],
    }
    logger.configure(**log_cfg)
    logger.enable("uvva")

    logger.info("Runing '"+__file__+"'")
    logger.debug("Config. file: '"+cfg["cfg_file"]+"'")
    logger.debug("Output log. file: '" + log_cfg["handlers"][1]["sink"]+"'")


if __name__ == '__main__':

    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--cfg', type=str, default=os.getcwd()+'/uvva_cfg.yaml')
    args = parser.parse_args()

    # Get UVVA configuration file
    cfg = set_config(args.cfg)

    # Setup logger
    set_logger(cfg)

    # Load region fields
    rg = Region.load_from_config(cfg["observations"])
    for field_id in rg.tt_fields["field_id"]:
        logger.info("Analysing field:"+str(field_id))

        # Load field
        gf = GALEXField.from_MAST(
            obs_id=field_id, obs_filter=cfg["observations"]["obs_filter"])

        # Run clustering
        gf.cluster_meanshift(cfg["cluster"]["add_upper_limits"],
                             **cfg["cluster"]["meanshift"])

        # Get light curve
        lcs = gf.get_light_curve([0, 12])
        print(lcs)

        gf.plot_light_curve(12)
        plt.show()

        # Write field out
        field_file_name = cfg["general"]["out_dir"] + \
            "/field_" + str(field_id)+".fits"
        gf.write_to_fits(file_name=field_file_name)

    rg.info()

    # Write out regions
    region_out_name = cfg["general"]["out_dir"] + \
        "/region_" + cfg["general"]["name"]+".fits"
    rg.write_to_fits(file_name=region_out_name)

    # Write used config file
    yaml_out_name = cfg["general"]["out_dir"] + "/uvva_ran_cfg.yaml"
    with open(yaml_out_name, "w") as yaml_file:
        yaml.dump(cfg, yaml_file)
