#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that runs the VASCA pipeline.
"""

import argparse
import os
import sys
from multiprocessing import Pool


import yaml
from loguru import logger

from vasca.region import Region
from vasca.utils import get_region_field_id
from vasca.tables_dict import dd_vasca_tables


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


def run_field(obs_nr, field, vasca_cfg):
    """
    Run analysis on a single field

    Parameters
    ----------
    obs_nr : int
        Observation number in the config file.
    field : vasca.field
        Field to run analysis on.
    vasca_cfg :
        VASCA configuration file

    Returns
    -------
    field : vasca.field
        Modified field with results

    """
    logger.info("Analysing field:" + str(field.field_id))

    # Create directory structure for fields
    field_out_dir = (
        vasca_cfg["general"]["out_dir_base"]
        + "/"
        + vasca_cfg["general"]["name"]
        + "/fields/"
    )

    # Create folder
    if not os.path.exists(field_out_dir):
        os.makedirs(field_out_dir)

    # Get configuration for this field
    obs_cfg = vasca_cfg["observations"][obs_nr]

    # Apply selections
    field.select_rows(obs_cfg["selection"]["det_quality"], remove_unselected=True)

    # Run clustering
    field.cluster_meanshift(
        **obs_cfg["cluster_det"]["meanshift"],
    )

    # # Calculate source variables from light curve
    field.set_src_stats()

    # Write out field
    field.write_to_fits(field_out_dir + "field_" + field.field_id + ".fits")

    # # Remove some items which are not further needed to free memory
    # # del field.__dict__["tt_detections"]
    # # field._table_names.remove("tt_detections")
    # field.select_rows(obs_cfg["selection"]["det_association"], remove_unselected=True)
    field.ref_img = None
    field.ref_wcs = None

    # Keep only base class columns
    for tt_name in field._table_names:
        cols = dd_vasca_tables["base_field"][tt_name]["names"]
        field.__dict__[tt_name] = field.__dict__[tt_name][cols]

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

    # Load region fields
    rg = Region.load_from_config(vasca_cfg)
    rg.add_table_from_fields("tt_visits")

    # Setup output directory
    region_dir = (
        vasca_cfg["general"]["out_dir_base"] + "/" + vasca_cfg["general"]["name"] + "/"
    )
    if not os.path.exists(region_dir):
        os.makedirs(region_dir)

    # Prepare fields to run on
    fd_pars = list()
    obs_nr = 0
    rg.tt_fields.add_index("field_id")
    for obs in vasca_cfg["observations"]:
        for field_id in obs["obs_field_ids"]:
            field_id = get_region_field_id(
                obs_field_id=field_id,
                observaory=obs["observatory"],
                obs_filter=obs["obs_filter"],
            )
            fd = rg.tt_fields.loc["field_id", field_id]
            rg_fd_id = fd["rg_fd_id"]
            fd_pars.append([obs_nr, rg.fields[rg_fd_id]])
        obs_nr += 1

    # Run each field in a separate process in parallel
    with Pool(vasca_cfg["general"]["nr_cpus"]) as pool:
        pool_return = pool.starmap(
            run_field,
            [(*fd_par, vasca_cfg) for fd_par in fd_pars],
        )
    pool.join()

    # update region fields
    for field in pool_return:
        fd = rg.tt_fields.loc["field_id", field.field_id]
        rg_fd_id = fd["rg_fd_id"]
        rg.fields[rg_fd_id] = field

    # Add field tables to region
    rg.add_table_from_fields("tt_sources")
    rg.add_table_from_fields("tt_detections", only_selected=False)

    rg.cluster_meanshift(**vasca_cfg["cluster_src"]["meanshift"])

    # Calculate source statistics
    rg.set_src_stats()

    # Remove sources with quality cuts
    rg.select_rows(vasca_cfg["selection"]["src_quality"], remove_unselected=True)

    # Select variable sources
    rg.select_rows(vasca_cfg["selection"]["src_variability"], remove_unselected=False)

    # Remove detections not associated to selected sources
    rg.select_rows(vasca_cfg["selection"]["det_association"], remove_unselected=True)

    # Store reference sources
    if vasca_cfg["general"]["save_ref_srcs"]:
        rg.add_table_from_fields("tt_ref_sources")

    # Write out regions
    rg.write_to_fits(
        file_name=region_dir + "/region_" + vasca_cfg["general"]["name"] + ".fits"
    )

    # Write used config file
    yaml_out_name = region_dir + "/cfg_ran_" + vasca_cfg["general"]["name"] + ".yaml"
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
