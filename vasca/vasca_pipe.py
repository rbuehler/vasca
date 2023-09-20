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
from astropy import units as uu
from astropy.table import unique
from loguru import logger
from yamlinclude import YamlIncludeConstructor

from vasca.region import Region
from vasca.tables_dict import dd_vasca_tables
from vasca.utils import dd_obs_id_add

import numpy as np

YamlIncludeConstructor.add_to_loader_class(
    loader_class=yaml.FullLoader
)  # , base_dir="."


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
        # vasca_cfg = yaml.safe_load(file)  # yaml.
        vasca_cfg = yaml.load(file, Loader=yaml.FullLoader)

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
    os.makedirs(log_dir, exist_ok=True)

    # create log file name
    log_file_name = "log_" + vasca_cfg["general"]["name"] + ".log"
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
            # {
            #     "sink": log_dir + log_file_name,
            #     "serialize": True,
            #     "backtrace": True,
            #     "diagnose": True,
            # },
        ],
    }
    logger.configure(**log_cfg)
    # logger.add(log_dir + log_file_name)
    logger.enable("vasca")

    logger.info("Runing '" + __file__ + "'")
    logger.debug("Config. file: '" + vasca_cfg["cfg_file"] + "'")
    # logger.debug("Output log. file: '" + log_dir + log_file_name + "'")


def keep_base_field(field):
    "Remove field data that is not needed for further analysis to save memory"
    # Remove detections which did not pass selections and where not used in clustering
    field.remove_unselected("tt_detections")
    if hasattr(field, "tt_coadd_detections"):
        field.remove_unselected("tt_coadd_detections")
    # Remove image data
    field.ref_img = None
    field.ref_wcs = None
    field.vis_img = None

    # Keep only base class columns
    for tt_name in field._table_names:
        cols = dd_vasca_tables["base_field"][tt_name]["names"]
        field.__dict__[tt_name] = field.__dict__[tt_name][cols]


def run_field(obs_nr, field_id, rg, vasca_cfg):
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
    logger.info("Analysing field:" + str(field_id))

    try:
        field = rg.get_field(
            field_id=field_id,
            load_method="VASCA",
            mast_products=vasca_cfg["resources"]["load_products"],
            field_kwargs=vasca_cfg["resources"]["field_kwargs"],
        )
    except Exception as e:
        logger.exception(
            f"Faild to run pipeline for field '{field_id}'. "
            f"Returning None. Exception:\n {e}"
        )
        return None

    # Create directory structure for fields
    field_out_dir = (
        vasca_cfg["general"]["out_dir_base"]
        + "/"
        + vasca_cfg["general"]["name"]
        + "/fields/"
    )

    # Create folder
    os.makedirs(field_out_dir, exist_ok=True)

    # Get configuration for this field
    obs_cfg = vasca_cfg["observations"][obs_nr]

    # Apply detections selections
    # Visit detections
    field.select_rows(obs_cfg["selection"]["det_quality"], remove_unselected=False)
    # Coadd detections
    if hasattr(field, "tt_coadd_detections"):
        field.select_rows(
            obs_cfg["selection"]["coadd_det_quality"], remove_unselected=False
        )

    # Cluster and do source stats
    if len(field.tt_detections) > 0:
        # Run clustering
        logger.info("Clustering field detections")
        field.cluster_meanshift(
            **obs_cfg["cluster_det"]["meanshift"],
        )

        # Calculate source variables from light curve
        # field.set_src_stats(src_id_name="fd_src_id")

    # Select detections used in clustering
    field.select_rows(obs_cfg["selection"]["det_association"], remove_unselected=False)

    # Write out field
    field.write_to_fits(field_out_dir + "field_" + field.field_id + ".fits")

    # Keep only data needed for further analysis
    keep_base_field(field)
    logger.info("Done with field:" + str(field_id))

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

    # Setup output directory
    os.makedirs(rg.region_path, exist_ok=True)

    # Set field selection to False and only back to Trua below if it contains detections
    rg.tt_fields["sel"] = np.zeros(len(rg.tt_fields), dtype=bool)
    rg.tt_fields.add_index("field_id")

    # Ifrunnign for the first time, or with new field settings, run fields
    if vasca_cfg["general"]["run_fields"]:
        # Prepare fields to run in parallel
        fd_pars = list()  # List of observations and field_ids in the config file
        vobs = vasca_cfg["observations"]
        for obs_nr in range(len(vobs)):
            for field_nr in vobs[obs_nr]["obs_field_ids"]:
                field_id = dd_obs_id_add[
                    vobs[obs_nr]["observatory"] + vobs[obs_nr]["obs_filter"]
                ] + str(field_nr)
                fd_pars.append([obs_nr, field_id, rg, vasca_cfg])

        # Run each field in a separate process in parallel

        nr_cpus = vasca_cfg["general"]["nr_cpus"]
        logger.info(f"Analyzing {len(fd_pars)} fields on {nr_cpus} parallel threads.")
        with Pool(processes=nr_cpus) as pool:
            pool_return = pool.starmap(run_field, fd_pars)
        pool.join()
        logger.info("Done analyzing individual fields.")

        # update region fields
        for field in pool_return:
            # Check if field was filled or is empty for any reason
            if hasattr(field, "tt_detections") and len(field.tt_detections) > 0:
                rg.fields[field.field_id] = field
                fd_idx = rg.tt_fields.loc_indices["field_id", field.field_id]
                rg.tt_fields["sel"][fd_idx] = True
                logger.info(f"Added field {field.field_id} from pool to region")
            else:
                logger.warning(
                    "Ignoring field, as it was empty or had no detections from the pool"
                )
    else:
        for rg_fd_id in rg.tt_fields["rg_fd_id"]:
            field = rg.get_field(rg_fd_id=rg_fd_id, load_method="FITS", add_field=False)
            if hasattr(field, "tt_detections") and len(field.tt_detections) > 0:
                # Keep only data needed for further analysis
                keep_base_field(field)
                rg.fields[field.field_id] = field
                fd_idx = rg.tt_fields.loc_indices["field_id", field.field_id]
                rg.tt_fields["sel"][fd_idx] = True
                logger.info(f"Added field: {field.field_id} from file")
            else:
                logger.warning(
                    "Ignoring field, as it was empty or had no detections from the file"
                )

    # Add field tables to region
    # For visits merge obs_filter_id and remove doubles
    rg.add_table_from_fields("tt_visits", sum_obs_filter=True)
    rg.tt_visits = unique(rg.tt_visits, keys="vis_id")
    rg.add_table_from_fields("tt_sources")
    rg.add_table_from_fields("tt_detections", only_selected=False)
    if vasca_cfg["resources"]["coadd_exists"]:
        rg.add_table_from_fields("tt_coadd_detections")

    del rg.fields  # All that was needed has been transferred to region tables

    # Cluster field sources and co-adds
    rg.cluster_meanshift(**vasca_cfg["cluster_src"]["meanshift"])
    if vasca_cfg["resources"]["coadd_exists"]:
        rg.cluster_meanshift(**vasca_cfg["cluster_coadd_dets"]["meanshift"])

    # Calculate source statistics
    rg.set_src_stats(src_id_name="rg_src_id")
    if vasca_cfg["resources"]["coadd_exists"]:
        rg.set_src_stats(src_id_name="coadd_src_id")

    # Add hardness ratio (default obs_filter_id 1 & 2 for GALEX)
    rg.set_hardness_ratio()

    # Match sources to coadd sources
    if vasca_cfg["resources"]["coadd_exists"]:
        rg.cross_match(
            dist_max=vasca_cfg["assoc_src_coadd"]["dist_max"] * uu.arcsec,
            dist_s2n_max=vasca_cfg["assoc_src_coadd"]["dist_s2n_max"],
        )

    # Select variable sources
    rg.select_rows(vasca_cfg["selection"]["src_variability"], remove_unselected=False)
    if vasca_cfg["resources"]["coadd_exists"]:
        rg.select_rows(
            vasca_cfg["selection"]["src_coadd_diff"], remove_unselected=False
        )

    # Remove detections not associated to selected sources
    if "det_association" in vasca_cfg["selection"].keys():
        rg.select_rows(
            vasca_cfg["selection"]["det_association"], remove_unselected=True
        )

    # Set source id table
    rg.set_src_id_info()

    # Write out region and reduced "catalog" region
    fname_base = rg.region_path + "/region_" + vasca_cfg["general"]["name"]
    rg.write_to_fits(file_name=fname_base + ".fits")

    # Write out catalog region file
    rc = rg.get_region_catalog()
    rc.write_to_fits(file_name=fname_base + "_cat.fits")

    # Write used config file
    yaml_out_name = (
        rg.region_path + "/cfg_ran_" + vasca_cfg["general"]["name"] + ".yaml"
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
