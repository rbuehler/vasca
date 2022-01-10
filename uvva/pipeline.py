#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that runs the UVVA pipeline.
"""

# %% import modules
import yaml
import logging
import sys
import datetime
import os

CLASS_DIR = os.path.dirname(os.path.abspath(__file__))
print(CLASS_DIR)
cfg_file = CLASS_DIR + "/pipeline_config.yml"


def set_logger(loglevel, logfile):
    """

    Parameters
    ----------
    loglevel : str
        Logging level, NOTSET, DEBUG, INFO, WARNING, ERROR or CRITICAL.
    logfile : str
        output log file.

    Returns
    -------
    None.

    """
    # Set up logger, log to file and stdout
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
    a_logger = logging.getLogger()
    a_logger.setLevel(cfg["run"]["loglevel"])
    file_handler = logging.FileHandler(cfg["run"]["logfile"])
    file_handler.setFormatter(formatter)
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setFormatter(formatter)
    a_logger.addHandler(file_handler)
    a_logger.addHandler(stdout_handler)
    logging.info("Start time: " + str(datetime.datetime.now()))


def run(cfg_file):
    """

    Parameters
    ----------
    cfg_file : str
        YAML configuration file of pipeline parameters

    Returns
    -------
    None.

    """

    # Load YAML configuration
    with open(cfg_file, "r") as ymlfile:
        cfg = yaml.safe_load(ymlfile)

    set_logger(cfg["run"]["loglevel"], cfg["run"]["logfile"])
