#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that runs the UVVA pipeline.
"""

# %% import modules
import toml
from loguru import logger
import sys
import datetime
import os
import argparse

# %% Argument parsing and confg file loadings
# Load parsed arguments
parser = argparse.ArgumentParser()
parser.add_argument('--cfg', type=str, default=os.getcwd()+'/uvva_cfg.toml')
args = parser.parse_args()


# Open config file & set output directory
with open(args.cfg) as file:
    uvva_cfg = toml.load(file)

if uvva_cfg["general"]["out_dir"] == "CWD":
    uvva_cfg["general"]["out_dir"] = os.getcwd()


# %% Setup logger
log_config = {
    "handlers": [
        {"sink": sys.stdout,
         "format": "<green>{time:YYYY-MM-DD HH:mm:ss.SSSSSS}</green> <blue>{level}</blue>: {message}",
         "level": uvva_cfg["general"]["log_level"],
         "colorize":True},
        {"sink": uvva_cfg["general"]["log_file"],
         "serialize": True,
         "backtrace":True,
         "diagnose":True},
    ]
}
logger.configure(**log_config)
logger.enable("uvva")

if __name__ == '__main__':
    logger.info("Runing '"+__file__+"'")
    logger.debug("Config. file: '"+args.cfg+"'")
    logger.debug("Output directory: '" + uvva_cfg["general"]["out_dir"]+"'")
    logger.debug("Log. file: '" + uvva_cfg["general"]["log_file"]+"'")
