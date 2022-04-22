#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that runs the UVVA pipeline.
"""

# %% Import modules
import toml
from loguru import logger
import sys
import datetime
import os
import argparse
from uvva.region import Region
from uvva.field import GALEXField

# %% Argument parsing and confg file loadings
parser = argparse.ArgumentParser()
parser.add_argument('--cfg', type=str, default=os.getcwd()+'/uvva_cfg.toml')
args = parser.parse_args()

with open(args.cfg) as file:
    uvva_cfg = toml.load(file)

if uvva_cfg["general"]["out_dir"] == "CWD":
    uvva_cfg["general"]["out_dir"] = os.getcwd()


# %% Setup logger
log_config = {
    "handlers": [
        {"sink": sys.stdout,
         "format": "<green>{time:YYYY-MM-DD HH:mm:ss.SSSS}</green> <cyan>{name}</cyan>:<cyan>{line}</cyan> | <level>{level}:</level> {message} ",
         "level": uvva_cfg["general"]["log_level"],
         "colorize":True,
         "backtrace":True,
         "diagnose":True},
        {"sink": uvva_cfg["general"]["log_file"],
         "serialize": True,
         "backtrace":True,
         "diagnose":True},
    ],
}
logger.configure(**log_config)
logger.enable("uvva")

logger.info("Runing '"+__file__+"'")
logger.debug("Config. file: '"+args.cfg+"'")
logger.debug("Output directory: '" + uvva_cfg["general"]["out_dir"]+"'")
logger.debug("Log. file: '" + uvva_cfg["general"]["log_file"]+"'")

# %% Create region
rg = Region()
rg.load_from_config(uvva_cfg["observations"])
for field_id in rg.tt_fields["field_id"]:
    logger.info("Analysing field:"+str(field_id))
    #gf = GALEXField(obs_id=field_id, filter=obs["obsfilter"])


# if __name__ == '__main__':
