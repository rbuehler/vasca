import json
import os
import requests
import matplotlib.pyplot as plt
from astropy.io import fits
from tqdm import tqdm
from pprint import pprint
import numpy as np
from vasca.region import Region
from vasca.utils import nb_fig
from vasca.acms_api import acms_xmatch_query, get_acms_cat_info
from vasca.xmatch import xmatch_ampel, xmatch
from loguru import logger

data_dir = "/afs/ifh.de/group/ultrasat/scratch/schliwij/data/vasca_data/xmatch_scan"
out = f"{data_dir}/acms_queries/cat_scan"

logger.enable("vasca")

file = f"{data_dir}/region_TDS_All.fits"
rg = Region()
rg.load_from_fits(file)

vasca_cat = rg.tt_sources["ra","dec","rg_src_id"][rg.tt_sources["sel"]]

acms_cat_list = list(get_acms_cat_info().keys())

for catalog in tqdm(acms_cat_list, total=len(acms_cat_list), desc="ACMS Catalogs"):
    try:
        api = "ampel"
        radius = 15.0
        coords = vasca_cat["ra", "dec"]
        ids = vasca_cat["rg_src_id"]
    
        tt_matched = xmatch(
            api=api,
            catalog=catalog,
            radius=radius,
            coords=coords,
            ids=ids,
            dropna=False,
            hide_progress=False,
        )
    
        out_file = f"{out}/tt_matched_{catalog}.csv"
        if not os.path.isfile(out_file):
            logger.info(f"Writing query results to '{out_file}'")
            tt_matched.write(f"{out}")
    except Exception as e:
        logger.exception(f"Failed to query catalog '{catalog}'. Exception:\n{e}")

