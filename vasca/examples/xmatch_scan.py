import os

from loguru import logger
from tqdm import tqdm

from vasca.acms_api import get_acms_cat_info
from vasca.region import Region
from vasca.xmatch import xmatch

# Data location
data_dir = "/afs/ifh.de/group/ultrasat/scratch/schliwij/data/vasca_data/xmatch_scan"
out = f"{data_dir}/acms_queries/cat_scan"

# Enable logging
logger.enable("vasca")

# Load VASCA region file
file = f"{data_dir}/region_TDS_All.fits"
rg = Region()
rg.load_from_fits(file)

# Select variable sources
vasca_cat = rg.tt_sources["ra", "dec", "rg_src_id"][rg.tt_sources["sel"]]

# Fetch ACMS catalog list
acms_cat_list = list(get_acms_cat_info().keys())

# Loops over ACMS catalog for x-matching
for catalog in tqdm(acms_cat_list, total=len(acms_cat_list), desc="ACMS Catalogs"):
    # Save x-matching catalog with this file name
    out_file = f"{out}/tt_matched_{catalog}.csv"

    # X-matching in case it hasn't been done (successfully) before
    if not os.path.isfile(out_file):
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

            logger.info(f"Writing query results to '{out_file}'")
            tt_matched.write(f"{out_file}", format="csv", overwrite=True)
        except Exception as e:
            logger.exception(f"Failed to query catalog '{catalog}'. Exception:\n{e}")
    else:
        logger.info(
            f"Skipping x-matching for catalog '{catalog}'. "
            f"File already exists at '{out_file}'"
        )
