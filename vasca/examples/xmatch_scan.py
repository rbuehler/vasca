import os

from loguru import logger
from tqdm import tqdm

from vasca.acms_api import get_acms_cat_info
from vasca.region import Region
from vasca.xmatch import xmatch

data_dir = "/afs/ifh.de/group/ultrasat/scratch/schliwij/data/vasca_data/xmatch_scan"
out = f"{data_dir}/acms_queries/cat_scan"

logger.enable("vasca")

file = f"{data_dir}/region_TDS_All.fits"
rg = Region()
rg.load_from_fits(file)

vasca_cat = rg.tt_sources["ra", "dec", "rg_src_id"][rg.tt_sources["sel"]]

acms_cat_list = list(get_acms_cat_info().keys())

for catalog in tqdm(acms_cat_list, total=len(acms_cat_list), desc="ACMS Catalogs"):
    try:
        out_file = f"{out}/tt_matched_{catalog}.csv"

        if not os.path.isfile(out_file):
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
            tt_matched.write(f"{out}", format="csv", overwrite=True)
        else:
            logger.info(
                f"Skipping x-matching for catalog '{catalog}'. "
                f"File already exists at '{out_file}'"
            )
    except Exception as e:
        logger.exception(f"Failed to query catalog '{catalog}'. Exception:\n{e}")
