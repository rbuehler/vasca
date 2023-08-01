"""
Create a complete list of all GALEX CAUSE Kepler survey visits.

This script reads all mcat files relative to the specified root data directory and
extracts relevant information from the FITS header. A combined list for all mcat files
is created. All required information is contained that VASCA field and visits tables can
be created.
"""

import os

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

import vasca.resource_manager as vascarm
import vasca.utils as vutils


def get_hdr_info(hdr):
    """
    Compile data from an mcat file header for general information about
    a drift scan observation.
    """
    # Extract header info
    hdr_info = {
        k: hdr[k]
        for k in [
            "TILENUM",
            "TILENAME",
            "OBJECT",
            "VISIT",
            "SUBVIS",
            "OBSDATIM",
            "NEXPSTAR",
            "NEXPTIME",
            "RA_CENT",
            "DEC_CENT",
            "GLONO",
            "GLATO",
        ]
    }

    # Create names and IDs for fields and visits
    field_name = (
        f'{hdr_info["TILENUM"]}-{hdr_info["TILENAME"]}_sv{hdr_info["SUBVIS"]:02}'
    )
    field_id = vutils.name2id(field_name, bits=64)
    vis_name = f'{field_name}_{hdr_info["VISIT"]:04}-img'
    vis_id = vutils.name2id(vis_name, bits=64)

    hdr_info.update(
        {
            "field_name": field_name,
            "field_id": field_id,
            "vis_name": vis_name,
            "vis_id": vis_id,
        }
    )

    # Time stamp
    time_bin_start = Time(hdr_info["NEXPSTAR"], format="unix").mjd
    hdr_info["time_bin_start"] = time_bin_start

    # Other info
    hdr_info.update(
        {
            "observatory": "GALEX_DS",  # GALEX drift scan
            "obs_filter": "NUV",
            "fov_diam": -999.99,  # FoV undefined for drift scan
            "sel": 0,
        }
    )

    return hdr_info


# Settings

# Input/output directories
with vascarm.ResourceManager() as rm:
    root_data_dir = rm.get_path("gal_ds_fields", "lustre")
    out_dir = (os.sep).join(
        rm.get_path("gal_ds_visits_list", "lustre").split(os.sep)[:-1]
    )

# Load visual image quality table
df_img_quality = pd.read_csv(
    f"{out_dir}/GALEX_DS_GCK_visits_img_quality.csv", index_col=0
)

# Load visual image quality table
df_img_quality = pd.read_csv(
    f"{out_dir}/GALEX_DS_GCK_visits_img_quality.csv", index_col=0
)
# List of visits with bad image quality
visit_is_bad = df_img_quality.query("quality in ['bad']").vis_name.tolist()

# Dry-run, don't export final list
dry_run = False

# Loops over mcat files and saves info
info = list()
for path, subdirs, files in os.walk(root_data_dir):
    for name in files:
        # Gets visit name
        if name.endswith("-xd-mcat.fits"):
            vis_name = os.path.join(path, name).split(os.sep)[-2]
        else:
            vis_name = None
        # Select mcat file for visits if not bad image quality
        if name.endswith("-xd-mcat.fits") and vis_name not in visit_is_bad:
            # Load mcat file and get relevant info
            mcat_path = os.path.join(path, name)
            with fits.open(mcat_path) as hdul:
                hdr_info = get_hdr_info(hdul[0].header)

            # Cross-checks
            # File name matches 'OBJECT' key
            if (
                mcat_path.split(os.sep)[-1].rstrip("-xd-mcat.fits")
                != hdr_info["OBJECT"]
            ):
                print(
                    "Warning: OBJECT key inconsistent "
                    f'(OBJECT: {hdr_info["OBJECT"]}, '
                    f'mcat: {mcat_path.split(os.sep)[-1].rstrip("-xd-mcat.fits")})'
                )
            # Visit directory name matches 'vis_name' key
            if mcat_path.split(os.sep)[-2] != hdr_info["vis_name"]:
                print(
                    "Warning: vis_name key inconsistent: "
                    f'(vis_name: {hdr_info["vis_name"]}, '
                    f"mcat directory: {mcat_path.split(os.sep)[-2]})"
                )

            info.append(hdr_info)

# Combines to astropy table via DataFrame
# because of problematic dtype handling of IDs
# tt_info = Table(info)  # this fails second cross-check
df_info = pd.DataFrame(info)
tt_info = Table.from_pandas(df_info)

# Cross-checks
# All visit IDs are unique
vis_ids = np.unique(tt_info["vis_id"])
if not len(vis_ids) == len(tt_info):
    raise ValueError("Non-unique visit IDs")
# All visit IDs have been consistently created from visit name
if not all(
    [
        int(vis_id) == vutils.name2id(vis_name, bits=64)
        for vis_id, vis_name in zip(tt_info["vis_id"], tt_info["vis_name"])
    ]
):
    raise ValueError("Inconsistent mapping vis_name to vis_id.")

if not dry_run:
    # Export
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    # FITS
    tt_info.write(f"{out_dir}/GALEX_DS_GCK_visits_list.fits", overwrite=True)
    # CSV
    df_info.to_csv(f"{out_dir}/GALEX_DS_GCK_visits_list.csv")
    # HTML
    df_info.to_html(f"{out_dir}/GALEX_DS_GCK_visits_list.html")

# Loads image quality table
visits_list_dir = (os.sep).join(
    rm.get_path("gal_ds_visits_list", "lustre").split(os.sep)[:-1]
)
df_img_quality = pd.read_csv(
    f"{visits_list_dir}/GALEX_DS_GCK_visits_img_quality.csv", index_col=0
)

# Passes if not a single bad visit is included in verify table
assert all(
    [
        name not in tt_info["vis_name"]
        for name in df_img_quality.query("quality == 'bad'").vis_name
    ]
)
