import os
from glob import glob

import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import wcs
from IPython.display import Image, clear_output, display
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm

import vasca.utils as vutils
from vasca.resource_manager import ResourceManager

# Settings

# Input/output directories
with ResourceManager() as rm:
    root_data_dir = rm.get_path("gal_ds_fields", "lustre")

# root_data_dir = "/Users/julianschliwinski/GALEX_DS/GALEX_DS_GCK_fields"

# Dry-run, don't export final list
dry_run = False

# Debugging switch
is_debug = False

hide_progress = False

refresh = False

coadd_cnt_paths = list()
coadd_rrhr_paths = list()

# Loops over drift scan directories
scan_names = [path.split(os.sep)[-1] for path in glob(f"{root_data_dir}/*")]
for idx_scan, scan_name in tqdm(
    enumerate(scan_names), total=len(scan_names), desc="Scans", disable=hide_progress
):
    # Debugging
    if idx_scan > 0 and is_debug:
        break

    # Loops over fields
    field_names = [
        path.split(os.sep)[-1] for path in glob(f"{root_data_dir}/{scan_name}/*")
    ]
    for idx_field, field_name in tqdm(
        enumerate(field_names),
        total=len(field_names),
        desc="Fields",
        disable=hide_progress,
    ):
        # Debugging
        if idx_field > 0 and is_debug:
            break

        # Collects file paths across all visit directories
        # Counts maps
        img_int_paths = glob(
            f"{root_data_dir}/{scan_name}/{field_name}/*/*-nd-cnt.fits.gz"
        )
        # Effective exposure time map
        img_rrhr_paths = glob(
            f"{root_data_dir}/{scan_name}/{field_name}/*/*-nd-rrhr.fits.gz"
        )

        # Stack image maps
        # Loops over image types
        for path_list, path_out_file_suffix, coadd_path_list in zip(
            [img_int_paths, img_rrhr_paths],
            ["-nd-cnt-coadd.fits", "-nd-rrhr-coadd.fits"],
            [coadd_cnt_paths, coadd_rrhr_paths],
        ):
            path_out_file = f"{root_data_dir}/{scan_name}/{field_name}/{field_name}{path_out_file_suffix}"
            if not refresh and not os.path.isfile(path_out_file):
                # Loops over individual files
                for idx_img, path in enumerate(path_list):

                    # Initialize stack image from first file
                    with fits.open(path) as hdul:
                        if idx_img == 0:
                            img = hdul[0].data
                            img_wcs = wcs.WCS(hdul[0].header)
                        else:
                            img += hdul[0].data

                # Export
                hdu = fits.CompImageHDU(img, header=img_wcs.to_header())
                hdu.writeto(path_out_file, overwrite=True)

            coadd_path_list.append(path_out_file)
