import os
from glob import glob
from pprint import pprint

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import wcs
from tqdm import tqdm

from vasca.resource_manager import ResourceManager

# Settings

# Input/output directories
with ResourceManager() as rm:
    root_data_dir = rm.get_path("gal_ds_fields", "lustre")
    visits_list_dir = (os.sep).join(
        rm.get_path("gal_ds_visits_list", "lustre").split(os.sep)[:-1]
    )

# Load visual image quality table
df_img_quality = pd.read_csv(
    f"{visits_list_dir}/GALEX_DS_GCK_visits_img_quality.csv", index_col=0
)
# List of visits with bad image quality
visit_is_bad = df_img_quality.query("quality in ['bad']").vis_name.tolist()

# Dry-run, don't export final list
dry_run = False

# Debugging switch
is_debug = False

# Don't show the progress bar
hide_progress = False

# Recreate the coadd image
refresh = True

# Delete an existing coadd
delete_existing = True

# Store paths of created files
coadd_cnt_paths = list()
coadd_rrhr_paths = list()

# Loops over drift scan directories
scan_names = sorted([path.split(os.sep)[-1] for path in glob(f"{root_data_dir}/*")])
n_scans = len(scan_names)
for idx_scan, scan_name in tqdm(
    enumerate(scan_names),
    total=n_scans,
    desc="Scans",
    disable=hide_progress,
):
    # Debugging
    if idx_scan > 0 and is_debug:
        break

    # Loops over fields
    field_names = sorted(
        [path.split(os.sep)[-1] for path in glob(f"{root_data_dir}/{scan_name}/*")]
    )
    n_fields = len(field_names)
    for idx_field, field_name in tqdm(
        enumerate(field_names),
        total=n_fields,
        desc=f"Fields ({scan_name})",
        disable=hide_progress,
    ):
        # Debugging
        if idx_field > 0 and is_debug:
            break

        # Selects visit names
        visit_names = [
            path.split(os.sep)[-1]
            for path in sorted(glob(f"{root_data_dir}/{scan_name}/{field_name}/*-img"))
        ]
        visit_names = [name for name in visit_names if name not in visit_is_bad]
        n_visits = len(visit_names)

        # Don't continue if all visits are bad
        if n_visits > 0:

            img_cnt_paths = list()
            img_rrhr_paths = list()
            for visit_name in visit_names:
                img_cnt_path = glob(
                    f"{root_data_dir}/{scan_name}/{field_name}/{visit_name}/"
                    f"*-nd-cnt.fits.gz"
                )[0]
                img_rrhr_path = glob(
                    f"{root_data_dir}/{scan_name}/{field_name}/{visit_name}/"
                    f"*-nd-rrhr.fits.gz"
                )[0]
                img_cnt_paths.append(img_cnt_path)
                img_rrhr_paths.append(img_rrhr_path)

            if len(img_cnt_paths) != len(img_rrhr_paths):
                raise ValueError("Expected same number of counts and exposure maps.")

            # Stack image maps
            # Loops over image types
            for path_list, path_out_file_suffix, coadd_path_list in zip(
                [img_cnt_paths, img_rrhr_paths],
                ["-nd-cnt-coadd.fits.gz", "-nd-rrhr-coadd.fits.gz"],
                [coadd_cnt_paths, coadd_rrhr_paths],
            ):
                # Path to the output coadd file
                path_out_file = (
                    f"{root_data_dir}/{scan_name}/{field_name}/"
                    f"{field_name}{path_out_file_suffix}"
                )

                # Delete existing coadd file
                if delete_existing and os.path.isfile(path_out_file):
                    os.remove(path_out_file)
                    alt_path = path_out_file.rstrip(".gz")
                    if os.path.isfile(alt_path):
                        os.remove(alt_path)

                # Do stacking if coadd file does not exist or refresh is specified
                if refresh or not os.path.isfile(path_out_file):
                    # Loops over individual visits/files
                    for idx_img, path in enumerate(sorted(path_list)):
                        # Initialize stack image from first file
                        with fits.open(path) as hdul:
                            if idx_img == 0:
                                img = hdul[0].data.clip(min=0)
                                img_wcs = wcs.WCS(hdul[0].header)
                            else:
                                img += hdul[0].data.clip(min=0)

                    # Export FITS file
                    hdu = fits.PrimaryHDU(img, header=img_wcs.to_header())
                    hdu.writeto(path_out_file, overwrite=True)

                # Save path of created file
                coadd_path_list.append(path_out_file)

            # Create intensity coadd

            # Coadd file paths
            coadd_cnt_path = (
                f"{root_data_dir}/{scan_name}/{field_name}/"
                f"{field_name}-nd-cnt-coadd.fits.gz"
            )
            coadd_rrhr_path = (
                f"{root_data_dir}/{scan_name}/{field_name}/"
                f"{field_name}-nd-rrhr-coadd.fits.gz"
            )
            coadd_int_path = (
                f"{root_data_dir}/{scan_name}/{field_name}/"
                f"{field_name}-nd-int-coadd.fits.gz"
            )

            # Delete existing coadd file
            if delete_existing and os.path.isfile(coadd_int_path):
                os.remove(coadd_int_path)
                alt_path = coadd_int_path.rstrip(".gz")
                if os.path.isfile(alt_path):
                    os.remove(alt_path)

            # Compute ratio: Counts over effective exposure = intensity

            # Load images
            with fits.open(coadd_cnt_path) as hdul:
                img_coadd_cnt = hdul[0].data
                img_coadd_wcs = wcs.WCS(hdul[0].header)
            with fits.open(coadd_rrhr_path) as hdul:
                img_coadd_rrhr = hdul[0].data

            # Compute ratio, avoid divide-by-zero problems by setting zeros to NaN
            # for the computation and setting the same elements back to zeros again.
            # Saving NaN to FITS and using CompImageHDUs in VASCA led to corrupt images
            img_coadd_int = img_coadd_cnt / np.where(
                img_coadd_rrhr == 0.0, np.nan, img_coadd_rrhr
            )
            np.nan_to_num(img_coadd_int, copy=False)

            # Export FITS file
            hdu = fits.PrimaryHDU(
                img_coadd_int,
                header=img_coadd_wcs.to_header(),
            )
            hdu.writeto(coadd_int_path, overwrite=True)
