import os
from glob import glob

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

coadd_cnt_paths = list()
coadd_rrhr_paths = list()

# Loops over drift scan directories
scan_names = [path.split(os.sep)[-1] for path in glob(f"{root_data_dir}/*")]
for idx_scan, scan_name in tqdm(
    enumerate(sorted(scan_names)),
    total=len(scan_names),
    desc="Scans",
    disable=hide_progress,
):
    # Debugging
    if idx_scan > 0 and is_debug:
        break

    # Loops over fields
    field_names = [
        path.split(os.sep)[-1] for path in glob(f"{root_data_dir}/{scan_name}/*")
    ]
    for idx_field, field_name in tqdm(
        enumerate(sorted(field_names)),
        total=len(field_names),
        desc=f"Fields ({scan_name})",
        disable=hide_progress,
    ):
        # Debugging
        if idx_field > 0 and is_debug:
            break

        # Collects file paths across all visit directories
        # Counts maps
        img_cnt_paths = glob(
            f"{root_data_dir}/{scan_name}/{field_name}/*/*-nd-cnt.fits.gz"
        )
        # Effective exposure time map
        img_rrhr_paths = glob(
            f"{root_data_dir}/{scan_name}/{field_name}/*/*-nd-rrhr.fits.gz"
        )

        # Stack image maps
        # Loops over image types
        for path_list, path_out_file_suffix, coadd_path_list in zip(
            [img_cnt_paths, img_rrhr_paths],
            ["-nd-cnt-coadd.fits.gz", "-nd-rrhr-coadd.fits.gz"],
            [coadd_cnt_paths, coadd_rrhr_paths],
        ):
            # Path to the coadd file
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

            if refresh or not os.path.isfile(path_out_file):
                # Loops over individual files
                img_counter = 0
                for idx_img, path in enumerate(sorted(path_list)):
                    # Gets visit name
                    vis_name = path.split(os.sep)[-2]
                    # Load image only if not bad quality
                    if vis_name not in visit_is_bad:
                        # Initialize stack image from first file
                        with fits.open(path) as hdul:
                            if img_counter == 0:
                                img = hdul[0].data
                                img_wcs = wcs.WCS(hdul[0].header)
                            else:
                                img += hdul[0].data
                        img_counter += 1

                # Export
                if img_counter > 0:
                    hdu = fits.PrimaryHDU(img, header=img_wcs.to_header())
                    hdu.writeto(path_out_file, overwrite=True)

            if img_counter > 0:
                coadd_path_list.append(path_out_file)

        # Create intensity coadd
        coadd_cnt_path = (
            f"{root_data_dir}/{scan_name}/{field_name}/"
            f"{field_name}-nd-cnt-coadd.fits.gz"
        )
        coadd_rrhr_path = (
            f"{root_data_dir}/{scan_name}/{field_name}/"
            f"{field_name}-nd-rrhr-coadd.fits.gz"
        )
        with fits.open(coadd_cnt_path) as hdul:
            img_coadd_cnt = hdul[0].data
            img_coadd_wcs = wcs.WCS(hdul[0].header)
        with fits.open(coadd_rrhr_path) as hdul:
            img_coadd_rrhr = hdul[0].data

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

        hdu = fits.PrimaryHDU(
            img_coadd_cnt / np.where(img_coadd_rrhr == 0, np.nan, img_coadd_rrhr),
            header=img_coadd_wcs.to_header(),
        )
        hdu.writeto(coadd_int_path, overwrite=True)
