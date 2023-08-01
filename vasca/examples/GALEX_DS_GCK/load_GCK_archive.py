"""
Script to download archival data from the GALEX-CAUSE Kepler survey (GCK).
The data is hosted at MAST and only accessible through https here:
https://archive.stsci.edu/pub/galex/KS/

The Script uses the a text file listing all URLs of files hosted on the webserver.

More info: DEEP GALEX UV SURVEY OF THE KEPLER FIELD. I. POINT SOURCE CATALOG
https://iopscience.iop.org/article/10.1088/0004-637X/813/2/100
"""

import os
import pathlib

import pandas as pd
import requests
from tqdm import tqdm

# Read in URLs form the file manifest as DataFrame. Each column corresponds to a
# directory level but URLs have different depth. To handle this, the maximum number of
# levels are inferred by a look-ahead counting the number of levels line by line.

# Dynamically generate column names

# Input
try:
    manifest_path = pathlib.Path(__file__).parent / "pipeDirectoryListing.txt"
except NameError:
    # Handle jupyter notebook problem with __file__
    manifest_path = rf"{os.path.abspath('')}/pipeDirectoryListing.txt"

# Delimiter
data_file_delimiter = "/"

# The max column count a line in the file could have
largest_column_count = 0

# Count number of lines
lcounter = 0

# Loop the data lines
with open(manifest_path, "r") as temp_f:
    # Read the lines
    lines = temp_f.readlines()

    for line in lines:
        # Count the column count for the current line
        column_count = len(line.split(data_file_delimiter)) + 1

        # Set the new most column count
        largest_column_count = (
            column_count
            if largest_column_count < column_count
            else largest_column_count
        )

        lcounter += 1

# Generate column names (will be L0, L1, L2, ..., L<largest_column_count - 1>)
column_names = [f"L{i}" for i in range(0, largest_column_count)]

print("URL manifest info:")
print(f"    - Number of lines: {lcounter}")
print(f"    - Maximum number of directory levels: {largest_column_count}")

# Combine individual level columns with full path
# Important columns:
# L8: Drift scan (15 continuous scans)
# L12: Visit (i.e, scan repetition, 20 visits on average)
# L13: image/catalog data ("svXX" corresponds to field number,
#      short scans (1−3, 13−15) consist of 9 fields,
#      long scans (4-12) consist of 14 fields each).
df_manifest = pd.concat(
    [
        pd.read_csv(
            manifest_path,
            header=None,
            delimiter=data_file_delimiter,
            names=column_names,
        ),
        pd.read_csv(manifest_path, delimiter=",", header=None, names=["path"]),
    ],
    axis=1,
)

# Download


def get_n_fields(id_scan):
    """
    Returns the number of sub-fields per drift scan:
    9 and 14 fields for the short scans (1−3 and 13−15) and long scans (4−12)
    """
    scan_num = int(id_scan.split("_")[-1])
    if scan_num <= 3 or scan_num >= 13:
        # short scans
        n_fields = 9
    elif scan_num >= 4 and scan_num <= 12:
        # long scans
        n_fields = 14
    else:
        raise ValueError(f"Unexpected scan ID '{id_scan}'.")

    return n_fields


# Settings

# Path to root data directory

# Debugging
# root_data_dir = "/Users/julianschliwinski/GALEX_DS/galex_kepler"

# On WGS
root_data_dir = (
    "/lustre/fs24/group/ultrasat/vasca_data/uc_science/uvvarcat/GALEX_DS_GCK_fields"
)

# Refresh download (possibly overwrites existing files)
refresh = False

# Supress progress bar
hide_progress = False

# Dry-run, no download
dry_run = True

# Group by scan and visit
df_grpd = df_manifest.groupby(["L8", "L11"])

# Counters
n_files = 0
n_vis = 0

# Loops over drift scans and visits
for i, (id_scan, id_visit) in tqdm(
    enumerate(df_grpd.groups.keys()), disable=hide_progress
):
    # Debugging
    # if i > 100: break
    # if id_scan != "29200-KEPLER_SCAN_001":
    #     break
    # if id_visit != "0012_img":
    #     continue

    # Loops over fields (sv = "sub-visit")
    for k, id_sv in tqdm(
        enumerate([f"sv{num+1:02}" for num in range(get_n_fields(id_scan))]),
        total=get_n_fields(id_scan),
        desc=f"Fields ({id_scan}, {id_visit})",
        disable=hide_progress,
    ):
        # Debugging
        # if id_sv != "sv01":
        #     break

        # List of file name endings
        files_select = [
            "nd-cnt.fits",
            "nd-rrhr.fits",
            "nd-int.fits",
            "nd-flags.fits",
            "xd-mcat.fits",
        ]
        n_files_select = len(files_select)

        # Selects files corresponding to field ID
        files_select_str = "|".join(files_select)
        query_str = (
            f"L13.str.contains('{files_select_str}') and "
            f"L13.str.contains('{id_sv}')"
        )
        df_select = df_grpd.get_group((id_scan, id_visit)).query(query_str)
        df_select.sort_values("L13", inplace=True)

        # Download data
        if len(df_select) > 0:
            if len(df_select) != n_files_select:
                print(
                    f"Warning: unexpected number of files ({n_files_select}), "
                    f"got {len(df_select)} "
                    f"({id_scan}, {id_sv}, {id_visit})"
                )

            # Create output directory path
            # Data is sorted by scan and visit (swap field and visit directories
            # compared to web server)
            out_dir = (
                f"{root_data_dir}/"
                f"{id_scan}/{id_scan}_{id_sv}/"
                f"{id_scan}_{id_sv}_{id_visit}"
            )

            # Create the output directory if it doesn't exist
            os.makedirs(out_dir, exist_ok=True)

            # Loops over data files
            for idx, file_row in df_select.iterrows():
                # File name
                file_name = file_row.L13

                file_url = file_row.path

                # Construct the output file path
                out_file_path = os.path.join(out_dir, file_name)

                # Download the file if it doesn't already exist or download is forced
                if not dry_run:
                    if not os.path.isfile(out_file_path) or refresh:
                        response = requests.get(file_url)
                        if response.status_code == 200:
                            with open(out_file_path, "wb") as f:
                                f.write(response.content)
                            print(f"Downloaded: {file_url}")
                        else:
                            print(f"Failed to download: {file_url}")
                    else:
                        print(f"File exists: {file_url}")

            # counts number of files
            n_files += len(df_select)

    # Counts number of visits
    n_vis += 1

print(
    f"Downloaded {n_files/n_files_select:1.1f} exposures "
    f"({n_files} files) for {n_vis} visits."
)
