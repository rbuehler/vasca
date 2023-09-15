#!/usr/bin/env python
# coding: utf-8

import gPhoton
import numpy as np
import os
from astropy import units as uu
from vasca.region import Region
from astropy.table import Table

region_name = "TDS"  # "WD" #"MDIS_10-800" # _ELAISN1
srcs_ids = [
    1891,
    3403,
    16149,
    45461,
]
outdir = "./resources/gPhoton_out/"
region_fname = (
    "./vasca_pipeline/" + region_name + "/region_" + region_name + "_cat.fits"
)

ron = (6 * uu.arcsec).to(uu.deg).value  # Radius of signal annulus
roff1 = (10 * uu.arcsec).to(uu.deg).value  # Inner radius of background ring
roff2 = (15 * uu.arcsec).to(uu.deg).value  # Outer radius of background ring
bands = ["FUV"]  # "NUV",


rg = Region()
rg.load_from_fits(region_fname)

# Subselect sources based on choice
if len(srcs_ids) > 0:
    rg.tt_sources.add_index("rg_src_id")
    idx_srcs = rg.tt_sources.loc_indices["rg_src_id", srcs_ids]
    tt_srcs = Table(rg.tt_sources[idx_srcs])
else:
    tt_srcs = rg.tt_sources


for band in bands:
    print("--- Running for band", band)
    for rg_src_id in srcs_ids:

        tt_src = tt_srcs[tt_srcs["rg_src_id"] == rg_src_id]

        # Load file
        fname_base = (
            "gPhoton_ra"
            + str(round(tt_src["ra"][0], 5))
            + "_dec"
            + str(round(tt_src["dec"][0], 5))
        )
        outfile_fin = outdir + fname_base + "_" + band.lower() + "_fin.npy"
        outfile_app = outdir + fname_base + "_" + band.lower() + "_app.npy"

        # Run or load gFind
        if os.path.isfile(outfile_fin):
            dd_gfind = np.load(outfile_fin, allow_pickle="TRUE").item()
            t_bins = list(zip(dd_gfind[band]["t0"], dd_gfind[band]["t1"]))
            print("Number of time bins:", len(t_bins))
        else:
            print("Running query with gFind for:\n", tt_src)
            dd_gfind = gPhoton.gFind(
                band=band, skypos=[tt_src["ra"][0], tt_src["dec"][0]]
            )  # ,maxgap=100.,minexp=100.
            t_bins = list(zip(dd_gfind[band]["t0"], dd_gfind[band]["t1"]))
            np.save(outfile_fin, dd_gfind)

        # Run or load gApperture
        if os.path.isfile(outfile_app):
            print("gAperture file already found, not running")
        else:
            print("Running lightcurve with gAperture..")
            dd_gaperture = gPhoton.gAperture(
                band=band,
                skypos=[tt_src["ra"][0], tt_src["dec"][0]],
                radius=ron,
                annulus=[roff1, roff2],
                tranges=t_bins,
            )
            dd_gph = {"gAperture": dd_gaperture, "gFind": dd_gfind}
            np.save(outfile_app, dd_gph)
            print("..done")
