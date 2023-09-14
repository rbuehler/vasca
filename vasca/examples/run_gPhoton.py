#!/usr/bin/env python
# coding: utf-8

import gPhoton
import numpy as np
from astropy.io import ascii
import os
from astropy import units as uu

src_str = """label	ra	dec
GALEX_J221409.8+005245 333.5413675869932 0.8794347888579439
PB_5130 334.61900740472345 -0.0033813554025151345
"""

tt_srcs = ascii.read(src_str)
outdir = "./resources/gPhoton_out/"

ron = (6 * uu.arcsec).to(uu.deg).value  # Radius of signal annulus
roff1 = (10 * uu.arcsec).to(uu.deg).value  # Inner radius of background ring
roff2 = (15 * uu.arcsec).to(uu.deg).value  # Outer radius of background ring

for run_src_idx in range(len(tt_srcs)):
    pos_ra = tt_srcs["ra"][run_src_idx]
    pos_dec = tt_srcs["dec"][run_src_idx]

    outfile_fin = outdir + "gPhoton_" + str(tt_srcs["label"][run_src_idx]) + "_fin.npy"
    outfile_app = outdir + "gPhoton_" + str(tt_srcs["label"][run_src_idx]) + "_app.npy"

    # Run or load gFind
    if os.path.isfile(outfile_fin):
        dd_gfind = np.load(outfile_fin, allow_pickle="TRUE").item()
        t_bins = list(zip(dd_gfind["NUV"]["t0"], dd_gfind["NUV"]["t1"]))
        print("Number of time bins:", len(t_bins))
    else:
        print("Running query with gFind for:\n", tt_srcs[run_src_idx])
        dd_gfind = gPhoton.gFind(
            band="NUV", skypos=[pos_ra, pos_dec]
        )  # ,maxgap=100.,minexp=100.
        t_bins = list(zip(dd_gfind["NUV"]["t0"], dd_gfind["NUV"]["t1"]))
        np.save(outfile_fin, dd_gfind)

    # Run or load gApperture
    if os.path.isfile(outfile_app):
        dd_gph = np.load(outfile_app, allow_pickle="TRUE").item()
    else:
        print("Running lightcurve with gAperture..")
        dd_gaperture = gPhoton.gAperture(
            band="NUV",
            skypos=[pos_ra, pos_dec],
            radius=ron,
            annulus=[roff1, roff2],
            tranges=t_bins,
        )
        dd_gph = {"gAperture": dd_gaperture, "gFind": dd_gfind}
        np.save(outfile_app, dd_gph)
        print("..done")
