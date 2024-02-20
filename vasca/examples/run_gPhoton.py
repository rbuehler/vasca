#!/usr/bin/env python
# coding: utf-8

import gPhoton
import numpy as np
import os
from astropy import units as uu
from vasca.region import Region
from astropy.table import Table
from vasca.resource_manager import ResourceManager

region_name = "ALL_10-800"  # "TDS"  # "WD" #"MDIS_10-800" # _ELAISN1
t_binning = 5 # Binning in seconds, if -1 "visit binning"
t_name = ""
t_name = "" if t_binning < 0 else "_"+str(t_binning)

# srcs_ids = [
#     193067,
#     432606,
#     535864,
#     451644,
#     1551422,
#     541266,
#     581995,
#     625693,
#     187856,
#     8215,
#     494782,
#     166179,
#     172775,
#     34658,
#     98746,
#     1521738,
#     2136829,
#     297278,
#     426363,
#     426330,
#     151796,
#     305192,
#     259271,
#     388172,
#     265150,
#     54184,
#     472623,
#     419001,
#     25273,
#     26195,
#     32448,
#     199832,
# ]  # WD ALL_10-800

srcs_ids = [357455] # Massive pulsating WD

rm = ResourceManager()
outdir = rm.get_path("gal_gphoton", "sas_cloud")
region_fname = (
    "./vasca_pipeline/" + region_name + "/region_" + region_name + "_cat.fits"
)

ron = (6 * uu.arcsec).to(uu.deg).value  # Radius of signal annulus
roff1 = (10 * uu.arcsec).to(uu.deg).value  # Inner radius of background ring
roff2 = (15 * uu.arcsec).to(uu.deg).value  # Outer radius of background ring
bands = ["NUV", "FUV"]  # "NUV",

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
            + str(round(tt_src["ra"][0], 3))
            + "_dec"
            + str(round(tt_src["dec"][0], 3))
        )
        outfile_fin = outdir + "/" + fname_base + "_" + band.lower() + "_fin.npy"
        outfile_app = outdir + "/" + fname_base + "_" + band.lower() + t_name + "_app.npy"

        # Run or load gFind
        if os.path.isfile(outfile_fin):
            print("Loading file", outfile_fin)
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

        #If binning with fixed time bins is requested, define bin times
        if t_binning >0:
            t_bins_start = []
            t_bins_end = []
            for bin in t_bins:
                bins_fine = np.arange(bin[0], bin[1], t_binning)
                t_bins_start.extend(bins_fine[:-1])
                t_bins_end.extend(bins_fine[1:])
            t_bins = list(zip(t_bins_start,t_bins_end))
            print(t_bins)


        # Run or load gApperture
        if os.path.isfile(outfile_app):
            print("gAperture file already found, not running", outfile_app)
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
