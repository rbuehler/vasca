#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Defines dictionary for the tables used by vasca.tables.TableCollection
"""
# %% dtype guidelines
#
# >np.finfo(np.float16)
# >finfo(resolution=0.001, min=-6.55040e+04, max=6.55040e+04, dtype=float16)
#
# >np.finfo(np.float32)
# >finfo(resolution=1e-06, min=-3.4028235e+38, max=3.4028235e+38, dtype=float32)
#
# >np.finfo(np.float64)
# >ffinfo(resolution=1e-15, min=-1.7976931348623157e+308, max=1.7976931348623157e+308, dtype=float64)
#
# >np.iinfo(np.int32)
# >info(min=-2147483648, max=2147483647, dtype=int32)
#
# >np.iinfo(np.int64)
# >iinfo(min=-9223372036854775808, max=9223372036854775807, dtype=int64)
#
# Based on the above, we need float64 for MJD, ra, dec and int64 for external ID numbers
# For everything else float32 and int32 should be sufficient
# (e.g. for all flux variables and errors)

# import numpy as np

# global dictionaries defining the table structures


# %% Column definitions

dd_vasca_columns = {
    # %%% field_id
    "field_id": {
        "name": "field_id",
        "dtype": "S32",
        "unit": "1",
        "default": "none",
        "description": "Field source ID number",
    },
    "src_name": {
        "name": "src_name",
        "dtype": "S24",
        "unit": "1",
        "default": "none",
        "description": "VASCA catalog source name",
    },
    # %%% field_name
    "field_name": {
        "name": "field_name",
        "dtype": "S32",
        "unit": "",
        "default": "none",
        "description": "Field name",
    },
    # %%% project
    "project": {
        "name": "project",
        "dtype": "S32",
        "unit": "",
        "default": "none",
        "description": "Field project, typically survey name",
    },
    # %%% ra
    "ra": {
        "name": "ra",
        "dtype": "float64",
        "unit": "degree",
        "default": -1.0,
        "description": "Sky coordinate Right Ascension (J2000)",
    },
    # %%% dec
    "dec": {
        "name": "dec",
        "dtype": "float64",
        "unit": "degree",
        "default": -1.0,
        "description": "Sky coordinate Declination (J2000)",
    },
    # %%% observatory
    "observatory": {
        "name": "observatory",
        "dtype": "S22",
        "unit": "",
        "default": "none",
        "description": "Telescope of the observation (e.g. GALEX)",
    },
    # %%% obs_filter
    "obs_filter": {
        "name": "obs_filter",
        "dtype": "S8",
        "unit": "",
        "default": "none",
        "description": "Filter of the observation (e.g. NUV)",
    },
    # %%% obs_filter_idx
    "obs_filter_idx": {
        "name": "obs_filter_idx",
        "dtype": "int32",
        "unit": "1",
        "default": -1,
        "description": "Filter index in filter dependent arrays",
    },
    # %%% obs_filter_id
    "obs_filter_id": {
        "name": "obs_filter_id",
        "dtype": "int32",
        "unit": "1",
        "default": 0,
        "description": "Observation filter ID number",
    },
    # %%% wavelength
    "wavelength": {
        "name": "wavelength",
        "dtype": "float32",
        "unit": "AA",
        "default": -1.0,
        "description": "Effective wavelength of the observation filter",
    },
    # %%% fov_diam
    "fov_diam": {
        "name": "fov_diam",
        "dtype": "float32",
        "unit": "degree",
        "default": -1.0,
        "description": "Field radius or box size (depending on the observatory)",
    },
    # %%% vis_id
    "vis_id": {
        "name": "vis_id",
        "dtype": "uint64",
        "unit": "1",
        "default": -1,
        "description": "Visit ID number",
    },
    # %%% time_bin_start
    "time_bin_start": {
        "name": "time_bin_start",
        "dtype": "float64",
        "unit": "d",
        "default": -1.0,
        "description": "Visit exposure start date and time in MJD",
    },
    # %%% time_bin_size
    "time_bin_size": {
        "name": "time_bin_size",
        "dtype": "float32",
        "unit": "s",
        "default": -1.0,
        "description": "Visit exposure time in s",
    },
    # %%% fd_src_id
    "fd_src_id": {
        "name": "fd_src_id",
        "dtype": "int32",
        "unit": "1",
        "default": -1,
        "description": "Source ID associated to the visit detection",
    },
    # %%% pos_err
    "pos_err": {
        "name": "pos_err",
        "dtype": "float32",
        "unit": "arcsec",
        "default": -1.0,
        "description": "Sky coordinate position error",
    },
    # %%% flux
    "flux": {
        "name": "flux",
        "dtype": "float32",
        "unit": "1e-6Jy",
        "default": -1.0,
        "description": "Flux density",
    },
    "flux_model": {
        "name": "flux",
        "dtype": "float32",
        "unit": "1e-6Jy",
        "default": -1.0,
        "description": "Flux from SDSS model spectrum.",
    },
    # %%% flux_app_ratio
    "flux_app_ratio": {
        "name": "flux_app_ratio",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Flux ratio measured at a different apertures (3.8 arcsec / 6 arcsec radius for GALEX)",
    },
    # %%% flux_err
    "flux_err": {
        "name": "flux_err",
        "dtype": "float32",
        "unit": "1e-6Jy",
        "default": -1.0,
        "description": "Flux density error",
    },
    # %%% s2n
    "s2n": {
        "name": "s2n",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Flux signal to noise",
    },
    # %%% det_id
    "det_id": {
        "name": "det_id",
        "dtype": "int64",
        "unit": "1",
        "default": -1,
        "description": "Reference source ID number",
    },
    # %%% mag
    "mag": {
        "name": "mag",
        "dtype": "float32",
        "unit": "mag",
        "default": -1.0,
        "description": "AB magnitude",
    },
    # %%% mag_err
    "mag_err": {
        "name": "mag_err",
        "dtype": "float32",
        "unit": "mag",
        "default": -1.0,
        "description": "AB magnitude error",
    },
    # %%% nr_det
    "nr_det": {
        "name": "nr_det",
        "dtype": "int32",
        "unit": "1",
        "default": -1,
        "description": "Number of detections",
    },
    # %%% pos_xv
    "pos_xv": {
        "name": "pos_xv",
        "dtype": "float32",
        "unit": "arcsec2",
        "default": -1.0,
        "description": "Sky position excess variance",
    },
    # %%% pos_var
    "pos_var": {
        "name": "pos_var",
        "dtype": "float32",
        "unit": "arcsec2",
        "default": -1.0,
        "description": "Sky position variance",
    },
    # %%% pos_cpval
    "pos_cpval": {
        "name": "pos_cpval",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Sky position quality",
    },
    # %%% pos_rchiq
    "pos_rchiq": {
        "name": "pos_rchiq",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Sky position reduced chisquared of the constant mean",
    },
    # %%% coadd_id
    "coadd_id": {
        "name": "coadd_id",
        "dtype": "int64",
        "unit": "1",
        "default": -1,
        "description": "Associated co-add source or detection ID",
    },
    # %%% coadd_dist
    "coadd_dist": {
        "name": "coadd_dist",
        "dtype": "float32",
        "unit": "arcsec",
        "default": -1.0,
        "description": "Angular distance to associated source",
    },
    "cat_dist": {
        "name": "cat_dist",
        "dtype": "float32",
        "unit": "arcsec",
        "default": -1.0,
        "description": "Angular distance to associated catalog source",
    },
    # %%% match_distance
    "match_distance": {
        "name": "match_distance",
        "dtype": "float32",
        "unit": "arcsec",
        "default": -1.0,
        "description": "Angular distance to associated source",
    },
    # %%% sel
    "sel": {
        "name": "sel",
        "dtype": "bool",
        "unit": "1",
        "default": True,
        "description": "Selection of rows for VASCA analysis",
    },
    # %%% flux_nxv
    "flux_nxv": {
        "name": "flux_nxv",
        "dtype": "float32",
        "unit": "1",
        "default": -100.0,
        "description": "Flux normalised excess variance",
    },
    # %%% flux_ne
    "flux_ne": {
        "name": "flux_ne",
        "dtype": "float32",
        "unit": "1",
        "default": -100.0,
        "description": "Flux square root of normalised excess variance",
    },
    # %%% flux_var
    "flux_var": {
        "name": "flux_var",
        "dtype": "float32",
        "unit": "1e-12Jy2",
        "default": -1.0,
        "description": "Flux variance",
    },
    # %%% flux_cpval
    "flux_cpval": {
        "name": "flux_cpval",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Probability value for a constant flux from the chisquare test",
    },
    # %%% flux_rchiq
    "flux_rchiq": {
        "name": "flux_rchiq",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Flux reduced chisquared of the constant mean",
    },
    # %%% coadd_ffactor
    "coadd_ffactor": {
        "name": "coadd_ffactor",
        "dtype": "float32",
        "unit": "1",
        "default": -100.0,
        "description": "Source flux divided by flux of the associated co-add source",
    },
    # %%% time_start
    "time_start": {
        "name": "time_start",
        "dtype": "float64",
        "unit": "d",
        "default": -1.0,
        "description": "Start date and time in MJD",
    },
    # %%% time
    "time": {
        "name": "time",
        "dtype": "float64",
        "unit": "d",
        "default": -1.0,
        "description": "Time in MJD",
    },
    # %%% ul
    "ul": {
        "name": "ul",
        "dtype": "float32",
        "unit": "1e-6Jy",
        "default": -1.0,
        "description": "Flux density upper limit",
    },
    # %%% time_bin_size_alt_filt
    "time_bin_size_alt_filt": {
        "name": "time_bin_size_alt_filt",
        "dtype": "float64",
        "unit": "s",
        "default": -1.0,
        "description": "Visit exposure time of the alternative filter in s",
    },
    # %%% r_fov
    "r_fov": {
        "name": "r_fov",
        "dtype": "float32",
        "unit": "degree",
        "default": -1.0,
        "description": "Distance from center of FOV in degrees",
    },
    # %%% artifacts
    "artifacts": {
        "name": "artifacts",
        "dtype": "int64",
        "unit": "1",
        "default": -1,
        "description": "Logical OR of artifact flags",
    },
    # %%% flags
    "flags": {
        "name": "flags",
        "dtype": "int64",
        "unit": "1",
        "default": -1,
        "description": "Logical OR of artifact flags from gphoton",
    },
    # %%% class_star
    "class_star": {
        "name": "class_star",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Point-source probability: 0.0 (resolved), 1.0 (unresolved, mcat file filter_CLASS_STAR variable)",
    },
    # %%% chkobj_type
    "chkobj_type": {
        "name": "chkobj_type",
        "dtype": "int32",
        "unit": "1",
        "default": -1,
        "description": "Detection matched to a known star (bright_match=1, mcat file chkobj_type variable)",
    },
    # %%% size_world
    "size_world": {
        "name": "size_world",
        "dtype": "float32",
        "unit": "arcsec",
        "default": -1,
        "description": "Mean RMS of the ellipse size: (major axis + minor axis) / 2",
    },
    # %%% ellip_world
    "ellip_world": {
        "name": "ellip_world",
        "dtype": "float32",
        "unit": "1",
        "default": -1,
        "description": "Ellipticity of the detection or source",
    },
    # %%% flux_auto
    "flux_auto": {
        "name": "flux_auto",
        "dtype": "float32",
        "unit": "Jy",
        "default": -1.0,
        "description": "Flux from SExtractor auto option",
    },
    # %%% flux_auto_err
    "flux_auto_err": {
        "name": "flux_auto_err",
        "dtype": "float32",
        "unit": "Jy",
        "default": -1.0,
        "description": "Flux error from a fixed circular aperture (3.8 arcsec radius for GALEX)",
    },
    # %%% E_bv
    "E_bv": {
        "name": "E_bv",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Galactic reddening expressed as E(B-V)",
    },
    # %%% nr_vis
    "nr_vis": {
        "name": "nr_vis",
        "dtype": "int32",
        "unit": "1",
        "default": -1,
        "description": "Total number of visits of the field",
    },
    # %%% time_bin_size_sum
    "time_bin_size_sum": {
        "name": "time_bin_size_sum",
        "dtype": "float32",
        "unit": "s",
        "default": -1.0,
        "description": "Total exposure time",
    },
    # %%% time_stop
    "time_stop": {
        "name": "time_stop",
        "dtype": "float64",
        "unit": "d",
        "default": -1.0,
        "description": "End time of last exposure",
    },
    # %%% pix_id
    "pix_id": {
        "name": "pix_id",
        "dtype": "uint32",
        "unit": "1",
        "default": 0,
        "description": "Healpix ID",
    },
    # %%% exp
    "exp": {
        "name": "exp",
        "dtype": "float32",
        "unit": "1",
        "default": -1,
        "description": "Total exposure",
    },
    # %%% coadd_fdiff_s2n
    "coadd_fdiff_s2n": {
        "name": "coadd_fdiff_s2n",
        "dtype": "float32",
        "unit": "1",
        "default": -10000.0,
        "description": "Signal to noise of the flux difference",
    },
    # %%% rg_src_id
    "rg_src_id": {
        "name": "rg_src_id",
        "dtype": "int32",
        "unit": "1",
        "default": -1,
        "description": "Region source ID number",
    },
    # %%% rg_fd_id
    "rg_fd_id": {
        "name": "rg_fd_id",
        "dtype": "int64",
        "unit": "1",
        "default": -1,
        "description": "Region field ID number",
    },
    # %%% nr_fds
    "nr_fds": {
        "name": "nr_fds",
        "dtype": "int32",
        "unit": "1",
        "default": -1,
        "description": "Number of fields",
    },
    # %%% nr_fd_srcs
    "nr_fd_srcs": {
        "name": "nr_fd_srcs",
        "dtype": "int32",
        "unit": "1",
        "default": -1,
        "description": "Number of field sources",
    },
    # %%% nr_fd_dets
    "nr_fd_dets": {
        "name": "nr_fd_dets",
        "dtype": "int32",
        "unit": "1",
        "default": -1,
        "description": "Number of field detections",
    },
    # %%% coadd_src_id
    "coadd_src_id": {
        "name": "coadd_src_id",
        "dtype": "int64",
        "unit": "1",
        "default": -1,
        "description": "Co-add source ID number",
    },
    "cat_src_id": {
        "name": "cat_src_id",
        "dtype": "int64",
        "unit": "1",
        "default": -1,
        "description": "Catalog source ID number",
    },
    # %%% otype
    "otype": {
        "name": "otype",
        "dtype": "S32",
        "unit": "1",
        "default": "none",
        "description": "SIMBAD source type",
    },
    # %%% ogrp
    "ogrp": {
        "name": "ogrp",
        "dtype": "S8",
        "unit": "1",
        "default": "none",
        "description": "SIMBAD source type group in VASCA",
    },
    # %%% hr
    "hr": {
        "name": "hr",
        "dtype": "float32",
        "unit": "1",
        "default": -1,
        "description": "Flux hardness ratio, only simultaneous detections considered",
    },
    # %%% hr_err
    "hr_err": {
        "name": "hr_err",
        "dtype": "float32",
        "unit": "1",
        "default": -1,
        "description": "Flux hardness ratio error",
    },
    # %%% match_id
    "match_id": {
        "name": "match_id",
        "dtype": "int32",
        "unit": "1",
        "default": -1,
        "description": "VASCA internal source ID number for associated sources",
    },
    "sp_type": {
        "name": "sp_type",
        "dtype": "S32",
        "unit": "",
        "default": "none",
        "description": "SIMBAD spectral type",
    },
    "main_id": {
        "name": "main_id",
        "dtype": "S32",
        "unit": "",
        "default": "none",
        "description": "SIMBAD main ID",
    },
    "Source": {
        "name": "Source",
        "dtype": "S32",
        "unit": "",
        "default": "none",
        "description": "GAIA DR3 source ID",
    },
    "gfcat_objid": {
        "name": "gfcat_objid",
        "dtype": "S32",
        "unit": "",
        "default": "none",
        "description": "GFCAT object ID",
    },
    "WDJname": {
        "name": "WDJname",
        "dtype": "S32",
        "unit": "",
        "default": "none",
        "description": "GAIA-EDR3-WD object name",
    },
    # %%% origin
    "origin": {
        "name": "origin",
        "dtype": "S22",
        "unit": "",
        "default": "none",
        "description": "Origin of the data",
    },
    "PQSO": {
        "name": "PQSO",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "Probability to be a Quasar from GAIA-DR3",
    },
    "PGal": {
        "name": "PGal",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "Probability to be a Galaxy from GAIA-DR3",
    },
    "PSS": {
        "name": "PSS",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "Probability to be a Single (non-WD) star from GAIA-DR3",
    },
    "RPlx": {
        "name": "RPlx",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "Parallax divided by its standard error from GAIA-DR3",
    },
    "RFRP": {
        "name": "RFRP",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "Red magnitude by its standard error from GAIA-DR3",
    },
    "RFBP": {
        "name": "RFBP",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "Blue magnitude by its standard error from GAIA-DR3",
    },
    "AG": {
        "name": "AG",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "G magnitude absorption from gsp-phot for GAIA-DR3",
    },
    "E_BP-RP": {
        "name": "E_BP-RP",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "Reddening from gsp-phot for GAIA-DR3",
    },
    "VarFlag": {
        "name": "VarFlag",
        "dtype": "S32",
        "unit": "1",
        "default": "none",
        "description": "Vizier GAIA-DR3 VarFlag",
    },
    "Plx": {
        "name": "Plx",
        "dtype": "float32",
        "unit": "1e-3 arcsec",
        "default": -1,
        "description": "Parallax from GAIA-DR3",
    },
    "e_Plx": {
        "name": "e_Plx",
        "dtype": "float32",
        "unit": "1e-3 arcsec",
        "default": -1,
        "description": "Parallax error from GAIA-DR3",
    },
    "Plx_dist": {
        "name": "Plx_dist",
        "dtype": "float32",
        "unit": "pc",
        "default": -1,
        "description": "Parallax distance from GAIA-DR3",
    },
    "Gmag": {
        "name": "Gmag",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "G-band mean magnitude from GAIA-DR3",
    },
    "Gmag_abs": {
        "name": "Gmag_abs",
        "dtype": "float32",
        "unit": "",
        "default": -100,
        "description": "G-band mean absolute magnitude calculated from GAIA-DR3",
    },
    "BP-RP": {
        "name": "BP-RP",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "BP-RP colour from GAIA-DR3",
    },
    "Pwd": {
        "name": "Pwd",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "The probability of being a white dwarf in GAIA-EDR3 WD catalog",
    },
    "ls_peak_power": {
        "name": "ls_peak_power",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "LombScargle power at peak frequency",
    },
    "ls_peak_freq": {
        "name": "ls_peak_freq",
        "dtype": "float32",
        "unit": "1/d",
        "default": -1,
        "description": "LombScargle peak frequency",
    },
    "ls_peak_pval": {
        "name": "ls_peak_pval",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "LombScargle power probability value",
    },
    "ls_pval_alt_flt": {
        "name": "ls_pval_alt_flt",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "LombScargle power probability value for alternate filter",
    },
    "ls_model_rchiq": {
        "name": "ls_model_rchiq",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "Reduced chisquare of the LombScargle peak frequency sine wave",
    },
    "ls_model_pval": {
        "name": "ls_model_pval",
        "dtype": "float32",
        "unit": "",
        "default": -1,
        "description": "Chisquare probability of the LombScargle peak frequency model ",
    },
    "maincat_match_id": {
        "name": "maincat_match_id",
        "dtype": "int32",
        "unit": "",
        "default": -1,
        "description": "Internal source ID from associated source in the GAIA-EDR3 White Dwarf catalog",
    },
    "simbad_match_id": {
        "name": "simbad_match_id",
        "dtype": "int32",
        "unit": "",
        "default": -1,
        "description": "Internal source ID from associated source in the SIMBAD data base",
    },
    "gaiadr3_match_id": {
        "name": "gaiadr3_match_id",
        "dtype": "int32",
        "unit": "",
        "default": -1,
        "description": "Internal source ID from associated source in the GAIA-DR3 catalog",
    },
    "gfcat_src_id": {
        "name": "gfcat_match_id",
        "dtype": "int32",
        "unit": "",
        "default": -1,
        "description": "Internal source ID from associated source in the GFCAT catalog",
    },
}

# %% Table definitions
# %%% base_field
base_field = {
    "tt_fields": {
        "names": [
            "field_id",
            "field_name",
            "project",
            "ra",
            "dec",
            "observatory",
            "obs_filter",
            "fov_diam",
            "sel",
        ],
        "meta": {"DATAPATH": "None", "INFO": "Field information table"},
    },
    "tt_visits": {
        "names": ["vis_id", "time_bin_start", "time_bin_size", "sel", "obs_filter_id"],
        "meta": {"INFO": "Visit information table"},
    },
    "tt_detections": {
        "names": [
            "vis_id",
            "fd_src_id",
            "ra",
            "dec",
            "pos_err",
            "flux",
            "flux_err",
            "flux_app_ratio",
            "s2n",
            "obs_filter_id",
            "sel",
            "r_fov",
            "artifacts",
            "class_star",
            "chkobj_type",
            "size_world",
            "ellip_world",
            "flux_auto",
            "flux_auto_err",
        ],
        "meta": {"INFO": "Visit detections table"},
    },
    "tt_coadd_detections": {
        "names": [
            "det_id",
            "ra",
            "dec",
            "pos_err",
            "flux",
            "flux_err",
            "flux_app_ratio",
            "s2n",
            "obs_filter_id",
            "sel",
        ],
        "meta": {"INFO": "Reference detections table"},
    },
    "tt_sources": {
        "names": [
            "fd_src_id",
            "nr_det",
            "ra",
            "dec",
            "pos_err",
            "pos_xv",
            "pos_var",
            "pos_cpval",
            "pos_rchiq",
            "coadd_src_id",
            "coadd_dist",
            "obs_filter_id",
            "sel",
            "flux",
            "flux_err",
            "flux_nxv",
            "flux_var",
            "flux_cpval",
            "flux_rchiq",
            "coadd_ffactor",
            "coadd_fdiff_s2n",
        ],
        "meta": {"INFO": "Source infomation table", "CLUSTALG": "None"},
    },
    "tt_source_lc": {
        "names": [
            "time",
            "time_bin_size",
            "flux",
            "flux_err",
            "ul",
            "sel",
            "obs_filter",
            "obs_filter_id",
            "r_fov",
            "artifacts",
            "class_star",
            "chkobj_type",
            "size_world",
            "ellip_world",
            "flux_auto",
            "flux_auto_err",
            "vis_id",
        ],
        "meta": {"INFO": "Light curve flux table for one source"},
    },
}
# %%% galex_field
galex_field = {
    "tt_visits": {
        "names": [
            *base_field["tt_visits"]["names"],
            "ra",
            "dec",
        ],
        "meta": {"INFO": "Visit information table"},
    },
    "tt_detections": {
        "names": [
            *base_field["tt_detections"]["names"],
            "det_id",
            "E_bv",
        ],
        "meta": {"INFO": "Visit detections table", "PRECUTS": "List of pre-cuts"},
    },
    "tt_coadd_detections": {
        "names": [
            *base_field["tt_coadd_detections"]["names"],
            "r_fov",
            "artifacts",
            "class_star",
            "chkobj_type",
            "size_world",
            "ellip_world",
            "flux_auto",
            "flux_auto_err",
            "E_bv",
        ],
        "meta": {"INFO": "Reference detections table", "PRECUTS": "List of pre-cuts"},
    },
}
# %%% region
region = {
    "tt_fields": {
        "names": [
            *base_field["tt_fields"]["names"],
            "nr_vis",
            "time_bin_size_sum",
            "time_start",
            "time_stop",
            "rg_fd_id",
        ],
        "meta": {"DATAPATH": "None", "INFO": "Field information table"},
    },
    "tt_visits": {
        "names": [*base_field["tt_visits"]["names"]],
        "meta": {"INFO": "Visit information table"},
    },
    "tt_coverage_hp": {
        "names": ["pix_id", "nr_vis", "exp", "nr_fds"],
        "meta": {
            "DATAPATH": "None",
            "INFO": "Region observations properties in healpix binning. RING ordering and equatorial coordinates",
            "NSIDE": "None",
        },
    },
    "tt_coadd_detections": {
        "names": [
            *base_field["tt_coadd_detections"]["names"],
            "rg_fd_id",
            "coadd_src_id",
        ],
        "meta": {"INFO": "Reference detections table"},
    },
    "tt_detections": {
        "names": [*base_field["tt_detections"]["names"], "rg_fd_id", "rg_src_id"],
        "meta": {"INFO": "Visit detections table"},
    },
    "tt_sources": {
        "names": [
            *base_field["tt_sources"]["names"],
            "rg_fd_id",
            "rg_src_id",
            "nr_fd_srcs",
        ],
        "meta": {"INFO": "Source infomation table", "CLUSTALG": "None"},
    },
    "tt_coadd_sources": {
        "names": [
            *base_field["tt_sources"]["names"],
            "rg_fd_id",
            "nr_fd_dets",
        ],
        "meta": {"INFO": "Source infomation table", "CLUSTALG": "None"},
    },
    "tt_src_id_map": {
        "names": ["rg_src_id", "rg_fd_id", "fd_src_id", "sel"],
        "meta": {"INFO": "Map between region and field source IDs"},
    },
    "tt_filters": {
        "names": ["obs_filter_id", "obs_filter", "obs_filter_idx"],
        "meta": {
            "INFO": "Filters, their IDs and index, the last is specific for this region."
        },
    },
    "tt_sed": {
        "names": [
            "flux",
            "flux_err",
            "observatory",
            "obs_filter",
            "origin",
            "wavelength",
        ],
        "meta": {"INFO": "Spectral Energy Distribution from VizieR database"},
    },
    "tt_gphoton_lc": {
        "names": [
            "time",
            "time_bin_size",
            "flux",
            "flux_err",
            "sel",
            "obs_filter",
            "obs_filter_id",
            "flags",
            "s2n",
        ],
        "meta": {"INFO": "Light curve from gPhoton.gAperture"},
    },
    "tt_spectrum": {
        "names": ["flux", "wavelength", "s2n", "flux_model", "sel"],
        "meta": {"INFO": "Spectrum"},
    },
    "tt_lombscargle": {
        "names": [
            "rg_src_id",
            "ls_peak_power",
            "ls_peak_freq",
            "ls_peak_pval",
            "ls_pval_alt_flt",
            "ls_model_rchiq",
            "ls_model_pval",
        ],
        "meta": {"INFO": "LombScargle results information"},
    },
}


# global, combined dictionary
class_keys = ["base_field", "galex_field", "region"]
class_dicts = [base_field, galex_field, region]
dd_vasca_tables = {c_key: c_dict for c_key, c_dict in zip(class_keys, class_dicts)}

# Add columns to tables dictionary
for tab_type, table_group in dd_vasca_tables.items():
    for tab_name, tab in table_group.items():
        tab["dtype"] = []
        tab["units"] = []
        tab["defaults"] = []
        tab["descriptions"] = []
        for col_name in tab["names"]:
            col = dd_vasca_columns[col_name]
            tab["dtype"].append(col["dtype"])
            tab["units"].append(col["unit"])
            tab["defaults"].append(col["default"])
            tab["descriptions"].append(col["description"])
