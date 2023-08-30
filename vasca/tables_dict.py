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
# For enverything else float32 and int32 should be sufficient
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
        "description": "Field source ID nr",
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
        "description": "Field project, typicaly survey name",
    },
    # %%% ra
    "ra": {
        "name": "ra",
        "dtype": "float64",
        "unit": "degree",
        "default": -1.0,
        "description": "Center RA of the field (J2000)",
    },
    # %%% dec
    "dec": {
        "name": "dec",
        "dtype": "float64",
        "unit": "degree",
        "default": -1.0,
        "description": "Center Dec of the field (J2000)",
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
        "description": "Filter ID number",
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
        "description": "Visit ID nr",
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
        "description": "Visit position error",
    },
    # %%% flux
    "flux": {
        "name": "flux",
        "dtype": "float32",
        "unit": "1e-6Jy",
        "default": -1.0,
        "description": "Visit detection flux density",
    },
    # %%% flux_err
    "flux_err": {
        "name": "flux_err",
        "dtype": "float32",
        "unit": "1e-6Jy",
        "default": -1.0,
        "description": "Visit detection flux density error",
    },
    # %%% s2n
    "s2n": {
        "name": "s2n",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Signal to noise",
    },
    # %%% det_id
    "det_id": {
        "name": "det_id",
        "dtype": "int64",
        "unit": "1",
        "default": -1,
        "description": "Reference source ID nr",
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
        "description": "Position excess variance, entries for different filters",
    },
    # %%% pos_var
    "pos_var": {
        "name": "pos_var",
        "dtype": "float32",
        "unit": "arcsec2",
        "default": -1.0,
        "description": "Position variance of detections, entries for different filters",
    },
    # %%% pos_cpval
    "pos_cpval": {
        "name": "pos_cpval",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Position probability value for a constant from the chisquare test",
    },
    # %%% pos_rchiq
    "pos_rchiq": {
        "name": "pos_rchiq",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Position reduced chisquared of the constant mean",
    },
    # %%% assoc_id
    "assoc_id": {
        "name": "assoc_id",
        "dtype": "int64",
        "unit": "1",
        "default": -1,
        "description": "Associated source or detection ID",
    },
    # %%% assoc_dist
    "assoc_dist": {
        "name": "assoc_dist",
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
        "description": "Flux density normalized excess variance, entries for different filters",
    },
    # %%% flux_var
    "flux_var": {
        "name": "flux_var",
        "dtype": "float32",
        "unit": "1e-12Jy2",
        "default": -1.0,
        "description": "Flux density variance of detections, entries for different filters",
    },
    # %%% flux_cpval
    "flux_cpval": {
        "name": "flux_cpval",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Flux probability value for a constant from the chisquare test, entries for different filters",
    },
    # %%% flux_rchiq
    "flux_rchiq": {
        "name": "flux_rchiq",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Flux reduced chisquared of the constant mean, entries for different filters",
    },
    # %%% assoc_ffactor
    "assoc_ffactor": {
        "name": "assoc_ffactor",
        "dtype": "float32",
        "unit": "1",
        "default": -100.0,
        "description": "Source flux divided by flux of the associated source, entries for different filters",
    },
    # %%% time_start
    "time_start": {
        "name": "time_start",
        "dtype": "float64",
        "unit": "d",
        "default": -1.0,
        "description": "Visit exposure start date and time in MJD",
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
    # %%% flux_f60
    "flux_f60": {
        "name": "flux_f60",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Flux in a fixed circular 6.0 arcsec radius aperture in cts/sec",
    },
    # %%% flux_f60_err
    "flux_f60_err": {
        "name": "flux_f60_err",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Flux error in a fixed circular 6.0 arcsec radius aperture in cts/sec",
    },
    # %%% flux_f38
    "flux_f38": {
        "name": "flux_f38",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Flux in a fixed circular 3.8 arcsec radius aperture in cts/sec",
    },
    # %%% flux_f38_err
    "flux_f38_err": {
        "name": "flux_f38_err",
        "dtype": "float32",
        "unit": "1",
        "default": -1.0,
        "description": "Flux error in a fixed circular 3.8 arcsec radius aperture in cts/sec",
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
    # %%% assoc_fdiff_s2n
    "assoc_fdiff_s2n": {
        "name": "assoc_fdiff_s2n",
        "dtype": "float32",
        "unit": "1",
        "default": -10000.0,
        "description": "Signal to noise of the flux difference, entries for different filters",
    },
    # %%% rg_src_id
    "rg_src_id": {
        "name": "rg_src_id",
        "dtype": "int32",
        "unit": "1",
        "default": -1,
        "description": "Region source ID nr",
    },
    # %%% rg_fd_id
    "rg_fd_id": {
        "name": "rg_fd_id",
        "dtype": "int64",
        "unit": "1",
        "default": -1,
        "description": "Region field ID nr",
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
        "description": "Region coadd source ID nr",
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
        "description": "SIMBAD source type group, defined in vasca.utils",
    },
    # %%% hr
    "hr": {
        "name": "hr",
        "dtype": "float32",
        "unit": "1",
        "default": -1,
        "description": "Flux hardness ratio, only strictly simulateneous detections considered",
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
            "s2n",
            "obs_filter_id",
            "sel",
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
            "assoc_id",
            "assoc_dist",
            "obs_filter_id",
            "sel",
            "flux",
            "flux_err",
            "flux_nxv",
            "flux_var",
            "flux_cpval",
            "flux_rchiq",
            "assoc_ffactor",
            "assoc_fdiff_s2n",
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
            "r_fov",
            "artifacts",
            "class_star",
            "chkobj_type",
            "size_world",
            "ellip_world",
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
            "flux_f60",
            "flux_f60_err",
            "flux_f38",
            "flux_f38_err",
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
            "coadd_src_id",
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
