#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Defines dictionary for the tables used by vasca.tables.TableCollection
"""
# **dtype guidelines**
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
# (e.g. for all magnitude/flux variables and errors)

import numpy as np

# global dictionaries defining the table structures

# %% Base field
base_field = {
    "tt_fields": {
        "names": [
            "field_id",
            "name",
            "ra",
            "dec",
            "observatory",
            "obs_filter",
            "fov_diam",
            "sel",
        ],
        "dtype": ["S22", "S22", "float64", "float64", "S22", "S22", "float32", "bool"],
        "units": ["1", "", "degree", "degree", "", "", "degree", "1"],
        "defaults": [-1, "none", -1.0, -1.0, "none", "none", -1.0, True],
        "descriptions": [
            "Field source ID nr.",
            "Field name",
            "Center RA of the field (J2000)",
            "Center Dec of the field (J2000)",
            "Telescope of the observation (e.g. GALEX)",
            "Filter of the observation (e.g. NUV)",
            "Field radius or box size (depending on the observatory)",
            "Selection of rows for VASCA analysis.",
        ],
        "meta": {"DATAPATH": "None", "INFO": "Field information table"},
    },
    "tt_visits": {
        "names": ["vis_id", "time_bin_start", "time_bin_size", "sel"],
        "dtype": ["int64", "float64", "float32", "bool"],
        "units": ["1", "d", "s", "1"],
        "defaults": [-1, -1.0, -1.0, True],
        "descriptions": [
            "Visit ID nr.",
            "Visit exposure start date and time in MJD",
            "Visit exposure time in s",
            "Selection of rows for VASCA analysis.",
        ],
        "meta": {"INFO": "Visit information table"},
    },
    "tt_detections": {
        "names": [
            "vis_id",
            "fd_src_id",
            "ra",
            "dec",
            "pos_err",
            "mag",
            "mag_err",
            "s2n",
            "sel",
        ],
        "dtype": [
            "int64",
            "int32",
            "float64",
            "float64",
            "float32",
            "float32",
            "float32",
            "float32",
            "bool",
        ],
        "units": [
            "1",
            "1",
            "degree",
            "degree",
            "degree",
            "mag",
            "mag",
            "1",
            "1",
        ],
        "defaults": [-1, -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, True],
        "descriptions": [
            "Visit ID associated to the visit detection",
            "Source ID associated to the visit detection",
            "Visit detection RA (J2000)",
            "Visit detection Dec (J2000)",
            "Visit position error",
            "Visit detection AB magnitude",
            "Visit detection magnitude error",
            "Signal to noise",
            "Selection of rows for VASCA analysis.",
        ],
        "meta": {"INFO": "Visit detections table"},
    },
    "tt_ref_sources": {
        "names": ["det_id", "ra", "dec", "pos_err", "mag", "mag_err", "s2n", "sel"],
        "dtype": [
            "int64",
            "float64",
            "float64",
            "float32",
            "float32",
            "float32",
            "float32",
            "bool",
        ],
        "units": ["1", "degree", "degree", "degree", "mag", "mag", "1", "1"],
        "defaults": [-1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, True],
        "descriptions": [
            "Reference source ID nr.",
            "Reference source RA (J2000)",
            "Reference source Dec (J2000)",
            "Reference position error",
            "Reference source AB magnitude",
            "Reference source magnitude error",
            "Signal to noise",
            "Selection of rows for VASCA analysis.",
        ],
        "meta": {"INFO": "Reference detections table"},
    },
    "tt_sources": {
        "names": [
            "fd_src_id",
            "nr_det",
            "nr_uls",
            "ra",
            "dec",
            "pos_err",
            "pos_var_ex",
            "pos_var",
            "mag",
            "mag_err",
            "mag_var_ex",
            "mag_var",
            "flux_cpval",
            "flux_rchiq",
            "mag_dmagmax_abs",
            "mag_sigmax",
            "mag_sigmax_dmag",
            "mag_skew",
            "sel",
        ],
        "dtype": [
            "int32",
            "int32",
            "int32",
            "float64",
            "float64",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "bool",
        ],
        "units": [
            "1",
            "1",
            "1",
            "degree",
            "degree",
            "degree",
            "degree2",
            "degree2",
            "mag",
            "mag",
            "mag2",
            "mag2",
            "1",
            "1",
            "mag",
            "1",
            "mag",
            "mag2",
            "1",
        ],
        "defaults": [
            -1,
            -1,
            -1,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -100.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -100.0,
            True,
        ],
        "descriptions": [
            "Source ID nr.",
            "Number of visit detections of the source",
            "Number of upper limits",
            "Source RA (J2000)",
            "Source Dec (J2000)",
            "Position error",
            "Position excess variance",
            "Position variance of detections",
            "Magnitude for a constant flux",
            "Magnitude error",
            "AB magnitude excess variance",
            "AB magnitude variance of detections",
            "Flux probability value for a constant from the chisquare test",
            "Flux reduced chisquared of the constant mean",
            "Maximum magnitude of magnitude variation in detections",
            "Maxmimum significance of magnitude variations",
            "Magnitude variation for maximum significance variation",
            "Magnitude skewness",
            "Selection of rows for VASCA analysis.",
        ],
        "meta": {"INFO": "Source infomation table", "CLUSTALG": "None"},
    },
    "tt_source_lc": {
        "names": ["time_start", "time_delta", "mag", "mag_err", "ul", "sel"],
        "dtype": [
            "float64",
            "float32",
            "float32",
            "float32",
            "float32",
            "bool",
        ],
        "units": ["d", "s", "mag", "mag", "mag", "1"],
        "defaults": [-1.0, -1.0, -1, -1, -1.0, True],
        "descriptions": [
            "Visit exposure start date and time in MJD",
            "Visit exposure stop date and time in MJD",
            "Flux AB magnitude",
            "Flux magnitude error",
            "Flux AB magnitude upper limit",
            "Selection of rows for VASCA analysis.",
        ],
        "meta": {"INFO": "Light curve magnitude flux table for one source"},
    },
    "ta_sources_lc": {
        "names": [
            "fd_src_id",
            "mag",
            "mag_err",
            "ul",
            "time_bin_start",
            "time_bin_size",
            "ra",
            "dec",
            "pos_err",
        ],
        "dtype": [
            "int64",
            np.object_,
            np.object_,
            np.object_,
            np.object_,
            np.object_,
            np.object_,
            np.object_,
            np.object_,
        ],
        "units": ["1", "mag", "mag", "mag", "d", "s", "deg", "deg", "deg"],
        "defaults": [  # For dtype accuracies when saving, see tables.write_to_fits()
            -1,
            np.array([-1.0], dtype=np.object_),
            np.array([-1.0], dtype=np.object_),
            np.array([-1.0], dtype=np.object_),
            np.array([-1.0], dtype=np.object_),
            np.array([-1.0], dtype=np.object_),
            np.array([-1.0], dtype=np.object_),
            np.array([-1.0], dtype=np.object_),
            np.array([-1.0], dtype=np.object_),
        ],
        "descriptions": [
            "Field source ID Nr.",
            "Flux AB magnitude, set to default when not meassured",
            "Flux magnitude error, set to default when not meassured",
            "Flux AB magnitude upper limit",
            "Visit exposure start date and time in MJD",
            "Visit exposure time in s",
            "Visit detection RA (J2000)",
            "Visit detection Dec (J2000)",
            "Visit position error",
        ],
        "meta": {
            "INFO": "Light curve table for many sources. Stored in a variable-length-array format."
        },
    },
}
# For GALEX data products mcat descriptions see
# http://www.galex.caltech.edu/wiki/Public:Documentation/Appendix_A.1
# %% GALEX field
galex_field = {
    "tt_visits": {
        "names": [
            *base_field["tt_visits"]["names"],
            "time_bin_size_alt_filt",
            "ra",
            "dec",
        ],
        "dtype": [*base_field["tt_visits"]["dtype"], "float64", "float64", "float64"],
        "units": [*base_field["tt_visits"]["units"], "s", "degree", "degree"],
        "defaults": [*base_field["tt_visits"]["defaults"], -1.0, -1.0, -1.0],
        "descriptions": [
            *base_field["tt_visits"]["descriptions"],
            "Visit exposure time of the alternative filter in s",
            "Center RA of the visit FoV (J2000)",
            "Center Dec of the visit FoV (J2000)",
        ],
        "meta": {**base_field["tt_visits"]["meta"]},
    },
    "tt_detections": {
        "names": [
            *base_field["tt_detections"]["names"],
            "det_id",
            "r_fov",
            "artifacts",
            "point_src_prob",
            "bright_match",
            "flux_f60",
            "flux_f60_err",
            "flux_f38",
            "flux_f38_err",
            "E_bv",
        ],
        "dtype": [
            *base_field["tt_detections"]["dtype"],
            "int64",
            "float32",
            "int64",
            "float32",
            "int32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
        ],
        "units": [
            *base_field["tt_detections"]["units"],
            "1",
            "degree",
            "1",
            "1",
            "1",
            "1",
            "1",
            "1",
            "1",
            "mag",
        ],
        "defaults": [
            *base_field["tt_detections"]["defaults"],
            -1,
            -1.0,
            -1,
            -1.0,
            -1,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
        ],
        "descriptions": [
            *base_field["tt_detections"]["descriptions"],
            "Visit detection ID",
            "Distance from center of FOV in degrees",
            "Logical OR of artifact flags",
            "Point-source probability: 0.0 (resolved), 1.0 (unresolved, mcat file filter_CLASS_STAR variable)",  # noqa E501
            "Detection matched to a known star (bright_match=1, mcat file chkobj_type variable)",  # noqa E501
            "Flux in a fixed circular 6.0 arcsec radius aperture in cts/sec",
            "Flux error in a fixed circular 6.0 arcsec radius aperture in cts/sec",
            "Flux in a fixed circular 3.8 arcsec radius aperture in cts/sec",
            "Flux error in a fixed circular 3.8 arcsec radius aperture in cts/sec",
            "Galactic reddening expressed as E(B-V)",
        ],
        "meta": {**base_field["tt_detections"]["meta"], "PRECUTS": "List of pre-cuts"},
    },
    "tt_ref_sources": {
        "names": [
            *base_field["tt_ref_sources"]["names"],
            "r_fov",
            "artifacts",
            "point_src_prob",
            "bright_match",
            "flux_f60",
            "flux_f60_err",
            "flux_f38",
            "flux_f38_err",
            "E_bv",
        ],
        "dtype": [
            *base_field["tt_ref_sources"]["dtype"],
            "float32",
            "int64",
            "float32",
            "int32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
        ],
        "units": [
            *base_field["tt_ref_sources"]["units"],
            "degree",
            "1",
            "1",
            "1",
            "1",
            "1",
            "1",
            "1",
            "mag",
        ],
        "defaults": [
            *base_field["tt_ref_sources"]["defaults"],
            -1.0,
            -1,
            -1.0,
            -1,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
            -1.0,
        ],
        "descriptions": [
            *base_field["tt_ref_sources"]["descriptions"],
            "Distance from center of FOV in degrees",
            "Logical OR of artifact flags",
            "Point-source probability: 0.0 (resolved), 1.0 (unresolved)",
            "Detection matched to a known star (bright_match=1)",
            "Flux in a fixed circular 6.0 arcsec radius aperture in cts/sec",
            "Flux error in a fixed circular 6.0 arcsec radius aperture in cts/sec",
            "Flux in a fixed circular 3.8 arcsec radius aperture in cts/sec",
            "Flux error in a fixed circular 3.8 arcsec radius aperture in cts/sec",
            "Galactic reddening expressed as E(B-V)",
        ],
        "meta": {**base_field["tt_ref_sources"]["meta"], "PRECUTS": "List of pre-cuts"},
    },
}

# %% Region
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
        "dtype": [
            *base_field["tt_fields"]["dtype"],
            "int32",
            "float32",
            "float64",
            "float64",
            "int32",
        ],
        "units": [*base_field["tt_fields"]["units"], "1", "s", "d", "d", "1"],
        "defaults": [
            *base_field["tt_fields"]["defaults"],
            -1,
            -1.0,
            -1.0,
            -1.0,
            -1,
        ],
        "descriptions": [
            *base_field["tt_fields"]["descriptions"],
            "Total number of visits of the field",
            "Total exposure time",
            "Start time of first exposure",
            "End time of last exposure",
            "Region field ID Nr.",
        ],
        "meta": {**base_field["tt_fields"]["meta"]},
    },
    "tt_visits": {
        "names": [*base_field["tt_visits"]["names"], "rg_fd_id"],
        "dtype": [*base_field["tt_visits"]["dtype"], "int32"],
        "units": [*base_field["tt_visits"]["units"], "1"],
        "defaults": [*base_field["tt_visits"]["defaults"], -1],
        "descriptions": [
            *base_field["tt_visits"]["descriptions"],
            "Field source ID nr.",
        ],
        "meta": {**base_field["tt_visits"]["meta"]},
    },
    "tt_coverage_hp": {
        "names": ["pix_id", "nr_vis", "exp", "nr_fds"],
        "dtype": ["uint32", "uint32", "float32", "uint32"],
        "units": ["1", "1", "1", "1"],
        "defaults": [0, 0, -1, 0],
        "descriptions": [
            "Healpix ID",
            "Nr. of visits",
            "Total exposure",
            "Nr. of fields",
        ],
        "meta": {
            "DATAPATH": "None",
            "INFO": "Region observations properties in healpix binning.\
                     RING ordering and equatorial coordinates",
            "NSIDE": "None",
        },
    },
    "tt_ref_sources": {
        "names": [*base_field["tt_ref_sources"]["names"], "rg_fd_id"],
        "dtype": [*base_field["tt_ref_sources"]["dtype"], "int32"],
        "units": [*base_field["tt_ref_sources"]["units"], "1"],
        "defaults": [*base_field["tt_ref_sources"]["defaults"], -1],
        "descriptions": [
            *base_field["tt_ref_sources"]["descriptions"],
            "Field source ID nr.",
        ],
        "meta": {**base_field["tt_ref_sources"]["meta"]},
    },
    "tt_detections": {
        "names": [*base_field["tt_detections"]["names"], "rg_fd_id", "rg_src_id"],
        "dtype": [*base_field["tt_detections"]["dtype"], "int32", "int32"],
        "units": [*base_field["tt_detections"]["units"], "1", "1"],
        "defaults": [*base_field["tt_detections"]["defaults"], -1, -1],
        "descriptions": [
            *base_field["tt_detections"]["descriptions"],
            "Field source ID nr.",
            "Region source ID nr.",
        ],
        "meta": {**base_field["tt_detections"]["meta"]},
    },
    "tt_sources": {
        "names": [
            *base_field["tt_sources"]["names"],
            "rg_fd_id",
            "rg_src_id",
            "nr_fd_srcs",
        ],
        "dtype": [*base_field["tt_sources"]["dtype"], "int32", "int32", "int32"],
        "units": [*base_field["tt_sources"]["units"], "1", "1", "1"],
        "defaults": [*base_field["tt_sources"]["defaults"], -1, -1, -1],
        "descriptions": [
            *base_field["tt_sources"]["descriptions"],
            "Region field ID nr.",
            "Region source ID nr.",
            "Nr. of fields sources in the cluster.",
        ],
        "meta": {**base_field["tt_sources"]["meta"]},
    },
    "ta_sources_lc": {
        "names": [*base_field["ta_sources_lc"]["names"], "rg_fd_id", "rg_src_id"],
        "dtype": [*base_field["ta_sources_lc"]["dtype"], "int32", "int32"],
        "units": [*base_field["ta_sources_lc"]["units"], "1", "1"],
        "defaults": [*base_field["tt_sources"]["defaults"], -1, -1],
        "descriptions": [
            *base_field["ta_sources_lc"]["descriptions"],
            "Field source ID nr.",
            "Region source ID nr.",
        ],
        "meta": {**base_field["ta_sources_lc"]["meta"]},
    },
    "tt_src_id_map": {
        "names": ["rg_src_id", "rg_fd_id", "fd_src_id", "sel"],
        "dtype": ["int32", "int32", "int32", "bool"],
        "units": ["1", "1", "1", "1"],
        "defaults": [-1, -1, -1, True],
        "descriptions": [
            "Region source ID nr.",
            "Region field ID nr.",
            "Field source ID nr.",
            "Selection of rows for VASCA analysis.",
        ],
        "meta": {"INFO": "Map between region and field source IDs."},
    },
}

# global, combined dictionary
class_keys = ["base_field", "galex_field", "region"]
class_dicts = [base_field, galex_field, region]
dd_vasca_tables = {c_key: c_dict for c_key, c_dict in zip(class_keys, class_dicts)}
