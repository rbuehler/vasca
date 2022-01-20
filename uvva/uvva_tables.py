import os
from pprint import pprint

import yaml

# global paths
FILE_DIR = os.path.dirname(os.path.abspath(__file__))  # path to the dir. of this file
ROOT_DIR = FILE_DIR + "/../"  # path to the root directory of the repository

base_field = {
    "tt_field": {
        "names": ["id", "name", "ra", "dec", "observatory", "obsfilter"],
        "dtype": ["uint32", "S32", "float16", "float16", "S32", "S32"],
        "units": ["1", "", "degree", "degree", "", ""],
        "descriptions": [
            "Field ID nr.",
            "Field name",
            "Center RA of the field (J2000)",
            "Center Dec of the field (J2000)",
            "Telescope of the observation (e.g. GALEX)",
            "Filter of the observation (e.g. NUV)",
        ],
        "meta": {"DATAPATH": "None", "INFO": "Field information table"},
    },
    "tt_visits": {
        "names": ["id", "time", "exposure"],
        "dtype": ["uint32", "float64", "float16"],
        "units": ["1", "1", "s"],
        "descriptions": [
            "Visit ID nr.",
            "Visit central time in MJD",
            "Visit exposure time in s",
        ],
        "meta": {"INFO": "Visit information table"},
    },
    "tt_visit_sources": {
        "names": ["id", "ra", "dec", "mag", "mag_err", "s2n", "flags", "src_id"],
        "dtype": [
            "uint32",
            "float16",
            "float16",
            "float16",
            "float16",
            "float16",
            "int32",
            "uint32",
        ],
        "units": ["1", "degree", "degree", "1", "1", "1", "1", "1"],
        "descriptions": [
            "Visit source ID nr.",
            "Visit source RA (J2000)",
            "Visit source Dec (J2000)",
            "Visit source magnitude",
            "Visit source magnitude error",
            "Visit source signal to noise",
            "Visit source flags",
            "Source ID associated to the visit source",
        ],
        "meta": {
            "INFO": "Visit detections table",
            "Name": "visit_sources",
        },
    },
    "tt_sources": {
        "names": ["src_id", "ra", "dec", "nr_vis_det", "flag"],
        "dtype": ["uint32", "float16", "float16", "uint32", "int32"],
        "units": ["1", "degree", "degree", "1", "1"],
        "descriptions": [
            "Source ID nr.",
            "Source RA (J2000)",
            "Source Dec (J2000)",
            "Number of visit detections of the source",
            "Source flags",
        ],
        "meta": {"INFO": "Source infomation table", "CLUSTALG": "None"},
    },
    "t_sources_mag": {
        "meta": {"INFO": "AB Magnitude flux table"},
    },
    "tt_sources_s2n": {
        "meta": {"INFO": "Signal to noise of the detection table"},
    },
    "tt_sources_ulmag": {
        "meta": {"INFO": "AB Magnitude upper limit table", "CONFLEVE": "0.95"}
    },
}
galex_field = {
    "tt_visits": {
        "names": [
            *base_field["tt_visits"]["names"],
            "exposure_alt_filt",
            "ra",
            "dec",
            "t_stop",
            "t_start",
        ],
        "dtype": [
            *base_field["tt_visits"]["dtype"],
            "float16",
            "float16",
            "float16",
            "S32",
            "S32",
        ],
        "units": [*base_field["tt_visits"]["units"], "s", "degree", "degree", "", ""],
        "descriptions": [
            *base_field["tt_visits"]["descriptions"],
            "Visit exposure time of the alternative filter in s",
            "Center RA of the visit FoV (J2000)",
            "Center Dec of the visit FoV (J2000)",
            "Visit start time in 'mm/dd/yyyy hh:mm:ss <AM/PM>' format",
            "Visit start time in 'mm/dd/yyyy hh:mm:ss <AM/PM>' format",
        ],
        "meta": {**base_field["tt_visits"]["meta"]},
    }
}
class_keys = ["base_field", "galex_field"]
class_dicts = [base_field, galex_field]
combined = {c_key: c_dict for c_key, c_dict in zip(class_keys, class_dicts)}

with open(f"{FILE_DIR}/uvva_tables_auto_gen.yml", "w") as yaml_file:
    yaml.dump(combined, yaml_file, default_flow_style=False)
