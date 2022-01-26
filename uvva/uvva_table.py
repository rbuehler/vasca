import os
from pprint import pprint

import numpy as np
import yaml
from astropy import units as uu
from astropy.table import Table
from loguru import logger

dimless = uu.dimensionless_unscaled

# global paths
FILE_DIR = os.path.dirname(os.path.abspath(__file__))  # path to the dir. of this file
ROOT_DIR = FILE_DIR + "/../"  # path to the root directory of the repository

# global dictionaries defining the table structures
base_field = {
    "tt_field": {
        "names": ["id", "name", "ra", "dec", "observatory", "obsfilter"],
        "dtype": ["uint64", "S64", "float64", "float64", "S64", "S64"],
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
        "names": ["id", "t_start", "t_stop", "t_exp"],
        "dtype": ["uint64", "float64", "float64", "float64"],
        "units": ["1", "1", "1", "s"],
        "descriptions": [
            "Visit ID nr.",
            "Visit exposure start date and time in MJD",
            "Visit exposure stop date and time in MJD",
            "Visit exposure time in s",
        ],
        "meta": {"INFO": "Visit information table"},
    },
    "tt_visit_sources": {
        "names": ["visit_id", "id", "ra", "dec", "mag", "mag_err", "s2n", "flags"],
        "dtype": [
            "uint64",
            "uint64",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
            "int32",
        ],
        "units": ["1", "1", "degree", "degree", "1", "1", "1", "1"],
        "descriptions": [
            "Visit ID",
            "Visit source ID nr.",
            "Visit source RA (J2000)",
            "Visit source Dec (J2000)",
            "Visit source magnitude",
            "Visit source magnitude error",
            "Visit source signal to noise",
            "Visit source flags",
        ],
        "meta": {
            "INFO": "Visit detections table",
            "Name": "visit_sources",
        },
    },
    "tt_reference_sources": {
        "names": ["id", "ra", "dec", "mag", "mag_err", "s2n", "flags"],
        "dtype": [
            "uint64",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
            "int32",
        ],
        "units": ["1", "degree", "degree", "1", "1", "1", "1"],
        "descriptions": [
            "Reference source ID nr.",
            "Reference source RA (J2000)",
            "Reference source Dec (J2000)",
            "Reference source magnitude",
            "Reference source magnitude error",
            "Reference source signal to noise",
            "Reference source flags",
        ],
        "meta": {
            "INFO": "Reference detections table",
            "Name": "reference_sources",
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
            "t_exp_alt_filt",
            "ra",
            "dec",
        ],
        "dtype": [
            *base_field["tt_visits"]["dtype"],
            "float64",
            "float64",
            "float64",
        ],
        "units": [*base_field["tt_visits"]["units"], "s", "degree", "degree"],
        "descriptions": [
            *base_field["tt_visits"]["descriptions"],
            "Visit exposure time of the alternative filter in s",
            "Center RA of the visit FoV (J2000)",
            "Center Dec of the visit FoV (J2000)",
        ],
        "meta": {**base_field["tt_visits"]["meta"]},
    }
}
# global, combined dictionary
class_keys = ["base_field", "galex_field"]
class_dicts = [base_field, galex_field]
dd_uvva_tables = {c_key: c_dict for c_key, c_dict in zip(class_keys, class_dicts)}


class UVVATable(Table):
    def __init__(self):
        # Configure logger
        # adds the class name as an extra key; accessible vie the handler format
        logger.configure(extra={"classname": self.__class__.__name__})

        self.combined_templates = dd_uvva_tables

    @staticmethod
    def from_template(data, template_name):
        """
        Creates a new astropy table.

        Parameters
        ----------
        data : list, array-like
            Data of the table with shape (n, n_cols), with ''n'' arbitrary number
            of rows and ''n_cols'' number of columns that must match the
            length of the table attributes ''colnames'', ''units'', ''dtype'',
            ''descriptions'' and ''meta''.
        template_name : str
            Identifier to select a table template. Templates are selected by
            setting the class key and a corresponding table key in one string
            separated by a colon, e.g. template_name=<class_key>:<table_key>.

        Returns
        -------
        astropy.table.Table
        """

        # Takes pre-defined template dictionary
        templates = dd_uvva_tables

        # Parse template identifier keys
        if template_name is not None:
            class_key = template_name.split(":")[0]
            table_key = template_name.split(":")[1]

        # Generate lists of available class and table keys
        class_keys = [clss for clss in templates.keys()]
        table_keys = [
            tbl for clss in templates.keys() for tbl in templates[clss].keys()
        ]

        # Check if template exists for template_name
        if class_key not in class_keys:
            raise KeyError(
                f"Unknown class key '{class_key}'. Choose one from {class_keys}"
            )
        if table_key not in table_keys:
            raise KeyError(
                f"Unknown table key '{table_key}'. Choose one from {table_keys}"
            )

        # Create table
        data = np.asarray(data)
        tt_out = Table(data=data, **templates[class_key][table_key])

        # logging
        logger.debug(f"Created new table from template '{template_name}'.")
        return tt_out

    @staticmethod
    def uvva_generic(
        data, colnames=None, units=None, dtype=None, descriptions=None, meta=None
    ):
        """
        Returns a new astropy table following the UVVA required format
        """

        # Check usage
        if any(v is None for v in [colnames, units, dtype, descriptions, meta]):
            raise ValueError(
                "Cannot add table from scratch without specifing all "
                "table parameters (colnames, units, dtype, descriptions, meta)"
            )

        # Create table
        data = np.asarray(data)
        # Check if table attributes exist for all columns
        if not all(
            len(attr) == data.shape[1]
            for attr in [colnames, units, dtype, descriptions, meta]
        ):
            raise ValueError(
                "Cannot create table from scratch. "
                "All table attributes and data must have the same length."
            )
        tt_out = Table(
            data=np.asarray(data),
            names=colnames,
            units=units,
            dtype=dtype,
            descriptions=descriptions,
            meta=meta,
        )

        # Logging
        logger.debug("Created new table from scratch")

        return tt_out
