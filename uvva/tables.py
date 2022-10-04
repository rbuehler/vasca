import os

import h5py
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from astropy import units as uu
from astropy.io import fits
from astropy.table import unique, Column, Table
from astropy.wcs import wcs
from astropy.nddata import bitmask
from loguru import logger

# import warnings
# from astropy.io.fits.verify import VerifyWarning

# deactivate warnings
# warnings.simplefilter('ignore', category=VerifyWarning)

dimless = uu.dimensionless_unscaled

# global paths
# path to the dir. of this file
FILE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = FILE_DIR + "/../"  # path to the root directory of the repository

# global dictionaries defining the table structures
base_field = {
    "tt_field": {
        "names": ["field_id", "name", "ra", "dec", "observatory", "obs_filter"],
        "dtype": ["int64", "S64", "float64", "float64", "S64", "S64"],
        "units": ["1", "", "degree", "degree", "", ""],
        "descriptions": [
            "Field ID nr.",
            "Field name",
            "Center RA of the field (J2000)",
            "Center Dec of the field (J2000)",
            "Telescope of the observation (e.g. GALEX)",
            "Filter of the observation (e.g. NUV)",
        ],
        "meta": {
            "DATAPATH": "None",
            "INFO": "Field information table",
        },
    },
    "tt_visits": {
        "names": [
            "vis_id",
            "time_bin_start",
            "time_bin_size",
        ],
        "dtype": [
            "int64",
            "float64",
            "float64",
        ],
        "units": ["1", "d", "s"],
        "descriptions": [
            "Visit ID nr.",
            "Visit exposure start date and time in MJD",
            "Visit exposure time in s",
        ],
        "meta": {"INFO": "Visit information table"},
    },
    "tt_detections": {
        "names": [
            "vis_id",
            "src_id",
            "det_id",
            "ra",
            "dec",
            "pos_err",
            "mag",
            "mag_err",
            "s2n",
        ],
        "dtype": [
            "int64",
            "int64",
            "int64",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
        ],
        "units": [
            "1",
            "1",
            "1",
            "degree",
            "degree",
            "degree",
            "1",
            "1",
            "1",
        ],
        "descriptions": [
            "Visit ID associated to the visit detection",
            "Source ID associated to the visit detection",
            "Visit detection ID",
            "Visit detection RA (J2000)",
            "Visit detection Dec (J2000)",
            "Visit position error",
            "Visit detection magnitude",
            "Visit detection magnitude error",
            "Signal to noise",
        ],
        "meta": {
            "INFO": "Visit detections table",
        },
    },
    "tt_ref_sources": {
        "names": [
            "det_id",
            "ra",
            "dec",
            "pos_err",
            "mag",
            "mag_err",
            "s2n",
        ],
        "dtype": [
            "int64",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
        ],
        "units": [
            "1",
            "degree",
            "degree",
            "degree",
            "1",
            "1",
            "1",
        ],
        "descriptions": [
            "Reference source ID nr.",
            "Reference source RA (J2000)",
            "Reference source Dec (J2000)",
            "Reference position error",
            "Reference source magnitude",
            "Reference source magnitude error",
            "Signal to noise",
        ],
        "meta": {
            "INFO": "Reference detections table",
        },
    },
    "tt_sources": {
        "names": [
            "src_id",
            "ra",
            "dec",
            "nr_det",
            "nr_uls",
            "mag_mean",
            "mag_var",
            "mag_rchiq",
            "mag_dmax",
            "mag_dmax_sig",
            "ul_weight",
        ],
        "dtype": [
            "int64",
            "float64",
            "float64",
            "int64",
            "int64",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
        ],
        "units": ["1", "degree", "degree", "1", "1", "1", "1", "1", "1", "1", "1"],
        "descriptions": [
            "Source ID nr.",
            "Source RA (J2000)",
            "Source Dec (J2000)",
            "Number of visit detections of the source",
            "Number of upper limits",
            "Mean flux magnitude",
            "Flux magnitude variace",
            "Flux magnitude reduced chisquared of the constant mean",
            "Maximum magnitude flux variation compared to average",
            "Maximum magnitude flux variation compared to average divided by magnitude error",
            "Nr of upper limits magnitude greater than mean magnitude divided by srt(nr_det)",
        ],
        "meta": {"INFO": "Source infomation table", "CLUSTALG": "None"},
    },
    "tt_source_lc": {
        "names": ["time_start", "time_delta", "mag", "mag_err", "ul"],
        "dtype": ["float64", "float64", "float64", "float64", "float64"],
        "units": ["d", "s", "1", "1", "1"],
        "descriptions": [
            "Visit exposure start date and time in MJD",
            "Visit exposure stop date and time in MJD",
            "Flux magnitude",
            "Flux magnitude error",
            "Flux magnitude upper limit",
        ],
        "meta": {"INFO": "Light curve magnitude flux table for one source"},
    },
    "ta_sources_lc": {
        "names": ["src_id", "mag", "mag_err", "ul"],
        "dtype": ["int64", np.object_, np.object_, np.object_],
        "units": ["1", "1", "1", "1"],
        "descriptions": [
            "Source ID Nr.",
            "Flux magnitude",
            "Flux magnitude error",
            "Flux magnitude upper limit",
        ],
        "meta": {
            "INFO": "Light curve table for many sources. Stored in a variable-length-array format."
        },
    },
}
# For GALEX data products mcat descriptions see
# http://www.galex.caltech.edu/wiki/Public:Documentation/Appendix_A.1
galex_field = {
    "tt_visits": {
        "names": [
            *base_field["tt_visits"]["names"],
            "time_bin_size_alt_filt",
            "ra",
            "dec",
        ],
        "dtype": [
            *base_field["tt_visits"]["dtype"],
            "float64",
            "float64",
            "float64",
        ],
        "units": [
            *base_field["tt_visits"]["units"],
            "s",
            "degree",
            "degree",
        ],
        "descriptions": [
            *base_field["tt_visits"]["descriptions"],
            "Visit exposure time of the alternative filter in s",
            "Center RA of the visit FoV (J2000)",
            "Center Dec of the visit FoV (J2000)",
        ],
        "meta": {
            **base_field["tt_visits"]["meta"],
        },
    },
    "tt_detections": {
        "names": [
            *base_field["tt_detections"]["names"],
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
            "float64",
            "int64",
            "float64",
            "int64",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
        ],
        "units": [
            *base_field["tt_detections"]["units"],
            "degree",
            "1",
            "1",
            "1",
            "1",
            "1",
            "1",
            "1",
            "1",
        ],
        "descriptions": [
            *base_field["tt_detections"]["descriptions"],
            "Distance from center of FOV in degrees",
            "Logical OR of artifact flags",
            "Point-source probability: 0.0 (resolved), 1.0 (unresolved, mcat file filter_CLASS_STAR variable)",
            "Detection matched to a known star (bright_match=1, mcat file chkobj_type variable)",
            "Flux in a fixed circular 6.0 arcsec radius aperture in cts/sec",
            "Flux error in a fixed circular 6.0 arcsec radius aperture in cts/sec",
            "Flux in a fixed circular 3.8 arcsec radius aperture in cts/sec",
            "Flux error in a fixed circular 3.8 arcsec radius aperture in cts/sec",
            "Galactic reddening expressed as E(B-V) [mag]",
        ],
        "meta": {
            **base_field["tt_detections"]["meta"],
            "PRECUTS": "List of pre-cuts",
        },
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
            "float64",
            "int64",
            "float64",
            "int64",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
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
            "1",
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
            "Galactic reddening expressed as E(B-V) [mag]",
        ],
        "meta": {
            **base_field["tt_ref_sources"]["meta"],
            "PRECUTS": "List of pre-cuts",
        },
    },
}
region = {
    "tt_fields": {
        "names": [
            *base_field["tt_field"]["names"],
            "fov_diam",  # TODO: Add sky angle for ULTRASAT squared FoV
            "nr_vis",
            "time_bin_size_sum",
            "time_start",
            "time_stop",
        ],
        "dtype": [
            *base_field["tt_field"]["dtype"],
            "float64",
            "int64",
            "float64",
            "float64",
            "float64",
        ],
        "units": [
            *base_field["tt_field"]["units"],
            "degree",
            "1",
            "s",
            "d",
            "d",
        ],
        "descriptions": [
            *base_field["tt_field"]["descriptions"],
            "Field radius or box size (depending on the observatory)",
            "Total number of visits of the field",
            "Total exposure time",
            "Start time of first exposure",
            "End time of last exposure",
        ],
        "meta": {
            **base_field["tt_field"]["meta"],
        },
    },
    "tt_visits": {
        "names": [
            *base_field["tt_visits"]["names"],
            "field_id",
        ],
        "dtype": [
            *base_field["tt_visits"]["dtype"],
            "int64",
        ],
        "units": [
            *base_field["tt_visits"]["units"],
            "1",
        ],
        "descriptions": [
            *base_field["tt_visits"]["descriptions"],
            "Field ID nr.",
        ],
        "meta": {
            **base_field["tt_visits"]["meta"],
        },
    },
    "tt_coverage_hp": {
        "names": ["pix_id", "nr_vis", "exp"],
        "dtype": ["uint32", "uint32", "float32"],
        "units": ["1", "1", "1"],
        "descriptions": [
            "Healpix ID",
            "Nr of visits",
            "Total exposure",
        ],
        "meta": {
            "DATAPATH": "None",
            "INFO": "Region observations properties in healpix binning.\
                     RING ordering and equatorial coordinates",
            "NSIDE": "None",
        },
    },
    "ta_sources_lc": {
        "names": [
            *base_field["ta_sources_lc"]["names"],
            "field_id",
        ],
        "dtype": [
            *base_field["ta_sources_lc"]["dtype"],
            "int64",
        ],
        "units": [
            *base_field["ta_sources_lc"]["units"],
            "1",
        ],
        "descriptions": [
            *base_field["ta_sources_lc"]["descriptions"],
            "Field ID nr.",
        ],
        "meta": {
            **base_field["ta_sources_lc"]["meta"],
        },
    },
}

# global, combined dictionary
class_keys = ["base_field", "galex_field", "region"]
class_dicts = [base_field, galex_field, region]
dd_uvva_tables = {c_key: c_dict for c_key, c_dict in zip(class_keys, class_dicts)}


class TableCollection(object):
    def __init__(self):
        # Configure logger
        # adds the class name as an extra key; accessible vie the handler format
        logger.configure(extra={"classname": self.__class__.__name__})

        self.combined_templates = dd_uvva_tables

        self._table_names = list()

    @staticmethod
    def table_from_template(data, template_name):
        """
        Creates a new astropy table.

        Parameters
        ----------
        data : list, array-like
            Data of the table with shape (n, n_cols) or as dictionaty with the
            key corresponding to the templates columns.
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
        tt_out = Table(data=data, **templates[class_key][table_key])
        # tt_out.meta["template"] = template_name

        # logging
        logger.debug(f"Created new table from template '{template_name}'.")
        return tt_out

    def add_table(self, data, template_name, add_sel_col=True):
        """
        Add a UVVA table to the field.

        Parameters
        ----------
        data : list, array-like
            Data of the table with shape (n, n_cols) or as dictionary with the
            key corresponding to the templates columns.
        template_name : str
            Identifier to select a table template. Templates are selected by
            setting the class key and a corresponding table key in one string
            separated by a colon, e.g. template_name=<class_key>:<table_key>.
        add_sel_col: bool
            Add a selection column to the table. Default is true.

        """
        logger.debug(f"Adding table '{template_name}'")

        table_key = template_name.split(":")[1]

        if table_key in self._table_names:
            logger.warning(f"Table '{table_key}' already exists, overwriting")
        self._table_names.append(table_key)

        tt = self.table_from_template(data, template_name)
        if add_sel_col:
            default_sel = np.ones(len(tt), dtype="int64")
            col = Column(
                data=default_sel,
                name="sel",
                dtype="bool",
                unit="1",
                description="Selection of rows for UVVA analysis.",
            )
            tt.add_column(col)

        setattr(self, table_key, tt)

    def write_to_fits(
        self, file_name="tables.fits", overwrite=True, fits_verify="warn"
    ):
        """
        Write tables and image of a field to a fits file.

        Parameters
        ----------
        file_name : str, optional
            File name. The default is "field_default.fits".
        overwrite : bool, optional
            Overwrite existing file. The default is True.
        fits_verfy: str, optional
            Verify if output is compatible with FITS format. Options are:
            'exception', 'ignore', 'fix', 'silentfix', 'warn'
            See https://docs.astropy.org/en/stable/io/fits/api/verification.html
            The default is 'warn'.

        Returns
        -------
        None.

        """
        logger.info(f"Writing file with name '{file_name}'")

        # Create HDU list and write
        hdup = fits.PrimaryHDU()

        # Check if image data is set and add to primary HDU
        if hasattr(self, "ref_img"):
            if self.ref_img is not None and self.ref_wcs is not None:
                logger.debug("Storing image data'")
                hdup = fits.PrimaryHDU(self.ref_img, header=self.ref_wcs.to_header())

        hdus = [hdup]
        new_hdul = fits.HDUList([hdup])
        new_hdul.writeto(file_name, overwrite=overwrite, output_verify=fits_verify)

        for key in self._table_names:
            if key in self.__dict__:
                logger.debug(f"Writing table '{key}'")
                if not key.startswith("ta_"):
                    self.__dict__[key].write(file_name, append=True)
                # Write table with vector entries, currently not supportet by astropy.Table
                else:
                    cols = list()
                    for colname in self.__dict__[key].colnames:
                        coldata = self.__dict__[key][colname].data

                        col_for = "PE()"
                        # TODO: Make this more general
                        if colname == "field_id" or colname == "src_id":
                            col_for = "K"
                        col = fits.Column(
                            name=colname,
                            format=col_for,
                            array=coldata,
                        )
                        cols.append(col)

                    hdu_ta = fits.BinTableHDU.from_columns(cols)
                    with fits.open(file_name, mode="append") as hdula:
                        hdula.append(hdu_ta)

        # Rename extensions to table names
        ext_nr = 0
        with fits.open(file_name, "update", output_verify=fits_verify) as ff:
            for key in self._table_names:
                if key in self.__dict__:
                    ext_nr += 1
                    ff[ext_nr].header["EXTNAME"] = key
                    ff.flush(output_verify=fits_verify)
        ff.close()

    def load_from_fits(self, file_name="tables.fits"):
        """
        Loads field from a fits file

        Parameters
        ----------
        file_name : str, optional
            File name. The default is "field_default.fits".

        Returns
        -------
        None.

        """
        logger.debug(f"Loading file with name '{file_name}'")
        with fits.open(file_name) as ff:
            # Load tables
            # get available table names
            tt_names = [
                hdu.header["EXTNAME"]
                for hdu in ff[1:]
                if hdu.header["EXTNAME"].startswith("tt_")
            ]
            # loop over tables
            for tt_name in tt_names:
                logger.debug(f"Loading table '{tt_name}'")
                # add to table manifest
                if tt_name in self._table_names:
                    logger.warning(f"Table '{tt_name}' already exists, overwriting.")
                else:
                    self._table_names.append(tt_name)
                setattr(self, tt_name, Table.read(file_name, hdu=tt_name))

            # Load image data
            if hasattr(self, "ref_img"):
                self.ref_img = ff[0].data
                if self.ref_img is not None:
                    self.ref_wcs = wcs.WCS(ff[0].header)
                else:
                    self.ref_wcs = None

            # Load tables with vectors
            ta_names = [
                hdu.header["EXTNAME"]
                for hdu in ff[1:]
                if hdu.header["EXTNAME"].startswith("ta_")
            ]
            for ta_name in ta_names:
                logger.debug(f"Loading table '{ta_name}'")
                # add to table manifest
                if ta_name in self._table_names:
                    logger.warning(f"Table '{ta_name}' already exists, overwriting.")
                else:
                    self._table_names.append(ta_name)

                # Load table data into dictionary
                ta_data = {}
                col_names = ff[ta_name].columns.names
                for col_name in col_names:
                    ta_data[col_name] = ff[ta_name].data[col_name]

                # TODO make loading more general.
                # Expect astropy handling of fits vecotors to simplify this in the future
                if "field_id" in col_names:
                    self.add_table(ta_data, "region:" + ta_name, add_sel_col=False)
                else:
                    self.add_table(ta_data, "base_field:" + ta_name, add_sel_col=False)

    def write_to_hdf5(self, file_name="tables.hdf5"):
        """
        Write tables of a field to a hdf5 file.

        Parameters
        ----------
        file_name : str, optional
            File name. The default is "field_default.fits".
        overwrite : bool, optional
            Overwrite existing file. The default is True.

        Returns
        -------
        None.

        """

        logger.info(f"Writing file with name '{file_name}'")

        ii = 0
        for key in self._table_names:
            # TODO: Make hdf5 work with numpy.dtype.object objects
            if key.startswith("ta_"):
                logger.warning(f"Not writing vector table '{key}'")
            elif key in self.__dict__:
                logger.debug(f"Writing table '{key}'")
                ii += 1
                if ii == 1:
                    self.__dict__[key].write(
                        file_name,
                        path="TABDATA/" + key,
                        overwrite=True,
                        serialize_meta=True,
                    )
                else:
                    self.__dict__[key].write(
                        file_name,
                        path="TABDATA/" + key,
                        append=True,
                        serialize_meta=True,
                    )

    def load_from_hdf5(self, file_name="tables.hdf5"):
        """
        Loads field from a hdf5 file

        Parameters
        ----------
        file_name : str, optional
            File name. The default is "field_default.hdf5".

        Returns
        -------
        None.

        """

        logger.info(f"Loading file with name '{file_name}'")
        in_file = h5py.File(file_name, "r")

        for table in in_file["TABDATA"].keys():
            if "meta" in str(table):
                continue
            logger.debug(f"Loading table '{table}'")
            self._table_names.append(str(table))
            setattr(
                self, str(table), Table.read(file_name, path="TABDATA/" + str(table))
            )

    def info(self):
        """
        Print out information about the field, its visits and sources.

        Returns
        -------
        None.

        """
        for key in self._table_names:
            if key in self.__dict__:
                print(f"\n {key}:")
                self.__dict__[key].info()
                print(self.__dict__[key].meta)

    def __str__(self):
        """
        Return string with information about the field, its visits and sources.

        Returns
        -------
        str

        """
        out_str = ""
        for key in self._table_names:
            if key in self.__dict__:
                out_str += "\n" + self.__dict__[key].__str__()

        return out_str
        return out_str

    def select_rows(self, selections, remove_unselected=False):
        """
        Apply selection to a passed table.

        Parameters
        ----------
        selections : dict
            Dictionary with selection table, variables and cut values
        remove_unselected: bool
            Remove table rows with entry False in 'sel' table.

        Returns
        -------
        None.

        """

        # Get table and check if selection column is available
        table_name = selections["table"]
        logger.info(f"Applying selection on table '{table_name}'")
        if table_name not in self.__dict__.keys():
            logger.error("Table does not exist, it need to be created beforehand.")
        tt = self.__dict__[table_name]
        if "sel" not in tt.colnames:
            logger.error("Table does not have selection column")

        sel = tt["sel"].data.astype("bool")
        nr_sel = sel.sum()
        sel = tt["sel"].data.astype("bool")
        nr_sel = sel.sum()

        if selections["sel_type"] == "and":

            # Apply min/max cuts
            if "range" in selections.keys():
                for var, vals in selections["range"].items():
                    sel = sel * (tt[var] >= vals[0]) * (tt[var] <= vals[1])
                    logger.debug(
                        f"AND selecting '{var}' {vals}, kept: {100*sel.sum()/nr_sel : .4f}%"
                    )

            # Apply bitmask cuts
            if "bitmask" in selections.keys():
                for var, vals in selections["bitmask"].items():
                    no_art = sel * (tt[var].data.astype("int") == 0)
                    bit = bitmask.bitfield_to_boolean_mask(tt[var], ignore_flags=vals)
                    sel = sel * (no_art + bit)
                    logger.debug(
                        f"AND selecting bitmask '{var}' keep {vals}, kept: {100*sel.sum()/nr_sel : .4f}%"
                    )
        elif selections["sel_type"] == "or":
            sel_or = np.zeros(len(sel))

            if "range" in selections.keys():
                for var, vals in selections["range"].items():
                    sel_or = sel_or + (tt[var] >= vals[0]) * (tt[var] <= vals[1])
                    logger.debug(
                        f"OR selecting '{var}' {vals}, kept: {100*sel_or.sum()/nr_sel : .4f}%"
                    )
                sel = sel * sel_or
        else:
            logger.error("Unkown selection type.")

        sel = sel.astype("bool")
        tt.replace_column("sel", sel)

        if remove_unselected:
            tt = tt[sel]

    def plot_hist(self, table_name, var, ax=None, logx=False, **hist_kwargs):
        """
        Plot histogram for passed table variable

        Parameters
        ----------
        table_name : str
            Table name.
        var : str, optional
            Variable name
        ax : matplotlib.axes, optional
            Axes to draw on. The default is None.
        logx : bool, optional
            Histoogram of log10(var) instead of var. The default is False.
        **hist_kwargs : dict
            Key word arguments passed tu plt.hist

        Returns
        -------
        ax : matplotlix.axes
            Axes that where used to draw.

        """
        logger.debug(f"Plotting histogram of variable '{var}' in table '{table_name}'")

        if ax is None:
            ax = plt.gca()

        # Set marker properties for sources
        plot_kwargs = {
            "bins": "auto",
            "histtype": "stepfilled",
            "log": True,
            "alpha": 0.5,
        }
        if hist_kwargs is not None:
            plot_kwargs.update(hist_kwargs)

        tt = self.__dict__[table_name]
        col = tt[var]
        sel = tt["sel"]
        str_nrsel = str(sel.sum())
        str_nrnotsel = str((~sel).sum())
        data = [col[~sel], col[sel]]
        xlabel = var + " [" + str(col.unit) + "]"
        if str(col.unit) == "None" or str(col.unit) == "":
            xlabel = var
        if logx:
            with np.errstate(divide="ignore", invalid="ignore"):
                data = [np.log10(data[0]), np.log10(data[1])]
            xlabel = "log10( " + xlabel + " )"

        ax.hist(
            data,
            label=["unselected_" + str_nrnotsel, "selected_" + str_nrsel],
            **plot_kwargs,
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Counts")

        return ax

    def plot_scatter(
        self,
        table_name,
        varx,
        vary,
        ax=None,
        xlim=None,
        ylim=None,
        invert_xaxis=None,
        invert_yaxis=None,
        xscale="linear",
        yscale="linear",
        **scatter_kwargs,
    ):
        """
        Plot scatter plot for passed table variable

        Parameters
        ----------
        table_name : str
            Table name.
        varx : str
            Variable name on X-axis
        vary : str
            Variable name on Y-axis
        ax : matplotlib.axes, optional
            Axes to draw on. The default is None.
        xlim : list, optional
            List with [xmin, xmax] axis value. Default is None.
        ylim : list, optional
            List with [ymin, ymax] axis value. Default is None.
        xscale : str, optional
            Type of x-scale ("log", "linear"). Default is "linear".
        yscale : str, optional
            Type of y-scale ("log", "linear"). Default is "linear".
        **scatter_kwargs : dict
            Key word arguments passed tu plt.scatter

        Returns
        -------
        ax : matplotlix.axes
            Axes that where used to draw.

        """
        logger.debug(
            f"Plotting of variables '{varx}' and '{vary}' in table '{table_name}'"
        )

        if ax is None:
            ax = plt.gca()

        # Set marker properties for sources
        plot_kwargs = {"s": 2.0, "alpha": 0.5}
        if scatter_kwargs is not None:
            plot_kwargs.update(scatter_kwargs)

        tt = self.__dict__[table_name]
        sel = tt["sel"].astype(bool)
        str_nrsel = str(sel.sum())
        str_nrnotsel = str((~sel).sum())
        ax.scatter(
            tt[varx][~sel],
            tt[vary][~sel],
            label="unselected_" + str_nrnotsel,
            **plot_kwargs,
        )
        ax.scatter(
            tt[varx][sel], tt[vary][sel], label="selected_" + str_nrsel, **plot_kwargs
        )

        # Set labels
        xlabel = varx + " [" + str(tt[varx].unit) + "]"
        if str(tt[varx].unit) == "None" or str(tt[varx].unit) == "":
            xlabel = varx
        ylabel = vary + " [" + str(tt[vary].unit) + "]"
        if str(tt[vary].unit) == "None" or str(tt[vary].unit) == "":
            ylabel = vary
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        # Set axis limits
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)

        ax.set_xscale(xscale)
        ax.set_yscale(yscale)

        if invert_xaxis:
            ax.invert_xaxis()
        if invert_yaxis:
            ax.invert_yaxis()

        return ax

    def get_light_curve(self, src_ids, field_ids=None):
        """
        Get a light curve for one source or a list of sources.

        Parameters
        ----------
        src_ids : list or int
            Source ID(s) to return the light curve.
        field_ids: list or int
            Field IDs. Array length has to match source IDs. If None
            assumes field_id from current field. Default is None.

        Returns
        -------
        lc_dict : list or table
            Light curve as an astropy Table compatible
            with astropy BinnedTimeSeries.

        """
        if not hasattr(src_ids, "__iter__"):
            src_ids = [src_ids]
        if not hasattr(field_ids, "__iter__"):
            field_ids = [field_ids]

        logger.debug(f"Getting lightcurve for Nr src: {len(src_ids)}")

        if "ta_sources_lc" not in self._table_names:
            logger.error(
                "Light curve table does not exist, run 'set_light_curve()' first."
            )

        # Dictionary to store light curve tables
        lc_dict = dict()

        # Loop over sources and get lc info
        # If called from a field no field ID is passed
        if field_ids[0] == None:
            for src_id in src_ids:
                sel = self.ta_sources_lc["src_id"] == src_id
                src_data = {
                    "time_start": self.tt_visits["time_bin_start"].data,
                    "time_delta": self.tt_visits["time_bin_size"].data,
                    "mag": self.ta_sources_lc[sel]["mag"].data[0].astype(np.float64),
                    "mag_err": self.ta_sources_lc[sel]["mag_err"]
                    .data[0]
                    .astype(np.float64),
                    "ul": self.ta_sources_lc[sel]["ul"].data[0].astype(np.float64),
                }
                # Create and store table
                tt_lc = self.table_from_template(src_data, "base_field:tt_source_lc")
                tt_lc.meta["src_id"] = src_id
                lc_dict[src_id] = tt_lc
        # If called from a region field_id needs to be matched too
        else:
            self.tt_visits.add_index("field_id")
            for src_id, field_id in zip(src_ids, field_ids):
                sel = (self.ta_sources_lc["src_id"] == src_id) * (
                    self.ta_sources_lc["field_id"] == field_id
                )

                # Get visit info and sort it
                vis_idx = self.tt_visits.loc_indices["field_id", field_id]
                tt_vis = unique(self.tt_visits[vis_idx], "time_bin_start")

                # Store data
                src_data = {
                    "time_start": tt_vis["time_bin_start"].data,
                    "time_delta": tt_vis["time_bin_size"].data,
                    "mag": self.ta_sources_lc[sel]["mag"].data[0],
                    "mag_err": self.ta_sources_lc[sel]["mag_err"].data[0],
                    "ul": self.ta_sources_lc[sel]["ul"].data[0],
                }

                # Create table
                tt_lc = self.table_from_template(src_data, "base_field:tt_source_lc")
                tt_lc.meta["src_id"] = src_id
                tt_lc.meta["field_id"] = src_id
                lc_dict[src_id] = tt_lc

        # If only one src_id was passed do not return as list
        if len(lc_dict) == 1:
            lc_dict = list(lc_dict.values())[0]

        return lc_dict

    def plot_light_curve(
        self,
        src_ids,
        field_ids=None,
        ax=None,
        ylim=None,
        legend_loc="upper center",
        plot_upper_limits=True,
        **errorbar_kwargs,
    ):
        """
        Plot the magnitude light curves of the passed sources.

        Parameters
        ----------
        src_ids : list or int
            List or single source IDs to plot.
        field_ids: list or int
            Field IDs. Array length has to match source IDs. If None
            assumes field_id from current field. Default is None.
        ax : axes, optional
                Matplotlib axes to plot on. The default is None.
        ylim : list, optional
            Limits of the y axis. Default is None
        legend_loc : string, optional
            Position of the legend in the figure. The default is "upper right".
        plot_upper_limits : bool
            Plot upper limits to the lightcurve. The default is True.
        **errorbar_kwargs : TYPE
            Key word arguments for pyplot.errorbars plotting.

        Returns
        -------
        ax : axes
            Used Matplotlib axes.

        """

        # If only one src_id/field_id was passed create a list
        if not hasattr(src_ids, "__iter__"):
            src_ids = [src_ids]
        if not hasattr(field_ids, "__iter__"):
            field_ids = [field_ids]
        if field_ids[0] == None:
            field_ids = [None] * len(src_ids)

        logger.info("Plotting light curves")

        # Setup plotting parameters
        if ax is None:
            ax = plt.gca()
        ax.invert_yaxis()
        if hasattr(ylim, "__iter__"):
            ax.set_ylim(ylim)

        plt_errorbar_kwargs = {
            "markersize": 4,
            "alpha": 0.6,
            "capsize": 0,
            "lw": 0.1,
            "linestyle": "dotted",
            "elinewidth": 0.6,
        }
        if errorbar_kwargs is not None:
            plt_errorbar_kwargs.update(errorbar_kwargs)

        # Loop over selected sources and plot
        colors = cycle("bgrcmykbgrcmykbgrcmykbgrcmyk")
        markers = cycle("osDd<>^v")
        ctr = 0
        for src_id, field_id, col, mar in zip(src_ids, field_ids, colors, markers):

            # Every 8 markers plot open symbols
            ctr += 1
            mfc = col
            if ctr > 8:
                mfc = "None"

            # Get light curve
            lc = self.get_light_curve(src_id, field_id)

            # Get arrays
            src_lab = str(src_id)
            uplims = np.zeros(len(lc))
            sel = lc["mag"] > 0
            mags = lc["mag"]
            mags_err = lc["mag_err"]
            ul = lc["ul"]

            # Modify arrays if upper limits are plotted
            if plot_upper_limits:
                uplims = lc["mag"] < 0
                sel = np.ones(len(lc), dtype=bool)
                mags = mags * ~uplims + ul * uplims
                mags_err = mags_err * ~uplims + 0.1 * uplims

            # Plot
            plt.errorbar(
                lc["time_start"][sel],  # TODO: Move this to the bin center
                mags[sel],
                yerr=mags_err[sel],
                lolims=uplims[sel],
                color=col,
                markeredgecolor=col,
                markerfacecolor=mfc,
                marker=mar,
                label=src_lab,
                **plt_errorbar_kwargs,
            )
        ax.legend(loc=legend_loc, ncol=7, fontsize="small", handletextpad=0.05)
        ax.set_xlabel("MJD")
        ax.set_ylabel("Magnitude")

        if field_ids[0] == None:
            uniq_field_ids = [self.field_id]
        else:
            uniq_field_ids = np.unique(field_ids)
        fig_title = "Field id:"
        for fid in uniq_field_ids:
            fig_title += " " + str(fid)
        ax.set_title(fig_title, fontsize="small", loc="left")

        return ax
