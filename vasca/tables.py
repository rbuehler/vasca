#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
astropy.Table collection class for VASCA
"""

import os

import h5py
import numpy as np
from astropy import units as uu
from astropy.io import fits
from astropy.nddata import bitmask
from astropy.table import Table
from astropy.wcs import wcs
from loguru import logger
from scipy.stats import chi2, skew
from sklearn.cluster import MeanShift, estimate_bandwidth

from vasca.tables_dict import dd_vasca_tables
from vasca.utils import add_rg_src_id, table_to_array

# import warnings
# from astropy.io.fits.verify import VerifyWarning

# deactivate warnings
# warnings.simplefilter('ignore', category=VerifyWarning)

dimless = uu.dimensionless_unscaled

# global paths
# path to the dir. of this file
FILE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = FILE_DIR + "/../"  # path to the root directory of the repository


class TableCollection(object):
    def __init__(self):
        # Configure logger
        # adds the class name as an extra key; accessible vie the handler format
        logger.configure(extra={"classname": self.__class__.__name__})

        self._table_names = list()

    @staticmethod
    def table_from_template(dd_data, template_name):
        """
        Creates a new astropy table.

        Parameters
        ----------
        dd_data : dict
            Data dictionaty with the key corresponding to the templates columns.
            Id data==None an empty template table is returned.
        template_name : str
            Identifier to select a table template. Templates are selected by
            setting the class key and a corresponding table key in one string
            separated by a colon, e.g. template_name=<class_key>:<table_key>.

        Returns
        -------
        astropy.table.Table
        """

        # Check dd_data type
        if not (isinstance(dd_data, dict) or (dd_data is None)):
            logger.error(f"Passed data format not supported: {type(dd_data)}")

        # Takes pre-defined template dictionary
        templates = dd_vasca_tables

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

        # Get template just for the required table
        template_copy = templates[class_key][table_key].copy()

        # Check if data misses columns and in case fill with default values
        if dd_data is not None:
            len_data_cols = len(dd_data[list(dd_data.keys())[0]])
            for col in template_copy["names"]:
                if col not in dd_data.keys():
                    idx = template_copy["names"].index(col)
                    dd_data[col] = [template_copy["defaults"][idx]] * len_data_cols

        # Create table, delete defaults entry first,
        # as astropy Table does not support this
        del template_copy["defaults"]
        tt_out = Table(data=dd_data, **template_copy)
        # tt_out.meta["template"] = template_name

        # logging
        # Remove "defaults" from templates dictionary,
        # as astropy.Table does not support this
        logger.debug(f"Created new table from template '{template_name}'.")
        return tt_out

    def add_table(self, data, template_name):
        """
        Add a VASCA table to the field.

        Parameters
        ----------
        data : list, array-like
            Data of the table with shape (n, n_cols) or as dictionary with the
            key corresponding to the templates columns.
        template_name : str
            Identifier to select a table template. Templates are selected by
            setting the class key and a corresponding table key in one string
            separated by a colon, e.g. template_name=<class_key>:<table_key>.

        """
        logger.debug(f"Adding table '{template_name}'")

        table_key = template_name.split(":")[1]

        if table_key in self._table_names:
            logger.warning(f"Table '{table_key}' already exists, overwriting")
        self._table_names.append(table_key)

        tt = self.table_from_template(data, template_name)

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

        new_hdul = fits.HDUList([hdup])
        new_hdul.writeto(file_name, overwrite=overwrite, output_verify=fits_verify)

        for key in self._table_names:
            if key in self.__dict__:
                logger.debug(f"Writing table '{key}'")
                if not key.startswith("ta_"):
                    self.__dict__[key].write(file_name, append=True)
                # Write table with vector entries,
                # currently not supported by astropy.Table
                else:
                    cols = list()
                    for colname in self.__dict__[key].colnames:
                        coldata = self.__dict__[key][colname].data
                        coltype = str(self.__dict__[key][colname].dtype)

                        # Setup column fits format see
                        # https://heasarc.gsfc.nasa.gov/docs/software/fitsio/quick/node10.html
                        # https://docs.astropy.org/en/stable/io/fits/usage/unfamiliar.html
                        # TODO: Make this more general and read type from tables_dict
                        col_for = "PD()"
                        if (
                            colname == "mag"
                            or colname == "mag_err"
                            or colname == "ul"
                            or colname == "pos_err"
                            or colname == "time_bin_size"
                        ):
                            col_for = "PE()"

                        if "int" in coltype:
                            col_for = "K"
                        elif "float" in coltype:
                            col_for = "D"
                        elif "|S" in coltype:
                            col_for = coltype[2:] + "A"
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

    def load_from_fits(self, file_name):
        """
        Loads field from a fits file

        Parameters
        ----------
        file_name : str, optional
            File name.

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
                # Expect astropy handling of fits vector
                # to simplify this in the future
                if "rg_fd_id" in col_names:
                    self.add_table(ta_data, "region:" + ta_name)
                else:
                    self.add_table(ta_data, "base_field:" + ta_name)

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
                out_str += "\n" + key + ":\n" + self.__dict__[key].__str__() + "\n"

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
                        f"AND selecting '{var}' {vals}, "
                        f"kept: {100*sel.sum()/nr_sel : .4f}%"
                    )

            # Apply bitmask cuts
            if "bitmask" in selections.keys():
                for var, vals in selections["bitmask"].items():
                    no_art = sel * (tt[var].data.astype("int") == 0)
                    bit = bitmask.bitfield_to_boolean_mask(tt[var], ignore_flags=vals)
                    sel = sel * (no_art + bit)
                    logger.debug(
                        f"AND selecting bitmask '{var}' keep {vals}, "
                        f"kept: {100*sel.sum()/nr_sel : .4f}%"
                    )
        elif selections["sel_type"] == "or":
            sel_or = np.zeros(len(sel))

            if "range" in selections.keys():
                for var, vals in selections["range"].items():
                    sel_or = sel_or + (tt[var] >= vals[0]) * (tt[var] <= vals[1])
                    logger.debug(
                        f"OR selecting '{var}' {vals}, "
                        f"kept: {100*sel_or.sum()/nr_sel : .4f}%"
                    )
                sel = sel * sel_or
        elif selections["sel_type"] == "is_in":
            tt_ref = self.__dict__[selections["ref_table"]]
            sel = np.in1d(tt[selections["var"]], tt_ref[selections["var"]])
        else:
            logger.error("Unkown selection type.")

        sel = sel.astype("bool")
        tt.replace_column("sel", sel)

        if remove_unselected:
            tt = tt[sel]

        self.__dict__[table_name] = tt

    def get_light_curve(self, fd_src_ids=None, rg_src_ids=None):
        """
        Get a light curves for one or list of sources, for regions or fields.

        Parameters
        ----------
        fd_src_ids : list or int
            List or single field source IDs to plot. Default is None.
        rg_src_ids : list or int
            List or single region source IDs to plot. Default is None.

        Returns
        -------
        lc_dict : dict
            Dictionary {src_id : light_curve). Light curve as an astropy Table
                compatible with astropy BinnedTimeSeries.

        """

        logger.debug(
            f"Getting light curve for fd_src_ids '{fd_src_ids}' and rg_src_ids '{rg_src_ids}'"
        )

        # Setup input values depending on field or region
        src_ids = None
        ids_type = "none"
        if hasattr(rg_src_ids, "__iter__"):
            src_ids = list(rg_src_ids)
            ids_type = "rg_src_id"
        elif rg_src_ids is not None:
            src_ids = [rg_src_ids]
            ids_type = "rg_src_id"

        if hasattr(fd_src_ids, "__iter__"):
            src_ids = list(fd_src_ids)
            ids_type = "fd_src_id"
        elif fd_src_ids is not None:
            src_ids = [fd_src_ids]
            ids_type = "fd_src_id"

        logger.debug(f"Getting lightcurve for {ids_type}: {src_ids}")

        # Dictionary to store light curve tables
        lc_dict = dict()

        # Index field or region for source ID
        self.tt_detections.add_index(ids_type)
        self.tt_visits.add_index("vis_id")

        # Loop over sources and get lc info
        for src_id in src_ids:
            tt_det_src = self.tt_detections.loc[[src_id]]
            vis_idx = self.tt_visits.loc_indices["vis_id", tt_det_src["vis_id"]]
            tt_det_src["time_bin_start"] = self.tt_visits[vis_idx]["time_bin_start"]
            tt_det_src["time_bin_size"] = self.tt_visits[vis_idx]["time_bin_size"]
            tt_det_src.sort("time_bin_start")
            src_data = {
                "time_start": np.array(tt_det_src["time_bin_start"]),
                "time_delta": np.array(tt_det_src["time_bin_size"]),
                "mag": np.array(tt_det_src["mag"]),
                "mag_err": np.array(tt_det_src["mag_err"]),
                #                "ul": np.array(src_lc["ul"]),
            }
            # Create and store table
            tt_lc = self.table_from_template(src_data, "base_field:tt_source_lc")
            tt_lc.meta[ids_type] = src_id
            lc_dict[src_id] = tt_lc

        return lc_dict

    def cluster_meanshift(self, **ms_kw):
        """
        Apply _MeanShift clustering algorithm using to derive sources.

        .. _MeanShift: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MeanShift.html

        Parameters
        ----------
        ms_kw : dict, optional
            Keywords passed to the scikit MeanShift function. Note that the
            bandwidth is assumed to be in units of arc seconds.

        Returns
        -------
        int
            Number of detected clusters.
        """

        # Select seed table, detections or sources
        clus_srcs = False
        tt_in = self.tt_detections
        if "Region" in self.__class__.__name__:
            clus_srcs = True
            tt_in = self.tt_sources
        logger.info(f"Clustering with MeanShift with clus_srcs: {clus_srcs}")

        # Selection
        sel = tt_in["sel"]

        # Get detection coordinates and run clustering
        coords = table_to_array(tt_in[sel]["ra", "dec"])

        # Do bandwidth determination "by hand" to print it out and convert
        # bandwidth unit from arc seconds into degerees
        dd_ms = ms_kw
        if "bandwidth" not in ms_kw or ms_kw["bandwidth"] is None:
            logger.debug("Estimating bandwidth")
            dd_ms["bandwidth"] = estimate_bandwidth(coords, quantile=0.2, n_samples=500)
        else:
            dd_ms["bandwidth"] = (ms_kw["bandwidth"] * uu.arcsec).to(uu.deg).value

        logger.debug(f"MeanShift with parameters (bandwith in degrees): '{dd_ms}' ")
        ms = MeanShift(**dd_ms)

        ms.fit(coords)

        # Fill in data into source tables
        src_ids, det_cts = np.unique(ms.labels_, return_counts=True)
        cluster_centers = ms.cluster_centers_
        nr_srcs = len(cluster_centers)

        # Store clusters in tt_source table
        if clus_srcs:
            srcs_data = {
                "rg_src_id": src_ids,
                "ra": cluster_centers[:, 0],
                "dec": cluster_centers[:, 1],
                "nr_fds": det_cts,
                "nr_uls": np.zeros(nr_srcs),
            }

            # Add rg_src_id to detections
            tt_in["rg_src_id"][sel] = ms.labels_
            add_rg_src_id(tt_in[sel], self.tt_detections)

            # Remove existing table and add new one
            del self.__dict__["tt_sources"]
            self._table_names.remove("tt_sources")
            self.add_table(srcs_data, "region:tt_sources")

            # Add total number of detections
            det_src_ids, src_nr_det = np.unique(
                self.tt_detections["rg_src_id"].data,
                return_counts=True,
            )
            self.tt_sources.add_index("rg_src_id")
            src_idx = self.tt_sources.loc_indices["rg_src_id", det_src_ids]
            self.tt_sources["nr_det"][src_idx] = src_nr_det

            # Log percentage of merged sources
            nr_merged = len(tt_in[sel]) - len(src_ids)
            perc_merged = np.round(100 * nr_merged / len(tt_in[sel]), 4)
            logger.debug(f"Merged sources: {nr_merged} ({perc_merged}%)")

        else:
            srcs_data = {
                "fd_src_id": src_ids,
                "ra": cluster_centers[:, 0],
                "dec": cluster_centers[:, 1],
                "nr_det": det_cts,
                "nr_uls": np.zeros(nr_srcs),
            }

            # Fill information into tables.
            self.add_table(srcs_data, "base_field:tt_sources")

            # Update fd_src_id entries
            self.tt_detections["fd_src_id"][sel] = ms.labels_

            # Fill light curve data into tables
            self.remove_double_visit_detections()

        self.tt_sources.meta["CLUSTALG"] = "MeanShift"

        return nr_srcs

    # Calculation is done on a field level for easy numpy paralleization,
    # as this calculation is computationally intensive.
    def set_src_stats(self):
        """
        Calculates source parameters from detections and stores them
        in the source table (tt_source).

        Returns
        -------
        None.

        """

        logger.debug("Calculating source statistics.")

        # Check if field or region
        id_name = "fd_src_id"
        if "Region" in self.__class__.__name__:
            id_name = "rg_src_id"

        # Prepare detection data
        tt_det = Table(self.tt_detections, copy=True)
        tt_det.sort(id_name)

        # Add flux and flux_err to table
        tt_det["flux"] = (np.array(tt_det["mag"]) * uu.ABmag).physical
        tt_det["flux_err"] = (
            (np.array(tt_det["mag"]) - np.array(tt_det["mag_err"])) * uu.ABmag
        ).physical - tt_det["flux"]

        src_ids, src_indices, src_nr_det = np.unique(
            tt_det[id_name].data, return_index=True, return_counts=True
        )

        # Buffer input data
        ll_bvar = ["mag", "mag_err", "pos_err", "ra", "dec", "flux", "flux_err"]
        dd_bvar = dict()
        for bvar in ll_bvar:
            dd_bvar[bvar] = tt_det[bvar].data

        # Prepare output data arrays
        ll_svar = [
            "ra",
            "dec",
            "pos_err_mean",
            "pos_err",
            "pos_var",
            "mag_mean",
            "mag_err_mean",
            "mag_var",
            "flux_cpval",
            "flux_rchiq",
            "mag_dmax",
            "mag_dmax_abs",
            "mag_dmax_sig",
            "mag_sig_int2",
            "mag_skew",
        ]

        dd_svar = dict()
        for svar in ll_svar:
            dd_svar[svar] = np.zeros(len(src_ids))

        # Do loop and calculate stats variables
        for isrc, srcid in enumerate(src_ids):
            idx1 = src_indices[isrc]
            idx2 = idx1 + src_nr_det[isrc]

            # Mean magnitude
            dd_svar["mag_mean"][isrc] = np.mean(dd_bvar["mag"][idx1:idx2])
            dd_svar["mag_err_mean"][isrc] = np.mean(dd_bvar["mag_err"][idx1:idx2])
            dd_svar["pos_err_mean"][isrc] = np.mean(dd_bvar["pos_err"][idx1:idx2])

            # Variances
            if src_nr_det[isrc] > 1:
                dd_svar["mag_var"][isrc] = np.var(dd_bvar["mag"][idx1:idx2], ddof=1)
                dd_svar["pos_var"][isrc] = np.var(
                    dd_bvar["ra"][idx1:idx2], ddof=1
                ) + np.var(dd_bvar["dec"][idx1:idx2], ddof=1)
            else:
                dd_svar["mag_var"][isrc] = -1
                dd_svar["pos_var"][isrc] = -1

            # Weighted averages
            pos_weight = 1.0 / dd_bvar["pos_err"][idx1:idx2] ** 2
            dd_svar["ra"][isrc] = np.average(
                dd_bvar["ra"][idx1:idx2], weights=pos_weight
            )
            dd_svar["dec"][isrc] = np.average(
                dd_bvar["dec"][idx1:idx2], weights=pos_weight
            )

            dd_svar["pos_err"][isrc] = np.sqrt(1.0 / np.sum(pos_weight))

            # Chi square
            flux_mean = np.mean(dd_bvar["flux"][idx1:idx2])
            chiq_el = np.power(dd_bvar["flux"][idx1:idx2] - flux_mean, 2) / np.power(
                dd_bvar["flux_err"][idx1:idx2], 2
            )
            chiq = np.sum(chiq_el)
            if src_nr_det[isrc] > 1:
                dd_svar["flux_rchiq"][isrc] = chiq / (src_nr_det[isrc] - 1)
                dd_svar["flux_cpval"][isrc] = chi2.sf(chiq, src_nr_det[isrc] - 1)

                # Intrinsic variability
                dd_svar["mag_sig_int2"][isrc] = (
                    np.var(dd_bvar["mag"][idx1:idx2], ddof=1)
                    - np.mean(dd_bvar["mag_err"][idx1:idx2]) ** 2
                )
            else:
                dd_svar["flux_rchiq"][isrc] = -1.0
                dd_svar["flux_cpval"][isrc] = -1.0
                dd_svar["mag_sig_int2"][isrc] = -1.0

            # Skewness
            if src_nr_det[isrc] > 2:
                dd_svar["mag_skew"][isrc] = skew(dd_bvar["mag"][idx1:idx2], bias=False)
            else:
                dd_svar["mag_skew"][isrc] = -100.0

            # Maximum variation significance
            dmag = np.abs(dd_bvar["mag"][idx1:idx2] - dd_svar["mag_mean"][isrc])
            dmag_sig = dmag / dd_bvar["mag_err"][idx1:idx2]
            dmax_sig_max_idx = np.argmax(dmag_sig)

            dd_svar["mag_dmax_sig"][isrc] = dmag_sig[dmax_sig_max_idx]
            dd_svar["mag_dmax"][isrc] = dmag[dmax_sig_max_idx]

            # Maximum absolute variation
            dmag_abs_max = dmag.max()
            dd_svar["mag_dmax_abs"][isrc] = dmag_abs_max

        # Write them into tt_sources
        self.tt_sources.add_index(id_name)
        src_idx = self.tt_sources.loc_indices[id_name, src_ids]  # src_ids]

        for svar in ll_svar:
            self.tt_sources[svar][src_idx] = dd_svar[svar]  # tt_mean["mag"]
