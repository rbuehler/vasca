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
from astropy.coordinates import SkyCoord

from vasca.tables_dict import dd_vasca_tables
from vasca.utils import add_rg_src_id, table_to_array, flux2mag

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
        hdulist = [fits.PrimaryHDU()]

        # Check if image data is set and add to primary HDU
        if hasattr(self, "ref_img"):
            if self.ref_img is not None and self.ref_wcs is not None:
                logger.debug(f"Storing image data of shape {self.ref_img.shape}")
                hdulist.append(
                    fits.CompImageHDU(
                        self.ref_img,
                        header=self.ref_wcs.to_header(),
                        name="ref_img",
                    )
                )
                if hasattr(self, "vis_img"):
                    hdulist.append(
                        fits.CompImageHDU(
                            data=self.vis_img,
                            header=self.ref_wcs.to_header(),
                            name="vis_img",
                        )
                    )

        new_hdul = fits.HDUList(hdulist)
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
                            colname == "flux"
                            or colname == "flux_err"
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
        ext_nr = len(hdulist)  # First extensions reserved for images
        with fits.open(file_name, "update", output_verify=fits_verify) as ff:
            for key in self._table_names:
                if key in self.__dict__:
                    ff[ext_nr].header["EXTNAME"] = key
                    ff.flush(output_verify=fits_verify)
                    ext_nr += 1
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
            img_names = [
                hdu.header["EXTNAME"]
                for hdu in ff[1:]
                if hdu.header["EXTNAME"].lower().endswith("_img")
            ]
            if hasattr(self, "ref_img"):
                for img_name in img_names:
                    setattr(self, img_name.lower(), ff[img_name].data)
                    if img_name.lower() == "ref_img":
                        setattr(self, "ref_wcs", wcs.WCS(ff[img_name].header))

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

        # Get selected events
        presel = tt["sel"].data.astype("bool")
        nr_presel = presel.sum()

        if selections["sel_type"] == "and":
            sel = np.ones(len(tt), dtype=bool)

            # Apply min/max cuts
            if "range" in selections.keys():
                for var, vals in selections["range"].items():
                    sel = sel * (tt[var] >= vals[0]) * (tt[var] <= vals[1])
                    logger.info(
                        f"AND selecting '{var}' {vals}, "
                        f"kept: {100*sel.sum()/nr_presel : .4f}%"
                    )

            # Apply bitmask cuts
            if "bitmask" in selections.keys():
                for var, vals in selections["bitmask"].items():
                    # See https://docs.astropy.org/en/stable/api/astropy.nddata.bitfield_to_boolean_mask.html
                    bit = bitmask.bitfield_to_boolean_mask(
                        tt[var],
                        flip_bits=True,
                        dtype=bool,
                        good_mask_value=0,
                        ignore_flags=vals,
                    )
                    sel = sel * ~bit
                    logger.info(
                        f"AND selecting bitmask '{var}' removing {vals}, "
                        f"kept: {100*sel.sum()/nr_presel : .4f}%"
                    )
        elif selections["sel_type"] == "or":
            sel = np.zeros(len(sel), dtype=bool)

            if "range" in selections.keys():
                for var, vals in selections["range"].items():
                    sel = sel + (tt[var] >= vals[0]) * (tt[var] <= vals[1])
                    logger.info(
                        f"OR selecting '{var}' {vals}, "
                        f"kept: {100*sel.sum()/nr_presel : .4f}%"
                    )
        elif selections["sel_type"] == "is_in":
            tt_ref = self.__dict__[selections["ref_table"]]
            sel = np.in1d(tt[selections["var"]], tt_ref[selections["var"]])
        else:
            logger.error("Unkown selection type.")

        # combine with preselection
        sel_tot = presel * sel
        if selections["presel_type"].casefold() == "or".casefold():
            sel_tot = presel + sel

        tt.replace_column("sel", sel_tot.astype("bool"))

        if remove_unselected:
            tt = tt[sel_tot]

        # Set minimum and maximum values if asked
        if "set_range" in selections.keys():
            for var, vals in selections["set_range"].items():
                logger.info(f"Setting range for '{var}' to '{vals}'")
                tt[var][:] = (tt[var] >= vals[0]) * tt[var] + (
                    tt[var] < vals[0]
                ) * vals[
                    0
                ]  # Set minimum
                tt[var][:] = (tt[var] <= vals[1]) * tt[var] + (
                    tt[var] > vals[1]
                ) * vals[
                    1
                ]  # Set maximum

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
            if not hasattr(vis_idx, "__iter__"):
                logger.warning(
                    f"Found only one light curve points, vis_idx: {vis_idx}, src_id: {src_id}"
                )
                vis_idx = np.array([vis_idx])
                tt_det_src = Table(tt_det_src)
            tt_det_src["time_bin_start"] = self.tt_visits[vis_idx]["time_bin_start"]
            tt_det_src["time_bin_size"] = self.tt_visits[vis_idx]["time_bin_size"]
            tt_det_src.sort("time_bin_start")
            src_data = {
                "time_start": np.array(tt_det_src["time_bin_start"]),
                "time_delta": np.array(tt_det_src["time_bin_size"]),
                "flux": np.array(tt_det_src["flux"]),
                "flux_err": np.array(tt_det_src["flux_err"]),
                #                "ul": np.array(src_lc["ul"]),
            }
            # Create and store table
            tt_lc = self.table_from_template(src_data, "base_field:tt_source_lc")
            tt_lc.meta[ids_type] = src_id
            lc_dict[src_id] = tt_lc

        return lc_dict

    def cluster_meanshift(self, **ms_kw):
        """
        Apply _MeanShift clustering algorithm using to derive sources. Runs only on selected
        detections or sources.

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
        table_name = ms_kw["table_name"]
        del ms_kw["table_name"]
        tt_in = self.__dict__[table_name]
        logger.info(f"Clustering with MeanShift with table: {table_name}")

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
        src_ids, clu_cts = np.unique(ms.labels_, return_counts=True)
        cluster_centers = ms.cluster_centers_
        nr_srcs = len(cluster_centers)

        # How to store clusters depends on table, check cases
        # Store clusters in tt_source table
        if table_name == "tt_sources":
            srcs_data = {
                "rg_src_id": src_ids,
                "ra": cluster_centers[:, 0],
                "dec": cluster_centers[:, 1],
                "nr_fd_srcs": clu_cts,
                "filter_id": list(),
            }

            # Add rg_src_id to detections
            tt_in["rg_src_id"][sel] = ms.labels_
            add_rg_src_id(tt_in[sel], self.tt_detections)

            # Add filter info
            # TODO: Write this vecotrized for speed increase
            tt_in.add_index("rg_src_id")
            for rg_src_id in src_ids:
                idxs = tt_in.loc_indices["rg_src_id", rg_src_id]
                filter_ids = np.unique(tt_in["filter_id"][idxs].data)
                srcs_data["filter_id"].append(filter_ids.sum())

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

        elif table_name == "tt_detections":

            # Assume all detections have the same filter_id
            # TODO: This might have to be extendet to the general case in the future.
            filter_id = self.tt_detections["filter_id"][0]

            srcs_data = {
                "fd_src_id": src_ids,
                "ra": cluster_centers[:, 0],
                "dec": cluster_centers[:, 1],
                "nr_det": clu_cts,
                "filter_id": np.zeros(nr_srcs) + filter_id,
            }

            # Fill information into tables.
            self.add_table(srcs_data, "base_field:tt_sources")

            # Update fd_src_id entries
            self.tt_detections["fd_src_id"][np.where(sel)] = ms.labels_

            # Remove two detections in the same visit (keep the closer one)
            self.remove_double_visit_detections()

        elif table_name == "tt_coadd_detections":
            coadd_data = {
                "coadd_src_id": src_ids,
                "ra": cluster_centers[:, 0],
                "dec": cluster_centers[:, 1],
                "nr_det": clu_cts,
                "filter_id": list(),
            }

            # Add source label to coadd detections
            self.tt_coadd_detections["coadd_src_id"][np.where(sel)] = ms.labels_

            # Add filter id
            self.tt_coadd_detections.add_index("coadd_src_id")
            for coadd_src_id in src_ids:
                idxs = self.tt_coadd_detections.loc_indices[
                    "coadd_src_id", coadd_src_id
                ]
                filter_ids = np.unique(self.tt_coadd_detections["filter_id"][idxs].data)
                coadd_data["filter_id"].append(filter_ids.sum())

            # Fill information into tables.
            self.add_table(coadd_data, "region:tt_coadd_sources")

        else:
            logger.error("Unkown table name")

        self.tt_sources.meta["CLUSTALG"] = "MeanShift"

        return nr_srcs

    # Calculation is done on a field level for easy numpy paralleization,
    # as this calculation is computationally intensive.
    def set_src_stats(self, src_id_name="fd_src_id"):
        """
        Calculates source parameters from detections and stores them
        in the source table (tt_source).

        Returns
        -------
        None.

        """

        logger.debug("Calculating source statistics.")

        def get_var_stat(vals, vals_err):
            "Helper function to calculate error weighted average mean and its error"
            rr = {}
            wght = 1.0 / vals_err**2
            rr["wght_mean"] = np.average(vals, weights=wght)
            rr["wght_mean_err"] = np.sqrt(1.0 / np.sum(wght))
            chiq_el = np.power(vals - rr["wght_mean"], 2) / np.power(vals_err, 2)
            chiq = np.sum(chiq_el)
            nr_vals = len(vals)

            if nr_vals > 1:
                rr["var"] = np.var(vals, ddof=1)
                rr["nxv"] = (rr["var"] - np.mean(vals_err**2)) / (
                    rr["wght_mean"] * rr["wght_mean"]
                )
                rr["rchiq"] = chiq / (nr_vals - 1)
                rr["cpval"] = chi2.sf(chiq, nr_vals - 1)
            else:
                rr["var"] = rr["nxv"] = -100
                rr["rchiq"] = rr["cpval"] = -1.0

            return rr

        # Setup table names
        if src_id_name == "fd_src_id" or src_id_name == "rg_src_id":
            tt_det_name = "tt_detections"
            tt_src_name = "tt_sources"
        else:
            tt_det_name = "tt_coadd_detections"
            tt_src_name = "tt_coadd_sources"

        # Prepare detection data
        sel_det = self.__dict__[tt_det_name]["sel"]
        tt_det = Table(self.__dict__[tt_det_name][sel_det], copy=True)
        tt_det.sort([src_id_name, "filter_id"])

        # Get src_ids, src_index and det_nr
        src_ids, src_indices, src_nr_det = np.unique(
            tt_det[src_id_name].data, return_index=True, return_counts=True
        )

        # Check which filter ids are present
        filter_ids = np.sort(np.unique(tt_det["filter_id"].data))
        nr_filters = len(filter_ids)

        # Buffer input data for speed and convert position errors to degrees
        dd_det_var = {"pos_err_deg": tt_det["pos_err"].data.astype(np.float64) / 3600.0}
        ll_det_var = ["flux", "flux_err", "ra", "dec", "filter_id"]
        for bvar in ll_det_var:
            dd_det_var[bvar] = tt_det[bvar].data

        ll_src_var = [
            "ra",
            "dec",
            "pos_err",
            "pos_nxv",
            "pos_var",
            "pos_cpval",
            "pos_rchiq",
            "flux",
            "flux_err",
            "flux_nxv",
            "flux_var",
            "flux_cpval",
            "flux_rchiq",
            "nr_det",
        ]

        # Dictionary to store calculated source parameters
        dd_src_var = dict()
        for svar in ll_src_var:
            dd_src_var[svar] = list()  # np.zeros(len(src_ids))

        # Do loop over all sources and calculate stats variables
        for isrc, srcid in enumerate(src_ids):
            idx1 = src_indices[isrc]
            idx2 = idx1 + src_nr_det[isrc]

            # Position variables
            rr_ra = get_var_stat(
                dd_det_var["ra"][idx1:idx2],
                dd_det_var["pos_err_deg"][idx1:idx2],
            )
            rr_dec = get_var_stat(
                dd_det_var["dec"][idx1:idx2],
                dd_det_var["pos_err_deg"][idx1:idx2],
            )
            dd_src_var["ra"].append(rr_ra["wght_mean"])
            dd_src_var["dec"].append(rr_dec["wght_mean"])
            dd_src_var["pos_err"].append(
                ((rr_ra["wght_mean_err"] + rr_dec["wght_mean_err"]) * 3600 / 2.0)
            )
            dd_src_var["pos_nxv"].append((rr_ra["nxv"] + rr_dec["nxv"]) / 2.0)
            dd_src_var["pos_var"].append(
                ((rr_ra["var"] + rr_dec["var"]) * (3600**2) / 2.0)
            )
            dd_src_var["pos_cpval"].append(rr_ra["cpval"] * rr_dec["cpval"])
            dd_src_var["pos_rchiq"].append((rr_ra["rchiq"] + rr_dec["rchiq"]) / 2.0)

            # Check at what index filter changes and analyse separatelly
            idxfs = np.where(
                np.diff(dd_det_var["filter_id"][idx1:idx2], prepend=np.nan)
            )[0]
            idxfs = np.append(idxfs, idx2 - idx1)

            # Create empty array
            # TODO: bring defaults, etc automatically in synch with tables_dict.py info
            aa_zero = np.zeros(nr_filters, dtype=np.float32)
            dd_src_var["flux"].append(np.copy(aa_zero) - 1.0)
            dd_src_var["flux_err"].append(np.copy(aa_zero) - 1.0)
            dd_src_var["flux_nxv"].append(np.copy(aa_zero) - 100.0)
            dd_src_var["flux_var"].append(np.copy(aa_zero) - 1.0)
            dd_src_var["flux_cpval"].append(np.copy(aa_zero) - 1.0)
            dd_src_var["flux_rchiq"].append(np.copy(aa_zero) - 1.0)
            dd_src_var["nr_det"].append(np.zeros(nr_filters, dtype=np.int32) - 1)

            # Loop over all filters
            for ii in range(len(idxfs) - 1):
                rr_flux = get_var_stat(
                    dd_det_var["flux"][idx1 + idxfs[ii] : idx2 + idxfs[ii + 1]],
                    dd_det_var["flux_err"][idx1 + idxfs[ii] : idx2 + idxfs[ii + 1]],
                )
                filter_id = dd_det_var["filter_id"][idx1 + idxfs[ii]]
                filter_nr = np.where(filter_ids == filter_id)
                dd_src_var["flux"][-1][filter_nr] = rr_flux["wght_mean"]
                dd_src_var["flux_err"][-1][filter_nr] = rr_flux["wght_mean_err"]
                dd_src_var["flux_nxv"][-1][filter_nr] = rr_flux["nxv"]
                dd_src_var["flux_var"][-1][filter_nr] = rr_flux["var"]
                dd_src_var["flux_cpval"][-1][filter_nr] = rr_flux["cpval"]
                dd_src_var["flux_rchiq"][-1][filter_nr] = rr_flux["rchiq"]
                dd_src_var["nr_det"][-1][filter_nr] = idxfs[ii + 1] - idxfs[ii]

        for svar in ll_src_var:
            self.__dict__[tt_src_name].replace_column(svar, dd_src_var[svar])

    def cross_match(
        self,
        tt_cat,
        cat_id_name="coadd_src_id",
        table_name="tt_sources",
        dist_max=1 * uu.arcsec,
        dist_s2n_max=3,
    ):
        """
        Cross match sources to a catalog.

        Parameters
        ----------
        tt_cat : astropy.Table
            Catalog table. Has to contain "ra","dec" (in deg), "flux" (in microJy)
            and "cat_id_name" columns.  Marks associated catalog sources
            in the "sel" column of the catalog table, if it exists.
        cat_id_name : str, optional
            Catalog ID Br. variable name. The default is "coadd_src_id".
        table_name : str, optional
            Table to crossmatch to catalog. The default is "tt_sources".
        dist_max astropy.Quantity
            Maximum angular distance under which all associations are done, independent
            of dist_s2n. The default is "1 arcsec".
        dist_s2n_max float
            Maximum distance in units of position error. All sources below this cut are
            associated, independently of the dist_max selection.

        Returns
        -------
        None.

        """

        logger.debug(f"Cross matching table {table_name}")

        # Get source positions
        tt_srcs = self.__dict__[table_name]
        pos_srcs = SkyCoord(
            ra=tt_srcs["ra"], dec=tt_srcs["dec"], unit="deg", frame="icrs"
        )
        pos_cat = SkyCoord(ra=tt_cat["ra"], dec=tt_cat["dec"], unit="deg", frame="icrs")

        idx_cat, dist_cat, _ = pos_srcs.match_to_catalog_sky(pos_cat)

        sel_dist = dist_cat.to("arcsec") < dist_max.to("arcsec")

        # Check how compatible positions are within errors
        sigma_dist = np.sqrt(
            tt_srcs["pos_err"] ** 2 + tt_cat["pos_err"][idx_cat] ** 2
        )  # Convert to 68% containment radius for 2D gaussian
        sel_dist_s2n = (dist_cat.to("arcsec") / sigma_dist.to("arcsec")) < dist_s2n_max

        sel = sel_dist + sel_dist_s2n

        tt_srcs["assoc_id"][sel] = tt_cat[idx_cat[sel]][cat_id_name]
        tt_srcs["assoc_dist"][sel] = dist_cat[sel].to("arcsec")
        tt_srcs["assoc_ffactor"][sel] = (
            tt_srcs["flux"][sel] / tt_cat[idx_cat[sel]]["flux"]
        )
        flux_diff = tt_srcs["flux"][sel] - tt_cat[idx_cat[sel]]["flux"]
        flux_diff_err = np.sqrt(
            tt_srcs["flux_err"][sel] ** 2 + tt_cat[idx_cat[sel]]["flux_err"] ** 2
        )
        tt_srcs["assoc_fdiff_s2n"][sel] = flux_diff / flux_diff_err

        # Mark as selected, if selection column exists
        if "sel" in tt_cat.colnames:
            tt_cat["sel"][:] = np.zeros(len(tt_cat), dtype=bool)
            tt_cat["sel"][idx_cat[sel]] = np.ones(len(idx_cat[sel]), dtype=bool)
