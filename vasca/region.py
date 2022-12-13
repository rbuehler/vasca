#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import healpy as hpy
import numpy as np
from astropy import units as uu
from loguru import logger
from sklearn.cluster import MeanShift, estimate_bandwidth

from vasca.field import BaseField, GALEXField
from vasca.tables import TableCollection
from vasca.tables_dict import dd_vasca_tables
from vasca.utils import table_to_array


class Region(TableCollection):
    """
    `~vasca.Region` defines a region in the sky as a
    list of vasca.field objects. It provides functionality to
    loop over fields to derive source lists, etc.
    """

    def __init__(self):
        """

        Notes
        -----
        Many class attributes are stored in astropy.table.Tables_. To see a
        description of each of their columns run :meth: `~vasca.Regions.info`.

        .. _astropy.table.Tables: https://docs.astropy.org/en/stable/api/astropy.table.Table.html

        Returns
        -------
        None.

        """
        # Sets skeleton
        super().__init__()

        # Setup empty tables to fill
        self.add_table(None, "region:tt_fields")

        self.fields = {}  # dictionary of field IDs and objects

    @classmethod
    def load_from_config(cls, vasca_cfg):
        """
        Loads region from configuration from dictionary

        Parameters
        ----------
        vasca_cfg : dict
            Dictionary with region parameters derived from the vasca pipeline
            YAML configuration file.

        Returns
        -------
        None.

        """

        rg = cls()

        logger.debug("Loading fields from config file")

        for obs in vasca_cfg["observations"]:

            if obs["observatory"] == "GALEX":

                # Loop over fields and store info
                rg_fd_id = 0
                for gfield_id in obs["obs_field_ids"]:

                    gf = GALEXField.load(
                        gfield_id,
                        obs_filter=obs["obs_filter"],
                        method=vasca_cfg["ressources"]["load_method"],
                        load_products=vasca_cfg["ressources"]["load_products"],
                        **vasca_cfg["ressources"]["field_kwargs"],
                    )
                    field_info = dict(gf.tt_field[0])
                    field_info["fov_diam"] = 1.10
                    field_info["nr_vis"] = gf.nr_vis
                    field_info["time_bin_size_sum"] = gf.time_bin_size_sum
                    field_info["time_start"] = gf.time_start.mjd
                    field_info["time_stop"] = gf.time_stop.mjd
                    field_info["rg_fd_id"] = rg_fd_id
                    rg.tt_fields.add_row(field_info)

                    if vasca_cfg["ressources"]["load_products"]:
                        rg.fields[rg_fd_id] = gf
                    else:
                        rg.fields[rg_fd_id] = None
                    rg_fd_id += 1
            else:
                logger.waring(
                    "Selected observatory `" + obs["observatory"] + "` not supported"
                )

        return rg

    def add_table_from_fields(self, table_name, only_selected=False):
        """
        Add tables from the fields to the region by stacking them,
        adding the field_id column.


        Parameters
        ----------
        table_name : str
            Table to be added
        only_selected : bool, optional
            Only selected rows from the field tables are copied over to the region.
            The default is False.

        Returns
        -------
        None.

        """

        logger.debug(f"Adding table from fields: {table_name}")
        if table_name in self._table_names:
            logger.warning(f"Table '{table_name}' already exists, overwriting")

        # Loop over fields and add field_id column and field id table
        ll_tt = []  # List of "table_name" tables for all fields
        for rg_fd_id, field in self.fields.items():
            tt = field.__dict__[table_name]

            # Apply row selection
            sel = np.ones(len(tt), dtype="bool")
            if only_selected:
                sel = tt["sel"]
            tt_sel = tt[sel]
            tt_sel["rg_fd_id"] = len(tt_sel) * [rg_fd_id]
            ll_tt.append(tt_sel)

        # colnames = dd_vasca_tables["region"][table_name]["names"]
        colnames = list(
            set(dd_vasca_tables["region"][table_name]["names"]) & set(ll_tt[0].colnames)
        )

        # Create empty data structure and then fill it with field tables
        dd_data = dict(zip(colnames, [list() for ii in range(len(colnames))]))
        for tt in ll_tt:
            for colname in colnames:
                dd_data[colname].extend(tt[colname].tolist())

        # For vector columns convert to numpy arrays of type object_
        # This is needed for correct writing to fits in Astropy v5.0.4
        for colname in colnames:
            if len(np.array(dd_data[colname], dtype=object).shape) > 1:
                dd_data[colname] = np.array(dd_data[colname], dtype=np.object_)

        self.add_table(dd_data, "region:" + table_name)

    def add_coverage_hp(self, nside=4096):

        npix = hpy.nside2npix(nside)
        pix_diam = hpy.nside2resol(nside, arcmin=True) / 60 * uu.deg

        logger.debug(
            f"Healpix NSIDE: {nside}, NPIX: {npix}, pixel diameter: {pix_diam}"
        )
        pix_nrs = np.arange(npix)
        hp_vis = np.zeros(npix, dtype="float32")
        hp_exp = np.zeros(npix, dtype="float32")
        for field in self.tt_fields:
            pos_vec = hpy.ang2vec(field["ra"], field["dec"], lonlat=True)
            rdisc = field["fov_diam"] / 2.0
            # TODO: Here a more general querry_polygon will have to be done for ULTRASAT
            ipix_disc = hpy.query_disc(
                nside=nside, vec=pos_vec, radius=np.radians(rdisc)
            )

            hp_vis[ipix_disc] += field["nr_vis"]
            hp_exp[ipix_disc] += field["time_bin_size_sum"]

        # Write to table
        sel_pix = hp_vis > 0
        keys_data = ["pix_id", "nr_vis", "exp"]
        ll_data = [pix_nrs[sel_pix], hp_vis[sel_pix], hp_exp[sel_pix]]
        dd_data = dict(zip(keys_data, ll_data))
        self.add_table(dd_data, "region:tt_coverage_hp")
        self.tt_coverage_hp.meta["NSIDE"] = nside

        hp_exp[hp_exp < 1e-6] = hpy.UNSEEN
        hp_vis[hp_vis < 1e-6] = hpy.UNSEEN

        return hp_vis, hp_exp

    def load_from_fits(self, file_name, load_fields=True):
        """
        Loads field from a fits file

        Parameters
        ----------
        file_name : str, optional
            Region file name. The default is "field_default.fits".
        load_fields : bool,
            Load the fields, which have to be located as fits in the subfolder "fields"
            of the region file in "file_name". Default is True.

        Returns
        -------
        None.

        """

        # Load file from TableCollection base class
        super().load_from_fits(file_name)
        region_path = os.path.dirname(file_name)

        # Load fields
        if load_fields:
            for ff in self.tt_fields:
                fd = BaseField()
                fd.load_from_fits(
                    region_path + "/fields/field_" + ff["field_id"] + ".fits"
                )
                self.fields[ff["field_id"]] = fd


# def cluster_srcs_meanshift(self, **ms_kw):
#     """
#     Apply _MeanShift clustering algorithm using to sources for overlapping fields.

#     .. _MeanShift: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MeanShift.html

#     Parameters
#     ----------
#     ms_kw : dict, optional
#         Keywords passed to the scikit MeanShift function. Note that the
#         bandwidth is assumed to be in units of arc seconds.

#     Returns
#     -------
#     int
#         Number of detected clusters.
#     """
#     logger.info("Clustering sources with MeanShift")

#     # Selection
#     sel = self.tt_sources["sel"]

#     # Get detection coordinates and run clustering
#     coords = table_to_array(self.tt_sources[sel]["ra", "dec"])

#     # Do bandwidth determination "by hand" to print it out and convert
#     # bandwidth unit from arc seconds into degerees
#     dd_ms = ms_kw
#     if "bandwidth" not in ms_kw or ms_kw["bandwidth"] is None:
#         logger.debug("Estimating bandwidth")
#         dd_ms["bandwidth"] = estimate_bandwidth(coords, quantile=0.2, n_samples=500)
#     else:
#         dd_ms["bandwidth"] = (ms_kw["bandwidth"] * uu.arcsec).to(uu.deg).value

#     logger.debug(f"MeanShift with parameters (bandwith in degrees): '{dd_ms}' ")
#     ms = MeanShift(**dd_ms)

#     ms.fit(coords)

#     # Fill in data into field tables
#     fd_src_ids, det_cts = np.unique(ms.labels_, return_counts=True)

#     cluster_centers = ms.cluster_centers_
#     nr_srcs = len(cluster_centers)
#     srcs_data = {
#         "fd_src_id": fd_src_ids,
#         "ra": cluster_centers[:, 0],
#         "dec": cluster_centers[:, 1],
#         "nr_det": det_cts,
#         "nr_uls": np.zeros(nr_srcs),
#     }

#     # Fill information into tables.
#     self.add_table(srcs_data, "base_field:tt_sources")
#     self.tt_sources.meta["CLUSTALG"] = "MeanShift"

#     # Update fd_src_id entries
#     self.tt_detections["fd_src_id"][sel] = ms.labels_

#     # Fill light curve data into tables
#     self.remove_double_visit_detections()

#     return nr_srcs
