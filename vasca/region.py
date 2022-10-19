#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from loguru import logger

import numpy as np
import matplotlib.pyplot as plt
import healpy as hpy
from astropy.table import Column, vstack
from astropy import units as uu

from vasca import tables
from vasca.field import GALEXField
from vasca.tables import TableCollection
from vasca.tables_dict import dd_vasca_tables


class Region(TableCollection):
    """
    :class: `~vasca.Region` defines a region in the sky as a
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

        if vasca_cfg["observations"]["observatory"] == "GALEX":

            # Loop over fields and store info
            for field_id in vasca_cfg["observations"]["field_ids"]:

                gf = GALEXField.load(
                    field_id,
                    obs_filter=vasca_cfg["observations"]["obs_filter"],
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
                rg.tt_fields.add_row(field_info)

                if vasca_cfg["ressources"]["load_products"]:
                    rg.fields[field_id] = gf
                else:
                    rg.fields[field_id] = None
        else:
            logger.waring(
                "Selected observatory `"
                + vasca_cfg["observations"]["observatory"]
                + "` not supported"
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
        for field_id, field in self.fields.items():
            tt = field.__dict__[table_name]

            # Apply row selection
            sel = np.ones(len(tt), dtype="bool")
            if only_selected:
                sel = tt["sel"]
            tt_sel = tt[sel]
            tt_sel["field_id"] = np.ones(len(tt_sel), dtype="int64") * field_id
            ll_tt.append(tt_sel)

        colnames = dd_vasca_tables["region"][table_name]["names"]

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

    def plot_sky_gnomeview(
        self, ra, dec, sel_srcs=True, gw_kwargs=None, ps_kwargs=None
    ):
        """
        Plot the nr visits and optionally (selected) sources (optinally) on the sky.

        Parameters
        ----------
        ra : float
            RA of the center of the sky figure.
        dec: float
            DEC of the center of the sky figure.
        sel_srcs : bool, optional
            Show only selected sources. The default is True.
        gw_kwargs : dict, optional
            Keyword arguments for healpy..gnomview. The default is None.
        ps_kwargs : dict, optional
            Keyword arguments for healpy.projscatter. The default is None.

        Returns
        -------
        ax : axes
            Used Matplotlib axes.

        """

        # Get healpix map of Nr of visits
        nside = self.tt_coverage_hp.meta["NSIDE"]
        npix = hpy.nside2npix(nside)
        hp_map = np.zeros(npix, dtype=np.float64)
        pix_ids = self.tt_coverage_hp["pix_id"].data.astype(np.int64)
        hp_map[pix_ids] = self.tt_coverage_hp["nr_vis"]

        # Plot background map
        plt_gw_kwargs = {
            "title": "Nr. of visits",
            "coord": "C",
            "reso": 5 / 60.0,
            "xsize": 2000,
            "ysize": 2000,
            "cmap": "gray",
        }
        if gw_kwargs is not None:
            plt_gw_kwargs.update(gw_kwargs)
        hpy.gnomview(hp_map, rot=[ra, dec], **plt_gw_kwargs)

        # Plots selected sources
        plt_ps_kwargs = {"lonlat": True, "marker": "o", "s": 0.2}
        if ps_kwargs is not None:
            plt_ps_kwargs.update(ps_kwargs)

        tt_srcs = self.tt_sources.group_by("field_id")
        for tt in tt_srcs.groups:
            sel = tt["sel"]
            if not sel_srcs:
                sel = np.ones(len(sel)).astype(bool)
            hpy.projscatter(tt[sel]["ra"], tt[sel]["dec"], **plt_ps_kwargs)

        # hpy.projscatter(
        #    [ra], [dec], lonlat=True, marker="o", s=4.0
        # )  # Mark center of gnomeview
        # hpy.mollview(hp_map, title="Nr. of visits in log10",nest=False,cmap="nipy_spectral",xsize=4800)
        # hpy.graticule(local=False,coord="C",dpar=1.0, color="white") # show graticules every 0.5 deg
        return plt.gca()
