#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from loguru import logger

import numpy as np
import healpy as hpy
from astropy.table import Column, vstack
from astropy import units as uu

from uvva import tables
from uvva.field import GALEXField
from uvva.tables import TableCollection


class Region(TableCollection):
    """
    :class: `~uvva.Region` defines a region in the sky as a
    list of uvva.field objects. It provides functionality to
    loop over fields to derive source lists, etc.
    """

    def __init__(self):
        """

        Notes
        -----
        Many class attributes are stored in astropy.table.Tables_. To see a
        description of each of their columns run :meth: `~uvva.Regions.info`.

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
    def load_from_config(cls, uvva_cfg):
        """
        Loads region from configuration from dictionary

        Parameters
        ----------
        uvva_cfg : dict
            Dictionary with region parameters derived from the uvva pipeline
            YAML configuration file.

        Returns
        -------
        None.

        """

        rg = cls()

        logger.debug("Loading fields from config file")

        if uvva_cfg["observations"]["observatory"] == "GALEX":

            # Loop over fields and store info
            for field_id in uvva_cfg["observations"]["field_ids"]:

                gf = GALEXField.load(
                    field_id,
                    obs_filter=uvva_cfg["observations"]["obs_filter"],
                    method=uvva_cfg["ressources"]["load_method"],
                    load_products=uvva_cfg["ressources"]["load_products"],
                    **uvva_cfg["ressources"]["field_kwargs"],
                )
                field_info = dict(gf.tt_field[0])
                field_info["fov_diam"] = 0.55
                field_info["nr_vis"] = gf.nr_vis
                field_info["time_bin_size_sum"] = gf.time_bin_size_sum
                field_info["time_start"] = gf.time_start.mjd
                field_info["time_stop"] = gf.time_stop.mjd
                rg.tt_fields.add_row(field_info)

                if uvva_cfg["ressources"]["load_products"]:
                    rg.fields[field_id] = gf
                else:
                    rg.fields[field_id] = None
        else:
            logger.waring(
                "Selected observatory `"
                + uvva_cfg["observations"]["observatory"]
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

        ll_tt = []

        for field_id, field in self.fields.items():

            tt = field.__dict__[table_name]

            field_id_col = np.ones(len(tt), dtype="int64") * field_id
            col = Column(
                data=field_id_col,
                name="field_id",
                dtype="int64",
                unit="1",
                description="Field ID nr.",
            )
            tt.add_column(col)

            # Apply row selection
            sel = np.ones(len(tt), dtype="bool")
            if only_selected:
                sel = tt["sel"]

            ll_tt.append(tt[sel])

        # Add stacked table, for tables with vecotr entries this needs to be done by hand
        if table_name.startswith("tt_"):
            tt_data = vstack(ll_tt)
            if table_name in self._table_names:
                logger.warning(f"Table '{table_name}' already exists, overwriting")
            self._table_names.append(table_name)
            setattr(self, table_name, tt_data)
        elif table_name.startswith("ta_"):
            colnames = ll_tt[0].colnames

            # Create empty data structure
            dd_data = dict(zip(colnames, [list() for ii in range(len(colnames))]))
            print(dd_data)

            for tt in ll_tt:
                for colname in colnames:
                    dd_data[colname].extend(tt[colname].tolist())

            # For vector columns convert to numpy arrays of type object
            for colname in colnames:
                if len(np.array(dd_data[colname]).shape) > 1:
                    dd_data[colname] = np.array(dd_data[colname], dtype=np.object_)

            self.add_table(dd_data, "region:" + table_name, add_sel_col=False)

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
        self.add_table(
            [pix_nrs[sel_pix], hp_vis[sel_pix], hp_exp[sel_pix]],
            "region:tt_coverage_hp",
            add_sel_col=False,
        )
        self.tt_coverage_hp.meta["NSIDE"] = nside

        hp_exp[hp_exp < 1e-6] = hpy.UNSEEN
        hp_vis[hp_vis < 1e-6] = hpy.UNSEEN

        return hp_vis, hp_exp
