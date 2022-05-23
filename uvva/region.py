#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from loguru import logger

import numpy as np
from uvva import tables
from uvva.field import GALEXField
from uvva.tables import TableCollection
from astropy.table import Column, vstack


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
        self.add_table(None, "region:tt_visits")

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
                field_info["size"] = 0.55
                field_info["n_visits"] = gf.n_visits
                field_info["time_bin_size_sum"] = gf.time_bin_size_sum
                field_info["time_start"] = gf.time_start.mjd
                field_info["time_stop"] = gf.time_stop.mjd
                rg.tt_fields.add_row(field_info)

                # Loop over visits and store info
                keys_store = tables.dd_uvva_tables["base_field"]["tt_visits"]["names"]
                for ii in range(0, len(gf.tt_visits)):
                    visits_info = dict(gf.tt_visits[keys_store][0])
                    visits_info["field_id"] = field_id
                    rg.tt_visits.add_row(visits_info)
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

        # Add stacked table
        tt_data = vstack(ll_tt)
        if table_name in self._table_names:
            logger.warning(f"Table '{table_name}' already exists, overwriting")
        self._table_names.append(table_name)
        setattr(self, table_name, tt_data)
