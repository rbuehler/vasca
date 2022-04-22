#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 14:33:38 2022

@author: buehler
"""
from uvva.field import BaseField, GALEXField
from uvva.tables import TableCollection
from uvva import tables
from loguru import logger


class Region(TableCollection):
    """
    :class: `~uvva.Region` defines a region in the sky as a
    list of uvva.field objects. It provides funtionality to
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

        # Add empty tables
        self.add_table(None, "region:tt_fields")

    def load_from_config(self, obs):
        """
        Loads region from configuration from dictionary

        Parameters
        ----------
        obs : dict
            Dictionary with region parameters derived from the uvva pipeline
            TOML configuration file.

        Returns
        -------
        None.

        """
        logger.debug("Loading fields from config file")
        if obs["observatory"] == "GALEX":
            for field_id in obs["field_ids"]:
                gf = GALEXField.from_MAST(
                    obs_id=field_id, obs_filter=obs["obs_filter"], load_products=False)
                field_info = dict(gf.tt_field[0])
                field_info["size"] = 0.55
                field_info["n_visits"] = gf.n_visits
                field_info["time_bin_size_sum"] = gf.time_bin_size_sum
                field_info["time_start"] = gf.time_start.mjd
                field_info["time_stop"] = gf.time_stop.mjd
                self.tt_fields.add_row(field_info)
        else:
            logger.waring("Selected observatory `"+obs["observatory"]+"` not supportet")
