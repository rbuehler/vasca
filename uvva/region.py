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
        logger.debug("Loading fields from config file")
        if obs["observatory"] == "GALEX":
            for field_id in obs["field_ids"]:
                gf = GALEXField(obs_id=field_id, obs_filter=obs["obs_filter"])
                gf._load_galex_field_info(obs_id=field_id, obs_filter=obs["obs_filter"])
                field_info = dict(gf.tt_field[0])
                field_info["size"] = 0.55
                self.tt_fields.add_row(field_info)
