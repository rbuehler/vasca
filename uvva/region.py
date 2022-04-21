#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 14:33:38 2022

@author: buehler
"""
from uvva.field import BaseField, GALEXField
from .uvva_table import UVVATable


class Region(object):
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
        #: Internal list of tables holding the central field data
        #: All of these are defined in ``uvva_table.py``
        self._table_names = list()

    def add_fields(self, field_id, name, "ra", "dec", "observatory", "obsfilter"):
