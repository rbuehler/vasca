#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 10:17:54 2023

@author: buehler
"""
from vasca.tables import TableCollection


class Source(TableCollection):
    """
    `~vasca.Source` is  class to store all vasca information for one particular
    source. This class is for conviniece, the same data is also found in the Field
    and Region classes containing this source."""

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
