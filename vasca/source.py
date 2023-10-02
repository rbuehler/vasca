#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 10:17:54 2023

@author: buehler
"""
import numpy as np
from astropy import constants as cc
from astropy import units as uu
from loguru import logger

from vasca.tables import TableCollection
from vasca.utils import query_vizier_sed


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

    def add_vizier_SED(self, vizier_radius=1 * uu.arcsec):
        """
        Add spectral energy distribution table (tt_vizier_sed) with all
        spectral points from VizieR within given radius

        Parameters
        ----------
        self vasca.TableCollection
            Table collection with all source information. SED table
            will be added to this collection.

        vizier_radius astropy.quantity
            Radius within which to add flux points from VizieR

        Returns
        -------

        """

        # Search for Vizier flux around source, or simbad associated source if present
        ra, dec = self.tt_sources["ra"][0], self.tt_sources["dec"][0]
        if "tt_simbad" in self._table_names:
            ra, dec = self.tt_simbad["ra"][0], self.tt_simbad["dec"][0]
        tt_vizier = query_vizier_sed(ra, dec, radius=vizier_radius.to(uu.arcsec).value)

        # Add columns in right formats for tt_sed later
        tt_vizier["wavelength"] = (cc.c / tt_vizier["sed_freq"]).to(uu.AA)
        tt_vizier["flux"] = tt_vizier["sed_flux"].quantity.to(uu.Unit("1e-6 Jy"))
        tt_vizier["flux_err"] = tt_vizier["sed_eflux"].quantity.to(uu.Unit("1e-6 Jy"))

        self.add_table(None, "region:tt_vizier_sed")
        for row in tt_vizier:
            if np.isnan(row["flux"]) or np.isnan(row["flux_err"]):
                logger.warning("Skipping row as flux contains nan")
            else:
                obs, flt = str(row["sed_filter"]).split(":")
                dd_vdat = {
                    "flux": row["flux"],
                    "flux_err": row["flux_err"],
                    "wavelength": row["wavelength"],
                    "observatory": obs,
                    "obs_filter": flt,
                    "origin": row["_tabname"],
                }

                self.tt_vizier_sed.add_row(dd_vdat)

        # for ii in range(len(flts)):
        #     if self.tt_sources["flux"].quantity[:, ii][0] > 0:
        #         tt_vasca.add_row(
        #             [
        #                 self.tt_sources["flux"].quantity[:, ii][0],
        #                 self.tt_sources["flux_err"].quantity[:, ii][0],
        #                 dd_filter2wavelength[flts[ii]],
        #             ]
        #         )
