#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 10:17:54 2023

@author: buehler
"""
import os
import numpy as np
from astropy import constants as cc
from astropy import units as uu
from astropy.table import Table, vstack
from loguru import logger

from vasca.tables import TableCollection
from vasca.utils import query_vizier_sed, dd_id2filter,dd_filter2id,dd_filter2wavelength, mag2flux, tgalex_to_astrotime
from vasca.resource_manager import ResourceManager
rm = ResourceManager()

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
        Add spectral energy distribution table (tt_sed) with all
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

        # Add Vizier info to SED table
        self.add_table(None, "region:tt_sed")
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

                self.tt_sed.add_row(dd_vdat)

        # Add mean VASCA flux for all filters
        flt_ids = np.array(self.tt_sources["obs_filter_id"]).flatten()
        for flt_idx in range(len(flt_ids)):
            flt_id = self.tt_sources["obs_filter_id"][:, flt_idx][0]
            if self.tt_sources["flux"].quantity[:, flt_idx][0] > 0:
                obs_flt = dd_id2filter[flt_id]
                dd_vdat = {
                    "flux": self.tt_sources["flux"][:, flt_idx][0],
                    "flux_err": self.tt_sources["flux_err"][:, flt_idx][0],
                    "wavelength": dd_filter2wavelength[obs_flt],
                    "observatory": "GALEX", #TODO: make this general
                    "obs_filter": obs_flt,
                    "origin": "VASCA",
                }
                self.tt_sed.add_row(dd_vdat)

        # Sort by wavelength
        self.tt_sed.sort("wavelength")

    def add_gphoton_lc(self, s2n_min=3.):
        """
        Add light curve from gPhoton. Only include poitns with no flags.
        Parameters
        ----------
        s2n_min: float, optional
            Minimum significance of points for selection in light curve.
        Returns
        -------
            None
        """
        def get_lc_from_gphoton_npfile(file_name, obs_filter):
            "Helper function to load gphoton results pickel with numpy"
            dd_gph = np.load(file_name, allow_pickle='TRUE').item()

            # Get gphoton lc
            keep_keys = ('t_mean', 'exptime', 'flux_bgsub', 'flux_bgsub_err', "flags",
                         "mag_mcatbgsub", "mag_mcatbgsub_err_2","flags")
            dd_gap = {x: dd_gph["gAperture"][x] for x in keep_keys if
                      x in dd_gph["gAperture"]}
            dd_gap["s2n"] = dd_gap['flux_bgsub'] / dd_gap['flux_bgsub_err']

            #Rename key and change units
            dd_gap["time_bin_size"] = dd_gap.pop('exptime')
            # Units of flux_bgsub are in erg sec^-1 cm^-2 Å^-1. . Get also Jy flux from AB magnitude
            dd_gap['flux'], dd_gap['flux_err'] = mag2flux(dd_gap["mag_mcatbgsub"],
                                                          dd_gap["mag_mcatbgsub_err_2"])
            # Units of time are in "GALEX Time" = "UNIX Time" - 315964800, change to MJD
            dd_gap['time'] = tgalex_to_astrotime(dd_gap['t_mean'], "mjd")

            dd_gap["obs_filter"] = [obs_filter] * len(dd_gap['flux'])
            dd_gap["obs_filter_id"] = [dd_filter2id[obs_filter]] * len(dd_gap['flux'])
            return Table(dd_gap)

        # Get location of gphoton files
        gphot_dir = rm.get_path("gal_gphoton", "sas_cloud")

        #Prepare info for file reading
        rg_src_id = self.tt_sources["rg_src_id"][0]
        ra_src = round(self.tt_sources["ra"][0], 5)
        dec_src = round(self.tt_sources["dec"][0], 5)

        #Check if NUV file is present and load it, this is requires
        fname_nuv = gphot_dir + "/gPhoton_ra" + str(ra_src) + "_dec" + str(dec_src) + "_nuv_app.npy"
        if os.path.exists(fname_nuv):
            tt_lc = get_lc_from_gphoton_npfile(fname_nuv,"NUV")
        else:
            logger.warning(f"gPhoton file not found {fname}")
            return

        #If FUV file present add it to table
        fname_fuv = fname_nuv.replace("_nuv","_fuv")
        if os.path.exists(fname_fuv):
            tt_lc_fuv = get_lc_from_gphoton_npfile(fname_fuv, "FUV")
            tt_lc = vstack([tt_lc, tt_lc_fuv])

        #Add light curve table
        self.add_table(tt_lc, "region:tt_gphoton_lc")

        #Modify selection
        sel = (self.tt_gphoton_lc["flags"]<0.5) * (self.tt_gphoton_lc["s2n"]>s2n_min)
        self.tt_gphoton_lc["sel"] = sel