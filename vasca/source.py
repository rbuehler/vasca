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
from astropy.coordinates import SkyCoord
from loguru import logger
from astroquery.sdss import SDSS

from vasca.tables import TableCollection
from vasca.utils import (
    query_vizier_sed,
    dd_id2filter,
    dd_filter2id,
    dd_filter2wavelength,
    mag2flux,
    tgalex_to_astrotime,
    get_var_stat,
    get_lc_from_gphoton_npfile,
)
from vasca.resource_manager import ResourceManager

rm = ResourceManager()


class Source(TableCollection):
    """
    Class to store all VASCA information for one particular
    source. This class is for convenience, the same data is also found in the Field
    and Region classes containing this source."""

    def __init__(self):
        """
        Many class attributes are stored in astropy.table.Table_.

        .. _astropy.table.Table: https://docs.astropy.org/en/stable/api/astropy.table.Table.html

        Returns
        -------
        None

        """
        # Sets skeleton
        super().__init__()

    def add_vizier_SED(self, vizier_radius=1 * uu.arcsec):
        """
        Add spectral energy distribution table (tt_sed) with all
        spectral points from VizieR within given radius. Uses
        ``vasca.utils.query_vizier_sed()``

        Parameters
        ----------
        vizier_radius astropy.Quantity
            Radius within which to add flux points from VizieR

        Returns
        -------
        None
        """

        # Search for Vizier flux around source, or simbad associated source if present
        ra, dec = self.tt_sources["ra"][0], self.tt_sources["dec"][0]
        if "tt_simbad" in self._table_names:
            ra, dec = self.tt_simbad["ra"][0], self.tt_simbad["dec"][0]

        self.add_table(None, "region:tt_sed")

        try:
            tt_vizier = query_vizier_sed(
                ra, dec, radius=vizier_radius.to(uu.arcsec).value
            )

            # Add columns in right formats for tt_sed later
            tt_vizier["wavelength"] = (cc.c / tt_vizier["sed_freq"]).to(uu.AA)
            tt_vizier["flux"] = tt_vizier["sed_flux"].quantity.to(uu.Unit("1e-6 Jy"))
            tt_vizier["flux_err"] = tt_vizier["sed_eflux"].quantity.to(
                uu.Unit("1e-6 Jy")
            )
            # tt_vizier.sort("wavelength")
            # tt_vizier.pprint_all()

            # Add Vizier info to SED table

            for row in tt_vizier:
                ll_flt = str(row["sed_filter"]).split(":")
                # Check if Vizier data ok
                if (
                    np.isnan(row["flux"])
                    or np.isnan(row["flux_err"])
                    or len(ll_flt) != 2
                    or len(ll_flt[0]) == 0
                    or len(ll_flt[1]) == 0
                ):
                    logger.warning("Skipping row as flux contains nan")
                else:
                    obs, flt = ll_flt
                    dd_vdat = {
                        "flux": row["flux"],
                        "flux_err": row["flux_err"],
                        "wavelength": row["wavelength"],
                        "observatory": obs,
                        "obs_filter": flt,
                        "origin": row["_tabname"],
                    }

                    self.tt_sed.add_row(dd_vdat)
        except Exception:
            print("No entries found in Vizier database")

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
                    "observatory": "GALEX",  # TODO: make this general
                    "obs_filter": obs_flt,
                    "origin": "VASCA",
                }
                self.tt_sed.add_row(dd_vdat)

        # Sort by wavelength
        self.tt_sed.sort("wavelength")

    def add_gphoton_lc(self, s2n_min=3.0, tbin=-1):
        """
        Add light curve from gPhoton. Only include points with no flags.
        Assumes gPhoton flux is given for a 6 arcsec aperture.

        Parameters
        ----------
        s2n_min: float, optional
            Minimum significance of points for selection in light curve.
        tbin_name int, optional
            Time binning in seconds used in the gphoton analysis. If negative assume visit time binning.
        Returns
        -------
        None
        """

        # Get location of gphoton files
        gphot_dir = rm.get_path("gal_gphoton", "sas_cloud")

        # Prepare info for file reading
        rg_src_id = self.tt_sources["rg_src_id"][0]
        ra_src = round(self.tt_sources["ra"][0], 3)
        dec_src = round(self.tt_sources["dec"][0], 3)

        tname = "" if tbin < 0 else "_" + str(tbin)

        # Check if NUV file is present and load it, this is requires
        fname_nuv = (
            gphot_dir
            + "/gPhoton_ra"
            + str(ra_src)
            + "_dec"
            + str(dec_src)
            + "_nuv"
            + tname
            + "_app.npy"
        )
        if os.path.exists(fname_nuv):
            tt_lc = get_lc_from_gphoton_npfile(fname_nuv, "NUV")
        else:
            logger.warning(f"gPhoton file not found {fname_nuv}")
            return

        # If FUV file present add it to table
        fname_fuv = fname_nuv.replace("_nuv", "_fuv")
        if os.path.exists(fname_fuv):
            tt_lc_fuv = get_lc_from_gphoton_npfile(fname_fuv, "FUV")
            tt_lc = vstack([tt_lc, tt_lc_fuv])

        # Add light curve table
        self.add_table(tt_lc, "region:tt_gphoton_lc")

        # Modify selection
        sel = (self.tt_gphoton_lc["flags"] < 0.5) * (
            self.tt_gphoton_lc["s2n"] > s2n_min
        )
        self.tt_gphoton_lc["sel"] = sel

        # Calculate variability statistic for each filter and add it to tt_source
        tt_lc = tt_lc[sel]
        filter_ids = np.sort(np.unique(tt_lc["obs_filter_id"].data))
        dd_gp_var = {"rg_src_id": [rg_src_id], "nr_det": [[]]}
        for flt_id in filter_ids:
            # Get stats variables and write them all to dictionary
            sel_flt = tt_lc["obs_filter_id"] == flt_id
            dd_var = get_var_stat(tt_lc["flux"][sel_flt], tt_lc["flux_err"][sel_flt])
            for var, val in dd_var.items():
                if var in dd_gp_var.keys():
                    dd_gp_var[var][0].append(val)
                else:
                    dd_gp_var[var] = [[val]]
            # Add number of bins
            dd_gp_var["nr_det"][0].append(sel_flt.sum())
        self.add_table(dd_gp_var, "tt_gphoton_stats")

    def add_spectrum(self, search_radius=2 * uu.arcsec):
        """
        Get spectrum from SDSS, if available
        Parameters
        ----------
        search_radius, astropy.Quantity, optional
            Search radius around multifrequency counterpart, or if this does not excist, around VASCA position.

        Returns
        -------
        None
        """

        # Prepare spectral query data
        ra, dec = self.tt_sources["ra"].quantity[0], self.tt_sources["dec"].quantity[0]
        if "tt_simbad" in self._table_names:
            ra, dec = (
                self.tt_simbad["ra"].quantity[0],
                self.tt_simbad["dec"].quantity[0],
            )
        pos = SkyCoord(ra, dec, frame="icrs")

        # Get spectrum from SDSS
        tt_xid = SDSS.query_region(pos, radius=search_radius, spectro=True)
        if type(tt_xid) != type(None):
            ll_sp = SDSS.get_spectra(matches=tt_xid)
            for ii in range(len(ll_sp)):
                hdu_spec = ll_sp[ii]
                tt_spec = Table(hdu_spec["COADD"].data)
                c_Aps = cc.c.to(uu.AA / uu.s)
                spec_flux = tt_spec["flux"] * 1e-17 * uu.erg / (uu.cm**2 * uu.s * uu.AA)
                model_flux = (
                    tt_spec["model"] * 1e-17 * uu.erg / (uu.cm**2 * uu.s * uu.AA)
                )
                spec_wave = np.power(10, tt_spec["loglam"]) * uu.AA
                dd_spec = {
                    "wavelength": spec_wave,
                    "flux": (spec_flux * spec_wave**2 / c_Aps).to(uu.Unit("1e-6 Jy")),
                    "flux_model": (model_flux * spec_wave**2 / c_Aps).to(
                        uu.Unit("1e-6 Jy")
                    ),
                    "s2n": tt_spec["flux"] * np.sqrt(tt_spec["ivar"]),
                }

                # Add data to table collection
                tt_ii = self.table_from_template(dd_spec, "region:tt_spectrum")
                # Select only significant flux pponts
                tt_ii["sel"] = tt_ii["s2n"] > 2
                self.add_table(tt_ii, "tt_spectrum_" + str(ii))
        else:
            logger.warning("No spectrum found")
