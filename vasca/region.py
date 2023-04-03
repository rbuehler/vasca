#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import healpy as hpy
import numpy as np
from astropy import units as uu
from astropy.coordinates import SkyCoord
from astropy.table import Table, unique
from loguru import logger

from vasca.field import BaseField, GALEXField
from vasca.source import Source
from vasca.tables import TableCollection
from vasca.tables_dict import dd_vasca_tables


class Region(TableCollection):
    """
    `~vasca.Region` defines a region in the sky as a
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

        for obs in vasca_cfg["observations"]:
            if obs["observatory"] == "GALEX":
                # Loop over fields and store info
                rg_fd_id = 0
                for gfield_id in obs["obs_field_ids"]:
                    gf = GALEXField.load(
                        gfield_id,
                        obs_filter=obs["obs_filter"],
                        method=vasca_cfg["ressources"]["load_method"],
                        load_products=vasca_cfg["ressources"]["load_products"],
                        **vasca_cfg["ressources"]["field_kwargs"],
                    )
                    field_info = dict(gf.tt_field[0])
                    field_info["fov_diam"] = gf.fov_diam
                    field_info["nr_vis"] = gf.nr_vis
                    field_info["time_bin_size_sum"] = gf.time_bin_size_sum
                    field_info["time_start"] = gf.time_start.mjd
                    field_info["time_stop"] = gf.time_stop.mjd
                    field_info["rg_fd_id"] = rg_fd_id
                    rg.tt_fields.add_row(field_info)

                    if vasca_cfg["ressources"]["load_products"]:
                        rg.fields[rg_fd_id] = gf
                    else:
                        rg.fields[rg_fd_id] = None
                    rg_fd_id += 1
            else:
                logger.waring(
                    "Selected observatory `" + obs["observatory"] + "` not supported"
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
        for rg_fd_id, field in self.fields.items():
            tt = field.__dict__[table_name]

            # Apply row selection
            sel = np.ones(len(tt), dtype="bool")
            if only_selected:
                sel = tt["sel"]
            tt_sel = tt[sel]
            tt_sel["rg_fd_id"] = len(tt_sel) * [rg_fd_id]
            ll_tt.append(tt_sel)

        # colnames = dd_vasca_tables["region"][table_name]["names"]
        colnames = list(
            set(dd_vasca_tables["region"][table_name]["names"]) & set(ll_tt[0].colnames)
        )

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

    def add_coverage_hp(self, nside=4096, coord_sys="galactic"):
        """
        Creates healpix arrays of Nr visits, fields and total exposure.

        Parameters
        ----------
        nside : int, optional
            NSIDE of healpix binning. The default is 4096.
        coord_sys : str, optional
            Coordinate system, "galactic" or "icrs"

        Returns
        -------
        hp_nr_vis : [int]
            Array with number of visits per pixel
        hp_exp : TYPE
            Array with exposure per pixel
        hp_nr_fds : TYPE
            Array with number of fields

        """

        npix = hpy.nside2npix(nside)
        pix_diam = hpy.nside2resol(nside, arcmin=True) / 60 * uu.deg

        logger.debug(
            f"Healpix NSIDE: {nside}, NPIX: {npix}, pixel diameter: {pix_diam}"
        )
        pix_nrs = np.arange(npix)
        hp_nr_vis = np.zeros(npix, dtype="float32")
        hp_exp = np.zeros(npix, dtype="float32")
        hp_nr_fds = np.zeros(npix, dtype="float32")
        for field in self.tt_fields:
            pos_vec = hpy.ang2vec(field["ra"], field["dec"], lonlat=True)
            if coord_sys == "galactic":
                cel = SkyCoord(field["ra"], field["dec"], frame="icrs", unit="deg")
                pos_vec = hpy.ang2vec(
                    cel.galactic.l.degree, cel.galactic.b.degree, lonlat=True
                )

            rdisc = field["fov_diam"] / 2.0
            # TODO: Here a more general querry_polygon will have to be done for ULTRASAT
            ipix_disc = hpy.query_disc(
                nside=nside, vec=pos_vec, radius=np.radians(rdisc)
            )

            hp_nr_vis[ipix_disc] += field["nr_vis"]
            hp_exp[ipix_disc] += field["time_bin_size_sum"]
            hp_nr_fds[ipix_disc] += 1

        # Write to table
        sel_pix = hp_nr_vis > 0
        keys_data = ["pix_id", "nr_vis", "exp", "nr_fds"]
        ll_data = [
            pix_nrs[sel_pix],
            hp_nr_vis[sel_pix],
            hp_exp[sel_pix],
            hp_nr_fds[sel_pix],
        ]
        dd_data = dict(zip(keys_data, ll_data))
        self.add_table(dd_data, "region:tt_coverage_hp")
        self.tt_coverage_hp.meta["NSIDE"] = nside
        self.tt_coverage_hp.meta["COOR_SYS"] = coord_sys

        hp_exp[hp_exp < 1e-6] = hpy.UNSEEN
        hp_nr_vis[hp_nr_vis < 1e-6] = hpy.UNSEEN
        hp_nr_vis[hp_nr_fds < 1e-6] = hpy.UNSEEN

        return hp_nr_vis, hp_exp, hp_nr_fds

    def load_from_fits(self, file_name, load_fields=True):
        """
        Loads field from a fits file

        Parameters
        ----------
        file_name : str, optional
            Region file name. The default is "field_default.fits".
        load_fields : bool,
            Load the fields, which have to be located as fits in the subfolder "fields"
            of the region file in "file_name". Default is True.

        Returns
        -------
        None.

        """

        # Load file from TableCollection base class
        super().load_from_fits(file_name)
        region_path = os.path.dirname(file_name)

        # Load fields
        if load_fields:
            for ff in self.tt_fields:
                fd = BaseField()
                fd.load_from_fits(
                    region_path + "/fields/field_" + ff["field_id"] + ".fits"
                )
                self.fields[ff["field_id"]] = fd

    def get_src_from_id(self, rg_src_id):
        """
        Get Source object containng all region table entries
        relevant for the passed rg_src_id

        Parameters
        ----------
        rg_src_id : int
            Region source ID.

        Returns
        -------
        src : vasca.table.Source
            VASCA source object.

        """

        logger.debug(f"Getting Source object for rg_src_id: {rg_src_id}")

        src = Source()

        # Adding source and detection table
        for tt_name in ["tt_sources", "tt_detections"]:
            self.__dict__[tt_name].add_index("rg_src_id")
            src._table_names.append(tt_name)
            setattr(
                src,
                tt_name,
                Table(self.__dict__[tt_name].loc["rg_src_id", [rg_src_id]]),
            )

        # Add visits
        self.tt_visits.add_index("vis_id")
        vis_ids = (unique(src.tt_detections, keys="vis_id"))["vis_id"]
        src._table_names.append("tt_visits")
        setattr(src, "tt_visits", Table(self.tt_visits.loc["vis_id", vis_ids]))

        # Add fields
        self.tt_fields.add_index("rg_fd_id")
        rg_fd_ids = (unique(src.tt_visits, keys="rg_fd_id"))["rg_fd_id"]
        src._table_names.append("tt_fields")
        setattr(src, "tt_fields", Table(self.tt_fields.loc["rg_fd_id", rg_fd_ids]))

        # Add light curve
        src._table_names.append("tt_source_lc")
        setattr(src, "tt_source_lc", src.get_light_curve(rg_src_ids=rg_src_id))

        # Add fd_src_ids to each field
        coord_src = SkyCoord(src.tt_sources["ra"], src.tt_sources["dec"], frame="icrs")
        fd_src_ids = list()
        for field_id in src.tt_fields["field_id"]:
            fd = self.fields[field_id]
            coord_fd_srcs = SkyCoord(
                fd.tt_sources["ra"], fd.tt_sources["dec"], frame="icrs"
            )
            idx, d2d, d3d = coord_src.match_to_catalog_sky(coord_fd_srcs)
            fd_src_ids.append(fd.tt_sources[idx]["fd_src_id"])
        # print(fd_src_ids)
        # print(src.tt_fields)
        src.tt_fields["fd_src_id"] = fd_src_ids

        return src

    def set_src_id_info(self):
        """
        Stores the mapping of rg_src_id to rg_fd_id and fd_src_id
        into tt_src_id_map table.

        Returns
        -------
        None.

        """
        logger.debug("Settign table tt_src_id_map.")

        dd_ids = {"rg_src_id": [], "rg_fd_id": [], "fd_src_id": [], "sel": []}
        tt_ids = unique(
            self.tt_detections["rg_src_id", "rg_fd_id", "fd_src_id", "sel"],
            keys=["rg_fd_id", "fd_src_id"],
        )
        tt_ids.add_index("rg_src_id")

        # TODO: This could be done much more effieciently
        for rg_src_id in self.tt_sources["rg_src_id"]:
            # tt_det = Table()
            ids_idx = np.array(tt_ids.loc_indices["rg_src_id", [rg_src_id]]).flatten()
            for key in dd_ids.keys():
                if len(ids_idx) == 1:
                    dd_ids[key].append(tt_ids[ids_idx[0]][key])
                else:
                    dd_ids[key].extend(tt_ids[ids_idx][key])

        self.add_table(dd_ids, "region:tt_src_id_map")

    def get_src_from_sky_pos(self, coordx, coordy, frame="icrs"):
        """
        Get Source object containng all region table entries
        relevant for the passed source position. The nearest source is matched.

        Parameters
        ----------
        coordx : str, float
            First coordinate component in any astropy.SkyCoord compatible format.
        coordx : str, float
            Second coordinate component in any astropy.SkyCoord compatible format.
        frame : str
            Coordinate system. Any astropy compatible format,
            see https://docs.astropy.org/en/stable/coordinates/skycoord.html

        Returns
        -------
        src : vasca.table.Source
            VASCA source object.
        dist: astropy.Quantity
            Distance to the nearest source
        """

        logger.debug(f"Getting Source object for coord: {coordx},{coordy}")

        coord_srcs = SkyCoord(
            self.tt_sources["ra"], self.tt_sources["dec"], frame="icrs"
        )
        coord_ref = SkyCoord(coordx, coordy, frame=frame)
        idx, d2d, d3d = coord_ref.match_to_catalog_sky(coord_srcs)

        rg_src_id = self.tt_sources[idx]["rg_src_id"]

        src = self.get_src_from_id(rg_src_id)

        return [src, *d2d]
