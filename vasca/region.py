#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import warnings

import healpy as hpy
import numpy as np
from astropy import units as uu
from astropy.coordinates import SkyCoord
from astropy.table import Table, unique, join
from astropy.utils.exceptions import AstropyWarning
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from loguru import logger

from vasca.field import BaseField, GALEXDSField, GALEXField
from vasca.source import Source
from vasca.tables import TableCollection
from vasca.tables_dict import dd_vasca_tables
from vasca.utils import dd_filter2id


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

        # Path where region and fields are stored
        self.region_path = None

        self.fields = {}  # dictionary of field_id and field objects

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
        rg.add_table(None, "region:tt_fields")

        logger.debug("Loading fields from config file")

        # Loop over all observations
        rg_fd_id = 0

        for obs in vasca_cfg["observations"]:
            if obs["observatory"] == "GALEX" or obs["observatory"] == "GALEX_DS":
                gfield_load_func = (
                    GALEXDSField.load
                    if obs["observatory"] == "GALEX_DS"
                    else GALEXField.load
                )
                # Loop over fields and store info
                for gfield_id in obs["obs_field_ids"]:
                    try:
                        # Data needs to be preloaded serially, to avoid many parallel
                        # requests to the MAST server later which results in errors
                        gf = gfield_load_func(
                            gfield_id,
                            obs_filter=obs["obs_filter"],
                            method=vasca_cfg["resources"]["load_method"],
                            load_products=vasca_cfg["resources"]["load_products"],
                            **vasca_cfg["resources"]["field_kwargs"],
                        )

                        rg_fd_id += 1
                        if type(gf) == type(None):
                            continue

                        field_info = dict(gf.tt_fields[0])
                        field_info["fov_diam"] = gf.fov_diam
                        field_info["nr_vis"] = gf.nr_vis
                        field_info["time_bin_size_sum"] = gf.time_bin_size_sum
                        field_info["time_start"] = gf.time_start.mjd
                        field_info["time_stop"] = gf.time_stop.mjd
                        field_info["rg_fd_id"] = rg_fd_id
                        rg.tt_fields.add_row(field_info)
                    except Exception as e:
                        logger.exception(
                            f"Faild to load field '{gfield_id}'. "
                            f"Moving on to next field. Exception:\n {e}"
                        )
                        continue
            else:
                logger.warning(
                    "Selected observatory `" + obs["observatory"] + "` not supported"
                )

        # Add region path
        rg.region_path = (
            vasca_cfg["general"]["out_dir_base"] + "/" + vasca_cfg["general"]["name"]
        )

        # Add filter table
        rg.add_table(None, "region:tt_filters")
        for flt, flt_id in dd_filter2id.items():
            rg.tt_filters.add_row(
                {"obs_filter": flt, "obs_filter_id": flt_id, "obs_filter_idx": -1}
            )

        return rg

    def add_table_from_fields(
        self, table_name, only_selected=False, sum_obs_filter=False
    ):
        """
        Add tables from the fields to the region by stacking them,
        adding the rg_fd_id column.


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
        self.tt_fields.add_index("field_id")

        ll_tt = []  # List of "table_name" tables for all fields
        for field_id, field in self.fields.items():
            rg_fd_id = self.tt_fields.loc["field_id", [field_id]]["rg_fd_id"]
            if table_name in field.__dict__.keys():
                tt = field.__dict__[table_name]
            else:
                logger.warning(f"No table {table_name} found for field {field_id}")
                continue

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
        # for colname in colnames:
        #    if len(np.array(dd_data[colname], dtype=object).shape) > 1:
        #        dd_data[colname] = np.array(dd_data[colname], dtype=np.object_)
        # tt_tmp = self.table_from_template(dd_data, "region:" + table_name)
        # print("l", len(tt_tmp))
        # if sum_obs_filter:
        #     tt_tmp.add_index("vis_id")
        #     vis_ids = np.unique(tt_tmp["vis_id"])
        #     for vis_id in vis_ids:
        #         idxs = tt_tmp.loc_indices["vis_id", vis_id]

        #         flt_id = np.unique(tt_tmp["obs_filter_id"][idxs].data).sum()
        #         # flt_id_vec = flt_id.sum() * np.ones(len(flt_id))

        #         tt_tmp[idxs]["vis_id"] = flt_id
        #     tt_tmp = unique(tt_tmp, keys="vis_id")

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
            # TODO: Here a more general query_polygon will have to be done for ULTRASAT
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

    def load_from_fits(self, file_name, load_fields=False):
        """
        Loads field from a fits file

        Parameters
        ----------
        file_name : str, optional
            Region file name.
        load_fields : bool,
            Load the fields, which have to be located as fits in the subfolder "./fields/"
            of the region file in "file_name". Default is False.

        Returns
        -------
        None.

        """

        # Load file from TableCollection base class
        super().load_from_fits(file_name)
        self.region_path = os.path.dirname(file_name)

        # Load fields
        if load_fields:
            for ff in self.tt_fields:
                self.get_field(
                    field_id=ff["field_id"], load_method="FITS", add_field=True
                )

    def get_src_from_id(
        self,
        rg_src_id,
        load_from_file=True,
        write_to_file=True,
        add_sed=True,
        add_gphoton=True,
        add_spectrum=True,
    ):
        """
        Get Source object containing all region table entries
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

        # Check if source shall be loaded for file and file exists
        fname_src = self.region_path + "/sources/src_" + str(rg_src_id) + ".fits"
        if os.path.exists(fname_src) and load_from_file:
            logger.debug(f"Loading source from file {fname_src}")
            src.load_from_fits(fname_src)
            return src
        # Otherwise derive from region
        else:
            # Adding source and detection table
            # Tables to copy entries from, if they exist in region
            ll_tables = [
                "tt_sources",
                "tt_detections",
                "tt_coadd_sources",
                "tt_simbad",
                "tt_gaiadr3",
            ]

            for tt_name in ll_tables:
                if hasattr(self, tt_name):
                    if rg_src_id in self.__dict__[tt_name]["rg_src_id"]:
                        self.__dict__[tt_name].add_index("rg_src_id")
                        src.add_table(
                            Table(self.__dict__[tt_name].loc["rg_src_id", [rg_src_id]]),
                            tt_name,
                        )

            # Add visits
            self.tt_visits.add_index("vis_id")
            vis_ids = (unique(src.tt_detections, keys="vis_id"))["vis_id"]
            src.add_table(Table(self.tt_visits.loc["vis_id", vis_ids]), "tt_visits")

            # Add fields
            self.tt_fields.add_index("rg_fd_id")
            rg_fd_ids = (unique(src.tt_detections, keys="rg_fd_id"))["rg_fd_id"]
            src.add_table(Table(self.tt_fields.loc["rg_fd_id", rg_fd_ids]), "tt_fields")

            # Add filter
            src.add_table(Table(self.tt_filters), "tt_filters")

            # Add light curve
            src.add_table(
                src.get_light_curve(rg_src_ids=rg_src_id)[rg_src_id], "tt_source_lc"
            )

            #Add additional info if asked
            if add_sed:
                src.add_vizier_SED()
            if add_gphoton:
                src.add_gphoton_lc()
            if add_spectrum:
                src.add_spectrum()

            #Write if asked
            if write_to_file:
                # create source directory if it does not exist
                os.makedirs(self.region_path + "/sources", exist_ok=True)
                src.write_to_fits(fname_src)

            return src

    def set_src_id_info(self):
        """
        Stores the mapping of rg_src_id to rg_fd_id and fd_src_id
        into tt_src_id_map table.

        Returns
        -------
        None.

        """
        logger.debug("Setting table tt_src_id_map.")

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
        Get Source object containing all region table entries
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

    def get_field(
        self,
        field_id=None,
        rg_fd_id=None,
        load_method="FITS",
        add_field=False,
        mast_products="TABLES",
        field_kwargs=dict(),
    ):
        """
        Load a field from a region, tt_fields table needs to include this field.

        Parameters
        ----------
        field_id : TYPE, str
            Field ID. The default is None.
        rg_fd_id : TYPE, int
            Region field ID. The default is None.
        load_method : TYPE, optional
            Method to load the field. Either from "FITS", from "MAST" or from "VASCA".
            The default is "FITS". Note that if the field is already in the
            region.fields dictionary, this will be ignored and the later be returned.
        add_field : TYPE, optional
            Add the field to the region.fields dictionary. The default is True.
        mast_products : str, optional
            load_products option of field.load_from_MAST. Sets to what level data
            products are downloaded.

        Returns
        -------
        field
            VASCA field

        """

        # Select field and set both field ids
        if type(field_id) is not type(None):
            sel_fd = self.tt_fields["field_id"] == field_id
            rg_fd_id = self.tt_fields[sel_fd]["rg_fd_id"][0]
        elif type(rg_fd_id) is not type(None):
            sel_fd = self.tt_fields["rg_fd_id"] == rg_fd_id
            field_id = self.tt_fields[sel_fd]["field_id"][0]
        else:
            logger.warning("Need to specify either field_id or rg_fd_id")

        fd_row = self.tt_fields[sel_fd]
        if len(fd_row) != 1:
            logger.warning(
                f"Could not identify single field for {field_id}, {rg_fd_id}"
            )
        # ** Load field accorading to passed method **

        # If already loaded into field dictionary return it
        if field_id in self.fields.keys():
            return self.fields[field_id]

        # If it shall be loaded from a fits file
        elif load_method == "FITS":
            fd = BaseField()
            fd.load_from_fits(self.region_path + "/fields/field_" + field_id + ".fits")
            if add_field:
                self.fields[field_id] = fd
            return fd

        # If it should be loaded from MAST or VASCA-MAST file database
        elif fd_row["observatory"] == "GALEX" or fd_row["observatory"] == "GALEX_DS":
            gfield_load_func = (
                GALEXDSField.load
                if fd_row["observatory"] == "GALEX_DS"
                else GALEXField.load
            )

            gf = gfield_load_func(
                str(fd_row["field_name"][0])
                if fd_row["observatory"] == "GALEX_DS"
                else int(str(fd_row["field_id"][0])[3:]),
                obs_filter=str(fd_row["obs_filter"][0]),
                method=load_method,
                load_products=mast_products,
                **field_kwargs,
            )
            if add_field:
                self.fields[field_id] = gf
            return gf
        else:
            logger.warning("Usupported loading method")

    # def set_lomb_scargle(self):

    def cross_match_cds(
        self,
        query_radius=1 * uu.arcsec,
        query_table="I/355/gaiadr3",
        vizier_columns=["*", "PQSO", "PGal", "PSS", "RPlx"],
        overwrite=False,
    ):
        """
        Match sources in region with SIMBAD-catalogs or Vizier database catalog. Runs only over selected sources.

        Parameters
        ----------
        query_radius : astropy.quantity, optional
            Query up to this distance from VASCA sources. The default is 1 * uu.arcsec.
        query_table : str, optional
            Vizier table to query, if "simbad" query SIMBAD instead. The default is the main GAIA-DR3 table.
        vizier_columns: list(str)
            Vizier catalog columns to get from the catalog. "*" is for default columns,
             "**" for all columns. The default is for selected GAIA columns, see:
            https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3
        overwrite: bool, optional
            Overwrite preexisting association for this source. The default is False.

        Returns
        -------
        None.

        """
        logger.debug(f"Query {query_table} table")

        # Define variables for later, depending on catalog source
        if query_table.lower() == "simbad":
            cat_name = "simbad"
        else:
            cat_name = query_table.split("/")[-1]
        tab_name = "tt_" + cat_name

        if tab_name in self._table_names:
            logger.warning(f"Region already contained {tab_name} info")
            if overwrite:
                self.remove_tables([tab_name])
            else:
                logger.warning(f"As overwrite is {overwrite}, query stopped")
                return

        # Get selected sources
        tt_src = self.tt_sources[self.tt_sources["sel"]]

        # Get coordinates for query
        coords = SkyCoord(
            tt_src["ra"].quantity,
            tt_src["dec"].quantity,
            frame="icrs",
        )#[:100]

        # ---- Run SIMBAD query and modify query results
        if query_table.lower() == "simbad":
            customSimbad = Simbad()
            customSimbad.TIMEOUT = 180

            # Get only this subset of SIMBAD variables
            vo_entries = [
                "otype(opt)",
                "otypes",
                "distance",
                "distance_result",
                "velocity",
                "z_value",
                "sptype",
            ]
            customSimbad.add_votable_fields(*vo_entries)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                logger.debug(f"Starting SIMBAD query for {len(coords)} sources..")
                tt_qr = customSimbad.query_region(coords, radius=query_radius)
                logger.debug("..query done.")

            # Add rg_src_id
            tt_qr["rg_src_id"] = tt_src["rg_src_id"][tt_qr["SCRIPT_NUMBER_ID"] - 1]

            # Change type for astropy
            vo_change_type = [
                "MAIN_ID",
                "COO_BIBCODE",
                "OTYPE_opt",
                "OTYPES",
                "RVZ_BIBCODE",
                "SP_TYPE",
                "SP_QUAL",
                "SP_BIBCODE",
            ]
            for vo in vo_change_type:
                tt_qr[vo] = tt_qr[vo].data.astype("S32")

            # Lower case column names
            for col in tt_qr.colnames:
                tt_qr.rename_column(col, col.lower())

            # Change unit to astropy format and change some names
            skycoords = SkyCoord(
                tt_qr["ra"].data,
                tt_qr["dec"].data,
                frame="icrs",
                unit=(uu.hourangle, uu.deg),
            )
            tt_qr["ra"] = skycoords.ra.degree * uu.deg
            tt_qr["dec"] = skycoords.dec.degree * uu.deg
            tt_qr.rename_column("otype_opt", "otype")
            tt_qr.rename_column("distance_result", "match_distance")

        # ---- Run Vizier query and modify query results
        else:
            logger.debug(f"Starting Vizier query for {len(coords)} sources..")

            # Define columns from GAIA to get.
            # https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=I/355/gaiadr3

            customVizier = Vizier(
                columns=vizier_columns,
                catalog=query_table,
                timeout=180,
            )
            tt_qr = (customVizier.query_region(coords, radius=query_radius))[0]
            logger.debug("..query done.")

            # Add rg_src_id
            tt_qr["rg_src_id"] = tt_src["rg_src_id"][tt_qr["_q"] - 1]

            # Rename columns
            tt_qr.rename_column("RA_ICRS", "ra")
            tt_qr.rename_column("DE_ICRS", "dec")

            # Add distance between Vizier and VASCA source
            match_coords = SkyCoord(
                tt_qr["ra"].quantity, tt_qr["dec"].quantity, frame="icrs"
            )
            idx_cat, dist_cat, _ = match_coords.match_to_catalog_sky(coords)
            tt_qr["match_distance"] = dist_cat.to(uu.arcsec)

        # Add vasca internal coutnerpart ID
        tt_qr[cat_name + "_match_id"] = np.array(range(0, len(tt_qr)), dtype=np.int32)

        # # Sort table
        tt_qr.sort(["rg_src_id", "match_distance"])
        #
        # Add match_id to tt_sources
        tu_qr = unique(tt_qr, keys="rg_src_id")
        self.tt_sources.add_index("rg_src_id")
        idx = self.tt_sources.loc_indices["rg_src_id", tu_qr["rg_src_id"]]
        self.tt_sources[cat_name + "_match_id"] = -1 * np.ones(
            len(self.tt_sources), dtype=np.int32
        )
        self.tt_sources[cat_name + "_match_id"][idx] = tu_qr[cat_name + "_match_id"]

        # Add table
        self.add_table(tt_qr, tab_name)

        # Add table explaining the different otypes
        if query_table.lower() == "simbad":
            self.add_simbad_otype_info()

    def add_simbad_otype_info(self):
        """
        Add table explaing SIMBAD ogroups

        Returns
        -------

        """

        tu_qr = unique(self.tt_simbad, keys="rg_src_id")

        # Get all associated otypes
        otypes_all, otype_cts_all = np.unique(tu_qr["otype"], return_counts=True)

        # Read files with otype description
        # TODO make the ressource manager handle this
        fpath = os.path.dirname(__file__)
        tt_nodes = Table.read(
            fpath + "/examples/resources/SIMBAD_otypes/otypes_nodes.csv"
        )

        # Get index for all associated otypes in deeescription table
        ids, ids_idx, _ = np.intersect1d(
            tt_nodes["Id"], np.array(otypes_all), return_indices=True
        )

        # Consider also unused otypes (marked by?)
        candidate = np.asarray(
            np.ma.masked_array(
                data=tt_nodes["Candidate"], mask=False, fill_value="none"
            )
        )
        can, can_idx, _ = np.intersect1d(
            candidate, np.array(otypes_all), return_indices=True
        )

        # Merge index of sure and unsure otypes
        all_idx = np.unique(np.append(ids_idx, can_idx))
        tt_nodes.rename_column("Id", "otype")

        # Add table and ogroup
        self.add_table(tt_nodes[all_idx], "tt_otypes")

    def synch_src_sel(self):
        """
        Synchronize selections among tables. All tables containing "rg_src_id" only keep rows
        for selected sources in tt_sources

        Returns
        -------
        None.

        """
        for tab_name in self._table_names:
            if "rg_src_id" in self.__dict__[tab_name].colnames:
                logger.debug(f"Synchronizing selection in table {tab_name}")
                sel = np.in1d(
                    self.__dict__[tab_name]["rg_src_id"],
                    self.tt_sources["rg_src_id"][self.tt_sources["sel"]],
                )
                self.__dict__[tab_name]["sel"] = sel

    def get_region_catalog(self):
        """
        Create a reduced region, which only contains info on selected sources

        Returns
        -------
        rg : vasca.Region
            Region with only the selected sources

        """

        # Create catalog region
        rg = Region()

        # keep only selected sources and detections, etc
        self.synch_src_sel()

        # Copy the tables listed below fully into catalog region
        for tab_name in self._table_names:
            if "sel" in self.__dict__[tab_name].colnames:
                rg.add_table(
                    self.__dict__[tab_name][self.__dict__[tab_name]["sel"]],
                    tab_name,
                )
            else:
                rg.add_table(self.__dict__[tab_name], tab_name)
        return rg
