#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Field classes for VASCA
"""

import os
import time
import warnings
from datetime import datetime
from glob import glob

import numpy as np
from astropy import units as uu
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, unique, vstack
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import wcs
from astroquery.mast import Observations
from loguru import logger
from regions import CircleSkyRegion
from requests.exceptions import HTTPError

from vasca.resource_manager import ResourceManager
from vasca.tables import TableCollection
from vasca.tables_dict import dd_vasca_columns
from vasca.utils import dd_filter2id, get_field_id, name2id

# global paths
# path to the dir. of this file
FILE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = FILE_DIR + "/../"  # path to the root directory of the repository

dimless = uu.dimensionless_unscaled
# conf.replace_warnings = ["always"]


class BaseField(TableCollection):
    """
    Class that defines the basic
    data structure for field-based analysis. One *field* is generally
    the area in the sky covered by a telescope in one observation.
    A field is generally composed of several *visits* of the telescope
    at different times.

    This class contains the main functionality for source
    detection and drawing. To be inherited by field analysis classes,
    which can then be tailored to the needs of the observatories supported
    by the VASCA pipeline.
    """

    def __init__(self):
        """
        Many class attributes are stored in astropy.table.Table_. To see a
        description of each of their columns run :meth: `~vasca.field.BaseField.info`.

        .. _astropy.table.Table: https://docs.astropy.org/en/stable/api/astropy.table.Table.html

        Returns
        -------
        None

        """
        # Sets skeleton
        super().__init__()

        #: Internal dictionary of important parameters
        #: with corresponding tables
        #: to be set as class attributes for convenience.
        self._dd_attr_names = {
            "tt_fields": [
                "field_id",
                "field_name",
                "ra",
                "dec",
                "fov_diam",
                "observatory",
                "obs_filter",
                "center",
            ],
            "tt_visits": [
                "nr_vis",
                "time_bin_size_sum",
                "time_start",
                "time_stop",
            ],
        }

        #: Reference image
        self.ref_img = None

        #: Reference WCS
        self.ref_wcs = None

        #: Visit images (assumes same WCS as reference image)
        self.vis_img = None

    def load_sky_map(self, file_name, img_attr="ref_img"):
        """
        Load reference image and WCS from passed fits file

        Parameters
        ----------
        filename : str
            File name of image FITS
        img_attr: str, optional
            Class attribute to store the image in, 'ref_img' or 'vis_img'.
            Default is 'ref_img'

        Returns
        -------
        None

        """
        logger.debug(f"Loading skypmap from file: '{file_name}'")
        with fits.open(file_name) as ff:
            if img_attr == "ref_img":
                self.ref_img = np.array(ff[0].data, dtype=np.float32)
                self.ref_wcs = wcs.WCS(ff[0].header)
            if img_attr == "vis_img":
                img_data = np.array(ff[0].data, dtype=np.float16)
                if type(self.vis_img) == type(None):
                    self.vis_img = img_data
                # If second or later image append
                else:
                    if self.vis_img.ndim == 2:
                        self.vis_img = [self.vis_img]
                    self.vis_img = np.append(self.vis_img, [img_data], axis=0)

    def get_upper_limits(self):
        """
        Calculates magnitude upper limits on non-detections to the tt_visits table.

        Returns
        -------
        float array
            Array of upper limits for each visit

        """

        observatory = (
            self.observatory
            if hasattr(self, "observatory")
            else self.get_field_par("observatory", "tt_fields")
        )
        upper_limit = None

        # Call upper limit function according to observatory.
        # String matching is case insensitive
        if observatory.casefold() == "GALEX".casefold():
            upper_limit = GALEXField.get_visit_upper_limits(self.tt_visits)

        return upper_limit

    def set_light_curve(self, add_upper_limits=True):
        """
        Adds detections information into ta_sources_lc for sources
        listed in self.tt_sources

        Parameters
        ----------
        add_upper_limits : bool, optional
            Add upper limits to the tt_sources_lc, for visits with no detection.
            Upper limits are stored in the flux_err columns for none detections.
            The default is True.

        Returns
        -------
        None

        """
        logger.info(
            "Creating light curve table. "
            f"Upper limits option set to {add_upper_limits}."
        )

        # Selection
        sel = self.tt_detections["sel"]

        # Prepare visit info
        nr_vis = len(self.tt_visits)
        self.tt_visits.add_index("vis_id")

        # Dictionary to collect lightcurve data
        tdata = {
            "fd_src_id": list(),
            "flux": list(),
            "flux_err": list(),
            "ul": list(),
            "time_bin_start": list(),
            "time_bin_size": list(),
            "ra": list(),
            "dec": list(),
            "pos_err": list(),
        }

        # Loop over sources and add them to dictionary
        tt_det_grp = self.tt_detections[sel].group_by(["fd_src_id"])
        for tt_det in tt_det_grp.groups:
            # Add src id
            fd_src_id = tt_det["fd_src_id"][0]
            tdata["fd_src_id"].append(fd_src_id.tolist())

            # Add flux and errors, array has length of Nr of visits.
            # Values is 0, except for detections.
            vis_idxs = self.tt_visits.loc_indices["vis_id", tt_det["vis_id"]]

            det_vars = ["flux", "flux_err", "ra", "dec", "pos_err"]

            for det_var in det_vars:
                np_var = np.zeros(nr_vis) - 1
                np_var[vis_idxs] = tt_det[det_var]
                tdata[det_var].append(np_var.tolist())

            tdata["time_bin_start"].append(
                np.array(self.tt_visits["time_bin_start"]).tolist()
            )
            tdata["time_bin_size"].append(
                np.array(self.tt_visits["time_bin_size"]).tolist()
            )

            # Store upper limits if no detection in a visit
            # TODO: make this more general and not GALEX specific
            if add_upper_limits:
                np_flux_ul = self.get_upper_limits()
                tdata["ul"].append(np_flux_ul.tolist())

        # Add light curve table, convert type to allow for variable length array
        tdata["flux"] = np.array(tdata["flux"], dtype=np.object_)
        tdata["flux_err"] = np.array(tdata["flux_err"], dtype=np.object_)
        tdata["ul"] = np.array(tdata["ul"], dtype=np.object_)
        tdata["time_bin_start"] = np.array(tdata["time_bin_start"], dtype=np.object_)
        tdata["time_bin_size"] = np.array(tdata["time_bin_size"], dtype=np.object_)

        self.add_table(tdata, "base_field:ta_sources_lc")

    def remove_double_visit_detections(self):
        """
        Remove multiple detections of one source in one visit from tt_detections.
        Keep only the closest detection.

        Returns
        -------
        list
            List of removed detection IDs

        """
        logger.debug("Scanning for doubled visit detections")

        # Selection
        sel = self.tt_detections["sel"]

        # Index tables
        self.tt_sources.add_index("fd_src_id")
        self.tt_detections.add_index("det_id")

        # Determine detection_id of srcs with more than one detection in one visit
        rm_det_ids = []
        tt_det_grp = self.tt_detections[sel].group_by(["vis_id", "fd_src_id"])

        for tt_det in tt_det_grp.groups:
            if len(tt_det) > 1:
                # Get source coordinate
                fd_src_id = tt_det["fd_src_id"].data[0]
                fd_src_idx = self.tt_sources.loc_indices["fd_src_id", fd_src_id]
                src_ra = self.tt_sources["ra"].quantity[fd_src_idx]
                src_dec = self.tt_sources["dec"].quantity[fd_src_idx]
                src_coord = SkyCoord(src_ra, src_dec, frame="icrs")

                # Determine distance
                det_coords = SkyCoord(
                    tt_det["ra"].quantity, tt_det["dec"].quantity, frame="icrs"
                )
                sep = det_coords.separation(src_coord).to(uu.arcsec)
                min_sep_idx = np.where(sep == np.amin(sep))[0][0]

                # Save all detection but the closest for deletion
                rm_det_ids.extend(tt_det["det_id"].data[:min_sep_idx])
                rm_det_ids.extend(tt_det["det_id"].data[min_sep_idx + 1 :])

        if len(rm_det_ids) > 0:
            nr_rm_det = len(rm_det_ids)
            perc_rm_det = 100 * nr_rm_det / len(self.tt_detections[sel])
            logger.warning(
                f"Removed double-visit detections: {nr_rm_det} ({perc_rm_det: .4f} %)"
            )

            # remove the doubled detections from tt_detections
            rm_idx = self.tt_detections.loc_indices["det_id", rm_det_ids]
            self.tt_detections.remove_rows(rm_idx)

            # Update detection count in tt_sources
            sel = self.tt_detections["sel"]
            fd_src_ids, det_cts = np.unique(
                self.tt_detections[sel]["fd_src_id"], return_counts=True
            )
            self.tt_sources.replace_column("nr_det", det_cts)

        return rm_det_ids

    def set_field_attr(self, dd_names=None):
        """
        Sets the most important field parameters as class attributes.

        Parameters
        ----------
        dd_names : dict, None
            List of strings containing the parameter names. If None (default),
            a set of pre-defined parameters is set.

        Returns
        -------
        None
        """

        if dd_names is None:
            dd_names = self._dd_attr_names
        for table in dd_names:
            for par in dd_names[table]:
                setattr(self, par, self.get_field_par(par, table))

        logger.debug(f"Set attributes: {self._dd_attr_names}")

    def get_field_par(self, par_name, table_name):
        """
        Returns metadata parameter ``par_name``
        derived from astropy table data ``table_name``.

        Parameter
        ---------
        par_name : str
            Name of the parameter to be read either directly or calculated
            from table data.
        table_name : str
            Name of table data. Must be initialized as class attribute.

        Returns
        -------
        str, None
            None is returned if either ``par_name`` or ``table_name`` are
            not available. In this case warnings are issued to logging.

        Raises
        ------
        AssertionError
            If :attr: `~vasca.field.BaseField.tt_fields` has not exactly one row

        """
        # Check data availability
        # Returns none if table is unknown
        if table_name not in self.__dict__:
            logger.warning(
                f"Data '{table_name}' not available. "
                f"Returning None for parameter '{par_name}'."
            )
            return None

        # Consistency checks
        # tt_fields must be single-rowed
        if table_name == "tt_fields" and len(self.tt_fields) != 1:
            raise ValueError(
                "Inconsistent data. Expeted single-rowed "
                f"field metadata table 'tt_fields', got {len(self.tt_fields)}."
            )

        # Indirectly derived parameters
        if par_name == "center":
            par = SkyCoord(
                self.get_field_par("ra", "tt_fields"),
                self.get_field_par("dec", "tt_fields"),
                frame="icrs",
            )
        elif par_name == "nr_vis":
            par = len(self.tt_visits)
        elif par_name == "time_bin_size_sum":
            par = self.tt_visits["time_bin_size"].sum() * uu.s
        elif par_name == "time_start":
            par = Time(self.tt_visits["time_bin_start"][0], format="mjd")
        elif par_name == "time_stop":
            par = Time(
                self.tt_visits["time_bin_start"][-1]
                + self.tt_visits["time_bin_size"][-1],
                format="mjd",
            )
        # Directly derived parameters
        # Return None if parameter is not in table
        elif par_name not in getattr(self, table_name).colnames:
            logger.warning(
                f"Unknown parameter in data set '{table_name}'. "
                f"Known parameters are {self.__dict__[table_name].colnames}. "
                f"Returning None for parameter '{par_name}'."
            )
            par = None
        else:
            # retrieve parameter
            par = self.tt_fields[par_name].data[0]
            # Maintains unit
            # Some combinations of unit and dtype need special handling
            unit = self.tt_fields[par_name].unit
            dtype_kind = self.tt_fields[par_name].dtype.kind
            logger.debug(f"Getting parameter '{par_name}': {par}, {unit}, {dtype_kind}")
            # check (byte-)string or unicode
            if dtype_kind in ["U", "S"]:
                par = par.decode("utf-8")
            # check uu.dimensionless_unscaled and integer
            elif (unit == "" or unit is None) and dtype_kind in ["i", "u"]:
                par = int(par)
            else:
                par = par * unit

        return par

    def get_sky_region(self):
        """
        Get region on the sky of the field.

        Returns
        -------
        regions.SkyRegion
            Region on the sky of the field.
        """
        if self.observatory.casefold() in ["GALEX".casefold(), "GALEX_DS".casefold()]:
            return CircleSkyRegion(center=self.center, radius=self.fov_diam / 2.0)
        else:
            logger.warning(f"No region known for observatory {self.observatory}")
            return None

    def load_from_fits(self, file_name="tables.fits"):
        """
        Loads field from a fits file and sets field attributes.

        Parameters
        ----------
        file_name : str, optional
            File name. The default is "field_default.fits".

        Returns
        -------
        None

        """
        super().load_from_fits(file_name)
        self.set_field_attr()


class GALEXField(BaseField):
    """
    Instance of one GALEX field
    """

    def __init__(self, obs_id, obs_filter=None, data_path=None, visits_data_path=None):
        """
        Initializes a new GALEXField instance with
        skeleton VASCA data structure.

        Parameters
        ----------
        obs_id : int
            GALEX field ID
        obs_filter : str, optional
            Selects the GALEX obs_filter for which the corresponding
            observation data is loaded. Needs to be either from:
            'FUV' -> 135-175 nm
            'NUV' -> 175-280 nm  (default)
        data_path : str, optional
            Path to root location of cloud-synced data associated with the GALEX field.
            Defaults to a path given by the resource manager.
        visits_data_path : str, optional
            Path to a pre-downloaded table holding the complete list of GALEX visits.
            Defaults to a path given by the resource manager.

        Attributes
        ----------
        data_path : str
            Path to root location of cloud-synced data associated with the GALEX field
        visits_data_path : str
            Path to a pre-downloaded table holding the complete list of GALEX visits
        vasca_file_prefix : str
            File name prefix following VASCA naming convention:
            'VASCA_<observatory>_<filed_id>_<obs_filter>'
        """

        # check obs_filter name, default is "NUV"
        if (
            obs_filter is None
        ):  # TODO: I think this is redundant, the default should then be "NUV"
            obs_filter = "NUV"
        elif not isinstance(obs_filter, str):
            raise TypeError(
                f"Expected argument of type 'str', got '{type(obs_filter).__name__}'."
            )
        elif obs_filter not in ["NUV", "FUV"]:
            raise ValueError(
                f"Unknown obs_filter '{obs_filter}'. "
                "Available filters for GALEX are 'NUV' or 'FUV'."
            )

        # Sets skeleton
        super().__init__()

        # Set root location of cloud-synced data associated with the field
        if data_path is not None:
            # Sets the data path following the pattern <root_path>/<obs_id>
            # If only <root_path> is passed, <obs_id> will be added
            self.data_path = (
                data_path
                if data_path.rstrip(os.sep).endswith(str(obs_id))
                else f"{data_path}/{obs_id}"
            )
        if visits_data_path is not None:
            self.visits_data_path = visits_data_path
        if data_path is None or visits_data_path is None:
            with ResourceManager() as rm:
                # Path to read and write data relevant to the pipeline
                if data_path is None:
                    self.data_path = (
                        rm.get_path("gal_fields", "sas_cloud") + "/" + str(obs_id)
                    )
                # Path to a pre-downloaded table holding
                # the complete list of GALEX visits
                if visits_data_path is None:
                    self.visits_data_path = rm.get_path("gal_visits_list", "sas_cloud")

        # Create and check existence of directory
        # that holds field data and VASCA outputs
        os.makedirs(self.data_path, exist_ok=True)

        # File name prefix for VASCA/GALEXField outputs
        self.vasca_file_prefix = f"VASCA_GALEX_{obs_id}_{obs_filter}"

        logger.debug(f"Field data path set to: '{self.data_path}'")
        logger.debug(f"Visits data path set to: '{self.visits_data_path}'")

    @classmethod
    def from_VASCA(cls, obs_id, obs_filter="NUV", fits_path=None, **kwargs):
        """
        Constructor to initialize a GALEXField instance
        from a VASCA-generated FITS file

        Parameters
        ----------
        obs_id : int
            GALEX field ID
        obs_filter : str, optional
            Selects the GALEX filter for which the corresponding
            observation data is loaded. Needs to be either from:
            'FUV' -> 135-175 nm
            'NUV' -> 175-280 nm  (default)
        fits_path : str, optional
            Path to the fits file. Defaults to a path handled by ResourceManager.
        **kwargs
            All additional keyword arguments are passed to `~GALEXField.__init__()`

        Returns
        -------
        vasca.field.GALEXField
        """
        # Bootstrap the initialization procedure using the base class
        gf = cls(obs_id, obs_filter, **kwargs)  # new GALEXField instance

        if fits_path is None:
            # Construct the file name from field ID and filter
            fits_path = f"{gf.data_path}/{gf.vasca_file_prefix}_field_data.fits"
        # Check if file exists
        if not os.path.isfile(fits_path):
            raise FileNotFoundError(
                "Wrong file or file path to VASCA data "
                f"for GALEX field '{obs_id}' with obs_filter '{obs_filter}'."
            )
        else:
            # Reads the VASCA-generated field data
            gf.load_from_fits(fits_path)
            # Sets convenience class attributes
            gf.set_field_attr()
        # field_id = get_field_id(
        #     obs_field_id=obs_id, observatory="GALEX", obs_filter=obs_filter
        # )

        # # Check consistency
        # if not gf.field_id == field_id:
        #     raise ValueError(
        #         "Inconsistent data: Missmatch for 'field_id'. "
        #         f"Expected '{obs_id}' but got '{gf.field_id}' "
        #         f"from file '{fits_path.split(os.sep)[-1]}."
        #     )
        # elif not gf.obs_filter == obs_filter:
        #     raise ValueError(
        #         "Inconsistent data: Missmatch for 'obs_filter'. "
        #         f"Expected '{obs_filter}' but got '{gf.obs_filter}' "
        #         f"from file '{fits_path.split(os.sep)[-1]}."
        #     )

        # logger.info(
        #     f"Loaded VASCA data for GALEX field '{obs_id}' "
        #     f"with obs_filter '{obs_filter}'."
        # )

        return gf

    @classmethod
    def from_MAST(
        cls,
        obs_id,
        obs_filter="NUV",
        refresh=False,
        load_products="TABLES",
        write=True,
        **kwargs,
    ):
        """
        Constructor to initialize a GALEXField instance either
        fresh from the MAST archive (refresh=True) or if available
        from cached raw data (refresh=False).

        The procedure uses vasca.resource_manager.ResourceManager
        to handle file locations.

        Parameters
        ----------
        obs_id : int
            GALEX field ID
        obs_filter : str, optional
            Selects the GALEX obs_filter for which the corresponding
            observation data is loaded. Needs to be either from:
            'FUV' -> 135-175 nm
            'NUV' -> 175-280 nm  (default)
        refresh : bool, optional
            Selects if data is freshly loaded from MAST (refresh=True) or
            from cashed data on disc (refresh=False, default).
        load_products : str, optional
            Selects if data products shall be loaded: "NONE", "TABLES" for tables,
            and reference image, "ALL" for tables and visit images.
        write: bool, optional
            If load_products is enabled, stores the data as VASCA tables in the cloud
            for faster loading in the future. Default is True.
        **kwargs
            All additional keyword arguments are passed to `~GALEXField.__init__()`

        Returns
        -------
        vasca.field.GALEXField
            Returns None if no entry is found
        """
        # Checks
        if not isinstance(refresh, bool):
            raise TypeError(f"Expected boolean argument, got {type(refresh).__name__}.")

        # Bootstrap the initialization procedure using the base class
        gf = cls(obs_id, obs_filter, **kwargs)  # new GALEXField instance

        # Sets ``gf.tt_fields``
        gf._load_galex_field_info(obs_id, obs_filter, refresh=refresh)

        # If no field entries returned return None
        if len(gf.tt_fields) == 0:
            logger.warning(f"No field found for ID {obs_id} and filter {obs_filter}")
            return None
        # Sets ``gf.tt_visits``
        gf._load_galex_visits_info(obs_id, obs_filter)

        # Sets convenience class attributes
        gf.set_field_attr()

        # Sets ``gf.tt_detections``, ``gf.tt_coadd_detections`` and loads the ref image
        if load_products != "NONE":
            if load_products == "ALL":
                gf._load_galex_archive_products(
                    obs_id, obs_filter, refresh=refresh, ref_maps_only=False
                )
            elif load_products == "TABLES":
                gf._load_galex_archive_products(
                    obs_id, obs_filter, refresh=refresh, ref_maps_only=True
                )
            else:
                logger.warning(f"load_product option not valid: {load_products}")
            if write:
                fits_name = f"{gf.data_path}/{gf.vasca_file_prefix}_field_data.fits"
                gf.write_to_fits(fits_name)

        meta_only = ", metadata-only." if load_products == "NONE" else "."
        logger.info(
            f"Loaded new GALEX field '{obs_id}' with obs_filter '{obs_filter}'"
            f"from MAST data {meta_only}"
        )

        return gf

    @staticmethod
    def load(
        gfield_id,
        obs_filter="NUV",
        method="MAST_LOCAL",
        load_products="TABLES",
        **field_kwargs,
    ):
        """
        Loads GALEX field data according to a given method and
        returns a GALEXField instance.

        Parameters
        ----------
        gfield_id : int
            GALEX field ID
        obs_filter : str, optional
            Selects the GALEX obs_filter for which the corresponding
            observation data is loaded. Needs to be either from:
            'FUV' -> 135-175 nm
            'NUV' -> 175-280 nm  (default)
        method : str, optional
            Specification of the load method. Four methods are implemented:
            MAST_REMOTE: Standard method to query data from MAST archives. This requires
            an internet connection to the MAST servers. Local data will be overwritten
            MAST_LOCAL: Dafault. Builds a new GALEXfield instance based on MAST
            archival data cached on local disc. If no data is found,
            the fallback is MAST_REMOTE.
            VASCA: Builds a new GALEXField based on a VASCA-generated field data file.
            AUTO: Attempts to load field data by using VASCA as method and
            falls back to MAST_LOCAL if no VASCA file is found on disc.

            The default directory where field data availability is checked is
            defined by the "data_path" attribute of GALEXField and can be passed
            explicitly via the "field_kwargs".
        load_products : str, optional
            Selects if data products shall be loaded: "NONE", "TABLES" for tables,
            and reference image, "ALL" for tables and visit images.
        field_kwargs
            All additional keyword arguments are passed to the load methods
            `~GALEXField.from_MAST()` and `~GALEXField.from_VASCA()`,
            as well as to `~GALEXField.__init__()`.

        Raises
        ------
        TypeError
            If the specified load method is not a string.
        ValueError
            If the specified load method is not one of
            '["mast_remote", "mast_local", "vasca", "auto"]'. String matching is
            case insensitive.

        Returns
        -------
        vasca.field.GALEXField

        """
        logger.info(
            f"Loading field '{gfield_id}' with method '{method}' for filter '{obs_filter}' and load_products '{load_products}'."
        )

        # Checks method argument
        # String matching is case insensitive (converts to all to lower case).
        if not isinstance(method, str):
            raise TypeError(
                "Expected string type for load method specification, "
                f"got '{type(method).__name__}'."
            )

        method = method.casefold()  # ensures all-lower-case string
        method_spec = ["mast_remote", "mast_local", "vasca", "auto"]
        if method not in method_spec:
            raise ValueError(
                f"Expected load method specification from {method_spec}, "
                f"got '{method}'."
            )

        # Sets refresh option for MAST methods
        if method == "mast_remote":
            refresh = field_kwargs.pop("refresh", True)
        else:
            refresh = field_kwargs.pop("refresh", False)

        # Loads field according to load method specification
        if method in ["mast_remote", "mast_local"]:
            gf = GALEXField.from_MAST(
                obs_id=gfield_id,
                obs_filter=obs_filter,
                refresh=refresh,
                load_products=load_products,
                **field_kwargs,
            )
        elif method == "vasca":
            # removes unused options
            field_kwargs.pop("refresh", None)
            field_kwargs.pop("write", None)

            try:
                gf = GALEXField.from_VASCA(
                    obs_id=gfield_id,
                    obs_filter=obs_filter,
                    **field_kwargs,
                )
            except FileNotFoundError as e:
                logger.warning(
                    f"Load method '{method}' failed due to FileNotFoundError ('{e}'). "
                    "Falling back to method 'mast_remote'."
                )
                gf = GALEXField.from_MAST(
                    obs_id=gfield_id,
                    obs_filter=obs_filter,
                    refresh=True,
                    load_products=load_products,
                    **field_kwargs,
                )
        elif method == "auto":
            raise NotImplementedError(f"method '{method}' not operational.")
            # TODO: Lookahead via rm to check data availability.
            # Then "VASCA" is preferred for performance reasons.
            # Fallback to "MAST" & refresh=True if "VASCA" fails for some reason
            # (e.g. not complete set of tables stored in the fits file).
            # Update (2023-01-11): Reevaluate if this option is still needed

        return gf

    @staticmethod
    def get_visit_upper_limits(tt_visits):
        """
        Calculates upper limits on non-detections to the tt_visits table
        following the article_ of Gezari et al. ApJ 766 (2013) p4

        .. _article: https://iopscience.iop.org/article/10.1088/0004-637X/766/1/60

        Returns
        -------
        float array
            Array of upper limits for each visit

        """

        B_sky = 3e-3
        N_pix = 16 * np.pi
        C_app = -0.23
        T_exp = tt_visits["time_bin_size"].data
        upper_limit = -2.5 * np.log(5 * (B_sky * N_pix / T_exp)) + C_app
        return upper_limit

    def _load_galex_field_info(self, obs_id, obs_filter, col_names=None, refresh=False):
        """
        Loads the archival metadata associated to a given field ID.

        Parameters
        ----------
        obs_id : int
            GALEX field ID
        obs_filter : str, optional
            Selects the GALEX filter for which the corresponding
            observation data is loaded. Needs to be either from:
            'FUV' -> 135-175 nm
            'NUV' -> 175-280 nm  (default)
        col_names : dict
            Dictionary with keys of MAST column names to load, and values the
            corresponding VASCA table column names. Default in None.
        refresh : bool, optional
            Selects if data is freshly loaded from MAST (refresh=True) or
            from cashed data on disc (refresh=False, default).

        """

        # Uses default columns if not otherwise specified
        # Already sets the order in which columns are added later on
        if col_names is None:
            # values represent variables in VASCA, keys the variable names
            # in the MAST database
            col_names = {
                "obs_id": "field_id",
                "target_name": "field_name",
                "s_ra": "ra",
                "s_dec": "dec",
                "instrument_name": "observatory",
                "filters": "obs_filter",
                "project": "project",
            }

        elif not isinstance(col_names, dict):
            raise TypeError(
                "Expected list type for argument 'col_names', "
                f"got '{type(col_names).__name__}'."
            )

        mast_col_names = list(col_names.keys())

        # Path to store raw data
        path_tt_coadd = f"{self.data_path}/MAST_{obs_id}_{obs_filter}_coadd.fits"

        # Reads cached file if found on disc
        if not os.path.isfile(path_tt_coadd) or refresh:
            # download raw table
            logger.debug(
                f"Downloading archive field info and saving to {path_tt_coadd}."
            )
            tt_coadd = Observations.query_criteria(
                obs_id=obs_id,
                dataRights="PUBLIC",
                instrument_name="GALEX",
                filters=obs_filter,
                dataproduct_type="image",
            )
            # save to disc
            tt_coadd.write(path_tt_coadd, overwrite=True)
        else:
            # read cached
            logger.debug(
                f"Reading archive field info from cashed file '{path_tt_coadd}'"
            )
            tt_coadd = Table.read(path_tt_coadd)

        # Constructs field info table
        logger.debug("Constructing 'tt_fields'.")

        # Selects columns and strip mask
        if tt_coadd.masked is True:
            tt_coadd_select = Table(tt_coadd[mast_col_names].as_array().data)
        else:
            tt_coadd_select = tt_coadd[mast_col_names]

        # Converts field id dtype from unicode ('U19') or S to bytestring ('S64')
        obsid_kind = tt_coadd_select["obs_id"].dtype.kind
        if obsid_kind == "U" or obsid_kind == "S":
            tt_coadd_select.replace_column(
                "obs_id",
                tt_coadd_select["obs_id"].astype(np.dtype("S64")),
            )

        # Converts obs_id colimn  into VASCA field_id
        for row in tt_coadd_select:
            row["obs_id"] = get_field_id(
                obs_field_id=row["obs_id"], observatory="GALEX", obs_filter=obs_filter
            )

        # Convert into dictionary with correct VASCA column names
        dd_coadd_select = {}
        for col in mast_col_names:
            dd_coadd_select[col_names[col]] = tt_coadd_select[col].data

        # Add fov size info
        dd_coadd_select["fov_diam"] = 1.2 * np.ones(len(tt_coadd_select["obs_id"]))

        self.add_table(dd_coadd_select, "base_field:tt_fields")

    def _load_galex_visits_info(self, obs_id, obs_filter, col_names=None):
        """
        Loads the archival metadata associated to a given visit ID.

        Parameters
        ----------
        obs_id : int
            GALEX field ID
        obs_filter : str, optional
            Selects the GALEX filter for which the corresponding
            observation data is loaded. Needs to be either from:
            'FUV' -> 135-175 nm
            'NUV' -> 175-280 nm  (default)
        col_names : dict
            Dictionary with keys of MAST column names to load, and values the
            corresponding VASCA table column names. Default in None.
        """
        logger.debug("Loading GALEX visit info :" f"{obs_id} {obs_filter} {col_names}")

        # Uses default columns if not otherwise specified
        # Already sets the order in which columns are added later on
        if col_names is None:
            mast_col_names = [
                "imgRunID",
                "minPhotoObsDate",
                "nexptime" if obs_filter == "NUV" else "fexptime",
                "fexptime" if obs_filter == "NUV" else "nexptime",
                "RATileCenter",
                "DECTileCenter",
            ]
            vasca_col_names = [
                "vis_id",
                "time_bin_start",
                "time_bin_size",
                "time_bin_size_alt_filt",
                "ra",
                "dec",
            ]
            col_names = dict(zip(mast_col_names, vasca_col_names))
        elif not isinstance(col_names, dict):
            raise TypeError(
                "Expected list type for argument 'col_names', "
                f"got '{type(col_names).__name__}'."
            )

        mast_col_names = list(col_names.keys())

        # Reads cached file if found on disc
        logger.debug(
            "Reading archive visit info from cashed file " f"'{self.visits_data_path}'"
        )
        tt_visits_raw = Table.read(self.visits_data_path)

        # Filters for visits corresponding to field id and selects specified columns
        tt_visits_raw_select = tt_visits_raw[tt_visits_raw["ParentImgRunID"] == obs_id][
            mast_col_names
        ]

        # Remove zero exposure visits
        exp_col_name = "nexptime" if obs_filter == "NUV" else "fexptime"
        sel_visits = tt_visits_raw_select[exp_col_name] > 1e-10
        tt_visits_raw_select = tt_visits_raw_select[sel_visits]

        if len(tt_visits_raw_select) == 0:
            logger.warning(f"No exposure for field {obs_id} and filter  {obs_filter}")

        # Check which filters where observing
        obs_filter_ids = np.zeros(len(tt_visits_raw_select))
        obs_filter_ids += (tt_visits_raw_select["nexptime"] > 1e-10) * dd_filter2id[
            "NUV"
        ]
        obs_filter_ids += (tt_visits_raw_select["fexptime"] > 1e-10) * dd_filter2id[
            "FUV"
        ]

        # Converts string time format to mjd
        ll_times = [
            datetime.strptime(date_str.decode("utf-8"), "%m/%d/%Y %I:%M:%S %p")
            for date_str in tt_visits_raw_select["minPhotoObsDate"].data
        ]
        tt_visits_raw_select.update({"minPhotoObsDate": Time(ll_times).mjd})

        # Convert into dictionary with correct VASCA column names
        dd_visits_raw_select = {}
        for col in mast_col_names:
            dd_visits_raw_select[col_names[col]] = tt_visits_raw_select[col].data
        dd_visits_raw_select["obs_filter_id"] = obs_filter_ids

        # Sets table as class attribute
        self.add_table(dd_visits_raw_select, "galex_field:tt_visits")

        # Add filter_id
        # flt_ids = dd_filter2id[obs_filter] * np.ones(len(self.tt_visits))
        # self.add_column("tt_visits", "obs_filter_id", flt_ids)

    def _load_galex_archive_products(
        self,
        obs_id,
        obs_filter,
        col_names=None,
        dd_products=None,
        ref_maps_only=False,
        refresh=False,
    ):
        """
        Loads the relevant data products from MAST servers and
        stores them using the ResourceManager.

        Parameters
        ----------
        obs_id : int
            GALEX field ID.
        obs_filter : str
            Selects the GALEX filter for which the corresponding
            observation data is loaded. Needs to be either from ['NUV', 'FUV'].
        col_names : list, optional
            List of columns to store from raw data.
        dd_products : dict, optional
            Dictionary of listing the data products to be loaded.
        ref_maps_only : bool, optional
            If True, only the reference/coadd maps are loaded (default).
            If False, also the visit maps are included.
        refresh : bool, optional
            Set True to refresh the raw data via a MAST query. On default
            cashed raw data is used.

        """
        # Path to MAST helper files
        path_tt_coadd = f"{self.data_path}/MAST_{obs_id}_{obs_filter}_coadd.fits"
        path_tt_data = f"{self.data_path}/MAST_{obs_id}_{obs_filter}_data.fits"
        path_tt_down = f"{self.data_path}/MAST_{obs_id}_{obs_filter}_down.ecsv"

        # tt_data: List of all archival data products
        # Reads cached file if found found on disc
        if not os.path.isfile(path_tt_data) or refresh:
            # Get cached field info for archive query input
            logger.debug(
                f"Reading archive field info from cashed file '{path_tt_coadd}'"
            )
            tt_coadd = Table.read(path_tt_coadd)
            # Download
            logger.debug(
                f"Downloading archive data products list and saving to {path_tt_data}."
            )
            # -> Todo: Create a generic decorator wrapper for Astroquery requests
            try:
                tt_data = Table(
                    Observations.get_product_list(tt_coadd).as_array().data
                )  # Returns unmasked table
            except HTTPError as e:
                # Need to check its an 404, 503, 500, 403 etc.
                status_code = e.response.status_code

                # Try again
                if status_code == 503:
                    logger.info(
                        "HTTPError Astroquery response "
                        "503 Server Error 'Service Unavailable'"
                    )
                    sleep_time = 5
                    logger.info(f"Retrying Astroquery request in {sleep_time} s.")
                    time.sleep(sleep_time)
                    tt_data = Table(
                        Observations.get_product_list(tt_coadd).as_array().data
                    )  # Returns unmasked table
            # Saves to disc
            tt_data.write(path_tt_data, overwrite=True)
        else:
            logger.debug(
                f"Reading archive data products list from cashed file '{path_tt_data}'"
            )
            tt_data = Table.read(path_tt_data)

        # Sets default products to download if a selection was not explicitly passed
        if dd_products is None:
            # Default: Catalog and NUV intensity map
            dd_products = {
                0: {
                    "name": "mcat",
                    "file_name": "xd-mcat.fits.gz",
                    "product_type": "catalog",
                },
                1: {
                    "name": "int_map",
                    "file_name": (
                        "nd-int.fits.gz" if obs_filter == "NUV" else "fd-int.fits.gz"
                    ),
                    "product_type": "map",
                },
            }
        # Verifies if products are listed in products list
        product_list_missing = [
            prod["file_name"]
            for _, prod in dd_products.items()
            if not any(prod["file_name"] in uri for uri in tt_data["dataURI"])
        ]
        if not len(product_list_missing) == 0:
            raise ValueError(
                f"Missing entry for data products {product_list_missing} "
                "in product list."
            )

        # Filters for products of interest
        aa_sel_prod = np.full(len(tt_data), False)  # bool array
        # Product files
        aa_data_uri = tt_data["dataURI"].data.astype(str)  # string array
        # obs_id
        aa_obs_id = tt_data["obs_id"].data.astype(int)  # int array

        for _, prod in dd_products.items():
            if prod["product_type"] == "map" and ref_maps_only:
                aa_sel_prod += np.logical_and(
                    np.char.endswith(aa_data_uri, prod["file_name"]),
                    aa_obs_id == obs_id,
                )
            else:
                aa_sel_prod += np.char.endswith(
                    aa_data_uri, prod["file_name"]
                )  # Logical AND operation

        # Download data products
        # tt_down: manifest of files downloaded
        if not os.path.isfile(path_tt_down) or refresh:
            # Download
            logger.debug(
                f"Downloading archive data products. Manifest saved to {path_tt_down}."
            )
            tt_down = Observations.download_products(
                tt_data,
                download_dir=self.data_path,
                cache=True,
                dataURI=tt_data[aa_sel_prod]["dataURI"],
            )["Local Path", "Status"]
            # Saves table to disc
            tt_down.replace_column(
                "Local Path",
                [str(s).split(self.data_path)[1] for s in tt_down["Local Path"]],
            )  # Keeps path only relative to self.data_path
            tt_down.write(path_tt_down, overwrite=True)
        else:
            # Reading cached manifest
            logger.debug(f"Reading archive data product manifest from {path_tt_down}.")
            tt_down = Table.read(path_tt_down)

        # Adds column with corresponding IDs
        tt_down["ID"] = [int(path.split(os.sep)[-2]) for path in tt_down["Local Path"]]

        # Verifies completed download
        aa_status_mask = np.where(
            tt_down["Status"].data.astype(str) == "COMPLETE", True, False
        )  # bool array
        if not aa_status_mask.all():
            raise ValueError(
                "Data products download incomplete. "
                "Missing products are:\n"
                f"{tt_down[~aa_status_mask]}"
            )

        # Collects data products

        # Source catalogs
        aa_sel_mcat = np.char.endswith(
            tt_down["Local Path"].data.astype(str), "xd-mcat.fits.gz"
        )
        # Selects only the coadd path
        path_tt_ref = (
            self.data_path
            + tt_down[np.logical_and(aa_sel_mcat, tt_down["ID"] == obs_id)][
                "Local Path"
            ][0]
        )  # String
        # Selects everything else but the coadd path
        path_tt_detections = [
            self.data_path + path
            for path in tt_down[np.logical_and(aa_sel_mcat, tt_down["ID"] != obs_id)][
                "Local Path"
            ].data.astype(str)
        ]  # list

        # Opens reference catalog
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            tt_coadd_detections_raw = Table.read(path_tt_ref)

        # Opens visit catalogs, add columns for visit and source IDs
        # and stack all catalogs
        tt_detections_raw = None
        # loop over visits
        for path, id in zip(path_tt_detections, self.tt_visits["vis_id"]):
            # read new visits catalog

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", AstropyWarning)
                tt_vis_mcat = Table.read(path)

            # set initial
            if tt_detections_raw is None:
                # column visit id
                tt_vis_mcat.add_column(
                    np.full(len(tt_vis_mcat), id), name="visit_id", index=0
                )
                # column source id
                tt_vis_mcat.add_column(
                    np.full(len(tt_vis_mcat), 999999999), name="source_id", index=1
                )
                tt_detections_raw = tt_vis_mcat
            # append
            else:
                # column visit id
                tt_vis_mcat.add_column(
                    np.full(len(tt_vis_mcat), id), name="visit_id", index=0
                )
                # column source id
                tt_vis_mcat.add_column(
                    np.full(len(tt_vis_mcat), 999999999), name="source_id", index=1
                )
                tt_detections_raw = vstack([tt_detections_raw, tt_vis_mcat])
            # clean up
            del tt_vis_mcat

        obs_filter_l = obs_filter.lower()  # lower case obs_filter name

        # Add to tables as class attributes
        if col_names is None:
            mast_col_names = [
                "visit_id",
                "ggoid_dec",
                f"{obs_filter}_ALPHA_J2000",  # take band-merged quantities?
                f"{obs_filter}_DELTA_J2000",  # take band-merged quantities?
                f"{obs_filter_l}_poserr",
                f"{obs_filter_l}_flux",
                f"{obs_filter_l}_fluxerr",
                f"{obs_filter_l}_s2n",
                "fov_radius",
                f"{obs_filter_l}_artifact",
                f"{obs_filter}_CLASS_STAR",
                "chkobj_type",
                f"{obs_filter}_ELLIPTICITY",
                f"{obs_filter}_FLUX_APER_4",
                f"{obs_filter}_FLUXERR_APER_4",
                f"{obs_filter}_FLUX_APER_3",
                f"{obs_filter}_FLUXERR_APER_3",
                f"{obs_filter}_FLUX_AUTO",
                "E_bv",
                f"{obs_filter_l}_mag",
                f"{obs_filter_l}_magerr",
            ]

            vasca_col_names = [
                "vis_id",
                "det_id",
                "ra",
                "dec",
                "pos_err",
                "flux_auto",
                "flux_auto_err",
                "s2n",
                "r_fov",
                "artifacts",
                "class_star",
                "chkobj_type",
                "ellip_world",
                "flux_r60",
                "flux_r60_err",
                "flux_r38",
                "flux_r38_err",
                "flux_auto_flt",
                "E_bv",
                "mag",
                "mag_err",
            ]

            col_names = dict(zip(mast_col_names, vasca_col_names))
        elif not isinstance(col_names, dict):
            raise TypeError(
                "Expected dict type for argument 'col_names', "
                f"got '{type(col_names).__name__}'."
            )
        mast_col_names = list(col_names.keys())

        # set data as class attributes

        # Correct for flux outside of aperture, as listed here:
        # http://www.galex.caltech.edu/researcher/techdoc-ch5.html
        acorr60 = 1.116863247 if obs_filter == "NUV" else 1.09647819

        # acorr38 = 1.21338885 if obs_filter == "NUV" else 1.202264
        def get_det_dict(tt_det):
            """Retrieve dictionary with detection/coadd_detections data,
            including some variables derived from other MAST variables"""

            # Convert into dictionary with correct VASCA column names
            # *** Detections
            dd_det = {}
            # Keep only entries with detections
            sel_s2n = tt_det[f"{obs_filter_l}_s2n"] > 0

            if len(tt_det[sel_s2n]) == 0:
                logger.warning(
                    f"No detection with s2n>0 in MAST for field {obs_id} filter {obs_filter}"
                )

            for col in mast_col_names:
                if col not in tt_det.colnames:
                    continue
                dd_det[col_names[col]] = tt_det[col][sel_s2n].data

            # Add size_world in arcsec
            dd_det["size_world"] = (
                3600
                * (
                    tt_det[f"{obs_filter}_A_WORLD"][sel_s2n].data
                    + tt_det[f"{obs_filter}_B_WORLD"][sel_s2n].data
                )
                / 2
            )

            # Add flux aperture ratio and covert counts to Jansky for r=3.8 arcsec aperture flux
            # Check to avoid divisions by zero
            idx_fauto = np.where(dd_det["flux_auto_flt"] > 0)
            cts2Jy = np.zeros(len(dd_det["flux_auto_flt"]))
            cts2Jy[idx_fauto] = (
                dd_det["flux_auto"][idx_fauto] / dd_det["flux_auto_flt"][idx_fauto]
            )

            # Use aperture flux instead of AUTO flux, as it is more reliable for variability studies
            dd_det["flux"] = dd_det["flux_r60"] * cts2Jy * acorr60
            dd_det["flux_err"] = dd_det["flux_r60_err"] * cts2Jy * acorr60

            # Add filter_id
            dd_det["obs_filter_id"] = np.array(
                dd_filter2id[obs_filter.upper()] + np.zeros(np.sum(sel_s2n))
            )

            # Ratio between the flux-counts, correcting for difference in aperture
            idx_fr60 = np.where(dd_det["flux_r60"] > 0)
            dd_det["flux_app_ratio"] = (
                np.zeros(
                    len(dd_det["flux_r60"]),
                    dtype=dd_vasca_columns["flux_app_ratio"]["dtype"],
                )
                + dd_vasca_columns["flux_app_ratio"]["default"]
            )
            dd_det["flux_app_ratio"][idx_fr60] = (
                dd_det["flux_r38"][idx_fr60] / dd_det["flux_r60"][idx_fr60]
            )
            return dd_det

        # *** detections
        dd_detections_raw = get_det_dict(tt_detections_raw)
        self.add_table(dd_detections_raw, "galex_field:tt_detections")
        logger.debug("Constructed 'tt_detections'.")

        # *** Coadds
        dd_coadd_detections_raw = get_det_dict(tt_coadd_detections_raw)
        self.add_table(dd_coadd_detections_raw, "galex_field:tt_coadd_detections")
        logger.debug("Constructed 'tt_coadd_detections'.")

        # *** Intensity maps
        aa_sel_int_map = np.char.endswith(
            tt_down["Local Path"].data.astype(str),
            "nd-int.fits.gz" if obs_filter == "NUV" else "fd-int.fits.gz",
        )
        # Selects the reference/coadd path
        path_int_map_ref = (
            self.data_path
            + tt_down[np.logical_and(aa_sel_int_map, tt_down["ID"] == obs_id)][
                "Local Path"
            ][0]
        )  # string

        self.load_sky_map(path_int_map_ref)

        if not ref_maps_only:
            # Selects everything else but the coadd path
            path_int_map_visits = [
                self.data_path + path
                for path in tt_down[
                    np.logical_and(aa_sel_int_map, tt_down["ID"] != obs_id)
                ]["Local Path"].data.astype(str)
            ]  # list
            for vis_img_name in path_int_map_visits:
                self.load_sky_map(vis_img_name, "vis_img")


class GALEXDSField(BaseField):
    """
    VASCA implementation of GALEX drift-scan mode observations

    The archival data is expected to be downloaded separately.
    """

    def __init__(self, field_name, data_path=None, visits_data_path=None, **kwargs):
        """
        Initializes a new GALEXDSField instance with
        skeleton VASCA data structure.

        Parameters
        ----------
        field_name : str
            GALEX drift-scan field name
        data_path : str, optional
            Path to root location of data associated with the GALEX field.
            Defaults to a path given by the resource manager.
        visits_data_path : str, optional
            Path to a pre-downloaded table holding the complete list of GALEX DS visits.
            Defaults to a path given by the resource manager.

        Attributes
        ----------
        data_path : str
            Path to root location of data associated with the GALEX field
        visits_data_path : str
            Path to a pre-downloaded table holding the complete list of GALEX DS visits
        vasca_file_prefix : str
            File name prefix following VASCA naming convention:
            'VASCA_<observatory>_<filed_id>_<obs_filter>'
        """

        # Sets skeleton
        super().__init__()

        # Parse field name (e.g. "29208-KEPLER_SCAN_009_sv03")
        scan_name = field_name.split("_sv")[0]
        field_id = name2id(field_name, bits=64)

        # Set root location of data associated with the field
        if data_path is None or visits_data_path is None:
            with ResourceManager() as rm:
                # Path to read and write data relevant to the pipeline
                if data_path is None:
                    self.data_path = (
                        f'{rm.get_path("gal_ds_fields", "lustre")}/'
                        f"{scan_name}/{field_name}"
                    )
                else:
                    self.data_path = data_path
                # Path to a pre-downloaded table holding
                # the complete list of GALEX visits
                if visits_data_path is None:
                    self.visits_data_path = rm.get_path("gal_ds_visits_list", "lustre")
                else:
                    self.visits_data_path = visits_data_path

        # Checks if field name exists
        if not field_name in Table.read(self.visits_data_path)["field_name"]:
            raise ValueError(
                f"Unkown field name '{field_name}'. "
                f"No matching entry in visits table ('{self.visits_data_path}') "
            )

        # Create and check existence of directory
        # that holds field data and VASCA outputs
        os.makedirs(self.data_path, exist_ok=True)

        # Checks if visits table is available
        if not os.path.isfile(self.visits_data_path):
            raise FileNotFoundError(
                "GALEX drift scan visits table not found at "
                f"'{self.visits_data_path}'."
            )

        # Create name-ID map for field and visits
        logger.debug(f"Creating name-ID map for field '{field_name}'.")
        tt_name_id = Table.read(self.visits_data_path)
        tt_name_id = tt_name_id[tt_name_id["field_name"] == field_name][
            "field_name", "field_id", "vis_name", "vis_id"
        ]
        self.tt_name_id = tt_name_id

        # File name prefix for VASCA/GALEXField outputs
        self.vasca_file_prefix = f"VASCA_GALEX-DS_{field_id}_NUV"

        logger.debug(f"Field data path set to: '{self.data_path}'")
        logger.debug(f"Visits data path set to: '{self.visits_data_path}'")

    def name_id_map(self, a):
        """
        Return the field/visit ID if the name is supplied and vice versa.

        Note its not necessary to specify if identifier corresponds to the
        field or visit since names/IDs are unique.

        Parameters
        ----------
        a : str, int
            Input name (if string) or ID (if integer)

        Returns
        -------
        str, int
            Output name or ID
        """

        # Gets field ID and name
        field_id = np.unique(self.tt_name_id["field_id"])[0]
        field_name = np.unique(self.tt_name_id["field_name"])[0]

        # Checks if input is name or ID
        is_name = False
        is_id = False
        if isinstance(a, str):
            is_name = True
        elif isinstance(a, int):
            is_id = True

        out = None
        # Name -> return ID
        if is_name:
            # Field
            if a == field_name:
                out = field_id
            # Visit
            else:
                try:
                    out = self.tt_name_id[self.tt_name_id["vis_name"] == a]["vis_id"][0]
                except IndexError as e:
                    logger.exception(
                        f"No match for name '{a}'. Returning None. Exception:\n"
                        f"{type(e).__name__}: {e}"
                    )
        # ID -> return name
        elif is_id:
            # Field
            if a == field_id:
                out = field_name
            # Visit
            else:
                try:
                    out = self.tt_name_id[self.tt_name_id["vis_id"] == a]["vis_name"][0]
                except IndexError as e:
                    logger.exception(
                        f"No match for ID '{a}'. Returning None. Exception:\n"
                        f"{type(e).__name__}: {e}"
                    )
        # Wrong input type
        else:
            logger.warning(
                f"Expected string or integer type for input {a},"
                f" got ({type(a).__name__}'). Returning None."
            )

        return out

    def _load_info(self, field_name):
        """
        Load field and associated visits information.

        Parameters
        ----------
        field_name : str
            GALEX field name.
        """

        # Central visits table
        tt_visits_full = Table.read(self.visits_data_path)

        # Constructs field info table
        logger.debug("Constructing 'tt_fields'.")

        # Create field information table
        # Map table keys to VASCA-table column names
        col_names_map = {
            "field_id": "field_id",
            "field_name": "field_name",
            "RA_CENT": "ra",
            "DEC_CENT": "dec",
            "observatory": "observatory",
            "obs_filter": "obs_filter",
        }
        col_names = list(col_names_map.keys())

        # Selects columns and filters for field ID
        tt_visits_select = tt_visits_full[tt_visits_full["field_name"] == field_name]
        tt_visits_select = tt_visits_select[col_names]
        tt_fields = unique(tt_visits_select)

        # Checks if resulting table has only one row
        if not len(tt_fields) == 1:
            raise ValueError(
                "Expected single-row table for tt_fields, "
                f"got {len(tt_fields)}, indicating inconsistent gloabl visits info."
            )

        # Converts field_id into VASCA field_id
        tt_fields.replace_column(
            "field_id",
            [
                get_field_id(
                    obs_field_id=tt_fields["field_id"][0],
                    observatory="GALEX_DS",
                    obs_filter="NUV",
                )
            ],
        )

        # Convert into dictionary with correct VASCA column names
        dd_fields = {}
        for col in col_names:
            dd_fields[col_names_map[col]] = tt_fields[col].data

        # Add FoV size info
        dd_fields["fov_diam"] = [3.0]  # Not defined for drift scan fields

        # Add table as class attribute
        self.add_table(dd_fields, "base_field:tt_fields")

        # Constructs field info table
        logger.debug("Constructing 'tt_visits'.")

        # Map table keys to VASCA-table column names
        col_names_map = {
            "vis_id": "vis_id",
            "time_bin_start": "time_bin_start",
            "NEXPTIME": "time_bin_size",
            "RA_CENT": "ra",
            "DEC_CENT": "dec",
        }

        col_names = list(col_names_map.keys())

        # Selects columns and filters for field ID
        tt_visits_select = tt_visits_full[tt_visits_full["field_name"] == field_name]
        tt_visits_select = tt_visits_select[col_names]

        # Convert into dictionary with correct VASCA column names
        dd_visits_select = {}
        for col in col_names:
            dd_visits_select[col_names_map[col]] = tt_visits_select[col].data

        # Add filter info
        n_visits = len(tt_visits_select)
        dd_visits_select["obs_filter_id"] = [dd_filter2id["NUV"]] * n_visits

        # Sets table as class attribute
        self.add_table(dd_visits_select, "galex_field:tt_visits")

        # Sets convenience class attributes
        self.set_field_attr()

        # Cross-check IDs
        # Field ID
        if not int(self.field_id[3:]) == name2id(field_name, bits=64):
            raise ValueError(
                "Inconsistent field ID. "
                f"Expected {name2id(field_name, bits=64)} for name '{field_name}', "
                f"got {self.field_id[3:]}."
            )
        # Visit IDs
        if not all(
            [
                vis_id0 == vis_id1
                for vis_id0, vis_id1 in zip(
                    self.tt_visits["vis_id"], tt_visits_select["vis_id"]
                )
            ]
        ):
            raise ValueError(
                "Inconsistent visit IDs. Missmatch between GALEXDSField.tt_visit "
                "and global DS visits table."
            )

    def _load_data_products(
        self,
        field_name,
        ref_maps_only=False,
    ):
        """
        Loads the relevant data products form storage.

        Parameters
        ----------
        field_name : str
            GALEX field name.
        ref_maps_only : bool, optional
            If True, only the reference/coadd maps are loaded (default).
            If False, also the visit maps are included.
        """

        # Specifies required columns for VASCA

        # Maps table columns to VASCA-table column names
        col_names_map = {
            "vis_id": "vis_id",
            "ggoid_dec": "det_id",
            "alpha_j2000": "ra",
            "delta_j2000": "dec",
            "nuv_poserr": "pos_err",
            "nuv_flux": "flux_auto",
            "nuv_fluxerr": "flux_auto_err",
            "nuv_s2n": "s2n",
            "fov_radius": "r_fov",
            "nuv_artifact": "artifacts",
            "NUV_CLASS_STAR": "class_star",
            "chkobj_type": "chkobj_type",
            "NUV_ELLIPTICITY": "ellip_world",
            "NUV_FLUX_APER_4": "flux_r60",
            "NUV_FLUXERR_APER_4": "flux_r60_err",
            "NUV_FLUX_APER_3": "flux_r38",
            "NUV_FLUXERR_APER_3": "flux_r38_err",
            "NUV_FLUX_AUTO": "flux_auto_flt",
            "E_bv": "E_bv",
            "nuv_mag": "mag",
            "nuv_magerr": "mag_err",
        }
        col_names = list(col_names_map.keys())

        # Visit detections

        # Opens and stacks all visit-detection catalogs, adds column for visit IDs,
        logger.debug(f"Loading visit-detection catalogs for field '{self.field_name}'.")

        # Gets paths to all mcat (visit-detection catalogs) files
        # via the global visits table. This ensures to exclude bad quality visits

        # Central visits info table
        tt_visits_full = Table.read(self.visits_data_path)
        tt_visits_select = tt_visits_full[
            tt_visits_full["field_name"] == self.field_name
        ]
        mcat_vis_paths = [
            glob(f"{self.data_path}{os.sep}{vis_name}/*xd-mcat.fits")[0]
            for vis_name in tt_visits_select["vis_name"]
        ]
        # Loops over visits to read mcat files
        for i, (mcat_path, vis_id) in enumerate(
            zip(mcat_vis_paths, self.tt_visits["vis_id"])
        ):
            # Reads visit detections catalog
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", AstropyWarning)
                tt_mcat = Table.read(mcat_path)

            # Initializes table in first iteration
            # Column selection is applied
            if i == 0:
                tt_mcat.add_column(
                    np.full(len(tt_mcat), vis_id), name="vis_id", index=0
                )
                tt_detections_raw = tt_mcat
            # Appends all catalogs that follow
            else:
                tt_mcat.add_column(
                    np.full(len(tt_mcat), vis_id), name="vis_id", index=0
                )
                tt_detections_raw = vstack([tt_detections_raw, tt_mcat])
            # Clean-up
            del tt_mcat

        # Modify and augment visit-detections table
        # and convert into dictionary with correct VASCA column names
        dd_det = {}

        # Keep only entries with detections, check that by looking at s2n column
        sel_s2n = tt_detections_raw["nuv_s2n"] > 0

        # Warn in cases where there is a filed entry in MAST without detections
        if len(tt_detections_raw[sel_s2n]) == 0:
            logger.warning(
                "No detection with s2n>0 in MAST for "
                f"field {self.tt_visits['vis_name']}."
            )

        # Collect all columns that don't need modifications
        for mast_col in col_names:
            if mast_col not in tt_detections_raw.colnames:
                continue
            else:
                vasca_col = col_names_map[mast_col]
                dd_det[vasca_col] = tt_detections_raw[mast_col][sel_s2n].data

        # Correct flux flux outside of aperture, as listed here:
        # http://www.galex.caltech.edu/researcher/techdoc-ch5.html
        acorr60 = 1.116863247  # 6.0 arcsec correction factor for NUV filter

        # Add flux aperture ratio and covert counts to Jansky
        # for r=3.8 arcsec aperture flux
        # Check to avoid divisions by zero
        idx_fauto = np.where(dd_det["flux_auto_flt"] > 0)
        cts2Jy = np.zeros(len(dd_det["flux_auto_flt"]))
        cts2Jy[idx_fauto] = (
            dd_det["flux_auto"][idx_fauto] / dd_det["flux_auto_flt"][idx_fauto]
        )

        # Use aperture flux instead of AUTO flux,
        # as it is more reliable for variability studies
        dd_det["flux"] = dd_det["flux_r60"] * cts2Jy * acorr60
        dd_det["flux_err"] = dd_det["flux_r60_err"] * cts2Jy * acorr60

        # Add filter_id
        dd_det["obs_filter_id"] = np.array(
            dd_filter2id["NUV"] + np.zeros(np.sum(sel_s2n))
        )

        # Ratio between the flux-counts, correcting for difference in aperture
        idx_fr60 = np.where(dd_det["flux_r60"] > 0)
        dd_det["flux_app_ratio"] = (
            np.zeros(
                len(dd_det["flux_r60"]),
                dtype=dd_vasca_columns["flux_app_ratio"]["dtype"],
            )
            + dd_vasca_columns["flux_app_ratio"]["default"]
        )
        dd_det["flux_app_ratio"][idx_fr60] = (
            dd_det["flux_r38"][idx_fr60] / dd_det["flux_r60"][idx_fr60]
        )

        # Compute geometric sizes of flux profile from average quantity in arcsec
        dd_det["size_world"] = (
            3600
            * (
                tt_detections_raw["NUV_A_WORLD"][sel_s2n].data
                + tt_detections_raw["NUV_B_WORLD"][sel_s2n].data
            )
            / 2
        )

        # Create VASCA visit-detections table
        self.add_table(dd_det, "galex_field:tt_detections")
        del tt_detections_raw

        # Co-add/reference image
        logger.debug(
            "Creating reference/coadd sky-map for field " f"'{self.field_name}'"
        )
        # Gets paths to intensity map FITS file
        coadd_int_file_name = glob(f"{self.data_path}{os.sep}*-nd-int-coadd.fits.gz")[0]

        # Loads the reference intensity map
        self.load_sky_map(coadd_int_file_name)

        # Loads visit intensity maps if requested
        logger.debug("Loading all visit-level sky maps.")
        if not ref_maps_only:
            visit_int_file_names = [
                glob(f"{self.data_path}{os.sep}{vis_name}/*nd-int.fits.gz")[0]
                for vis_name in tt_visits_select["vis_name"]
            ]
            for path in visit_int_file_names:
                self.load_sky_map(path, "vis_img")

    @classmethod
    def from_VASCA(cls, field_name, fits_path=None, **kwargs):
        """
        Constructor to initialize a GALEXField instance
        from a VASCA-generated FITS file

        Parameters
        ----------
        field_name : int
            GALEX field ID
        fits_path : str, optional
            Path to the fits file. Defaults to a path handled by ResourceManager.
        **kwargs
            All additional keyword arguments are passed to `~GALEXField.__init__()`

        Returns
        -------
        vasca.field.GALEXField
        """
        # Bootstrap the initialization procedure using the base class
        fd = cls(field_name, **kwargs)  # new GALEXField instance

        if fits_path is None:
            # Construct the file name from field ID and filter
            fits_path = f"{fd.data_path}/{fd.vasca_file_prefix}_field_data.fits"
        # Check if file exists
        if not os.path.isfile(fits_path):
            raise FileNotFoundError(
                "Wrong file or file path to VASCA data "
                f"for GALEX drift scan field '{field_name}'."
            )
        else:
            # Reads the VASCA-generated field data
            fd.load_from_fits(fits_path)
            # Sets convenience class attributes
            fd.set_field_attr()

        return fd

    @classmethod
    def from_MAST(
        cls,
        field_name,
        load_products="TABLES",
        write=True,
        **kwargs,
    ):
        """
        Constructor to initialize a GALEXDSField instance from pre-downloaded MAST
        archival data.

        The procedure uses vasca.resource_manager.ResourceManager
        to handle file locations.

        Parameters
        ----------
        field_name : str
            GALEX drift scan field name
        load_products : str, optional
            Selects if data products shall be loaded: "NONE", "TABLES" for tables,
            and reference image, "ALL" for tables and visit images.
        write: bool, optional
            If load_products is enabled, stores the data as VASCA tables in the cloud
            for faster loading in the future. Default is True.
        **kwargs
            All additional keyword arguments are passed to `~GALEXField.__init__()`

        Returns
        -------
        vasca.field.GALEXDSField
            Returns None if no entry is found
        """

        # Initialize field instance
        fd = cls(field_name, **kwargs)

        # Gets basic information of field and associated visits
        try:
            fd._load_info(field_name)
        except Exception as e:
            logger.exception(
                f"Could not load GALEX drift scan field: {field_name}. Exception:\n"
                f"{type(e).__name__}: {e}"
            )
            return None

        # Gets detections list and intensity maps
        if load_products != "NONE":
            if load_products == "ALL":
                fd._load_data_products(field_name, ref_maps_only=False)
            elif load_products == "TABLES":
                fd._load_data_products(field_name, ref_maps_only=True)
            else:
                logger.warning(f"load_product option not valid: {load_products}")

            if write:
                fits_name = f"{fd.data_path}/{fd.vasca_file_prefix}_field_data.fits"
                fd.write_to_fits(fits_name)

        meta_only = " (metadata only)" if load_products == "NONE" else ""
        logger.info(
            f"Finished loading new GALEX drift scan field '{field_name}'{meta_only} "
            f"from storage ({fd.data_path})."
        )

        return fd

    @staticmethod
    def load(
        field_name,
        method="MAST_LOCAL",
        load_products="TABLES",
        **field_kwargs,
    ):
        """
        Loads GALEX field data according to a given method and
        returns a GALEXField instance.

        Parameters
        ----------
        field_name : int
            GALEX drift scan field name
        method : str, optional
            Specification of the load method. Four methods are implemented:
            MAST_REMOTE: Standard method to query data from MAST archives. This requires
            an internet connection to the MAST servers. Local data will be overwritten
            MAST_LOCAL: Dafault. Builds a new GALEXfield instance based on MAST
            archival data cached on local disc. If no data is found,
            the fallback is MAST_REMOTE.
            VASCA: Builds a new GALEXField based on a VASCA-generated field data file.
            AUTO: Attempts to load field data by using VASCA as method and
            falls back to MAST_LOCAL if no VASCA file is found on disc.

            The default directory where field data availability is checked is
            defined by the "data_path" attribute of GALEXField and can be passed
            explicitly via the "field_kwargs".
        load_products : str, optional
            Selects if data products shall be loaded: "NONE", "TABLES" for tables,
            and reference image, "ALL" for tables and visit images.
        field_kwargs
            All additional keyword arguments are passed to the load methods
            `~GALEXField.from_MAST()` and `~GALEXField.from_VASCA()`,
            as well as to `~GALEXField.__init__()`.

        Raises
        ------
        TypeError
            If the specified load method is not a string.
        ValueError
            If the specified load method is not one of
            '["mast_remote", "mast_local", "vasca", "auto"]'. String matching is
            case insensitive.

        Returns
        -------
        vasca.field.GALEXField

        """
        logger.info(
            f"Loading field '{field_name}' with method '{method}' and "
            f"load_products '{load_products}'."
        )

        # Checks method argument
        # String matching is case insensitive (converts to all to lower case).
        if not isinstance(method, str):
            raise TypeError(
                "Expected string type for load method specification, "
                f"got '{type(method).__name__}'."
            )

        method = method.casefold()  # ensures all-lower-case string
        method_spec = ["mast_remote", "mast_local", "vasca", "auto"]
        if method not in method_spec:
            raise ValueError(
                f"Expected load method specification from {method_spec}, "
                f"got '{method}'."
            )

        # Removes parameters undefined by GALEXDSField
        for param in ["obs_filter", "refresh"]:
            field_kwargs.pop(param, None)

        # Loads field according to load method specification
        if method in ["mast_remote", "mast_local"]:
            gf = GALEXDSField.from_MAST(
                field_name=field_name,
                load_products=load_products,
                **field_kwargs,
            )
        elif method == "vasca":
            # removes unused options
            field_kwargs.pop("write", None)

            try:
                gf = GALEXDSField.from_VASCA(
                    field_name=field_name,
                    **field_kwargs,
                )
            except FileNotFoundError as e:
                logger.warning(
                    f"Load method '{method}' failed due to FileNotFoundError ('{e}'). "
                    "Falling back to method 'mast_remote'."
                )
                gf = GALEXField.from_MAST(
                    field_name=field_name,
                    load_products=load_products,
                    **field_kwargs,
                )
        elif method == "MAST-REMOTE":
            raise NotImplementedError(
                f"Method '{method}' not available for GALEX drift scan data."
            )
        elif method == "auto":
            raise NotImplementedError(f"Method '{method}' not operational.")

        return gf
