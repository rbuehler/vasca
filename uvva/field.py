import collections
import inspect
import itertools
import os
import sys
from collections import OrderedDict
from copy import copy, deepcopy
from datetime import datetime
from pprint import pprint

import healpy as hpy
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import yaml
from astropy import units as uu
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits.hdu.base import _BaseHDU
from astropy.table import Column, QTable, Table, conf, hstack, unique, vstack
from astropy.time import Time
from astropy.wcs import wcs
from astroquery.mast import Observations
from loguru import logger
from matplotlib import cm, colorbar, colors
from matplotlib.colors import LogNorm
from sklearn import cluster
from sklearn.cluster import MeanShift

from .resource_manager import ResourceManager
from .utils import get_time_delta, get_time_delta_mean, sky_sep2d
from .uvva_table import UVVATable, dd_uvva_tables

import matplotlib.pyplot as plt
from itertools import cycle

# global paths
FILE_DIR = os.path.dirname(os.path.abspath(__file__))  # path to the dir. of this file
ROOT_DIR = FILE_DIR + "/../"  # path to the root directory of the repository

dimless = uu.dimensionless_unscaled
conf.replace_warnings = ["always"]


class BaseField(object):
    """
    :class: `~uvva.field.BaseField` provides class that defines the basic data structure
    for field-based analysis. One *field* is generally the area in the sky covered
    by a telescope in one observation. A field is generally composed of several
    *visits* of the telescope at different times.

    This class contains the main functionality for source
    detection and drawing. To be inherited by field analysis classes,
    which can then be tailored to the needs of the observatories supported
    by the UVVA pipeline.
    """

    def __init__(self):
        """

        Notes
        -----
        Many class attributes are stored in astropy.table.Tables_. To see a
        description of each of their columns run :meth: `~uvva.field.BaseField.info`.

        .. _astropy.table.Tables: https://docs.astropy.org/en/stable/api/astropy.table.Table.html

        Returns
        -------
        None.

        """
        # Configure logger
        # adds the class name as an extra key; accessible vie the handler format
        logger.configure(extra={"classname": self.__class__.__name__})

        # Convenience class attributes

        #: Internal list of tables holding the central field data
        #: All of these are defined in ``uvva_table.py``
        self._table_names = list()

        #: Internal dictionary of important parameters
        #: with corresponding tables
        #: to be set as class attributes for convenience.
        self._dd_attr_names = {
            "tt_field": [
                "field_id",
                "name",
                "ra",
                "dec",
                "observatory",
                "obsfilter",
                "center",
            ],
            "tt_visits": [
                "n_visits",
                "t_exp_sum",
                "t_start",
                "t_stop",
            ],
        }

        #: List of 2D sky images of each visit,
        #: ordered as :attr: `~uvva.field.BaseField.tt_visits` table
        self.vis_imgs = None

        #: Astropy WCS object for one 2D images.
        #: The WCS has to be the same for all visit.
        self.vis_iwcs = None

    def plot_sky_sources(
        self, ax=None, plot_detections=True, src_kwargs=None, det_kwargs=None
    ):
        """
        Plot the sources and (optinally) the visit detections on the sky.

        Parameters
        ----------
        ax : axes, optional
            Matplotlib axes to plot on. The default is None.
        plot_detections : bool, optional
            Plot the visit detections below the sources. The default is True.
        src_kwargs : dict, optional
            Keyword arguments for pyplot.plot_ of the sources. The default is None.
        det_kwargs : dict, optional
            Keyword arguments for pyplot.plot_ of the detections. The default is None.

        Returns
        -------
        ax : axes
            DESCRIPTION.

        .. _pyplot.plot: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html

        """

        if ax is None:
            ax = plt.gca()

        # Set marker properties for sources
        plt_src_kwargs = {
            "marker": "o",
            "markerfacecolor": "None",
            "markersize": 12,
        }
        if src_kwargs is not None:
            plt_src_kwargs.update(src_kwargs)

        # Set marker properties for detections
        plt_det_kwargs = {
            "marker": "o",
            "markerfacecolor": "None",
            "markersize": 6,
            "alpha": 0.5,
            "lw": 0,
        }
        if det_kwargs is not None:
            plt_det_kwargs.update(det_kwargs)

        # Prepare data
        det_src_ids = self.tt_detections["src_id"].data
        src_ids = self.tt_sources["src_id"].data
        det_poss = np.array(
            list(zip(self.tt_detections["ra"].data, self.tt_detections["dec"].data))
        )
        src_poss = np.dstack((self.tt_sources["ra"].data, self.tt_sources["dec"].data))[
            0
        ]
        nr_srcs = len(src_ids)

        # Loop over all srcs and plot
        colors = cycle("bgrcmykbgrcmykbgrcmykbgrcmyk")
        for kk, col in zip(range(nr_srcs), colors):
            sel_dets = det_src_ids == kk
            src_pos = src_poss[kk]
            if plot_detections:
                ax.plot(
                    det_poss[sel_dets, 0],
                    det_poss[sel_dets, 1],
                    markeredgecolor=col,
                    **plt_det_kwargs,
                )
            ax.plot(src_pos[0], src_pos[1], markeredgecolor=col, **plt_src_kwargs)
        return ax

    # TODO: def plot_sky_map():

    def plot_sky(self, plot_detections=True, plot_map=False):
        """
        Plot all field sources and/or a background reference image in the sky.

        Parameters
        ----------
        plot_detections : bool, optional
            Plot sources. The default is True.
        plot_map : bool, optional
            Plot reference image in the background. The default is False.

        Returns
        -------
        fig : figure
            Matplotlib figure used to plot.

        """

        fig = plt.figure(1)
        plt.clf()
        if plot_detections:
            self.plot_sky_sources(plot_detections=plot_detections)
        plt.show()
        return fig

    def add_table(self, data, template_name):
        """
        Add a UVVA table to the field.

        Parameters
        ----------
        data : list, array-like
            Data of the table with shape (n, n_cols) or as dictionary with the
            key corresponding to the templates columns.
        template_name : str
            Identifier to select a table template. Templates are selected by
            setting the class key and a corresponding table key in one string
            separated by a colon, e.g. template_name=<class_key>:<table_key>.

        """
        logger.info(f"Adding table '{template_name}'")

        table_key = template_name.split(":")[1]

        if table_key in self._table_names:
            logger.warning(f"Table '{table_key}' already exists, overwriting")
        else:
            self._table_names.append(table_key)

        setattr(self, table_key, UVVATable.from_template(data, template_name))

    def cluster_meanshift(self, bandwidth=1.0, cluster_all=True):
        """
        Apply _MeanShift clustering algorithm using to derive sources.

        .. _MeanShift: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MeanShift.html

        Parameters
        ----------
        bandwidth : float, optional
            Bandwidth used in the RBF kernel. The default is 1.0.
        cluster_all : bool, optional
            If true, then all points are clustered, even those orphans that are not
            within any kernel. Orphans are assigned to the nearest kernel.
            The default is True.

        Returns
        -------
        int
            Number of detected clusters.
        """
        logger.info("Clustering sources")

        # Get detection coordinates and run clustering
        coords = np.array(
            list(zip(self.tt_detections["ra"].data, self.tt_detections["dec"].data))
        )
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=cluster_all)
        ms.fit(coords)

        # Fill in data into field tables
        self.tt_detections["src_id"] = ms.labels_
        src_ids, det_cts = np.unique(ms.labels_, return_counts=True)

        cluster_centers = ms.cluster_centers_
        nr_srcs = len(cluster_centers)
        srcs_data = {
            "src_id": src_ids,
            "ra": cluster_centers[:, 0],
            "dec": cluster_centers[:, 1],
            "nr_det": det_cts,
            "flag": np.zeros(nr_srcs),
        }

        self.add_table(srcs_data, "base_field:tt_sources")

        return nr_srcs

    def set_field_attr(self, dd_names=None):
        """
        Sets the most important field parameters as class attributes.

        Parameters
        ----------
        dd_names : dict, None
            List of strings containing the parameter names. If None (default),
            a set of pre-defined parameters is set.
        """

        if dd_names is None:
            dd_names = self._dd_attr_names
        for table in dd_names:
            for par in dd_names[table]:
                setattr(self, par, self.get_field_par(par, table))

        logger.info(f"Set attributes: {self._dd_attr_names}")

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
            If :attr: `~uvva.field.BaseField.tt_field` has not exactly one row

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
        # tt_field must be single-rowed
        if table_name == "tt_field" and len(self.tt_field) != 1:
            raise ValueError(
                "Inconsistent data. Expeted single-rowed "
                f"field metadata table 'tt_field', got {len(self.tt_field)}."
            )

        # Indirectly derived parameters
        if par_name == "center":
            par = SkyCoord(
                self.get_field_par("ra", "tt_field"),
                self.get_field_par("dec", "tt_field"),
                frame="icrs",
            )
        elif par_name == "n_visits":
            par = len(self.tt_visits)
        elif par_name == "t_exp_sum":
            par = self.tt_visits["t_exp"].sum() * uu.s
        elif par_name == "t_start":
            par = Time(self.tt_visits["t_start"][0], format="mjd")
        elif par_name == "t_stop":
            par = Time(self.tt_visits["t_stop"][-1], format="mjd")
        # Directly derived parameters
        # Return None if parameter is not in table
        elif par_name not in self.__dict__[table_name].colnames:
            logger.warning(
                f"Unknown parameter in data set '{table_name}'. "
                f"Known parameters are {self.__dict__[table_name].colnames}. "
                f"Returning None for parameter '{par_name}'."
            )
            par = None
        else:
            # retrieve parameter
            par = self.tt_field[par_name].data[0]
            # Maintains unit
            # Some combinations of unit and dtype need special handling
            unit = self.tt_field[par_name].unit
            dtype_kind = self.tt_field[par_name].dtype.kind
            # check (byte-)string or unicode
            if dtype_kind in ["U", "S"]:
                par = par.decode("utf-8")
            # check uu.dimensionless_unscaled and integer
            elif unit == "" and dtype_kind in ["i", "u"]:
                par = int(par)
            else:
                par = par * unit

        return par

    def write_to_fits(self, file_name="field_default.fits", overwrite=True):
        """

        Parameters
        ----------
        file_name : str, optional
            File name. The default is "field_default.fits".
        overwrite : bool, optional
            Overwrite existing file. The default is True.

        Returns
        -------
        None.

        """
        logger.info(f"Writing file with name '{file_name}'")

        # Create HDU list and write
        hdup = fits.PrimaryHDU()

        # Check if image data is set and add to primary HDU
        if self.vis_imgs is not None and self.vis_wcs is not None:
            logger.debug("Storing image data'")
            hdup = fits.PrimaryHDU(self.vis_imgs, header=self.vis_wcs.to_header())

        hdus = [hdup]
        new_hdul = fits.HDUList([hdup])
        new_hdul.writeto(file_name, overwrite=overwrite)

        for key in self._table_names:
            if key in self.__dict__:
                logger.debug(f"Writing table '{key}'")
                self.__dict__[key].write(file_name, append=True)

        # Rename extensions to table names
        ext_nr = 0
        ff = fits.open(file_name, "update")
        for key in self._table_names:
            if key in self.__dict__:
                ext_nr += 1
                ff[ext_nr].header["EXTNAME"] = key
        ff.close()

    def load_from_fits(self, file_name="field_default.fits"):
        """

        Parameters
        ----------
        file_name : str, optional
            File name. The default is "field_default.fits".

        Returns
        -------
        None.

        """
        # Load tables
        logger.info(f"Loading file with name '{file_name}'")
        for key in self._table_names:
            if key in self.__dict__:
                logger.debug(f"Loading table '{key}'")
                self.__dict__[key] = Table.read(file_name, hdu=key)

        # Load image data
        ff = fits.open(file_name)
        self.vis_imgs = ff[0].data
        if self.vis_imgs is not None:
            self.vis_wcs = wcs.WCS(ff[0].header)
        else:
            self.vis_wcs = None

        ff.close()

    def info(self):
        """
        Print out information about the field, its visits and sources.

        Returns
        -------
        None.

        """
        for key in self._table_names:
            if key in self.__dict__:
                print(f"\n {key}:")
                self.__dict__[key].info()
                print(self.__dict__[key].meta)

    def __str__(self):
        """
        Return string with information about the field, its visits and sources.

        Returns
        -------
        str

        """
        out_str = ""
        for key in self._table_names:
            if key in self.__dict__:
                out_str += "\n" + self.__dict__[key].__str__()

        return out_str


class GALEXField(BaseField):
    """
    Instance of one GALEX field (coadd + visits)

    Parameters
    ----------

    parobs_id : int

    include_coadd : bool

    Notes
    -----

    Note in the MAST data base "obs_id" is equal
    to the "ParentImgRunID" of the visits table.
    """

    def __init__(self, obs_id, filter=None, refresh=False):

        # check filter name, default is "NUV"
        if filter is None:
            filter = "NUV"
        elif not isinstance(filter, str):
            raise TypeError(
                f"Expected argument of type 'str', got '{type(filter).__name__}'."
            )
        elif filter not in ["NUV", "FUV"]:
            raise ValueError(
                f"Unknown filter '{filter}'. "
                "Available filters for GALEX are 'NUV' or 'FUV'."
            )

        # sets skeleton field
        super().__init__()

        # Collects coadd and visit metadata tables in order
        # to bootstrap the initialization procedure using the base class
        with ResourceManager() as rm:
            #: Path to read and write data relevant to the pipeline
            self.data_path = rm.get_path("gal_fields", "sas_cloud") + "/" + str(obs_id)
        # Sets ``self.tt_field``
        self._load_galex_field_info(obs_id, filter=filter, refresh=refresh)
        # Sets ``self.tt_visits``
        self._load_galex_visits_info(obs_id, filter=filter)
        # Sets ``self.tt_visit_sources``
        self._load_galex_archive_products(obs_id, filter=filter, refresh=refresh)

    def _load_galex_field_info(
        self, obs_id, col_names=None, filter=None, refresh=False
    ):
        """
        Loads the archival metadata associated to a given field ID.
        """
        # Use default columns if not otherwise specified
        # Already sets the order in which columns are added later on
        if col_names is None:
            col_names = [
                "obs_id",
                "target_name",
                "s_ra",
                "s_dec",
                "instrument_name",
                "filters",
            ]

        # path to store raw data
        path_tt_coadd = f"{self.data_path}/MAST_{obs_id}_{filter}_coadd.fits"

        # reads cached file if found found on disc
        if not os.path.isfile(path_tt_coadd) or refresh:
            # download raw table
            logger.debug(
                f"Downloading archive field info and saving to {path_tt_coadd}."
            )
            tt_coadd = Observations.query_criteria(
                obs_id=obs_id,
                dataRights="PUBLIC",
                instrument_name="GALEX",
                filters=filter,
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

        # construct field info table
        # Fill default columns from first row of the archive field info data table
        tt_coadd_select = tt_coadd[col_names]
        self.tt_field = UVVATable.from_template(tt_coadd_select, "base_field:tt_field")
        logger.info("Constructed 'tt_field'.")

    def _load_galex_visits_info(self, obs_id, col_names=None, filter=None):

        # Use default columns if not otherwise specified
        # Already sets the order in which columns are added later on
        if col_names is None:
            col_names = [
                "imgRunID",
                "minPhotoObsDate",
                "maxPhotoObsDate",
                "nexptime" if filter == "NUV" else "fexptime",
                "fexptime" if filter == "NUV" else "nexptime",
                "RATileCenter",
                "DECTileCenter",
            ]
        with ResourceManager() as rm:
            # read cached
            logger.debug(
                "Reading archive visit info from cashed file "
                f"'{rm.get_path('gal_visits_list', 'sas_cloud')}'"
            )
            tt_visits_raw = Table.read(rm.get_path("gal_visits_list", "sas_cloud"))

        # Filters for visits corresponding to field id and selects specified columns
        tt_visits_raw_select = tt_visits_raw[tt_visits_raw["ParentImgRunID"] == obs_id][
            col_names
        ]
        # Convert string time format to mjd
        for key in ["minPhotoObsDate", "maxPhotoObsDate"]:
            tt_visits_raw_select.update(
                {
                    key: Time(
                        [
                            datetime.strptime(
                                date_str.decode("utf-8"), "%m/%d/%Y %I:%M:%S %p"
                            )
                            for date_str in tt_visits_raw_select[key].data
                        ]
                    ).mjd
                }
            )

        # Set as class attribute
        self.tt_visits = UVVATable.from_template(
            tt_visits_raw_select, "galex_field:tt_visits"
        )
        logger.info("Constructed 'tt_visits'.")

    def _load_galex_archive_products(
        self, obs_id, col_names=None, filter=None, product_list=None, refresh=False
    ):
        """
        Loads the relevant data products from MAST servers and
        stores them using the ResourceManager
        """
        # Path to MAST helper files
        path_tt_coadd = f"{self.data_path}/MAST_{obs_id}_{filter}_coadd.fits"
        path_tt_data = f"{self.data_path}/MAST_{obs_id}_{filter}_data.fits"
        path_tt_down = f"{self.data_path}/MAST_{obs_id}_{filter}_down.ecsv"

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
            tt_data = Observations.get_product_list(tt_coadd)
            # Saves to disc
            tt_data.write(path_tt_data, overwrite=True)
        else:
            logger.debug(
                f"Reading archive data products list from cashed file '{path_tt_data}'"
            )
            tt_data = Table.read(path_tt_data)

        # Sets default products to download if a selection was not explicitly passed
        if product_list is None:
            # Default: Catalog and NUV intensity map
            product_list = [
                "xd-mcat.fits.gz",
                "nd-int.fits.gz" if filter == "NUV" else "fd-int.fits.gz",
            ]
        # Verifies if products are listed in products list
        product_list_missing = [
            prod
            for prod in product_list
            if not any(prod in uri for uri in tt_data["dataURI"])
        ]
        if not len(product_list_missing) == 0:
            raise ValueError(
                f"Missing entry for data products {product_list_missing} "
                "in product list."
            )

        # Filters for products of interest
        aa_sel_prod = np.full(len(tt_data), False)  # bool array
        aa_data_uri = tt_data["dataURI"].data.astype(str)  # string array
        for prod in product_list:
            aa_sel_prod += np.char.endswith(aa_data_uri, prod)  # Logical AND operation

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
            # Save to disc
            tt_down["Local Path"] = [
                str(s).split(self.data_path)[1] for s in tt_down["Local Path"]
            ]
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
            tt_down["Local Path"].data.astype(str), product_list[0]
        )
        # Selects only the coadd path
        path_tt_ref = (
            self.data_path
            + tt_down[np.logical_and(aa_sel_mcat, tt_down["ID"] == obs_id)][
                "Local Path"
            ][0]
        )  # string
        # Selects everything else but the coadd path
        path_tt_detections = [
            self.data_path + path
            for path in tt_down[np.logical_and(aa_sel_mcat, tt_down["ID"] != obs_id)][
                "Local Path"
            ].data.astype(str)
        ]  # list

        # Open reference catalog
        tt_ref_sources_raw = Table.read(path_tt_ref)
        # Open visit catalogs, add a columns for visit and source IDs
        # and stack all catalogs
        tt_detections_raw = None
        # loop over visits
        for path, id in zip(path_tt_detections, self.tt_visits["vis_id"]):
            # read new visits catalog
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

        # Add to tables as class attributes
        if col_names is None:
            col_names = [
                "visit_id",
                "source_id",
                "ggoid_dec",
                "alpha_j2000",  # take band-merged quantities?
                "delta_j2000",  # take band-merged quantities?
                "nuv_poserr" if filter == "NUV" else "fuv_poserr",
                "nuv_mag" if filter == "NUV" else "fuv_mag",
                "nuv_magerr" if filter == "NUV" else "fuv_magerr",
                # "nuv_s2n" if filter == "NUV" else "fuv_s2n",
                "fov_radius",
                "nuv_artifact" if filter == "NUV" else "fuv_artifact",
                "NUV_CLASS_STAR" if filter == "NUV" else "FUV_CLASS_STAR",
                "chkobj_type",
            ]
        # set data as class attributes
        self.add_table(tt_detections_raw[col_names], "galex_field:tt_detections")
        logger.info("Constructed 'tt_detections'.")
        self.add_table(tt_ref_sources_raw[col_names[2:]], "galex_field:tt_ref_sources")
        logger.info("Constructed 'tt_ref_sources'.")

        # Todo: Intensity maps


class Field:
    """
    Instance of one GALEX field (coadd + visits)

    Note in the MAST data base "obs_id" is equal
    to the "ParentImgRunID" of the visits table.
    """

    def __init__(self, parobs_id, cuts=None, refresh=False, verbose=False):

        # logging prefix
        log_prefix = str.format(
            "[{}.{}] ",
            self.__class__.__name__,
            inspect.currentframe().f_code.co_name,
        )

        # parse arguments --------------------------------------------------------------
        # <-- verify data types of input arguments
        self.parobs_id = parobs_id  # field of interest
        self.dd_cuts = self._parse_cuts(cuts)  # cuts
        self.flags = dict()  # dictionary holding all settings flags
        self.flags["refresh"] = refresh  # force MAST query
        self.flags["verbose"] = verbose  # verbose printing/logging

        # collect data -------------------------------------------------------------
        # directory to store data
        rm = ResourceManager()  # <-- use context manager in the future
        self.data_path = (
            rm.get_path("gal_fields", "sas_cloud") + "/" + str(self.parobs_id)
        )
        os.makedirs(self.data_path, exist_ok=True)
        # <-- assert data availability and warn about new downloads or refreshes
        # GALEX visits metadata
        self.tt_visits = self._load_visits(self.parobs_id)
        # GALEX coadd metadata
        self.tt_coadd = self._load_coadd(self.parobs_id)  # <-- check consistency
        # dictionary listing coadd and visit IDs
        self.dd_ids = OrderedDict()
        self.dd_ids["coadd"] = self.parobs_id
        for i, id in enumerate(self.tt_visits["imgRunID"]):
            self.dd_ids["visit" + str(i)] = id
        # GALEX archival data products
        # product_list = ["xd-mcat.fits.gz", "nd-int.fits.gz"]
        self.dd_mcats, self.dd_maps = self._load_archive_products(self.parobs_id)

        # source match + variability detection -----------------------------------------
        self.tt_match = self._source_match()

        # field metadata ----------------------------------------------------
        self.dd_meta = self._collect_metadata()

        # GALEX catalog and image data
        # self.dd_mcats, self.dd_imgs = self._load_archive_products(
        #     self.parobs_id
        # )

    def _parse_cuts(self, cuts):
        """
        Sets cuts
        """
        dd_cuts_default = {
            "s2n_coadd": 6,
            "s2n_visit": 2,
            "r_fov_coadd": 0.5,
            "r_fov_visit": 0.503,
            "mag_upper_lim_coadd": 16,
            "separation": 2.5,
            "ds2n": 4,
            "dmag": 0.2,
        }
        if cuts is not None:
            return {**dd_cuts_default, **cuts}
        else:
            return dd_cuts_default

    def _load_visits(self, parobs_id, hpy_nside=2 ** 10, cframe="galactic", raw=False):
        """
        Loads the table containing a list of GALEX visits.
        Reduces the list contents to all visits corresponding the the field ID
        Adds additional columns:
            - galactic coordinates derived from equatorial coordinates
            - healpix map containing field center coordinates
            - Time elapsed between consecutive visits
        If raw is requested, the unmodified table is returned
        """
        if self.flags["verbose"]:
            log_prefix = str.format(
                "[{}.{}] ",
                self.__class__.__name__,
                inspect.currentframe().f_code.co_name,
            )
            msg = "Loading GALEX visits table."
            print(log_prefix, msg)
        # <-- handle refresh cases
        rm = ResourceManager()  # <-- use context manager in the future
        tt_visits = Table.read(rm.get_path("gal_visits_list", "sas_cloud"))

        if raw:  # if only the raw table is requested, loading is completed here
            return tt_visits

        else:  # otherwise more columns are added and filter for field ID
            tt_visits = tt_visits[tt_visits["ParentImgRunID"] == parobs_id]
            # set coordinate system
            if cframe == "galactic":
                xx = "gall"
                yy = "galb"
            elif cframe == "icrs":
                xx = "RATileCenter"
                yy = "DECTileCenter"

            # add galactic coordinates
            cc_gvis = SkyCoord(
                ra=tt_visits["RATileCenter"],
                dec=tt_visits["DECTileCenter"],
                frame="icrs",
                unit=uu.degree,
            )
            tt_visits["gall"] = cc_gvis.galactic.l.degree
            tt_visits["galb"] = cc_gvis.galactic.b.degree

            # add units
            tt_visits["gall"].unit = uu.degree
            tt_visits["galb"].unit = uu.degree
            tt_visits["RATileCenter"].unit = uu.degree
            tt_visits["DECTileCenter"].unit = uu.degree

            # healpy settings
            hpy_npix = hpy.nside2npix(hpy_nside)
            pix_diam = hpy.nside2resol(hpy_nside, arcmin=True) / 60 * uu.deg

            if self.flags["verbose"]:
                msg = str.format(
                    "Generating healpix map. "
                    "Pixel diameter: {:1.3f} deg / {:1.1f} arcsec. "
                    "N_pix: {:d}, N_side: {:d}.",
                    pix_diam,
                    pix_diam.to(uu.arcsec),
                    hpy_npix,
                    hpy_nside,
                )
                print(log_prefix, msg)

            # convert field center coordinates to healpix map
            tt_visits["hpix"] = hpy.ang2pix(
                hpy_nside,
                tt_visits[xx],
                tt_visits[yy],
                lonlat=True if cframe == "galactic" else False,
                nest=False,
            )

            # time between visits
            tt_visits["PhotoObsDate_MJD_delta"] = np.insert(
                get_time_delta(tt_visits["PhotoObsDate_MJD"], uu.s), 0, 0
            )
            tt_visits["PhotoObsDate_MJD_delta"].unit = uu.s

            # store table
            tt_visits_path = self.data_path + "/MAST_" + str(parobs_id) + "_visits.fits"
            if not os.path.isfile(tt_visits_path):
                tt_visits.write(tt_visits_path, overwrite=True)

            return tt_visits

    def _load_coadd(self, parobs_id, filters=["NUV"]):
        """
        Loads the archival metadata associated to a given field/coadd ID.
        """
        if self.flags["verbose"]:
            log_prefix = str.format(
                "[{}.{}] ",
                self.__class__.__name__,
                inspect.currentframe().f_code.co_name,
            )
            msg = "Loading GALEX coadd table."
            print(log_prefix, msg)

        tt_coadd_path = self.data_path + "/MAST_" + str(parobs_id) + "_coadds.fits"
        # reads cached file if found found on disc
        if not os.path.isfile(tt_coadd_path) or self.flags["refresh"]:
            # download
            tt_coadd = Observations.query_criteria(
                obs_id=parobs_id,
                dataRights="PUBLIC",
                instrument_name="GALEX",
                filters=filters,
                dataproduct_type="image",
            )
            # save on disc
            tt_coadd.write(tt_coadd_path, overwrite=True)
        else:
            tt_coadd = Table.read(
                self.data_path + "/MAST_" + str(parobs_id) + "_coadds.fits"
            )
        return tt_coadd

    def _load_archive_products(self, parobs_id, product_list=None):
        """
        Loads archive data products.
        """

        if self.flags["verbose"]:
            log_prefix = str.format(
                "[{}.{}]",
                self.__class__.__name__,
                inspect.currentframe().f_code.co_name,
            )
            msg = "Loading GALEX data products."
            print(log_prefix, msg)

        tt_data_path = self.data_path + "/MAST_" + str(parobs_id) + "_data.fits"
        tt_down_path = self.data_path + "/MAST_" + str(parobs_id) + "_down.ecsv"

        # tt_data: list of all archival data products ----------------------------------
        # reads cached file if found found on disc
        if not os.path.isfile(tt_data_path) or self.flags["refresh"]:
            # download
            tt_data = Observations.get_product_list(self.tt_coadd)
            # save on disc
            tt_data.write(tt_data_path, overwrite=True)
        else:
            tt_data = Table.read(tt_data_path)
        # set default products to download if a selection was not explicitly passed
        if product_list is None:
            # default: catalog and NUV intensity map
            product_list = ["xd-mcat.fits.gz", "nd-int.fits.gz"]
        # verify if products are listed in products list
        product_list_missing = [
            prod
            for prod in product_list
            if not any(prod in uri for uri in tt_data["dataURI"])
        ]
        assert (
            len(product_list_missing) == 0
        ), "{} Missing entry for data products '[{}]' in product list .".format(
            log_prefix,
            ", ".join(["{}"] * len(product_list_missing)).format(*product_list_missing),
        )
        # filter for products of interest
        sel_prod = np.full(len(tt_data), False)  # bool array
        for prod in product_list:
            sel_prod += np.core.defchararray.find(tt_data["dataURI"], prod) != -1

        # tt_down: manifest of files downloaded ----------------------------------------
        if not os.path.isfile(tt_down_path) or self.flags["refresh"]:
            # download
            tt_down = Observations.download_products(
                tt_data,
                download_dir=self.data_path,
                cache=True,
                dataURI=tt_data[sel_prod]["dataURI"],
            )["Local Path", "Status"]
            # save to disc
            tt_down["Local Path"] = [
                str(s).split(self.data_path)[1] for s in tt_down["Local Path"]
            ]
            tt_down.write(tt_down_path, overwrite=True)
        else:
            tt_down = Table.read(tt_down_path)
        # add column with corresponding IDs
        tt_down["ID"] = [int(path.split(os.sep)[-2]) for path in tt_down["Local Path"]]

        # collect data_products --------------------------------------------------------

        # <-- assert tt_down["Local Path"] == "COMPLETE"

        dd_mcats = dict.fromkeys(self.dd_ids.keys())
        dd_maps = dict.fromkeys(self.dd_ids.keys())
        # select products
        sel_mcat = (
            np.core.defchararray.find(tt_down["Local Path"], "xd-mcat.fits.gz") != -1
        )
        sel_ndint = (
            np.core.defchararray.find(tt_down["Local Path"], "nd-int.fits.gz") != -1
        )

        # open <-- use less memory hungry method in the future (context manager!)
        if self.flags["verbose"]:
            msg = "Opening catalog files and NUV intesity maps."
            print(log_prefix, msg)
        for field, id in self.dd_ids.items():
            # catalog tabels
            sel = np.logical_and(sel_mcat, tt_down["ID"] == id)
            dd_mcats[field] = Table.read(self.data_path + tt_down[sel]["Local Path"][0])
            dd_mcats[field]["glon"].unit = uu.degree
            dd_mcats[field]["glat"].unit = uu.degree

            # NUV intensity maps
            sel = np.logical_and(sel_ndint, tt_down["ID"] == id)
            dd_maps[field] = fits.open(self.data_path + tt_down[sel]["Local Path"][0])[
                0
            ]

        return (dd_mcats, dd_maps)

    def _data_availability(self, parobs_id, data_names_list):
        """
        Verifies data availability.
        """

        if data_names_list is None:
            # default: catalog and NUV intensity map
            data_names_list = ["visits", "coadd", "data", "down"]

        for name in data_names_list:
            assert os.path.isfile(name)

    def _collect_metadata(self):

        if self.flags["verbose"]:
            log_prefix = str.format(
                "[{}.{}] ",
                self.__class__.__name__,
                inspect.currentframe().f_code.co_name,
            )

        n_visits = len(self.tt_visits)
        n_surveys = len(unique(self.tt_visits, "survey"))
        survey_names = list(unique(self.tt_visits, ["survey"])["survey"])
        survey_n_visits = [
            np.sum(self.tt_visits["survey"] == sur) for sur in survey_names
        ]
        nexptime_mean = self.tt_visits["nexptime"].mean() * uu.second
        nexptime_std = self.tt_visits["nexptime"].std() * uu.second
        obs_time_delta_mean, obs_time_delta_std = (
            get_time_delta_mean(self.tt_visits["PhotoObsDate_MJD"], uu.d) * uu.day
        )
        pos_delta_mean = (
            np.mean(
                sky_sep2d(
                    SkyCoord(
                        self.tt_visits["gall"],
                        self.tt_visits["galb"],
                        frame="galactic",
                    )
                )
            )
            * uu.degree
        )
        dd_meta = dict()
        meta_vars = dir()[2:]
        meta_vars.remove("self")
        for var in meta_vars:
            dd_meta[var] = eval(var)
        if self.flags["verbose"]:
            print(
                str.format(
                    (
                        "{:s}Report: GALEX field data. ID: {:d}. "
                        "Total number of visits {:d} from {:d} survey(s): {:s}. "
                        "Average NUV exposure time: {:1.1f} ± {:1.1f}. "
                        "Average time between visits: {:1.1f} ± {:1.1f}. "
                        "Averge field center variation: {:1.2f}."
                    ),
                    log_prefix,
                    self.parobs_id,
                    n_visits,
                    n_surveys,
                    ", ".join(
                        [
                            "{}({})".format(a, b)
                            for a, b in zip(survey_names, survey_n_visits)
                        ]
                    ),
                    nexptime_mean,
                    nexptime_std,
                    obs_time_delta_mean,
                    obs_time_delta_std,
                    pos_delta_mean,
                )
            )

        return dd_meta

    # Source matching & variability detection ------------------------------------------

    def _get_mcat_sel(self, mcat_name, s2n_cut, fov_cut, mag_upper_lim_cut=5):
        """
        Masking mcat tables.
        """
        tt = self.dd_mcats[mcat_name]
        sel_mag = (
            tt["nuv_mag"] > mag_upper_lim_cut
        )  # <-- replace this by a proper saturation cut
        sel_s2n = tt["nuv_s2n"] > s2n_cut
        sel_fov = tt["fov_radius"] < fov_cut
        sel_artifact = (
            (tt["nuv_artifact"] != 2)
            * (tt["nuv_artifact"] != 4)
            * (tt["nuv_artifact"] < 255)
        )  # For artifact ID see http://www.galex.caltech.edu/wiki/GCAT_Manual
        sel = sel_mag * sel_s2n * sel_fov * sel_artifact
        return sel

    def get_variability_mask(self, level=3, project_to_coadd=False, info=False):
        """
        Computes bool-array masks for the catalog table
        according to the cut level.

        l1: spatial separation of coadd and visit sources
        l2: signal-to-noise of the magnitude variation * l1
        l3: absolute cut on the magnitude variation * l1 * l2
        """
        cut_sep = self.dd_cuts["separation"] * uu.arcsec  # 2.5 arcsec separation
        cut_ds2n = self.dd_cuts["ds2n"]  # 4 signal-to-noise variation
        cut_dmag = self.dd_cuts["dmag"]  # 0.2 magnitude variation

        masks = [
            self.tt_match["visit_dist"] < cut_sep,
            self.tt_match["visit_nuv_ds2n"] > cut_ds2n,
            np.abs(self.tt_match["visit_nuv_dmag"]) > cut_dmag,
        ]

        # selects visit sources if all cut criteria are satisfied
        var_mask = np.prod(np.array(masks[:level]), axis=0).astype(
            bool
        )  # AND operation

        # generate statistics info about the source mask
        if info:
            pass  # <-- compute mask statistics

        if project_to_coadd:
            # selects coadd sources if at least one visit satisfies cut criterion
            return var_mask.any(axis=1)  # OR operation
        else:
            return var_mask

    def _source_match(self, masked=False):
        """
        Creates catalog of positionally matched sources
        where coadd sources are matched to nearest visit source.

        GALEX field catalog column description:
        http://www.galex.caltech.edu/wiki/Public:Documentation/Appendix_A.1

        Source classification:
        https://sextractor.readthedocs.io/en/latest/Position.html#class-star-def
        """

        if self.flags["verbose"]:
            log_prefix = str.format(
                "[{}.{}]",
                self.__class__.__name__,
                inspect.currentframe().f_code.co_name,
            )
            msg = "Matching visits to coadd sources."
            print(log_prefix, msg)

        # add coadd info ---------------------------------------------------------------
        coadd_keys = [
            "ggoid_dec",
            "glon",
            "glat",
            "nuv_mag",
            "nuv_magerr",
            "nuv_s2n",
            "fov_radius",
            "chkobj_type",
            "nuv_artifact",
            "nuv_poserr",
            "NUV_CLASS_STAR",
        ]
        # coadd source catalog
        tt_mcat_coadd = Table(
            self.dd_mcats["coadd"][
                self._get_mcat_sel(
                    "coadd",
                    self.dd_cuts["s2n_coadd"],
                    self.dd_cuts["r_fov_coadd"],
                    self.dd_cuts["mag_upper_lim_coadd"],
                )
            ][coadd_keys],
            masked=masked,
        )
        pos_coadd = SkyCoord(
            tt_mcat_coadd["glon"], tt_mcat_coadd["glat"], frame="galactic"
        )

        # add visit info ---------------------------------------------------------------
        n_sources_coadd = len(tt_mcat_coadd)  # number of sources satisfying coadd cuts
        visits_list = list(self.dd_ids.keys())[1:]  # names of the visits fields
        n_visits = len(visits_list)  # number of visit fields
        visit_keys = [
            "visit_" + key
            for key in [
                "dist",
                "ggoid_dec",
                "glat",
                "glon",
                "nuv_poserr",
                "nuv_mag",
                "nuv_magerr",
                "nuv_s2n",
                "nuv_artifact",
            ]
        ]
        # dictionary to store multi-dimensional visits info
        dd_visit_data = {
            key: np.zeros((n_sources_coadd, n_visits)) for key in visit_keys
        }
        # loop over visits
        for i, visit in enumerate(visits_list):
            # visit source catalog
            tt_mcat_visit = self.dd_mcats[visit][
                self._get_mcat_sel(
                    visit, self.dd_cuts["s2n_visit"], self.dd_cuts["r_fov_visit"]
                )
            ]

            pos_visit = SkyCoord(
                tt_mcat_visit["glon"], tt_mcat_visit["glat"], frame="galactic"
            )
            # distance to nearest source
            idx, d2d, _ = pos_coadd.match_to_catalog_sky(pos_visit)
            # fill data
            dd_visit_data["visit_dist"][:, i] = d2d.to(uu.arcsec).value
            for key in visit_keys[1:]:
                dd_visit_data[key][:, i] = np.array(
                    tt_mcat_visit[idx][key.lstrip("visit_")]
                )

        # add flux variability info
        dd_visit_data["visit_nuv_dmag"] = (
            dd_visit_data["visit_nuv_mag"] - tt_mcat_coadd["nuv_mag"][:, None]
        )
        dd_visit_data["visit_nuv_ds2n"] = np.abs(
            dd_visit_data["visit_nuv_dmag"] / dd_visit_data["visit_nuv_magerr"]
        )

        # combine ----------------------------------------------------------------------
        tt_match = hstack([tt_mcat_coadd, Table(dd_visit_data, masked=masked)])

        # fix units
        tt_match["visit_dist"].unit = uu.arcsec
        tt_match["visit_nuv_ds2n"].unit = tt_match["nuv_s2n"].unit
        for key in visit_keys[1:]:
            tt_match[key].unit = tt_match[key.lstrip("visit_")].unit

        return tt_match

    # Plotting functions ---------------------------------------------------------------

    def plot_style_cycler(self, name, shuffle=False):
        np.random.seed(42)
        if name == "color":
            # rgb colors
            sns_bright_rgb = [
                (0.00784313725490196, 0.24313725490196078, 1.0),
                (1.0, 0.48627450980392156, 0.0),
                (0.10196078431372549, 0.788235294117647, 0.2196078431372549),
                (0.9098039215686274, 0.0, 0.043137254901960784),
                (0.5450980392156862, 0.16862745098039217, 0.8862745098039215),
                (0.6235294117647059, 0.2823529411764706, 0.0),
                (0.9450980392156862, 0.2980392156862745, 0.7568627450980392),
                (0.6392156862745098, 0.6392156862745098, 0.6392156862745098),
                (1.0, 0.7686274509803922, 0.0),
                (0.0, 0.8431372549019608, 1.0),
            ]

            # rgba colors
            fc_alpha = 0.1
            sns_bright_rgba = [(*c, fc_alpha) for c in sns_bright_rgb]

            if shuffle:
                np.random.suffle(sns_bright_rgb)

            return itertools.cycle(sns_bright_rgb)

        elif name == "marker":
            # list of all matplotlib markers
            markers = [
                ".",
                ",",
                "o",
                "v",
                "^",
                "<",
                ">",
                "1",
                "2",
                "3",
                "4",
                "8",
                "s",
                "p",
                "*",
                "h",
                "H",
                "+",
                "x",
                "D",
                "d",
                "|",
                "_",
                "P",
                "X",
                # 0,
                # 1,
                # 2,
                # 3,
                # 4,
                # 5,
                # 6,
                # 7,
                # 8,
                # 9,
                # 10,
                # 11,
            ]
            if shuffle:
                # shuffle in-place
                np.random.shuffle(markers)
            return itertools.cycle(markers)

        elif name == "line_style":
            # from matplotlib docu:
            # https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
            # line_styles = ["-", "--", ".-", ":"]
            line_syteles = [
                (0, ()),  # "solid"
                (0, (1, 10)),  # "loosely dotted"
                (0, (1, 5)),  # "dotted"
                (0, (1, 1)),  # "densely dotted"
                (0, (5, 10)),  # "loosely dashed"
                (0, (5, 5)),  # "dashed"
                (0, (5, 1)),  # "densely dashed"
                (0, (3, 10, 1, 10)),  # "loosely dashdotted"
                (0, (3, 5, 1, 5)),  # "dashdotted"
                (0, (3, 1, 1, 1)),  # "densely dashdotted"
                (0, (3, 10, 1, 10, 1, 10)),  # "loosely dashdotdotted"
                (0, (3, 5, 1, 5, 1, 5)),  # "dashdotdotted"
                (0, (3, 1, 1, 1, 1, 1)),  # "densely dashdotdotted"
            ]
            if shuffle:
                np.random.suffle(line_syteles)
            return itertools.cycle(line_syteles)

    def show_sky_map(self, field_name=None, id=None, source_markers=True):
        # use coadd by default
        if id is None and field_name is None:
            id = self.parobs_id
            field_name = "coadd"
        # get field identifier key
        elif field_name is None:
            field_name = list(self.dd_ids.keys())[list(self.dd_ids.values()).index(id)]
        # get field id
        elif id is None:
            id = self.dd_ids[field_name]

        # fits image
        coadd_nint_map = self.dd_maps[field_name]
        # world coordinates
        iwcs = wcs.WCS(coadd_nint_map.header)
        # figure
        fig = plt.figure(figsize=(8.5, 7.2))

        fig.suptitle(
            "GALEX field {} - {}".format(
                self.parobs_id,
                field_name
                if "visit" in field_name
                else "{}: {}".format(
                    field_name,
                    ", ".join(
                        [
                            "{}({})".format(a, b)
                            for a, b in zip(
                                self.dd_meta["survey_names"],
                                self.dd_meta["survey_n_visits"],
                            )
                        ]
                    ),
                ),
            )
        )
        ax = plt.subplot(projection=iwcs)

        # plot data
        cmap = copy(cm.get_cmap("turbo"))
        cmap.set_bad([0.08, 0.06, 0.11])
        label_fontsize = 8

        img = ax.imshow(
            coadd_nint_map.data,
            cmap=cmap,
            norm=LogNorm(),
            origin="lower",
            interpolation="none",
        )
        ax.set_xlabel("Ra", fontsize=label_fontsize)
        ax.set_ylabel("Dec", fontsize=label_fontsize)
        ax.xaxis.set_tick_params(labelsize=label_fontsize)
        ax.yaxis.set_tick_params(labelsize=label_fontsize)
        ax.tick_params(axis="x", labelsize=label_fontsize, bottom=True, top=False)
        ax.tick_params(axis="y", labelsize=label_fontsize, left=True, right=False)

        # colorbar
        cbaxes = fig.add_axes([0.875, 0.1, 0.03, 0.75])
        cb = colorbar.Colorbar(
            ax=cbaxes,
            mappable=img,
            orientation="vertical",
        )
        cb.ax.tick_params(labelsize=label_fontsize)
        cb.set_label("NUV Flux [a.u.]", size=label_fontsize)

        if source_markers:
            # source positions
            for i, (name, id) in enumerate(self.dd_ids.items()):
                sel_mcat = self._get_mcat_sel(
                    name,
                    self.dd_cuts["s2n_visit"]
                    if "visit" in name
                    else self.dd_cuts["s2n_coadd"],
                    self.dd_cuts["r_fov_visit"]
                    if "visit" in name
                    else self.dd_cuts["r_fov_coadd"],
                    self.dd_cuts["mag_upper_lim_coadd"],
                )
                pos = self.dd_mcats[name][sel_mcat]["glon", "glat"]
                ax.scatter(
                    pos["glon"],
                    pos["glat"],
                    transform=ax.get_transform("galactic"),
                    s=4,
                    edgecolor="white",
                    facecolor="none",
                    lw=0.5,
                    label=(i // 2) * "_"
                    + "{}".format("visits" if "visit" in name else "coadd"),
                    marker="v" if "visit" in name else "s",
                    alpha=0.5,
                )
            lgnd = ax.legend(loc="upper right")
            for handle in lgnd.legendHandles:
                handle.set_sizes([20.0])

        # grid
        overlay = ax.get_coords_overlay("galactic")
        overlay.grid(color="white", ls="dotted")
        labels = ["Lon", "Lat"]
        ticks_pos = ["t", "r"]
        for i in range(2):
            overlay[i].set_axislabel(labels[i], fontsize=label_fontsize)
            overlay[i].set_ticklabel(color="k", size=label_fontsize)
            overlay[i].set_ticks(color="white", direction="in")
            overlay[i].set_ticks_position(ticks_pos[i])

        plt.subplots_adjust(
            left=-0.05, bottom=0.082, right=0.88, top=0.85, wspace=0.2, hspace=0.20
        )

        fig.canvas.toolbar_visible = True
        fig.canvas.header_visible = True
        fig.canvas.resizable = True

    def show_visits_timeline(self):
        """
        Used matplotlib examples:
            https://matplotlib.org/stable/gallery/lines_bars_and_markers/timeline.html
            https://matplotlib.org/stable/gallery/ticks_and_spines/date_concise_formatter.html
        """
        names = list(self.dd_ids.keys())
        names.remove("coadd")
        names = [str(i) for i in range(1, len(self.dd_ids.keys()))]
        dates = Time(self.tt_visits["PhotoObsDate_MJD"], format="mjd").to_value(
            "datetime"
        )

        nexpt = self.tt_visits["nexptime"].data

        years = np.unique([d.year for d in dates])
        n_years = len(years)

        fig, axs = plt.subplots(
            nrows=n_years, figsize=(8.8, 4 * n_years), constrained_layout=True
        )
        fig.suptitle(
            "GALEX field {}, {} - visits time line".format(
                self.parobs_id,
                ", ".join(
                    [
                        "{}({})".format(a, b)
                        for a, b in zip(
                            self.dd_meta["survey_names"],
                            self.dd_meta["survey_n_visits"],
                        )
                    ]
                ),
            )
        )

        for i, year in enumerate(years):

            if len(years) > 1:
                ax = axs[i]
            else:
                ax = axs

            # select data
            dates_within = list()
            names_within = list()
            nexpt_within = list()
            for j, (name, date, texp) in enumerate(zip(names, dates, nexpt)):
                if date.year == year:
                    dates_within.append(date)
                    names_within.append(name)
                    nexpt_within.append(texp)

            # Choose some nice levels
            # corresponding to NUV exposure time
            alternate = np.tile([1, -1], int(np.ceil(len(dates_within) / 2)))[
                : len(dates_within)
            ]
            levels = nexpt_within / nexpt.max() * alternate

            # Create figure and plot a stem plot with the date
            ax.set(title=year)

            ax.vlines(dates_within, 0, levels, color="tab:red")  # The vertical stems.
            ax.plot(
                dates_within,
                np.zeros_like(dates_within),
                "-o",
                color="k",
                markerfacecolor="w",
            )  # Baseline and markers on it.

            # annotate lines
            for d, l, r, t in zip(dates_within, levels, names_within, nexpt_within):
                ax.annotate(
                    r,
                    xy=(d, l),
                    xytext=(-0, np.sign(l) * 3),
                    textcoords="offset points",
                    horizontalalignment="center",
                    verticalalignment="bottom" if l > 0 else "top",
                )
                ax.annotate(
                    "{:1.0f}".format(t),
                    fontsize=6,
                    xy=(d, l),
                    xytext=(+0, np.sign(l) * 15),
                    textcoords="offset points",
                    horizontalalignment="center",
                    verticalalignment="bottom" if l > 0 else "top",
                )

            # format xaxis
            locator = mdates.AutoDateLocator(minticks=3, maxticks=15)
            formatter = mdates.ConciseDateFormatter(locator)
            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_major_formatter(formatter)

            # remove y axis and spines
            ax.yaxis.set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            # increase y axis limits
            yspace = ax.get_ylim()
            ax.set_ylim((yspace[0] - 0.3, yspace[1] + 0.3))

            ax.margins(y=0.1)

    def show_diff_map(self, field_name=None, id=None):
        pass

    def show_visit_statistics(
        self,
        key,
        coadd=True,
        bin_size=None,
        density=False,
        xscale=None,
        yscale=None,
        xlim_offset=None,
        figsize=(9, 3.5),
        label_fontsize=10,
    ):
        """
        Show histograms.

        Valid keys are:
            "s2n": Signal-to-Noise
            "mag": NUV magnitude
            "mag_err": NUV magnitude uncertainty
            "dist": Visit source distance to coadd sources
            "pos_err": Photometric positional uncertainty
            "artifact": GALEX catalog artifact ID
        """

        # use custom style sheet
        plt.style.use(ROOT_DIR + "/lib/mpl_style_sheets/spie_scout_testing.mplstyle")

        log_prefix = str.format(
            "[{}.{}]",
            self.__class__.__name__,
            inspect.currentframe().f_code.co_name,
        )

        # validate input key
        valid_keys = ["s2n", "mag", "magerr", "dist", "poserr", "artifact"]
        assert key in valid_keys, "{}: Unknown parameter key. Choose from [{}]".format(
            log_prefix, ", ".join(valid_keys)
        )

        # construct catalog key
        visit_key = "visit_nuv_" + key if not key == "dist" else "visit_" + key
        coadd_key = None
        if coadd:
            coadd_key = "nuv_" + key if not key == "dist" else None

        dd_plot_settings = {
            "title": {
                "s2n": "Signal-to-noise",
                "mag": "AB magnitude",
                "magerr": "AB magnitude uncertainty",
                "dist": "Coadd-visit separation",
                "poserr": "Positional uncertainty",
                "artifact": "Catalog artifacts",
            },
            "xaxis_label": {
                "s2n": "S/N",
                "mag": "AB mag",
                "magerr": "AB mag uncertainty",
                "dist": "Separation [arcsec]",
                "poserr": "Uncertainty [arcsec]",
                "artifact": "ID",
            },
            "yaxis_label": {
                "s2n": "Number of sources",
                "mag": "Number of sources",
                "magerr": "Number of sources",
                "dist": "Number of sources",
                "poserr": "Number of sources",
                "artifact": "Number of sources",
            },
            "bin_size": {
                "s2n": 1.0,
                "mag": 0.1,
                "magerr": 0.1,
                "dist": 0.75,
                "poserr": 0.1,
                "artifact": 1,
            },
            "density": {
                "s2n": False,
                "mag": False,
                "magerr": False,
                "dist": False,
                "poserr": False,
                "artifact": False,
            },
            "xscale": {
                "s2n": "log",
                "mag": "linear",
                "magerr": "linear",
                "dist": "log",
                "poserr": "linear",
                "artifact": "log",
            },
            "yscale": {
                "s2n": "linear",
                "mag": "linear",
                "magerr": "linear",
                "dist": "log",
                "poserr": "log",
                "artifact": "log",
            },
            "invert_xaxis": {
                "s2n": False,
                "mag": True,
                "magerr": False,
                "dist": False,
                "poserr": False,
                "artifact": False,
            },
            "xlim_offset": {
                "s2n": [5, 10],
                "mag": [0.75, 0.75],
                "magerr": None,
                "dist": None,
                "poserr": None,
                "artifact": None,
            },
        }

        # Selection of visits data
        sel_all = self.tt_match["visit_dist"] > -1
        sel_pos = self.get_variability_mask(level=1)
        sel_var = self.get_variability_mask(level=2)
        sel_flare = self.get_variability_mask(level=3)

        selection_list = [sel_all, sel_pos, sel_var, sel_flare]
        selection_label = [
            "All",
            "visit_dist < {}".format(self.dd_cuts["separation"] * uu.arcsec),
            "visit_nuv_ds2n > {}".format(self.dd_cuts["ds2n"]),
            "visit_nuv_dmag > {}".format(self.dd_cuts["dmag"]),
        ]
        n_sel = len(selection_list)

        # figure
        fig, axs = plt.subplots(
            ncols=n_sel, nrows=1, sharex=True, sharey=True, figsize=figsize
        )
        fig.suptitle(
            "Source matches statistics: {}".format(dd_plot_settings["title"][key])
        )

        # compute bins
        if coadd and key == "dist":
            par_key = visit_key
        elif coadd:
            par_key = coadd_key
        else:
            par_key = visit_key

        x_min = np.floor(self.tt_match[par_key].min())
        x_max = np.ceil(self.tt_match[par_key].max())
        if bin_size is None:
            bins = int((x_max - x_min) / dd_plot_settings["bin_size"][key])
        else:
            bins = bin_size

        for i, ax in enumerate(axs.flat):
            # x-axis range
            if xlim_offset is None and dd_plot_settings["xlim_offset"][key] is not None:
                off_lo = dd_plot_settings["xlim_offset"][key][0]
                off_up = dd_plot_settings["xlim_offset"][key][1]
                ax.set_xlim(x_min - off_lo, x_max + off_up)
            elif xlim_offset is not None:
                off_lo = xlim_offset[0]
                off_up = xlim_offset[1]
                ax.set_xlim(x_min - off_lo, x_max + off_up)

            # invert x-axis direction
            if dd_plot_settings["invert_xaxis"][key]:
                ax.invert_xaxis()

            # axis scale
            if yscale is None:
                yscale = dd_plot_settings["yscale"][key]
            if xscale is None:
                xscale = dd_plot_settings["xscale"][key]
            ax.set_yscale(yscale)
            ax.set_xscale(xscale)

            # normalization
            if density is None:
                density = dd_plot_settings["density"][key]

            # coadd hist
            if coadd and key != "dist":
                coadd_data = self.tt_match[coadd_key][selection_list[i].any(axis=1)]
                n_star_coadd = len(coadd_data)
                ax.hist(
                    coadd_data,
                    bins=bins,
                    density=density,
                    histtype="stepfilled",
                    edgecolor="b",
                    facecolor="b",
                    alpha=0.4,
                    label="coadd: $N = {}$".format(n_star_coadd),
                )
            # visits hist
            visit_data = np.where(selection_list[i], self.tt_match[visit_key], np.nan)
            n_star_visit = selection_list[i].sum()
            ax.hist(
                visit_data,
                bins=bins,
                density=density,
                histtype="step",
                edgecolor="k",
                facecolor="k",
                alpha=0.2,
                label=["visit: $N = {}$".format(n_star_visit)],
            )
            ax.legend(fontsize=6)
            ax.set_ylabel(dd_plot_settings["yaxis_label"][key], fontsize=label_fontsize)
            if i > 0:
                y_axis = ax.axes.get_yaxis()
                y_label = y_axis.get_label()
                y_label.set_visible(False)
            ax.set_xlabel(dd_plot_settings["xaxis_label"][key], fontsize=label_fontsize)
            ax.set_title(selection_label[i], fontsize=label_fontsize)

        fig.tight_layout()

        # reset style to default style sheet
        plt.style.use("default")

    def show_light_curves(self):
        # Select variable sources
        cut_sep = self.dd_cuts["separation"] * uu.arcsec  # 2.5 arcsec separation
        cut_ds2n = self.dd_cuts["ds2n"]  # 4 signal-to-noise variation
        cut_dmag = self.dd_cuts["dmag"]  # 0.2 magnitude variation

        sel_vis_pos = self.tt_match["visit_dist"] < cut_sep
        sel_vis_ds2n = self.tt_match["visit_nuv_ds2n"] > cut_ds2n
        sel_vis_dmag = np.abs(self.tt_match["visit_nuv_dmag"]) > cut_dmag

        sel_pos = self.get_variability_mask(level=1)
        sel_vis = self.get_variability_mask(level=3)
        sel_flare = self.get_variability_mask(level=3, project_to_coadd=True)

        n_flares = sel_flare.sum()
        n_matched = len(self.tt_match)

        # sort by coadd brightness
        nuv_mag_coadd = self.tt_match[sel_flare]["nuv_mag"]
        mag_idx = np.argsort(nuv_mag_coadd)
        # print(nuv_mag_coadd[mag_idx])

        # multiple selections to plot
        selection_list = [
            sel_pos[sel_flare][mag_idx],
            sel_vis[sel_flare][mag_idx],
        ]

        # catalog of flares sorted by coadd brightness
        tt_flares = self.tt_match[sel_flare][mag_idx]

        # computes some statistics of the flare catalog

        # # compute y-axis offset to disentangle light curves
        # # bins
        # x_min = np.floor(nuv_mag_coadd.min())
        # x_max = np.ceil(nuv_mag_coadd.max())
        # bin_size = 1.0
        # bins = np.arange(
        #     x_min + 1, x_max + 1 + bin_size, bin_size
        # )  # +1 needed otherwise np.digitize yields an empty bin
        #
        # # array of corresponding bin indices
        # digi = np.digitize(nuv_mag_coadd[mag_idx], bins, right=True)
        # print("idx", digi)
        # # count occurrence of bin index
        # digi_bins = np.bincount(digi)
        # print("count", digi_bins)
        # digi_density = digi_bins / len(nuv_mag_coadd)
        # print("density", digi_density)
        # # array with cumulative offset corresponding to number density
        # y_offset = np.cumsum(15 * np.take(digi_density, digi))
        # # y_offset = np.zeros(n_flares)
        #
        # print("y_offset", np.take(digi_density, digi))
        # print("cum. y_offset", y_offset)

        # color bar indicates date
        dates_str = [
            datetime.fromisoformat(d.to_value("iso"))
            for d in Time(
                self.tt_visits["PhotoObsDate_MJD"],
                format="mjd",
            )
        ]
        dates_num = mdates.date2num(dates_str)
        vmin = dates_num[0]
        vmax = dates_num[-1]
        cmap = copy(cm.get_cmap("viridis"))
        cmap.set_bad([0.08, 0.06, 0.11])
        norm = colors.Normalize(vmin=vmin, vmax=vmax)

        loc = mdates.AutoDateLocator()

        if n_flares > 2:

            n_lc_per_sub = 7
            n_sub_per_fig = 20
            # dice into subplots of "n_lc_per_sub" light curves each
            # total number of subfigures
            n_sub = int(np.ceil(n_flares / n_lc_per_sub))
            # list of index arrays for every chunk of light curves
            idx_list_sub = np.split(
                np.arange(n_flares),  # flare index
                np.arange(0, n_flares, n_lc_per_sub)[1:],  # split index
            )
            # create a new figure every "n_sub_per_fig" subplots
            n_fig = int(np.ceil(n_sub / n_sub_per_fig))

            # loop over figures
            for i in range(0, n_fig):

                idx_lo = i * n_sub_per_fig
                idx_hi = n_sub_per_fig + i * n_sub_per_fig
                idx_list_sub_current_fig = idx_list_sub[idx_lo:idx_hi]
                n_sub_current_fig = len(idx_list_sub_current_fig)
                # figure
                fig, axs = plt.subplots(
                    nrows=n_sub_current_fig,
                    figsize=(7, 5 * n_sub_current_fig),
                    constrained_layout=True,
                )
                fig.suptitle("Light curves")
                # style cycles
                cycle_marker = self.plot_style_cycler("marker", shuffle=True)
                cycle_color = self.plot_style_cycler("color")
                cycle_line_style = self.plot_style_cycler("line_style")

                # loop over subplots
                for j, (ax, idx) in enumerate(zip(axs, idx_list_sub_current_fig)):

                    # colorbar
                    fig.colorbar(
                        cm.ScalarMappable(norm=norm, cmap=cmap),
                        ax=ax,
                        orientation="horizontal",
                        ticks=loc,
                        format=mdates.AutoDateFormatter(loc),
                        shrink=0.75,
                        aspect=85,
                    )

                    # loop over light curves
                    for k, lc_idx in enumerate(idx):
                        marker = next(cycle_marker)
                        marker_size = 6
                        color = next(cycle_color)
                        color_list = ["gray", color]
                        alpha_list = [0.33, 1.0]
                        marker_size_list = [marker_size, marker_size]
                        line_style = next(cycle_line_style)
                        line_style_list = [line_style, ""]
                        x_offset_coadd = -1
                        # y_offset = 0.5 * k

                        # catalog source data (single row)
                        tt_source = tt_flares[lc_idx]

                        # loop over selections (cut levels)
                        for m, selection in enumerate(selection_list):
                            # select visits corresponding to the coadd flare source
                            sel_source = selection[lc_idx]
                            # brightness
                            vis_mag = tt_source["visit_nuv_mag"][sel_source]
                            vis_magerr = tt_source["visit_nuv_magerr"][sel_source]
                            # date
                            vis_mjd = self.tt_visits["PhotoObsDate_MJD"][sel_source]
                            vis_dates = dates_num[sel_source]
                            # visit number
                            vis_num = np.arange(selection.shape[1])[sel_source]

                            x_values = vis_num
                            # line with error bars
                            ax.errorbar(
                                x_values,
                                vis_mag,
                                vis_magerr,
                                ls=line_style_list[m],
                                c="gray",
                                marker=None,
                                label=tt_source["ggoid_dec"],
                                lw=1.0,
                                zorder=1,
                            )
                            # markers colored according to observation date
                            ax.scatter(
                                x_values,
                                vis_mag,
                                marker=marker,
                                s=marker_size_list[m] ** 2,
                                c=vis_dates,
                                alpha=alpha_list[m],
                                cmap=cmap,
                                norm=norm,
                                zorder=2,
                            )

                        ax.scatter(
                            x_offset_coadd,
                            tt_source["nuv_mag"],
                            color="k",
                            marker=marker,
                            s=marker_size ** 2,
                        )
                        ax.set_xticks(np.arange(selection.shape[1]))

                    ax.set_xlabel("Visit Number (first point coadd)")
                    ax.set_ylabel("nuv_mag")
                    ax.invert_yaxis()
