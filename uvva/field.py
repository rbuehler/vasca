import os
from datetime import datetime
from itertools import cycle

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from astropy import units as uu
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Column, Table, conf, vstack
from astropy.time import Time
from astropy.wcs import wcs
from astroquery.mast import Observations
from loguru import logger
from sklearn.cluster import MeanShift, estimate_bandwidth

from uvva.resource_manager import ResourceManager
from uvva.tables import TableCollection
from uvva.utils import get_time_delta, get_time_delta_mean, sky_sep2d, table_to_array

# global paths
# path to the dir. of this file
FILE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = FILE_DIR + "/../"  # path to the root directory of the repository

dimless = uu.dimensionless_unscaled
conf.replace_warnings = ["always"]


class BaseField(TableCollection):
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
        # Sets skeleton
        super().__init__()

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

        #: Reference wcs
        self.ref_wcs = None

    def plot_sky_sources(
        self,
        ax=None,
        plot_detections=True,
        plot_src_id=True,
        src_kwargs=None,
        det_kwargs=None,
    ):
        """
        Plot the selected sources and (optinally) the visit detections on the sky.

        Parameters
        ----------
        ax : axes, optional
            Matplotlib axes to plot on. The default is None.
        plot_detections : bool, optional
            Plot the visit detections below the sources. The default is True.
        plot_src_ids: bool, optional
            Write the source ID next to its marker. Default is True.
        src_kwargs : dict, optional
            Keyword arguments for pyplot.plot of the sources. The default is None.
        det_kwargs : dict, optional
            Keyword arguments for pyplot.plot of the detections. The default is None.

        Returns
        -------
        ax : axes
            Used Matplotlib axes.

        """
        logger.debug("Plotting sky sources")

        if ax is None:
            ax = plt.gca()

        # Set marker properties for sources
        plt_src_kwargs = {
            "marker": "o",
            "markersize": 1.5,
            "alpha": 0.5,
            "lw": 0,
            "markeredgewidth": 0.3,
            "fillstyle": "none",
        }
        if src_kwargs is not None:
            plt_src_kwargs.update(src_kwargs)

        # Set marker properties for detections
        plt_det_kwargs = {
            "marker": ".",
            "markersize": 0.2,
            "alpha": 0.5,
            "markeredgewidth": 0.0,
            "lw": 0,
        }
        if det_kwargs is not None:
            plt_det_kwargs.update(det_kwargs)

        # Prepare data
        sel = self.tt_sources["sel"]
        tt_src = self.tt_sources[sel]
        tt_det = self.tt_detections
        tt_det.add_index("src_id")

        # Loop over all srcs and plot
        colors = cycle("bgrcmykbgrcmykbgrcmykbgrcmyk")
        for src, col in zip(tt_src, colors):
            if plot_detections:
                det_idx = tt_det.loc_indices["src_id", src["src_id"]]
                ax.plot(
                    tt_det[det_idx]["ra"].data,
                    tt_det[det_idx]["dec"].data,
                    color=col,
                    **plt_det_kwargs,
                )
            ax.plot(src["ra"], src["dec"], color=col, **plt_src_kwargs)

            ax.text(
                src["ra"] + 0.005,
                src["dec"] + 0.004,
                str(src["src_id"]),
                transform=plt_src_kwargs["transform"],
                fontsize=1,
                color=col,
                alpha=0.5,
            )

        return ax

    def plot_sky_map(self, ax=None, **img_kwargs):
        """
        Plot the reference sky map.

        Parameters
        ----------
        ax : axes, optional
            Matplotlib axes to plot on. The default is None.
        **img_kwargs : dict
            Key word arguments for pyplot.imshow plotting.

        Returns
        -------
        graph : AxesImage
            Matplotlib axes of 2D image.

        """

        logger.debug("Plotting sky map'")

        if self.ref_img is None:
            logger.error("No map to draw")

        if ax is None:
            ax = plt.gca()

        plt_img_kwargs = {
            "interpolation": "None",
            "cmap": "gist_yarg",
            "origin": "lower",
            "norm": LogNorm(),
        }

        if img_kwargs is not None:
            plt_img_kwargs.update(img_kwargs)

        graph = ax.imshow(self.ref_img, **plt_img_kwargs)

        return graph

    def plot_sky(self, plot_detections=True, plot_map=True):
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

        logger.debug("Plotting sky map and/or sources'")

        fig = plt.figure(figsize=(8, 7))
        plt.clf()

        if self.ref_wcs is not None:
            ax = plt.subplot(projection=self.ref_wcs)  #
            ax.coords["ra"].set_major_formatter("d.dd")
            ax.coords["dec"].set_major_formatter("d.dd")
            ax.set_xlabel("Ra")
            ax.set_ylabel("Dec")

        if plot_map:
            graph = self.plot_sky_map(ax)
            fig.colorbar(graph, label="Intensity [a.u.]")

        if plot_map:
            plot_arg = {"transform": ax.get_transform("world")}
            self.plot_sky_sources(
                ax,
                plot_detections=plot_detections,
                src_kwargs=plot_arg,
                det_kwargs=plot_arg,
            )
        else:
            self.plot_sky_sources(plot_detections=plot_detections)

        return fig

    def plot_light_curve(
        self,
        src_ids,
        ax=None,
        ylim=[23.5, 16.5],
        legend_loc="upper right",
        plot_upper_limits=True,
        **errorbar_kwargs,
    ):
        """
        Plot the magnitude light curves of the passed sources.

        Parameters
        ----------
        src_ids : list or int
            List or single source IDs to plot.
        ax : axes, optional
                Matplotlib axes to plot on. The default is None.
        legend_loc : string, optional
            Position of the legend in the figure. The default is "upper right".
        plot_upper_limits : bool
            Plot upper limits to the lightcurve. The default is True.
        **errorbar_kwargs : TYPE
            Key word arguments for pyplot.errorbars plotting.

        Returns
        -------
        ax : axes
            Used Matplotlib axes.

        """

        # If only one src_id was passed create a list
        if not hasattr(src_ids, "__iter__"):
            src_ids = [src_ids]

        logger.debug("Plotting light curves'")

        # Setup plotting parameters
        if ax is None:
            ax = plt.gca()
        ax.invert_yaxis()
        ax.set_ylim(ylim)

        plt_errorbar_kwargs = {
            "markersize": 6,
            "capsize": 2,
            "lw": 0.2,
            "linestyle": "dotted",
            "elinewidth": 1,
        }
        if errorbar_kwargs is not None:
            plt_errorbar_kwargs.update(errorbar_kwargs)

        # Loop over selected sources and plot
        colors = cycle("bgrcmykbgrcmykbgrcmykbgrcmyk")
        markers = cycle("osDd.<>^vpP*")
        for src_id, col, mar in zip(src_ids, colors, markers):

            # Get light curve
            lc = self.get_light_curve(src_id)

            # Get arrays
            src_lab = "src_" + str(src_id)
            uplims = np.zeros(len(lc))
            sel = lc["mag"] > 0
            mags = lc["mag"]
            mags_err = lc["mag_err"]

            # Modify arrays if upper limits are plotted
            if plot_upper_limits:
                uplims = lc["ul"] < 0
                sel = np.ones(len(lc), dtype=bool)
                mags = mags * ~uplims + mags_err * uplims
                mags_err = mags_err * ~uplims + 1 * uplims

            # Plot
            plt.errorbar(
                lc["time_start"][sel],  # TODO: Move this to the bin center
                mags[sel],
                yerr=mags_err[sel],
                lolims=uplims[sel],
                color=col,
                marker=mar,
                label=src_lab,
                **plt_errorbar_kwargs,
            )
        ax.legend(loc=legend_loc)
        ax.set_xlabel("MJD")
        ax.set_ylabel("Magnitude")

        return ax

    def load_sky_map(self, file_name):
        """
        Load reference image and WCS from passed fits file

        Parameters
        ----------
        filename : str
            File name of image FITS

        Returns
        -------
        None.

        """
        logger.info(f"Loading skypmap from file: '{file_name}'")
        with fits.open(file_name) as ff:
            self.ref_wcs = wcs.WCS(ff[0].header)
            self.ref_img = ff[0].data

    def cluster_meanshift(self, add_upper_limits=True, **ms_kw):
        """
        Apply _MeanShift clustering algorithm using to derive sources.

        .. _MeanShift: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MeanShift.html

        Parameters
        ----------
        add_upper_limits : bool, optional
            Add upper limits to the tt_sources_lc, for visits with no detection.
            Upper limits are stored in the mag_err columns for none detections.
            The default is True.
        ms_kw : dict, optional
            Keywords passed to the scikit MeanShift function. Note that the
            bandwidth is assumed to be in units of arc seconds.

        Returns
        -------
        int
            Number of detected clusters.
        """
        logger.info("Clustering sources with MeanShift")

        # Selection
        sel = self.tt_detections["sel"]

        # Get detection coordinates and run clustering
        coords = table_to_array(self.tt_detections[sel]["ra", "dec"])

        # Do bandwidth determination "by hand" to print it out and convert
        # bandwidth unit from arc seconds into degerees
        dd_ms = ms_kw
        if not "bandwidth" in ms_kw or ms_kw["bandwidth"] == None:
            logger.debug(f"Estimating bandwidth")
            dd_ms["bandwidth"] = estimate_bandwidth(coords, quantile=0.2, n_samples=500)
        else:
            dd_ms["bandwidth"] = (ms_kw["bandwidth"] * uu.arcsec).to(uu.deg).value

        logger.debug(f"MeanShift with parameters (bandwith in degrees): '{dd_ms}' ")
        ms = MeanShift(**dd_ms)

        ms.fit(coords)

        # Fill in data into field tables
        src_ids, det_cts = np.unique(ms.labels_, return_counts=True)

        cluster_centers = ms.cluster_centers_
        nr_srcs = len(cluster_centers)
        srcs_data = {
            "src_id": src_ids,
            "ra": cluster_centers[:, 0],
            "dec": cluster_centers[:, 1],
            "nr_det": det_cts,
            "nr_det_meas": np.zeros(nr_srcs),
            "mag_mean": np.zeros(nr_srcs) - 1,
            "mag_var": np.zeros(nr_srcs) - 1,
            "mag_rchiq": np.zeros(nr_srcs) - 1,
            "mag_dmax": np.zeros(nr_srcs) - 1,
            "mag_dmax_sig": np.zeros(nr_srcs) - 1,
            "perc_ul_mean": np.zeros(nr_srcs) - 1,
        }

        # Fill information into tables.
        self.add_table(srcs_data, "base_field:tt_sources")
        self.tt_sources.meta["CLUSTALG"] = "MeanShift"

        # Update src_id entries
        self.tt_detections["src_id"][sel] = ms.labels_

        # Fill light curve data into tables
        self.remove_double_visit_detections()
        self.set_light_curve(add_upper_limits=add_upper_limits)

        return nr_srcs

    def get_upper_limits(self):
        """
        Calculates magnitude upper limits on non-detections to the tt_visits table.

        Returns
        -------
        upper_limits : float array
            Array of upper limits for each visit

        """

        observatory = (
            self.observatory
            if hasattr(self, "observatory")
            else self.get_field_par("observatory", "tt_field")
        )
        upper_limit = None

        # Call upper limit function according to observatory.
        # String matching is case insensitive
        if observatory.casefold() == "GALEX".casefold():
            upper_limit = GALEXField.get_visit_upper_limits(self.tt_visits)

        return upper_limit

    def set_light_curve(self, add_upper_limits=True):
        """
        Helper function of cluster_meanshift().
        Adds detections information into ta_sources_lc.

        Parameters
        ----------
        add_upper_limits : bool, optional
            Add upper limits to the tt_sources_lc, for visits with no detection.
            Upper limits are stored in the mag_err columns for none detections.
            The default is True.

        Returns
        -------
        None.

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
            "src_id": list(),
            "mag": list(),
            "mag_err": list(),
            "ul": list(),
        }

        # Loop over sources and add them to dictionary
        tt_det_grp = self.tt_detections[sel].group_by(["src_id"])
        for tt_det in tt_det_grp.groups:

            # Add src id
            src_id = tt_det["src_id"][0]
            tdata["src_id"].append(src_id.tolist())
            vis_idxs = self.tt_visits.loc_indices["vis_id", tt_det["vis_id"]]

            np_mag = np.zeros(nr_vis) - 1
            np_mag[vis_idxs] = tt_det["mag"]
            tdata["mag"].append(np_mag.tolist())

            np_mag_err = np.zeros(nr_vis) - 1
            np_mag_err[vis_idxs] = tt_det["mag_err"]
            tdata["mag_err"].append(np_mag_err.tolist())

            # Store upper limits if no detection in a visit
            # TODO: make this more general and not GALEX specific
            if add_upper_limits:
                np_mag_ul = self.get_upper_limits()
                tdata["ul"].append(np_mag_ul.tolist())

        # Add light curve table
        tdata["mag"] = np.array(tdata["mag"], dtype=np.object_)
        tdata["mag_err"] = np.array(tdata["mag_err"], dtype=np.object_)
        tdata["ul"] = np.array(tdata["ul"], dtype=np.object_)
        self.add_table(tdata, "base_field:ta_sources_lc", add_sel_col=False)
        self.set_var_stats()

    def set_var_stats(self):
        """
        Calculates source variability parameters and stores them
        in the source table (tt_source).

        Returns
        -------
        None.

        """

        logger.debug(f"Calculating source variability statistics.")

        if "ta_sources_lc" not in self._table_names:
            logger.error(
                "Light curve table does not exist, run 'set_light_curve()' first."
            )

        # Get lightcurve as numpy arrays to calculate stats
        mag = np.array((self.ta_sources_lc["mag"].data).tolist())
        mag_err = np.array((self.ta_sources_lc["mag_err"].data).tolist())
        mag_ul = np.array((self.ta_sources_lc["ul"].data).tolist())

        # Ignore entries with no magniture or valid error
        mask_mag = mag < 1e-6
        mask_mag_err = mag_err < 1e-6
        mask = mask_mag + mask_mag_err

        mag[mask] = np.nan
        mag_err[mask] = np.nan

        # Calculate variability parameters
        nr_mags = (~np.isnan(mag)).sum(axis=1)
        mag_mean = np.nanmean(mag, axis=1)

        mag_err_mean2 = np.nanmean(mag_err * mag_err, axis=1)
        mag_var = np.nanvar(mag, axis=1)
        rchiq_const = mag_var / mag_err_mean2

        # Get the maximum flux variation from the mean
        dmag_max = np.abs(np.nanmax(mag, axis=1) - mag_mean)
        dmag_min = np.abs(np.nanmin(mag, axis=1) - mag_mean)
        dmag = (dmag_max >= dmag_min) * dmag_max + (dmag_max < dmag_min) * dmag_min

        # Get maximum significance of flux variation compared to mean
        dmag_max_sig = np.nanmax(np.abs((mag - mag_mean[:, None]) / mag_err), axis=1)

        # Nr of upper limits below the mean flux (in magnitudes greater)
        nr_ulmean = (mag_mean[:, None] < mag_ul).sum(axis=1)

        # Write them into tt_sources
        src_ids = self.ta_sources_lc["src_id"]
        self.tt_sources.add_index("src_id")
        src_idx = self.tt_sources.loc_indices["src_id", src_ids]

        self.tt_sources["nr_det_meas"][src_idx] = nr_mags
        self.tt_sources["mag_mean"][src_idx] = mag_mean
        self.tt_sources["mag_var"][src_idx] = mag_var
        self.tt_sources["mag_rchiq"][src_idx] = rchiq_const
        self.tt_sources["mag_dmax"][src_idx] = dmag
        self.tt_sources["mag_dmax_sig"][src_idx] = dmag_max_sig
        self.tt_sources["perc_ul_mean"][src_idx] = (
            100.0 * nr_ulmean / len(self.tt_visits)
        )

    def remove_double_visit_detections(self):
        """
        Helper function of cluster_meanshift().
        Remove multiple detections of one source in one visit from tt_detections.
        Keep only the closest detection.

        Returns
        -------
        list
            List of removed det_ids

        """
        logger.debug("Scanning for doubled visit detections")

        # Selection
        sel = self.tt_detections["sel"]

        # Index tables
        self.tt_sources.add_index("src_id")
        self.tt_detections.add_index("det_id")

        # Determine detection_id of srcs with more than one detection in one visit
        rm_det_ids = []
        tt_det_grp = self.tt_detections[sel].group_by(["vis_id", "src_id"])

        for tt_det in tt_det_grp.groups:
            if len(tt_det) > 1:

                # Get source coordinate
                src_id = tt_det["src_id"].data[0]
                src_idx = self.tt_sources.loc_indices["src_id", src_id]
                src_ra = self.tt_sources["ra"].quantity[src_idx]
                src_dec = self.tt_sources["dec"].quantity[src_idx]
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
                f"Removed double-visit detections: {nr_rm_det} ({perc_rm_det: .2f} %)"
            )

            # remove the doubled detections from tt_detections
            rm_idx = self.tt_detections.loc_indices["det_id", rm_det_ids]
            self.tt_detections.remove_rows(rm_idx)

            # Update detection count in tt_sources
            sel = self.tt_detections["sel"]
            src_ids, det_cts = np.unique(
                self.tt_detections[sel]["src_id"], return_counts=True
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
            par = self.tt_field[par_name].data[0]
            # Maintains unit
            # Some combinations of unit and dtype need special handling
            unit = self.tt_field[par_name].unit
            dtype_kind = self.tt_field[par_name].dtype.kind
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


class GALEXField(BaseField):
    """
    Instance of one GALEX field
    """

    def __init__(self, obs_id, obs_filter=None, data_path=None, visits_data_path=None):
        """
        Initializes a new GALEXField instance with
        skeleton UVVA data structure.

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
        uvva_file_prefix : str
            File name prefix following UVVA naming convention:
            'UVVA_<observatory>_<filed_id>_<obs_filter>'
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

        # Create and check existence of directory that holds field data and UVVA outputs
        if not os.path.isdir(self.data_path):
            os.makedirs(self.data_path)

        # File name prefix for UVVA/GALEXField outputs
        self.uvva_file_prefix = f"UVVA_GALEX_{obs_id}_{obs_filter}"

        logger.debug(f"Field data path set to: '{self.data_path}'")
        logger.debug(f"Visits data path set to: '{self.visits_data_path}'")

    @classmethod
    def from_UVVA(cls, obs_id, obs_filter="NUV", fits_path=None, **kwargs):
        """
        Constructor to initialize a GALEXField instance
        from a UVVA-generated FITS file

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
        uvva.field.GALEXField
        """
        # Bootstrap the initialization procedure using the base class
        gf = cls(obs_id, obs_filter, **kwargs)  # new GALEXField instance

        if fits_path is None:
            # Construct the file name from field ID and filter
            fits_path = f"{gf.data_path}/{gf.uvva_file_prefix}_field_data.fits"
        # Check if file exists
        if not os.path.isfile(fits_path):
            raise FileNotFoundError(
                "Wrong file or file path to UVVA data "
                f"for GALEX field '{obs_id}' with obs_filter '{obs_filter}'."
            )
        # Reads the UVVA-generated field data
        gf.load_from_fits(fits_path)
        # Sets convenience class attributes
        gf.set_field_attr()
        # Check consistency
        if not gf.field_id == obs_id:
            raise ValueError(
                "Inconsistent data: Missmatch for 'field_id'. "
                f"Expected '{obs_id}' but got '{gf.field_id}' "
                f"from file '{fits_path.split(os.sep)[-1]}."
            )
        elif not gf.obs_filter == obs_filter:
            raise ValueError(
                "Inconsistent data: Missmatch for 'obs_filter'. "
                f"Expected '{obs_filter}' but got '{gf.obs_filter}' "
                f"from file '{fits_path.split(os.sep)[-1]}."
            )

        logger.info(
            f"Loaded UVVA data for GALEX field '{obs_id}' with obs_filter '{obs_filter}'."
        )

        return gf

    @classmethod
    def from_MAST(
        cls,
        obs_id,
        obs_filter="NUV",
        refresh=False,
        load_products=True,
        write=True,
        **kwargs,
    ):
        """
        Constructor to initialize a GALEXField instance either
        fresh from the MAST archive (refresh=True) or if available
        from cached raw data (refresh=False).

        The procedure uses uvva.resource_manager.ResourceManager
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
        load_products : bool, optional
            Selects if data products shall be loaded. Othervise only field and
            visit information is loaded.
        write: bool, optional
            If load_products is enabled, stores the data as UVVA tables in the cloud
            for faster loading in the future. Default is True.
        **kwargs
            All additional keyword arguments are passed to `~GALEXField.__init__()`

        Returns
        -------
        uvva.field.GALEXField

        """
        # Checks
        if not isinstance(refresh, bool):
            raise TypeError(f"Expected boolean argument, got {type(refresh).__name__}.")

        # Bootstrap the initialization procedure using the base class
        gf = cls(obs_id, obs_filter, **kwargs)  # new GALEXField instance

        # Sets ``gf.tt_field``
        gf._load_galex_field_info(obs_id, obs_filter, refresh=refresh)
        # Sets ``gf.tt_visits``
        gf._load_galex_visits_info(obs_id, obs_filter)

        # Sets convenience class attributes
        gf.set_field_attr()

        # Sets ``gf.tt_detections``, ``gf.tt_ref_sources`` and loads the ref image
        if load_products:
            gf._load_galex_archive_products(obs_id, obs_filter, refresh=refresh)
            meta_only = "."
            if write:
                fits_name = f"{gf.data_path}/{gf.uvva_file_prefix}_field_data.fits"
                gf.write_to_fits(fits_name)
        else:
            meta_only = ", metadata-only."

        logger.info(
            f"Loaded new GALEX field '{obs_id}' with obs_filter '{obs_filter}'"
            f"from MAST data {meta_only}"
        )

        return gf

    @staticmethod
    def load(
        field_id,
        obs_filter="NUV",
        method="MAST_LOCAL",
        load_products=True,
        **field_kwargs,
    ):
        """
        Loads GALEX field data according to a given method and
        returns a GALEXField instance.

        Parameters
        ----------
        field_id : int
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
            UVVA: Builds a new GALEXField based on a UVVA-generated field data file.
            AUTO: Attempts to load field data by using UVVA as method and
            falls back to MAST_LOCAL if no UVVA file is found on disc.

            The default directory where field data availability is checked is
            defined by the "data_path" attribute of GALEXField and can be passed
            explicitly via the "field_kwargs".
        load_products : bool, optional
            Specifies if GALEXField should be loaded completely with the
            full set of data products (default) or just containing metadata (False).
        field_kwargs
            All additional keyword arguments are passed to the load methods
            `~GALEXField.from_MAST()` and `~GALEXField.from_UVVA()`,
            as well as to `~GALEXField.__init__()`.

        Raises
        ------
        TypeError
            If the specified load method is not a string.
        ValueError
            If the specified load method is not one of
            '["mast_remote", "mast_local", "uvva", "auto"]'. String matching is
            case insensitive.

        Returns
        -------
        uvva.field.GALEXField

        """
        logger.info(f"Loading data for field '{field_id}' with method '{method}'.")

        # Checks method argument
        # String matching is case insensitive (converts to all to lower case).
        if not isinstance(method, str):
            raise TypeError(
                "Expected string type for load method specification, "
                f"got '{type(method).__name__}'."
            )

        method = method.casefold()  # ensures all-lower-case string
        method_spec = ["mast_remote", "mast_local", "uvva", "auto"]
        if method not in method_spec:
            raise ValueError(
                "Expected load method specification from {method_spec}, "
                "got '{method}'."
            )

        # Sets refresh option for MAST methods
        if method == "mast_remote":
            refresh = field_kwargs.pop("refresh", True)
        else:
            refresh = field_kwargs.pop("refresh", False)

        # Loads field according to load method specification
        if method in ["mast_remote", "mast_local"]:
            gf = GALEXField.from_MAST(
                obs_id=field_id,
                obs_filter=obs_filter,
                refresh=refresh,
                load_products=load_products,
                **field_kwargs,
            )
        elif method == "uvva":
            # removes unused options
            field_kwargs.pop("refresh", None)
            field_kwargs.pop("write", None)

            gf = GALEXField.from_UVVA(
                obs_id=field_id,
                obs_filter=obs_filter,
                **field_kwargs,
            )
        elif method == "auto":
            pass
            # TODO: Lookahead via rm to check data availability.
            # Then "UVVA" is preferred for performance reasons.
            # Fallback to "MAST" & refresh=True if "UVVA" fails for some reason
            # (e.g. not complete set of tables stored in the fits file).

        return gf

    @staticmethod
    def get_visit_upper_limits(tt_visits):
        """
        Calculates upper limits on non-detections to the tt_visits table
        following the article_ of Gezari et al. ApJ 766 (2013) p4

        .. _article: https://iopscience.iop.org/article/10.1088/0004-637X/766/1/60

        Returns
        -------
        upper_limits : float array
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
        """
        # Uses default columns if not otherwise specified
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
        elif not isinstance(col_names, list):
            raise TypeError(
                "Expected list type for argument 'col_names', "
                f"got '{type(col_names).__name__}'."
            )

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
        logger.debug("Constructing 'tt_field'.")
        # Selects columns and strip mask
        if tt_coadd.masked is True:
            tt_coadd_select = Table(tt_coadd[col_names].as_array().data)
        else:
            tt_coadd_select = tt_coadd[col_names]
        # Converts field id dtype from unicode ('U19') to bytestring ('S64')
        if tt_coadd_select["obs_id"].dtype.kind == "U":
            tt_coadd_select.replace_column(
                "obs_id",
                tt_coadd_select["obs_id"].astype(np.dtype("S64")),
            )
        # Sets table as class attribute
        self.add_table(tt_coadd_select, "base_field:tt_field")

    def _load_galex_visits_info(self, obs_id, obs_filter, col_names=None):

        # Uses default columns if not otherwise specified
        # Already sets the order in which columns are added later on
        if col_names is None:
            col_names = [
                "imgRunID",
                "minPhotoObsDate",
                "nexptime" if obs_filter == "NUV" else "fexptime",
                "fexptime" if obs_filter == "NUV" else "nexptime",
                "RATileCenter",
                "DECTileCenter",
            ]
        elif not isinstance(col_names, list):
            raise TypeError(
                "Expected list type for argument 'col_names', "
                f"got '{type(col_names).__name__}'."
            )
        # Reads cached file if found on disc
        logger.debug(
            "Reading archive visit info from cashed file " f"'{self.visits_data_path}'"
        )
        tt_visits_raw = Table.read(self.visits_data_path)

        logger.debug("Constructing 'tt_visits'.")
        # Filters for visits corresponding to field id and selects specified columns
        tt_visits_raw_select = tt_visits_raw[tt_visits_raw["ParentImgRunID"] == obs_id][
            col_names
        ]
        # Converts string time format to mjd
        for key in ["minPhotoObsDate"]:
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

        # Sets table as class attribute
        self.add_table(tt_visits_raw_select, "galex_field:tt_visits")

    def _load_galex_archive_products(
        self,
        obs_id,
        obs_filter,
        col_names=None,
        dd_products=None,
        ref_maps_only=True,
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
        re_maps_only : bool, optional
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
            tt_data = Table(
                Observations.get_product_list(tt_coadd).as_array().data
            )  # Returns unmasked table
            # Saves to disc
            # tt_data.write(path_tt_data, overwrite=True)
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
                    "file_name": "nd-int.fits.gz"
                    if obs_filter == "NUV"
                    else "fd-int.fits.gz",
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
        tt_ref_sources_raw = Table.read(path_tt_ref)
        # Opens visit catalogs, add columns for visit and source IDs
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

        obs_filter_l = obs_filter.lower()  # lower case obs_filter name

        # Convert positional error to degree
        # See GALEX docs mor details:
        # http://www.galex.caltech.edu/wiki/GCAT_Manual#Catalog_Column_Description
        for tbl in [tt_detections_raw, tt_ref_sources_raw]:
            tbl[f"{obs_filter_l}_poserr"].unit = uu.arcsec
            tbl.replace_column(
                f"{obs_filter_l}_poserr",
                Column(tbl[f"{obs_filter_l}_poserr"].to(uu.degree), dtype="float64"),
            )

        # Add to tables as class attributes
        if col_names is None:
            col_names = [
                "visit_id",
                "source_id",
                "ggoid_dec",
                "alpha_j2000",  # take band-merged quantities?
                "delta_j2000",  # take band-merged quantities?
                f"{obs_filter_l}_poserr",
                f"{obs_filter_l}_mag",
                f"{obs_filter_l}_magerr",
                f"{obs_filter_l}_s2n",
                "fov_radius",
                f"{obs_filter_l}_artifact",
                f"{obs_filter}_CLASS_STAR",
                "chkobj_type",
                f"{obs_filter}_FLUX_APER_4",
                f"{obs_filter}_FLUXERR_APER_4",
                f"{obs_filter}_FLUX_APER_3",
                f"{obs_filter}_FLUXERR_APER_3",
                "E_bv",
            ]
        elif not isinstance(col_names, list):
            raise TypeError(
                "Expected list type for argument 'col_names', "
                f"got '{type(col_names).__name__}'."
            )
        # set data as class attributes
        self.add_table(tt_detections_raw[col_names], "galex_field:tt_detections")
        logger.debug("Constructed 'tt_detections'.")
        self.add_table(tt_ref_sources_raw[col_names[2:]], "galex_field:tt_ref_sources")
        logger.debug("Constructed 'tt_ref_sources'.")

        # Intensity maps
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
        if not ref_maps_only:
            # Selects everything else but the coadd path
            path_int_map_visits = [
                self.data_path + path
                for path in tt_down[
                    np.logical_and(aa_sel_int_map, tt_down["ID"] != obs_id)
                ]["Local Path"].data.astype(str)
            ]  # list

        self.load_sky_map(path_int_map_ref)
        # Future: Optional loading of visit sky maps may happen here.
