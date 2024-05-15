#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualization related methods for VASCA
"""

from collections import OrderedDict
from itertools import cycle
import warnings

import astropy.units as uu
from astropy import constants as cc
import healpy as hpy
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.time import Time
from astropy.visualization.wcsaxes import SphericalCircle
from loguru import logger
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter
import matplotlib.cm as cm
from astropy.modeling import models, fitting
from astropy.utils.exceptions import AstropyUserWarning
import astropy.constants as const

from vasca.utils import (
    dd_id2filter,
    flux2mag,
    mag2flux,
    flux2mag_np,
    mag2flux_np,
    select_obs_filter,
    run_LombScargle,
    freq2period,
    period2freq,
    dd_spec_lines,
    mjd2yr,
    yr2mjd,
)

# %% sky plotting


def plot_sky_sources(
    tt_src,
    tt_det=None,
    only_selected=True,
    ax=None,
    src_id="rg_src_id",
    sky_region_wcs=None,
    draw_labels=True,
    src_kwargs=None,
    det_kwargs=None,
):
    """
    Plot the selected sources and (optinally) the visit detections on the sky.

    Parameters
    ----------
    tt_src : astropy.table.Table
        Source list to plot. Has to contain "ra", "dec" and srd id columns.
    tt_det : astropy.table.Table, optional
        Detection list to plot. Has to contain "ra", "dec" and src_id columns.
        Default is None.
    only_selected: bool, optional
        Show only selected sources. tt_src must have a "sel" column. default is True.
    ax : axes, optional
        Matplotlib axes to plot on. The default is None.
    src_id: str, optional
        Write the source ID next to its marker, saved in the passed column name,
        typically "rg_src_id", "fd_src_id" or None. Default is ""rg_src_id""
    sky_region_wcs: (regions.SkyRegion, WCS) , optional
        Plot only sources within the sky region. A WCS has to be passed along ina tuple.
        Default is None.
    draw_labels : bool, optional
        Draw labels next to sources showing "src_iud" entry. The default is True.
    src_kwargs : dict, optional
        Keyword arguments for pyplot.plot of the sources. The default is None.
    det_kwargs : dict, optional
        Keyword arguments for pyplot.patches.Polygon of the detections. The default is None.

    Returns
    -------
    matplotlib.axes
        Used Matplotlib axes.
    astropy.table.Table
        Sources plotted from tt_src
    """

    logger.debug("Plotting sky sources")

    if ax is None:
        ax = plt.gca()

    # Show only selected sources
    sel = np.ones(len(tt_src), dtype=bool)
    if only_selected and "sel" in tt_src.colnames:
        sel = tt_src["sel"]
    tt_src = tt_src[sel]

    # select sources within the region
    sel_reg = np.ones(len(tt_src), dtype=bool)
    if type(sky_region_wcs) is not type(None):
        src_coords = SkyCoord(tt_src["ra"], tt_src["dec"], frame="icrs")
        sel_reg = sky_region_wcs[0].contains(src_coords, sky_region_wcs[1])
        logger.debug(f"Limiting to region {sky_region_wcs[0]}")

    # Set marker properties for sources
    plt_src_kwargs = {
        "marker": "+",
        "markersize": 6.5,
        "alpha": 0.5,
        "lw": 0,
        "markeredgewidth": 2.0,
        "fillstyle": "none",
        "transform": ax.get_transform("world"),
    }
    if src_kwargs is not None:
        plt_src_kwargs.update(src_kwargs)

    # Set marker properties for detections
    plt_det_kwargs = {
        "fill": False,
        "lw": 0.8,
        "transform": ax.get_transform("world"),
    }
    if det_kwargs is not None:
        plt_det_kwargs.update(det_kwargs)

    plt_txt_kwargs = {
        "xycoords": ax.get_transform("world"),
        "fontsize": 12,
    }

    if type(tt_det) is not type(None):
        tt_det.add_index(src_id)

    # Loop over all srcs and plot
    colors = cycle("rbgcmrbgcmrbgcmrbgcm")
    lss = [":", "-", "--", "_."]
    for src, col in zip(tt_src[sel_reg], colors):
        # Set colors in tandem for srcs, det and label, unless specific color passed
        if src_kwargs is not None and "color" in src_kwargs.keys():
            col = plt_src_kwargs["color"]
        plt_src_kwargs["color"] = col
        plt_det_kwargs["color"] = col
        plt_txt_kwargs["color"] = col

        if type(tt_det) is not type(None) and src[src_id] in tt_det[src_id]:
            det_idxs = np.array(tt_det.loc_indices[src_id, src[src_id]]).flatten()
            for det_idx in det_idxs:
                coord_det = SkyCoord(
                    tt_det[det_idx]["ra"] * uu.deg,
                    tt_det[det_idx]["dec"] * uu.deg,
                    frame="icrs",
                )
                s_det = SphericalCircle(
                    coord_det,
                    tt_det[det_idx]["pos_err"] * uu.arcsec,
                    ls=lss[int(tt_det[det_idx]["obs_filter_id"] % 4)],
                    **plt_det_kwargs,
                )
                ax.add_patch(s_det)

        ax.plot(src["ra"], src["dec"], **plt_src_kwargs)

        # Add labels if src_id was passed
        if type(src_id) is not type(None) and draw_labels:
            ax.annotate(
                str(src[src_id]),
                (src["ra"], src["dec"]),
                **plt_txt_kwargs,
            )

    return ax, tt_src[sel_reg]


def plot_field_sky_map(
    field, fig=None, ax=None, img_idx=-1, sky_region=None, **img_kwargs
):
    """
    Plot the reference sky map.

    Parameters
    ----------
    field: vasca.BaseField
        VASCA field to be plotted.
    fig: figure, optional
        Matplotlib figure to draw on, if None a new figure is created. The default is None.
    ax : axes, optional
        Matplotlib axes to plot on. The default is None.
    img_idx : int, optional
        Index nr of the visit in the tt_visits table. If -1 is passed
        the reference image is shown. Default is -1.
    sky_region: regions.SkyRegion , optional
        Plot only within the sky region. The field WCS will be used
        Default is None.
    **img_kwargs : dict
        Key word arguments for pyplot.imshow plotting.

    Returns
    -------
    matplotlib.figure
        Matplotlib figure used to draw
    matplotlib.graph.AxesImage
        Matplotlib axes of 2D image.

    """

    logger.debug("Plotting sky map'")

    if field.ref_img is None and field.vis_img is None:
        logger.error("No map to draw")

    # Setup imshow parameters
    plt_img_kwargs = {
        "interpolation": "None",
        "cmap": "gray_r",
        "origin": "lower",
        "norm": LogNorm(),
    }
    if img_kwargs is not None:
        plt_img_kwargs.update(img_kwargs)

    # Check if reference or visit images should be plotted
    fig_title = f"Field {field.field_id} reference image"
    plot_img = field.ref_img
    if img_idx > -1:
        # Check if visit images are there and get the right image from the cube
        if img_idx > len(field.tt_visits) - 2:
            logger.error(f"Requested image index {img_idx} is out of range")
        elif type(field.vis_img) == type(None):
            raise TypeError(
                f"Visit images not set, got image type '{type(field.vis_img)}'."
            )
        elif field.vis_img.ndim == 3:
            plot_img = field.vis_img[img_idx, :, :]
        else:
            plot_img = field.vis_img
        fig_title = f"Visit {field.tt_visits[img_idx]['vis_id']} image"

    # Check if plotting should be restricted to one region only
    wcs = field.ref_wcs
    if type(sky_region) is not type(None):
        box = sky_region.to_pixel(field.ref_wcs).bounding_box
        cutout = Cutout2D(
            plot_img,
            position=[box.center[1], box.center[0]],
            size=[box.shape[1], box.shape[0]],
            wcs=field.ref_wcs,
            copy=True,
        )
        plot_img = cutout.data
        wcs = cutout.wcs

    # Start drawing
    # Check if figure was passed
    if type(fig) is type(None):
        fig = plt.figure(figsize=(8, 7))  # , constrained_layout=True
    if type(ax) is type(None):
        ax = plt.subplot(projection=wcs)
        ax.coords["ra"].set_major_formatter("d.dd")
        ax.coords["dec"].set_major_formatter("d.dd")
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
    else:
        plt.gcf()

    # Check if axis was passed and setup axis
    if ax is None:
        ax = plt.gca()

    # Add coordinate grid
    ax.coords.grid(True, color="grey", ls="-", lw=0.5)

    ax.imshow(plot_img, **plt_img_kwargs)
    # if type(pix_patch) is not type(None):
    # graph.set_clip_path(pix_patch)
    ax.set_title(fig_title)

    return ax, wcs


# TODO: Add coordinate system match check between tt_coverage_hp and this function
def plot_region_sky_gnomeview(
    region, ra, dec, sel_srcs=True, gw_kwargs=None, ps_kwargs=None
):
    """
    Plot the nr visits and optionally (selected) sources (optinally) on the sky.

    Parameters
    ----------
    region : vasca.region.Region
        Region for which is plotted.
    ra : float
        RA of the center of the sky figure.
    dec: float
        DEC of the center of the sky figure.
    sel_srcs : bool, optional
        Show only selected sources. The default is True.
    gw_kwargs : dict, optional
        Keyword arguments for healpy.gnomview. The default is None.
    ps_kwargs : dict, optional
        Keyword arguments for healpy.projscatter. The default is None.

    Returns
    -------
    matplotlib.axes
        Used Matplotlib axes.

    """
    if not hasattr(region, "tt_coverage_hp"):
        logger.warning(
            "Create healpix coverage table first with region.add_coverage_hp. "
            "Creating it now with default values."
        )
        region.add_coverage_hp()

    # Get healpix map of Nr of visits
    nside = region.tt_coverage_hp.meta["NSIDE"]
    npix = hpy.nside2npix(nside)
    hp_map = np.zeros(npix, dtype=np.float64)
    pix_ids = region.tt_coverage_hp["pix_id"].data.astype(np.int64)
    hp_map[pix_ids] = region.tt_coverage_hp["nr_vis"]

    # Plot background map
    plt_gw_kwargs = {
        "title": "Nr. of visits",
        "coord": "C",
        "reso": 10 / 60.0,
        "xsize": 1400,
        "ysize": 1400,
        "cmap": "gray",
    }
    if gw_kwargs is not None:
        plt_gw_kwargs.update(gw_kwargs)
    hpy.gnomview(hp_map, rot=[ra, dec], **plt_gw_kwargs)

    # Plots selected sources
    plt_ps_kwargs = {"lonlat": True, "marker": "o", "s": 0.2}
    if ps_kwargs is not None:
        plt_ps_kwargs.update(ps_kwargs)

    tt_srcs = region.tt_sources.group_by("rg_fd_id")
    for tt in tt_srcs.groups:
        sel = tt["sel"]
        if not sel_srcs:
            sel = np.ones(len(sel)).astype(bool)
        hpy.projscatter(tt[sel]["ra"], tt[sel]["dec"], **plt_ps_kwargs)

    # hpy.projscatter(
    #    [ra], [dec], lonlat=True, marker="o", s=4.0
    # )  # Mark center of gnomeview
    return plt.gca()


def plot_region_sky_mollview(region, var="nr_vis", mw_kwargs=None):
    """
    Plot the nr visits, fields or exposure on the sky. Coverage table has to be
    previously created with region.add_coverage_hp.

    Parameters
    ----------
    region : vasca.region.Region
        Region for which is plotted.
    coord : str
        Coordinate system, Galactic or ICKS
    var: str
        Variable to plot, exposure "exp", visits "nr_vis" or or fields "nr_fds"
    mw_kwargs : dict, optional
        Keyword arguments for healpy.mollview. The default is None.

    Returns
    -------
    matplotlib.axes
        Used Matplotlib axes.

    """
    if not hasattr(region, "tt_coverage_hp"):
        logger.warning(
            "Create healpix coverage table first with region.add_coverage_hp. "
            "Creating it now with default values."
        )
        region.add_coverage_hp()

    # Get healpix map of Nr of visits
    nside = region.tt_coverage_hp.meta["NSIDE"]
    npix = hpy.nside2npix(nside)
    hp_map = np.zeros(npix, dtype=np.float32)
    pix_ids = region.tt_coverage_hp["pix_id"].data
    hp_map[pix_ids] = region.tt_coverage_hp[var]

    coord_sys = str(region.tt_coverage_hp.meta["COOR_SYS"])
    dd_title = {
        "nr_vis": "Number of vists in " + coord_sys + " coord.",
        "nr_fds": "Number of fields in " + coord_sys + " coord.",
        "exp": "Exposure in " + coord_sys + " coord.",
    }

    # Plot background map
    plt_mw_kwargs = {
        "title": dd_title[var],
        "nest": False,
        "xsize": 4800,
        "cmap": "gist_stern",
    }
    if mw_kwargs is not None:
        plt_mw_kwargs.update(mw_kwargs)

    hpy.mollview(hp_map, **plt_mw_kwargs)

    return plt.gca()


# %% table variable plotting


def plot_table_hist(tt, var, ax=None, logx=False, obs_filter_id=None, **hist_kwargs):
    """
    Plot histogram for passed astropy.table.Table and variable

    Parameters
    ----------
    tt : astropy.table.Table
        Table containing a column with the plotted variable.
    var : str, optional
        Variable name
    ax : matplotlib.axes, optional
        Axes to draw on. The default is None.
    logx : bool, optional
        Histogram of log10(var) instead of var. The default is False.
    obs_filter_id: int, optional
        Observation filter ID Nr., if None all filters are shown. THe default is None.
    **hist_kwargs : dict
        Key word arguments passed tu plt.hist

    Returns
    -------
    matplotlib.axes
        Axes that where used to draw.
    list of float
        Histogram bin values.
    list of float
        Histogram bins, default is [selected, all] events.

    """
    logger.debug(f"Plotting histogram of variable '{var}'")

    if ax is None:
        ax = plt.gca()

    # Set marker properties for sources
    plot_kwargs = {
        "bins": "auto",
        "histtype": "step",
        "density": True,
        "log": True,
        "alpha": 0.5,
    }
    if hist_kwargs is not None:
        plot_kwargs.update(hist_kwargs)

    if obs_filter_id is not None:
        tt = select_obs_filter(tt, obs_filter_id)
    sel = tt["sel"]
    col = tt[var]

    str_nrsel = str(sel.sum())
    str_nrnotsel = str((~sel).sum())
    data = [np.array(col[~sel]).flatten(), np.array(col[sel]).flatten()]
    xlabel = var + " [" + str(col.unit) + "]"
    if str(col.unit) == "None" or str(col.unit) == "":
        xlabel = var
    if logx:
        with np.errstate(divide="ignore", invalid="ignore"):
            data = [np.log10(data[0]), np.log10(data[1])]
        xlabel = "log10( " + xlabel + " )"

    vals, bins, *_ = ax.hist(
        data,
        label=["unselected_" + str_nrnotsel, "selected_" + str_nrsel],
        **plot_kwargs,
    )

    if obs_filter_id is not None:
        xlabel = xlabel + " - " + dd_id2filter[obs_filter_id]

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Counts")

    return ax, vals, bins


def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)

    # now determine nice limits by hand:
    binwidth = 0.25
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax / binwidth) + 1) * binwidth

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation="horizontal")


def plot_table_scatter(
    tt,
    varx,
    vary,
    ax=None,
    xlim=None,
    ylim=None,
    invert_xaxis=None,
    invert_yaxis=None,
    xscale="linear",
    yscale="linear",
    obs_filter_id=None,
    grp_var="sel",
    grp_vals=None,
    add_projection=False,
    **scatter_kwargs,
):
    """
    Plot scatter plot for passed astropy.table.Table and variables

    Parameters
    ----------
    tt : astropy.table.Table
        Table containing columns with the plotted variables.
    varx : str
        Variable name on X-axis
    vary : str
        Variable name on Y-axis
    ax : matplotlib.axes, optional
        Axes to draw on. The default is None.
    xlim : list, optional
        List with [xmin, xmax] axis value. Default is None.
    ylim : list, optional
        List with [ymin, ymax] axis value. Default is None.
    xscale : str, optional
        Type of x-scale ("log", "linear"). Default is "linear".
    yscale : str, optional
        Type of y-scale ("log", "linear"). Default is "linear".
    obs_filter_id: int, optional
        Observation filter ID Nr., if None all filters are shown. The default is None.
    grp_var: str, optional
        Group scatter plot by colors based on this table variable. The default is "sel".
        If None is passed no groups will be done.
    add_projection: bool, optional
        Add histogram with projection on the sides
    **plot_kwargs : dict
        Key word arguments passed tu plt.plot

    Returns
    -------
    ax : matplotlix.axes
        Axes that where used to draw.

    """
    logger.debug(f"Plotting of variables '{varx}' and '{vary}'")

    if ax is None:
        ax = plt.gca()

    if add_projection:
        ax_histx = ax.inset_axes([0, 1.05, 1, 0.25], sharex=ax)
        ax_histy = ax.inset_axes([1.05, 0, 0.25, 1], sharey=ax)
        ax_histx.tick_params(axis="x", labelbottom=False)
        ax_histy.tick_params(axis="y", labelleft=False)

    if obs_filter_id is not None:
        tt = select_obs_filter(tt, obs_filter_id)

    # Set marker properties for sources
    plot_kwargs = {
        "markersize": 2.0,
        "linewidth": 0,
        "marker": "o",
        "markeredgewidth": 0,
    }
    if scatter_kwargs is not None:
        plot_kwargs.update(scatter_kwargs)

    if type(grp_var) == type(None):
        ax.plot(
            tt[varx],
            tt[vary],
            **plot_kwargs,
        )
    else:
        tt_grp = tt.group_by(grp_var)
        if type(grp_vals) == type(None):
            grp_vals = tt_grp.groups.keys[grp_var]
        for grp in grp_vals:
            mask = (
                tt_grp.groups.keys[grp_var] == grp
            )  # obs_by_name.groups.keys['name'] == 'M101'
            t_grp = tt_grp.groups[mask]
            pp = ax.plot(
                t_grp[varx],
                t_grp[vary],
                label=grp,
                **plot_kwargs,
            )
            if add_projection:
                ccol = pp[0].get_color()
                binsx = 10 ** np.linspace(np.log10(xlim[0]), np.log10(xlim[1]), 30)
                ax_histx.hist(
                    t_grp[varx], bins=binsx, log=True, histtype="step", color=ccol
                )  #
                binsy = 10 ** np.linspace(np.log10(ylim[0]), np.log10(ylim[1]), 30)
                ax_histy.hist(
                    t_grp[vary],
                    bins=binsy,
                    orientation="horizontal",
                    log=True,
                    histtype="step",
                    color=ccol,
                )

    # Set labels
    xlabel = varx + " [" + str(tt[varx].unit) + "]"
    if str(tt[varx].unit) == "None" or str(tt[varx].unit) == "":
        xlabel = varx
    ylabel = vary + " [" + str(tt[vary].unit) + "]"
    if str(tt[vary].unit) == "None" or str(tt[vary].unit) == "":
        ylabel = vary

    if obs_filter_id is not None:
        xlabel = xlabel + " - " + dd_id2filter[obs_filter_id]
        ylabel = ylabel + " - " + dd_id2filter[obs_filter_id]

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Set axis limits
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    if invert_xaxis:
        ax.invert_xaxis()
    if invert_yaxis:
        ax.invert_yaxis()

    return ax


# %% pipeline diagnostic plotting


# TODO: there are some GALEX specific variables
# will need to be adapted to other missions
def plot_pipe_diagnostic(
    tc, table_name, plot_type, fig_size=(12, 8), obs_filter_id=None
):
    """
    Diagnotic plots for VASCA pipe

    Parameters
    ----------
    tc : vasca.TableCollection
        Table collection, either field or region
    table_name : str
        Name of the table, tt_detections or tt_sources
    plot_type : str
        Type of plot, either "hist" or "scatter"
    fig_size: (float,float)
        Matplotlib figure (x,y) size in inches.
    obs_filter_id: int, optional
        Observation filter ID. Default is None.

    Returns
    -------
    matplotlib.figure
        Figure used for plotting
    dict
        Dictionary with plot variables keys and plot settings values

    """

    flt_name = ""
    if obs_filter_id is not None:
        flt_name = dd_id2filter[obs_filter_id]

    var_plt = OrderedDict()
    if plot_type == "hist":
        # Detections diagnostic
        if table_name == "tt_detections":
            var_plt["s2n"] = {"logx": True}
            var_plt["flux"] = {"logx": True}
            var_plt["r_fov"] = {"range": [0.0, 0.7]}
            var_plt["class_star"] = {}
            var_plt["artifacts"] = {"histtype": "step"}
            var_plt["chkobj_type"] = {}
            var_plt["ellip_world"] = {}  # "range": [0.0, 5]
            var_plt["size_world"] = {}  # "logx": True
            var_plt["flux_app_ratio"] = {"range": [0.1, 5]}  # "logx": True
            fig, axs = plt.subplots(
                3,
                3,
                figsize=(9, 9),
                squeeze=False,
                num="Detections Histograms " + flt_name,
            )
        elif table_name == "tt_sources":
            var_plt["nr_det"] = {}
            var_plt["flux_cpval"] = {"range": [0, 1]}
            var_plt["flux_nxv"] = {"range": [-0.01, 0.01]}  # "logx": True
            var_plt["coadd_fdiff_s2n"] = {"range": [-10, 25]}
            var_plt["coadd_ffactor"] = {"range": [-2, 5]}
            var_plt["pos_cpval"] = {}
            var_plt["flux"] = {"logx": True}
            var_plt["nr_fd_srcs"] = {}
            fig, axs = plt.subplots(
                2,
                4,
                figsize=(12, 6),
                squeeze=False,
                num="Sources Histograms " + flt_name,
            )
        elif table_name == "tt_coadd_sources":
            var_plt["nr_det"] = {}
            var_plt["flux_cpval"] = {}
            var_plt["pos_cpval"] = {}
            var_plt["flux"] = {"logx": True}
            fig, axs = plt.subplots(
                1,
                4,
                figsize=(12, 3),
                squeeze=False,
                num="Co-add Sources Histograms " + flt_name,
            )
        else:
            logger.warning("Diagnostic for table '{table_name}' not defined")

    elif plot_type == "scatter":
        # Detections diagnostic
        if table_name == "tt_detections":
            var_plt[("flux_app_ratio", "flux")] = {
                "yscale": "log",
                "xlim": [0, 5],
            }
            var_plt[("artifacts", "flux")] = {
                "yscale": "log",
            }
            var_plt[("class_star", "flux")] = {"yscale": "log"}
            var_plt[("size_world", "flux")] = {"yscale": "log"}  # {"yscale": "log"}
            var_plt[("ellip_world", "flux")] = {"yscale": "log"}
            var_plt[("class_star", "flux_app_ratio")] = {"ylim": [0, 5]}
            fig, axs = plt.subplots(
                2,
                3,
                figsize=(9, 6),
                squeeze=False,
                num="Detections Scatter " + flt_name,
            )
        elif table_name == "tt_sources":
            var_plt[("flux_cpval", "flux")] = {
                "xscale": "log",
                "xlim": [1e-23, 1.0],
                "yscale": "log",
            }
            var_plt[("pos_cpval", "flux")] = {
                "yscale": "log",
            }
            var_plt[("coadd_ffactor", "flux")] = {
                "yscale": "log",
                "xlim": [0.1, 100],
                "xscale": "log",
            }
            var_plt[("nr_det", "flux")] = {
                "yscale": "log",
            }
            var_plt[("coadd_fdiff_s2n", "coadd_ffactor")] = {
                "xlim": [-10, 25],
                "yscale": "log",
                "ylim": [0.1, 100],
            }
            var_plt[("flux_cpval", "nr_det")] = {
                "xscale": "log",
                "xlim": [1e-23, 1.0],
            }

            fig, axs = plt.subplots(
                2,
                3,
                figsize=(9, 6),
                squeeze=False,
                num="Sources Scatter " + flt_name,
            )
        elif table_name == "tt_coadd_sources":
            var_plt[("flux_cpval", "flux")] = {
                "xscale": "log",
                "yscale": "log",
            }
            var_plt[("pos_cpval", "flux")] = {
                "yscale": "log",
            }
            var_plt[("nr_det", "flux")] = {
                "yscale": "log",
            }

            fig, axs = plt.subplots(
                1,
                3,
                figsize=(9, 3),
                squeeze=False,
                num="Co-add Sources Scatter " + flt_name,
            )
        else:
            logger.warning("Diegnostic for table '{table_name}' not defined")
    else:
        logger.warning("Plot type '{plot_type}' unknown")

    # Select filter
    tt = tc.__dict__[table_name]
    if obs_filter_id is not None:
        tt = select_obs_filter(tt, obs_filter_id)

    # Start plotting
    axs = axs.flatten()
    ax_ctr = 0
    for var, plot_arg in var_plt.items():
        if plot_type == "hist":
            plot_table_hist(tt, var, axs[ax_ctr], **plot_arg)
        elif plot_type == "scatter":
            plot_table_scatter(tt, var[0], var[1], axs[ax_ctr], **plot_arg)
        ax_ctr += 1

    plt.tight_layout()
    plt.legend()
    return fig, var_plt


# %% light curve plotting


def plot_light_curves(
    tc,
    fd_src_ids=None,
    rg_src_ids=None,
    fig=None,
    ax=None,
    ylim=None,
    plot_upper_limits=True,
    flux_var="flux",
    **errorbar_kwargs,
):
    """
    Plot the light curves of the passed sources.

    Parameters
    ----------
    tc : VASCA.table.TableCollection
        Either field, region or source that contains the light curves.
    fd_src_ids : list or int
        List or single field source IDs to plot. Default is None.
    rg_src_ids : list or int
        List or single region source IDs to plot. Default is None.
    fig: figure, optional
        Matplotlib figure to draw on, if None a new figure is created. The default is None.
    ax : axes, optional
        Matplotlib axes to plot on. The default is None.
    ylim : list, optional
        Limits of the y axis. Default is None
    plot_upper_limits : bool
        Plot upper limits to the lightcurve. The default is True.
    flux_var: str, optional
        Variable in table to be used to get flux Jy
    **errorbar_kwargs : dict
        Key word arguments for pyplot.errorbars plotting.

    Returns
    -------
    matplotlib.figure
        Matplotlib figure used to draw
    matplotlib.axes
        Used Matplotlib axes.

    """

    logger.debug("Plotting lightcurves ")

    # Check if figure was passed
    if type(fig) is type(None) and type(ax) is type(None):
        fig = plt.figure(figsize=(6, 6))  # , constrained_layout=True
    else:
        fig = plt.gcf()

    # Check if axis was passed
    if ax is None:
        ax = plt.gca()

    # ax.invert_yaxis()
    ax.set_yscale("log")

    if hasattr(ylim, "__iter__"):
        ax.set_ylim(ylim)

    # Setup plotting parameters
    plt_errorbar_kwargs = {
        "markersize": 4,
        "alpha": 0.6,
        "capsize": 0,
        "lw": 0.1,
        "linestyle": "dotted",
        "elinewidth": 0.7,
    }
    if errorbar_kwargs is not None:
        plt_errorbar_kwargs.update(errorbar_kwargs)

    # Loop over selected sources and plot
    colors = cycle("bgrcmykbgrcmykbgrcmykbgrcmyk")
    markers = cycle("osDd<>^v")
    ctr = 0

    # Get light curves dictionary
    dd_lcs = tc.get_light_curve(fd_src_ids, rg_src_ids, flux_var=flux_var)
    src_ids = list(dd_lcs.keys())

    for src_id, col, mar in zip(src_ids, colors, markers):
        ctr += 1

        # Arrays to plot filters differently
        mfc = ["None", col]
        ls = ["--", ":"]

        # Get light curve
        lc = dd_lcs[src_id]
        fluxs = lc["flux"]
        fluxs_err = lc["flux_err"]

        # Check if one one filter present
        filter_ids = np.sort(np.unique(lc["obs_filter_id"].data))

        for flt_id in filter_ids:
            flt_plot = flt_id % 2
            #  Select filter to show
            sel = (lc["flux"] > 0) * (lc["obs_filter_id"] == flt_id)
            src_lab = str(src_id) + " " + str(lc[sel]["obs_filter"][0])

            # Draw mean value
            t_mean = [np.min(lc["time"][sel]), np.max(lc["time"][sel])]
            flux_weight = 1.0 / lc["flux_err"][sel] ** 2
            flux_mean = np.average(lc["flux"][sel], weights=flux_weight)
            ax.plot(
                t_mean,
                [flux_mean, flux_mean],
                ls=ls[flt_plot],
                color=col,
                linewidth=0.5,
            )

            # Plot
            ax.errorbar(
                lc["time"][sel],
                fluxs[sel],
                yerr=fluxs_err[sel],
                color=col,
                markeredgecolor=col,
                markerfacecolor=mfc[flt_plot],
                marker=mar,
                label=src_lab,
                **plt_errorbar_kwargs,
            )

    ax.legend(fontsize="small")  # bbox_to_anchor=(1.04, 1),
    ax.set_xlabel("MJD")
    ax.set_ylabel(r"Flux [$\mu$Jy]")

    # Add a second time axis on top showing years
    def mjd2yr(mjd):
        return Time(mjd, format="mjd").jyear

    def yr2mjd(jyr):
        return Time(jyr, format="jyear").mjd

    secax = ax.secondary_xaxis("top", functions=(mjd2yr, yr2mjd))
    secax.set_xlabel("Year")

    secay = ax.secondary_yaxis("right", functions=(flux2mag_np, mag2flux_np))

    # Avoid scientific notation for magnitudes

    # formatter.set_scientific(False)
    # secay.yaxis.set_minor_formatter(formatter)

    secay.set_ylabel("AB magnitude")

    return fig, ax


def plot_light_curve(
    tc_src, fig=None, ax=None, show_gphoton=True, add_axes=True, **errorbar_kwargs
):
    """
    Plots light curve

    Parameters
    ----------
    tc_src: vasca.TableCollection
        Table collection containing tt_sed table.
    fig: figure, optional
        Matplotlib figure to draw on, if None a new figure is created. The default is None.
    ax : axes, optional
        Matplotlib axes to plot on. The default is None.
    show_gphoton: bool, optional
        Show gphoton light curve too, if present in table collection?. The default is True.
    add_axes: bool, optional
        Add additional x- and y-axis on the empty sides
    **errorbar_kwargs : dict
        Key word arguments for pyplot.errorbars plotting.

    Returns
    -------
    matplotlib.figure
        Matplotlib figure used to draw
    matplotlib.axes
        Used Matplotlib axes.
    """

    logger.debug("Plotting spectral energy distribution ")

    # Check if lightcurve table exists
    if "tt_source_lc" not in tc_src._table_names:
        logger.warning("No light curve table found")
        return
    # Consider only selected points
    else:
        tt_lc = tc_src.tt_source_lc[tc_src.tt_source_lc["sel"]]

    # Get also gPhoton light curve if it excists
    if "tt_gphoton_lc" in tc_src._table_names:
        tt_gp_lc = tc_src.tt_gphoton_lc[tc_src.tt_gphoton_lc["sel"]]

    # Check if figure was passed
    if type(fig) is type(None) and type(ax) is type(None):
        fig = plt.figure(figsize=(12, 6))  # , constrained_layout=True
    else:
        plt.gcf()

    # Check if axis was passed
    if ax is None:
        ax = plt.gca()

    ax.set_yscale("log")

    # Setup plotting parameters
    plt_errorbar_kwargs = {
        "markersize": 6,
        "alpha": 0.6,
        "capsize": 0,
        "lw": 0.2,
        "linestyle": "None",  # "dotted",
        "elinewidth": 0.7,
    }
    if errorbar_kwargs is not None:
        plt_errorbar_kwargs.update(errorbar_kwargs)

    filter_ids = np.sort(np.unique(tt_lc["obs_filter_id"].data))

    colors = cycle("brgcmykbgrcmykbgrcmykbgrcmyk")
    markers = cycle("osDd<>^v")

    for flt_id, mar, col in zip(filter_ids, markers, colors):
        # Plot gPhoton light curve first, if present
        if show_gphoton and "tt_gphoton_lc" in tc_src._table_names:
            sel_gp = tt_gp_lc["obs_filter_id"] == flt_id
            if sel_gp.sum() > 0:
                ax.errorbar(
                    tt_gp_lc["time"][sel_gp],
                    tt_gp_lc["flux"][sel_gp],
                    yerr=tt_gp_lc["flux_err"][sel_gp],
                    marker=mar,
                    color="0.5",
                    **plt_errorbar_kwargs,
                )

        #  Select filter to show
        sel = tt_lc["obs_filter_id"] == flt_id
        src_lab = str(tt_lc[sel]["obs_filter"][0])

        # Plot
        if sel.sum() > 0:
            t_mean = [np.min(tt_lc["time"][sel]), np.max(tt_lc["time"][sel])]
            flux_weight = 1.0 / tt_lc["flux_err"][sel] ** 2
            flux_mean = np.average(tt_lc["flux"][sel], weights=flux_weight)
            ax.plot(
                t_mean,
                [flux_mean, flux_mean],
                ls="--",
                color=col,
                linewidth=0.5,
            )

            ax.errorbar(
                tt_lc["time"][sel],
                tt_lc["flux"][sel],
                yerr=tt_lc["flux_err"][sel],
                marker=mar,
                label=src_lab,
                color=col,
                **plt_errorbar_kwargs,
            )

    ax.legend(fontsize=12, frameon=False)  # bbox_to_anchor=(1.04, 1),
    ax.set_xlabel("MJD", fontsize=18)
    ax.set_ylabel(r"Flux [$\mu$Jy]", fontsize=18)
    ax.text(
        0.3,
        0.02,
        tc_src.tt_sources["src_name"][0],
        size=22,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
    )  # , color='purple'

    # Avoid scientific notation
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    ax.yaxis.set_minor_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)

    if add_axes:
        # Add a second time axis on top showing years

        secax = ax.secondary_xaxis("top", functions=(mjd2yr, yr2mjd))
        secax.ticklabel_format(useOffset=False, style="plain")
        secax.set_xlabel("Year", fontsize=18)

        secay = ax.secondary_yaxis("right", functions=(flux2mag_np, mag2flux_np))

        # Avoid scientific notation for magnitudes
        secay.yaxis.set_minor_formatter(formatter)

        secay.set_ylabel("AB magnitude", fontsize=18)

    return fig, ax


# %% LombScargle
def plot_lombscargle(
    tt_lc,
    fig=None,
    ax=None,
    ax_phase=None,
    ax_lc=None,
    obs_filter="NUV",
    nbins_min=10,
    logy=False,
    freq_range=[0.03, 2] / uu.d,
    plot_dtbins=True,
):
    """
    Runs and plots Lomb Scargle diagram

    Parameters
    ----------
    tc_src: vasca.TableCollection
        Table collection containing tt_sed table.
    fig: figure, optional
        Matplotlib figure to draw on, if None a new figure is created. The default is None.
    ax : axes, optional
        Matplotlib axes to plot LombScargle diagram on. The default is None.
    ax_lc : axes, optional
        Matplotlib axes to plot with the light curve to plot peak frequency model
        on. The default is None.
    ax_phase : axes, optional
        Matplotlib axes to plot phase diagram. The default is None.
    obs_filter : str, optional
        Observational filter to perform LombScargle on. The default is "NUV"
    nbins_min : int, optional
        Minimum number of time bins to perform LombScargle. The default is 20.
    logy : bool, optional
        Plot LombScargle diagram in log(Power). The default is True.
    freq_range : list
        Minimum and maximum Frequency. If None calculated automatically.

    Returns
    -------
    matplotlib.figure
        Matplotlib figure used to draw
    matplotlib.axes
        Used Matplotlib axes.
    dict
        Lomb Scargle results
    """

    logger.debug("Plotting Lomb-Scargle diagram ")
    label_size = 18

    # Prepare light curve
    sel = np.array(
        (tt_lc["obs_filter"] == obs_filter) * (tt_lc["sel"] == True), dtype=bool
    )
    if sel.sum() < nbins_min:
        return
    tt_lc = tt_lc[sel]
    tt_lc.sort("time")

    # Check if figure was passed
    if type(fig) is type(None) and type(ax) is type(None):
        fig = plt.figure(figsize=(6, 6))  # , constrained_layout=True
    else:
        plt.gcf()

    # Check if axis was passed
    if ax is None:
        ax = plt.gca()

    # Setup plotting parameters
    plt_plot_kwargs = {"alpha": 0.5}

    dd_ls_results = run_LombScargle(tt_lc, nbins_min=nbins_min, freq_range=freq_range)

    # Set labels
    ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel("Frequency [1/day]", fontsize=label_size)
    ax.set_ylabel("LS power", fontsize=label_size)
    secax = ax.secondary_xaxis("top", functions=(freq2period, period2freq))
    secax.set_xlabel("Period [day]", fontsize=label_size)

    # Plot frequency distribution of time bins
    if plot_dtbins:
        fbinsize = dd_ls_results["ls_freq"][1:] - dd_ls_results["ls_freq"][:-1]
        fedges = dd_ls_results["ls_freq"][:-1] + fbinsize / 2.0
        dtimes = tt_lc["time"][1:].data[None, :] - tt_lc["time"][:-1].data[:, None]
        dfs = 1 / dtimes.flatten()
        hist_f, _ = np.histogram(dfs, bins=fedges.data)
        hist_f = 0.5 * hist_f / np.max(hist_f)  # normalize maximum to 0.5
        ax.plot(
            dd_ls_results["ls_freq"][1:-1],
            hist_f,
            label="Frequency between visits",
            **plt_plot_kwargs,
        )

    # Plot LS
    pl = ax.plot(
        dd_ls_results["ls_freq"],
        dd_ls_results["ls_power"],
        **plt_plot_kwargs,
    )  # label="data",
    col = pl[0].get_color()

    probabilities = [
        0.95449973610364,
        0.002699796063,
        0.000063342484,
        0.000000573303,
    ]  # 2-5 sigma

    # Get confidence interval
    conf = dd_ls_results["ls"].false_alarm_level(
        probabilities, method="baluev"
    )  # "bootstrap"
    ax.axhline(
        conf[1],
        linewidth=0.5,
        ls="--",
        color=col,
        label=r"$3 \sigma$ - $5 \sigma$ confidence levels",
    )  # col
    ax.axhline(conf[2], linewidth=0.5, ls="--", color=col)  # col
    ax.axhline(conf[3], linewidth=0.5, ls="--", color=col)

    # ax.legend()

    # Plot model on lc
    if type(ax_lc) != type(None):
        t_min = np.min(tt_lc["time"].quantity)
        t_max = np.max(tt_lc["time"].quantity)
        nbins = int(
            16 * (t_max.value - t_min.value) * dd_ls_results["ls_peak_freq"].value
        )
        t_fit = np.linspace(t_min, t_max, nbins)
        flux_fit = dd_ls_results["ls"].model(t_fit, dd_ls_results["ls_peak_freq"])
        ax_lc.plot(t_fit, flux_fit, alpha=0.5)

    # Plot phase diagram
    if type(ax_phase) != type(None):
        period_peak = float(1 / dd_ls_results["ls_peak_freq"].value)
        # period_peak = 220.44/(24*60*60)
        times_phased = tt_lc["time"] % period_peak
        t_fit = np.linspace(0, period_peak, 40)
        flux_fit = dd_ls_results["ls"].model(
            t_fit * uu.d, dd_ls_results["ls_peak_freq"]
        )
        ax_phase.errorbar(
            times_phased / period_peak,
            tt_lc["flux"],
            yerr=tt_lc["flux_err"],
            linestyle="none",
            marker="o",
            color=col,
        )

        ax_phase.plot(t_fit / period_peak, flux_fit, **plt_plot_kwargs)
        ax_phase.set_xlabel("Phase ", fontsize=label_size)
        ax_phase.set_ylabel("Flux [Jy]", fontsize=label_size)
    return fig, ax, dd_ls_results


def plot_sed(
    tc_src, fig=None, ax=None, plot_spec_lines=False, plot_spec=False, **errorbar_kwargs
):
    """
    Plots spectral energy distribution

    Parameters
    ----------
    tc_src: vasca.TableCollection
        Table collection containing tt_sed table.
    fig: figure, optional
        Matplotlib figure to draw on, if None a new figure is created. The default is None.
    ax : axes, optional
        Matplotlib axes to plot on. The default is None.
    plot_spec_lines: bool
        Plot typical spectral lines in White Dwarfs. The default is False.
    plot_spec: bool
        Plot typical spectral lines in White Dwarfs. The default is False.
    **errorbar_kwargs : dict
        Key word arguments for pyplot.errorbars plotting.

    Returns
    -------
    matplotlib.figure
        Matplotlib figure used to draw
    matplotlib.axes
        Used Matplotlib axes.
    dict
        Dictionary with fit information
    """

    logger.debug("Plotting spectral energy distribution ")

    # Check if SED table exists
    if "tt_sed" not in tc_src._table_names:
        logger.warning("No SED table found")
        return
    else:
        tt_sed = tc_src.tt_sed

    # Remove SkyMapper data, as it seems off
    sel_obs = np.ones(len(tt_sed), dtype=bool)
    for obs in ["SkyMapper/SkyMapper", "GAIA/GAIA2", "Gaia", "GALEX", "SDSS"]:
        sel_obs = sel_obs * ~(tt_sed["observatory"] == obs)
    sel_obs = sel_obs + (tt_sed["origin"] == "VASCA")
    sel_obs = sel_obs + (tt_sed["origin"] == "J/MNRAS/486/2169/table3")

    # sel_obs = sel_obs + (tt_sed["origin"] == "V/154/sdss16")
    tt_sed = tt_sed[sel_obs]

    # Check if figure was passed
    if type(fig) is type(None) and type(ax) is type(None):
        fig = plt.figure(figsize=(6, 6))  # , constrained_layout=True
    else:
        plt.gcf()

    # Check if axis was passed
    if ax is None:
        ax = plt.gca()

    ax.set_xscale("log")
    ax.set_yscale("log")

    # Setup plotting parameters
    plt_errorbar_kwargs = {
        "capsize": 0,
        "lw": 0.1,
        "linestyle": "dotted",
        "elinewidth": 0.7,
    }
    if errorbar_kwargs is not None:
        plt_errorbar_kwargs.update(errorbar_kwargs)

    # Loop over selected sources and plot
    colors = cycle("bgcmykbgcmykbgcmykbgcmyk")
    markers = cycle("sDd<>^v")

    # Prepare tables, to plot vasca point separatelly
    sel = tt_sed["origin"] == "VASCA"

    # If ony VASCA point present, do not plot SED
    if sel.sum() == len(tt_sed):
        print("Only VASCA points in SED, not plotting")
        return
    tt_grp = tt_sed[~sel].group_by("observatory")

    # **********  Plot spectral lines in background first
    if "tt_spectrum_0" in tc_src._table_names and plot_spec_lines:
        for lwavs, llabel, lls, lcol in zip(
            dd_spec_lines.values(),
            dd_spec_lines.keys(),
            [":", "--", "-"],
            ["0.7", "0.8", "0.9"],
        ):
            valabel = llabel
            for lwav in lwavs:
                ax.axvline(
                    lwav.value,
                    0.3,
                    0.9,
                    linestyle=lls,
                    color=lcol,
                    linewidth=1.0,
                    label=valabel,
                )
                valabel = None

    # ********** Plot spectra, if present. Up to 5
    if plot_spec:
        for ii, col in zip(range(0, 5), colors):
            spec_name = "tt_spectrum_" + str(ii)
            lines_draw = True
            if spec_name in tc_src._table_names:
                tt_spec_ii = tc_src.__dict__[spec_name]
                sel_spec = tt_spec_ii["sel"]
                if sel_spec.sum() > 0:
                    ax.plot(
                        tt_spec_ii["wavelength"][sel_spec],
                        tt_spec_ii["flux"][sel_spec],
                        label="SDSS spectrum",  # + str(ii),
                        alpha=0.5,
                        color=col,
                    )
                    ax.plot(
                        tt_spec_ii["wavelength"][sel_spec],
                        tt_spec_ii["flux_model"][sel_spec],
                        label="SDSS model",  # + str(ii)
                        # alpha=0.5,
                        color=col,
                    )
            else:
                break

    # *********** Plot all none-VASCA points
    for tt, grp, col, mar in zip(tt_grp.groups, tt_grp.groups.keys, colors, markers):

        # Plot
        ax.errorbar(
            tt["wavelength"],
            tt["flux"],
            yerr=tt["flux_err"],
            color=col,
            markeredgecolor=col,
            marker=mar,
            label=str(grp[0]),
            markersize=6,
            alpha=0.4,
            **plt_errorbar_kwargs,
        )

    # ************* Plot VASCA points
    ax.errorbar(
        tt_sed[sel]["wavelength"],
        tt_sed[sel]["flux"],
        yerr=tt_sed[sel]["flux_err"],
        color="r",
        markeredgecolor="r",
        marker="o",
        markersize=6,
        alpha=1.0,
        label="GALEX/VASCA",
        **plt_errorbar_kwargs,
    )

    # ************* Plot Black Body
    # Fit Black Body spectrum
    uscale = uu.Unit("1e-6 Jy/sr")
    BB = models.BlackBody(temperature=2e4 * uu.K, scale=1e-24 * uscale)  #
    # fit = fitting.LMLSQFitter()
    fit = fitting.LevMarLSQFitter()

    # Restrict fit range of BB fit to reasonable range
    bb_flux = (tt_sed["flux"].quantity / uu.Unit("sr")).to(uscale)
    selfit = (
        (tt_sed["wavelength"] > 0)
        * (tt_sed["wavelength"] < 1e20)
        * (bb_flux.value > 1e-10)
        * (bb_flux.value < 1e10)
    )

    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.filterwarnings(
            "ignore",
            message="Model is linear in parameters",
            category=AstropyUserWarning,
        )

        bb_bins = (
            np.linspace(
                np.min(tt_sed["wavelength"][selfit]),
                np.max(tt_sed["wavelength"][selfit]),
                100,
            )
            * tt_sed["wavelength"].unit
        )
        fitted_bb = fit(
            BB,
            tt_sed["wavelength"].quantity[selfit],
            bb_flux[selfit],
            weights=None,  # 1 / tt_sed["flux_err"][selfit]
            maxiter=100,
            acc=1e-07,
            epsilon=1.4901161193847656e-08,
            estimate_jacobian=False,
            filter_non_finite=True,
        )
        if fit.fit_info["ierr"] <= 4 and fit.fit_info["ierr"] >= 1:
            bb_bins = (
                np.linspace(
                    np.min(tt_sed["wavelength"][selfit]),
                    np.max(tt_sed["wavelength"][selfit]),
                    100,
                )
                * tt_sed["wavelength"].unit
            )
            fit_flux = (fitted_bb(bb_bins) * uu.Unit("sr")).to(uu.Unit("1e-6 Jy"))
            sel_fit = tt_sed

            fit_temp = (
                np.round(fitted_bb.temperature.value, 0) * fitted_bb.temperature.unit
            )
            ax.plot(
                bb_bins,
                fit_flux,
                color="0.5",
                ls="-",
                # label="BB " + str(fit_temp),
            )
            print("BB fit temperature:", fit_temp)
        else:
            print("Black body fit did not converge.")

    # Plot star black body if distance is available

    if "tt_gaiadr3" in tc_src._table_names:
        gaia_match_id = tc_src.tt_sources["gaiadr3_match_id"][0]
        tc_src.tt_gaiadr3.add_index("gaiadr3_match_id")
        gaia_idx = tc_src.tt_gaiadr3.loc_indices["gaiadr3_match_id", gaia_match_id]
        if tc_src.tt_gaiadr3["Plx_dist"][gaia_idx] > 1:
            TS = 2700 * uu.K
            rS = (0.1 * const.R_sun).to(uu.cm)
            dS = (tc_src.tt_gaiadr3["Plx_dist"].quantity[gaia_idx]).to(uu.cm)
            BBS_scale = np.pi * (rS / dS) ** 2

            BBS = models.BlackBody(temperature=TS, scale=BBS_scale)

            wavS = np.linspace(7000, 110000, 100) * uu.AA
            fluxS = (BBS(wavS) * uu.Unit("sr")).to(uu.Unit("1e-6 Jy"))
            ax.plot(
                wavS,
                fluxS,
                color="0.5",
                ls="--",
                # label="BB star " + str(TS),
            )

    # Helper functions to define second axis

    def AA2ev_np(wave):
        return (cc.h * cc.c / (wave * uu.AA)).to(uu.eV).value

    def ev2AA_np(ener):
        return (cc.h * cc.c / (ener * uu.eV)).to(uu.AA).value

    # In units of Jy the Wien-displacement constant is given by the factor below
    # https://de.wikipedia.org/wiki/Wiensches_Verschiebungsgesetz
    w_kb_ev = cc.k_B.to(uu.eV / uu.K).value * 2.82

    def AA2K_np(wave):
        return AA2ev_np(wave) / w_kb_ev

    def K2AA_np(temp):
        return ev2AA_np(temp * w_kb_ev)

    # Axis and labels
    # secax = ax.secondary_xaxis("top", functions=(AA2ev_np, ev2AA_np))
    # secax.set_xlabel("eV")

    secax = ax.secondary_xaxis("top", functions=(AA2K_np, K2AA_np))
    secax.set_xlabel("Temperature [K]", fontsize=16)

    secay = ax.secondary_yaxis("right", functions=(flux2mag_np, mag2flux_np))

    # Avoid scientific notation for magnitudes
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    secay.yaxis.set_minor_formatter(formatter)
    secay.set_ylabel("AB magnitude", fontsize=16)

    ax.set_ylabel("Flux [$\mu$Jy]", fontsize=16)
    ax.set_xlabel("Wavelength [Angstom]", fontsize=16)

    ax.legend(frameon=False, loc="upper right")

    return fig, ax, fit.fit_info
