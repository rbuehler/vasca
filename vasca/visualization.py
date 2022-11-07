#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualization related methods for VASCA
"""

from collections import OrderedDict
from itertools import cycle

import healpy as hpy
import matplotlib.pyplot as plt
import numpy as np
from loguru import logger
from matplotlib.colors import LogNorm

# %% sky plotting


def plot_field_sky_sources(
    field,
    ax=None,
    plot_detections=True,
    plot_fd_src_id=True,
    src_kwargs=None,
    det_kwargs=None,
):
    """
    Plot the selected sources and (optinally) the visit detections on the sky.

    Parameters
    ----------
    field: vasca.BaseField
        VASCA field to be plotted.
    ax : axes, optional
        Matplotlib axes to plot on. The default is None.
    plot_detections : bool, optional
        Plot the visit detections below the sources. The default is True.
    plot_fd_src_ids: bool, optional
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
        "markersize": 4.5,
        "alpha": 0.5,
        "lw": 0,
        "markeredgewidth": 2.0,
        "fillstyle": "none",
    }
    if src_kwargs is not None:
        plt_src_kwargs.update(src_kwargs)

    # Set marker properties for detections
    plt_det_kwargs = {
        "marker": ".",
        "markersize": 3.0,
        "alpha": 0.5,
        "markeredgewidth": 1.0,
        "lw": 0.0,
    }
    if det_kwargs is not None:
        plt_det_kwargs.update(det_kwargs)

    # Prepare data
    sel = field.tt_sources["sel"]
    tt_src = field.tt_sources[sel]
    tt_det = field.tt_detections
    tt_det.add_index("fd_src_id")

    # Loop over all srcs and plot
    colors = cycle("bgrcmykbgrcmykbgrcmykbgrcmyk")
    for src, col in zip(tt_src, colors):
        if plot_detections:
            det_idx = tt_det.loc_indices["fd_src_id", src["fd_src_id"]]
            ax.plot(
                tt_det[det_idx]["ra"].data,
                tt_det[det_idx]["dec"].data,
                color=col,
                **plt_det_kwargs,
            )
        ax.plot(src["ra"], src["dec"], color=col, **plt_src_kwargs)

        ax.text(
            src["ra"] + 0.007,
            src["dec"] + 0.006,
            str(src["fd_src_id"]),
            transform=plt_src_kwargs["transform"],
            fontsize=7,
            color=col,
            alpha=0.7,
        )

    return ax


def plot_field_sky_map(field, ax=None, **img_kwargs):
    """
    Plot the reference sky map.

    Parameters
    ----------
    field: vasca.BaseField
        VASCA field to be plotted.
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

    if field.ref_img is None:
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

    graph = ax.imshow(field.ref_img, **plt_img_kwargs)

    return graph


def plot_field_sky(field, plot_detections=True, plot_map=True):
    """
    Plot all field sources and/or a background reference image in the sky.

    Parameters
    ----------
    field: vasca.BaseField
        VASCA field to be plotted.
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

    if field.ref_wcs is not None:
        ax = plt.subplot(projection=field.ref_wcs)  #
        ax.coords["ra"].set_major_formatter("d.dd")
        ax.coords["dec"].set_major_formatter("d.dd")
        ax.set_xlabel("Ra")
        ax.set_ylabel("Dec")

    if plot_map:
        graph = plot_field_sky_map(field, ax)
        fig.colorbar(graph, label="Intensity [a.u.]")

    if plot_map:
        plot_arg = {"transform": ax.get_transform("world")}
        plot_field_sky_sources(
            field,
            ax,
            plot_detections=plot_detections,
            src_kwargs=plot_arg,
            det_kwargs=plot_arg,
        )
    else:
        plot_field_sky_sources(field, plot_detections=plot_detections)

    return fig


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
        Keyword arguments for healpy..gnomview. The default is None.
    ps_kwargs : dict, optional
        Keyword arguments for healpy.projscatter. The default is None.

    Returns
    -------
    ax : axes
        Used Matplotlib axes.

    """

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
        "reso": 5 / 60.0,
        "xsize": 2000,
        "ysize": 2000,
        "cmap": "gray",
    }
    if gw_kwargs is not None:
        plt_gw_kwargs.update(gw_kwargs)
    hpy.gnomview(hp_map, rot=[ra, dec], **plt_gw_kwargs)

    # Plots selected sources
    plt_ps_kwargs = {"lonlat": True, "marker": "o", "s": 0.2}
    if ps_kwargs is not None:
        plt_ps_kwargs.update(ps_kwargs)

    tt_srcs = region.tt_sources.group_by("field_id")
    for tt in tt_srcs.groups:
        sel = tt["sel"]
        if not sel_srcs:
            sel = np.ones(len(sel)).astype(bool)
        hpy.projscatter(tt[sel]["ra"], tt[sel]["dec"], **plt_ps_kwargs)

    # hpy.projscatter(
    #    [ra], [dec], lonlat=True, marker="o", s=4.0
    # )  # Mark center of gnomeview
    # hpy.mollview(hp_map, title="Nr. of visits in log10",nest=False,cmap="nipy_spectral",xsize=4800)
    # hpy.graticule(local=False,coord="C",dpar=1.0, color="white") # show graticules every 0.5 deg
    return plt.gca()


# %% table variable plotting


def plot_table_hist(tt, var, ax=None, logx=False, **hist_kwargs):
    """
    Plot histogram for passed astropy.Table and variable

    Parameters
    ----------
    tt : astropy.Table
        Table containing a column with the plotted variable.
    var : str, optional
        Variable name
    ax : matplotlib.axes, optional
        Axes to draw on. The default is None.
    logx : bool, optional
        Histoogram of log10(var) instead of var. The default is False.
    **hist_kwargs : dict
        Key word arguments passed tu plt.hist

    Returns
    -------
    ax : matplotlix.axes
        Axes that where used to draw.

    """
    logger.debug(f"Plotting histogram of variable '{var}'")

    if ax is None:
        ax = plt.gca()

    # Set marker properties for sources
    plot_kwargs = {
        "bins": "auto",
        "histtype": "stepfilled",
        "log": True,
        "alpha": 0.5,
    }
    if hist_kwargs is not None:
        plot_kwargs.update(hist_kwargs)

    col = tt[var]
    sel = tt["sel"]
    str_nrsel = str(sel.sum())
    str_nrnotsel = str((~sel).sum())
    data = [col[~sel], col[sel]]
    xlabel = var + " [" + str(col.unit) + "]"
    if str(col.unit) == "None" or str(col.unit) == "":
        xlabel = var
    if logx:
        with np.errstate(divide="ignore", invalid="ignore"):
            data = [np.log10(data[0]), np.log10(data[1])]
        xlabel = "log10( " + xlabel + " )"

    ax.hist(
        data,
        label=["unselected_" + str_nrnotsel, "selected_" + str_nrsel],
        **plot_kwargs,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Counts")

    return ax


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
    **scatter_kwargs,
):
    """
    Plot scatter plot for passed astropy.Table and variables

    Parameters
    ----------
    tt : astropy.Table
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
    **scatter_kwargs : dict
        Key word arguments passed tu plt.scatter

    Returns
    -------
    ax : matplotlix.axes
        Axes that where used to draw.

    """
    logger.debug(f"Plotting of variables '{varx}' and '{vary}'")

    if ax is None:
        ax = plt.gca()

    # Set marker properties for sources
    plot_kwargs = {"s": 2.0, "alpha": 0.5}
    if scatter_kwargs is not None:
        plot_kwargs.update(scatter_kwargs)

    sel = tt["sel"].astype(bool)
    str_nrsel = str(sel.sum())
    str_nrnotsel = str((~sel).sum())
    ax.scatter(
        tt[varx][~sel],
        tt[vary][~sel],
        label="unselected_" + str_nrnotsel,
        **plot_kwargs,
    )
    ax.scatter(
        tt[varx][sel], tt[vary][sel], label="selected_" + str_nrsel, **plot_kwargs
    )

    # Set labels
    xlabel = varx + " [" + str(tt[varx].unit) + "]"
    if str(tt[varx].unit) == "None" or str(tt[varx].unit) == "":
        xlabel = varx
    ylabel = vary + " [" + str(tt[vary].unit) + "]"
    if str(tt[vary].unit) == "None" or str(tt[vary].unit) == "":
        ylabel = vary
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


# TODO: there are some GALEX specific variables
# will need to be adapted to other missions
def plot_pipe_diagnostic(tc, table_name, plot_type):
    """
    Diagnotic plots for VASCA pipe

    Parameters
    ----------
    tc : vasca.TableCollection
        Table collection
    table_name : str
        Name of the table
    plot_type : str
        Type of plot, either "hist" or "scatter"

    Returns
    -------
    fig : matplotlib.figure
        Figure used for plotting

    """

    det_vars = OrderedDict()
    if plot_type == "hist":
        # Detections diagnostic
        if table_name == "tt_detections":
            det_vars["s2n"] = {"logx": True}
            det_vars["mag"] = {"range": [14.5, 26.5]}
            det_vars["mag_err"] = {"range": [0.0, 3.0]}
            det_vars["r_fov"] = {"range": [0.0, 0.7]}
            det_vars["point_src_prob"] = {}
            det_vars["artifacts"] = {"histtype": "step"}
            fig, axs = plt.subplots(3, 2, figsize=(18, 12), squeeze=False)
        elif table_name == "tt_sources":
            det_vars["nr_det"] = {}
            det_vars["nr_uls"] = {}
            det_vars["mag_mean"] = {}
            det_vars["mag_rchiq"] = {"logx": True, "range": [-3, 3]}
            det_vars["mag_dmax"] = {}
            det_vars["mag_dmax_sig"] = {"logx": True, "range": [-3, 2]}
            det_vars["mag_var"] = {"logx": True, "range": [-3, 0]}
            det_vars["ul_weight"] = {"bins": 100}
            fig, axs = plt.subplots(2, 4, figsize=(22, 12), squeeze=False)
        else:
            logger.warning("Diegnostic for table '{table_name}' not defined")

    elif plot_type == "scatter":
        # Detections diagnostic
        if table_name == "tt_detections":
            det_vars[("s2n", "mag")] = {
                "invert_yaxis": True,
                "xlim": [1, 100],
                "ylim": [14.5, 26.5],
            }
            det_vars[("r_fov", "mag")] = {"invert_yaxis": True, "ylim": [15.5, 26.5]}
            det_vars[("point_src_prob", "mag")] = {
                "invert_yaxis": True,
                "ylim": [14.5, 26.5],
            }
            det_vars[("artifacts", "mag")] = {
                "invert_yaxis": True,
                "ylim": [14.5, 26.5],
            }
            det_vars[("r_fov", "artifacts")] = {}
            det_vars[("mag_err", "mag")] = {
                "xlim": [0.01, 3],
                "invert_yaxis": True,
                "ylim": [14.5, 26.5],
            }
            fig, axs = plt.subplots(3, 2, figsize=(14, 12), squeeze=False)
        elif table_name == "tt_sources":
            det_vars[("mag_rchiq", "mag_dmax_sig")] = {"xscale": "log"}
            det_vars[("mag_rchiq", "ul_weight")] = {"xscale": "log"}
            det_vars[("mag_dmax_sig", "ul_weight")] = {}
            det_vars[("mag_dmax_sig", "mag_mean")] = {
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            det_vars[("ul_weight", "mag_mean")] = {
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            det_vars[("mag_rchiq", "mag_mean")] = {
                "xscale": "log",
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            det_vars[("mag_var", "mag_mean")] = {
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            det_vars[("nr_uls", "mag_mean")] = {
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            fig, axs = plt.subplots(2, 4, figsize=(22, 12), squeeze=False)
        else:
            logger.warning("Diegnostic for table '{table_name}' not defined")
    else:
        logger.warning("Plot type '{plot_type}' unknown")

    axs = axs.flatten()
    ax_ctr = 0
    for var, plot_arg in det_vars.items():
        if plot_type == "hist":
            plot_table_hist(tc.__dict__[table_name], var, axs[ax_ctr], **plot_arg)
        elif plot_type == "scatter":
            plot_table_scatter(
                tc.__dict__[table_name], var[0], var[1], axs[ax_ctr], **plot_arg
            )
        ax_ctr += 1

    plt.tight_layout()
    plt.legend()
    return fig


# %% light curve plotting


def plot_light_curve(
    tc,
    fd_src_ids=None,
    rg_src_ids=None,
    ax=None,
    ylim=None,
    legend_loc="upper center",
    plot_upper_limits=True,
    **errorbar_kwargs,
):
    """
    Plot the magnitude light curves of the passed sources.

    Parameters
    ----------
    tc : VASCA.table.TableCollection
        Either field or region that contains the light curves.
    fd_src_ids : list or int
        List or single field source IDs to plot. Default is None.
    rg_src_ids : list or int
        List or single region source IDs to plot. Default is None.
    ax : axes, optional
        Matplotlib axes to plot on. The default is None.
    ylim : list, optional
        Limits of the y axis. Default is None
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

    logger.debug("Plotting lightcurves ")

    # Setup plotting parameters
    if ax is None:
        ax = plt.gca()
    ax.invert_yaxis()
    if hasattr(ylim, "__iter__"):
        ax.set_ylim(ylim)

    plt_errorbar_kwargs = {
        "markersize": 4,
        "alpha": 0.6,
        "capsize": 0,
        "lw": 0.1,
        "linestyle": "dotted",
        "elinewidth": 0.6,
    }
    if errorbar_kwargs is not None:
        plt_errorbar_kwargs.update(errorbar_kwargs)

    # Loop over selected sources and plot
    colors = cycle("bgrcmykbgrcmykbgrcmykbgrcmyk")
    markers = cycle("osDd<>^v")
    ctr = 0

    # Get light curves dictionary
    dd_lcs = tc.get_light_curve(fd_src_ids, rg_src_ids)
    src_ids = list(dd_lcs.keys())

    for src_id, col, mar in zip(src_ids, colors, markers):

        # Every 8 markers plot open symbols
        ctr += 1
        mfc = col
        if ctr > 8:
            mfc = "None"

        # Get light curve
        lc = dd_lcs[src_id]

        # Get arrays
        src_lab = str(src_id)
        uplims = np.zeros(len(lc))
        sel = lc["mag"] > 0
        mags = lc["mag"]
        mags_err = lc["mag_err"]
        ul = lc["ul"]

        # Modify arrays if upper limits are plotted
        if plot_upper_limits:
            uplims = lc["mag"] < 0
            sel = np.ones(len(lc), dtype=bool)
            mags = mags * ~uplims + ul * uplims
            mags_err = mags_err * ~uplims + 0.1 * uplims

        # Plot
        plt.errorbar(
            lc["time_start"][sel],  # TODO: Move this to the bin center
            mags[sel],
            yerr=mags_err[sel],
            lolims=uplims[sel],
            color=col,
            markeredgecolor=col,
            markerfacecolor=mfc,
            marker=mar,
            label=src_lab,
            **plt_errorbar_kwargs,
        )
    ax.legend(loc=legend_loc, ncol=7, fontsize="small", handletextpad=0.05)
    ax.set_xlabel("MJD")
    ax.set_ylabel("Magnitude")

    return ax
