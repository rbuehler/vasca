#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualization related methods for VASCA
"""

from collections import OrderedDict
from itertools import cycle

import astropy.units as u
import healpy as hpy
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.wcs import WCS
from loguru import logger
from matplotlib.colors import LogNorm

import vasca.utils as vutils

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
        "transform": ax.get_transform("world"),
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
        "transform": ax.get_transform("world"),
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
        "cmap": "gray_r",
        "origin": "lower",
        "norm": LogNorm(),
    }

    if img_kwargs is not None:
        plt_img_kwargs.update(img_kwargs)

    graph = ax.imshow(field.ref_img, **plt_img_kwargs)

    return graph


def plot_field_sky(
    field,
    plot_sources=True,
    plot_detections=True,
    plot_map=True,
    kwargs_fig=None,
    kwargs_map=None,
    kwargs_sources=None,
    kwargs_detections=None,
):
    """
    Optionally, plot field sources, with or without individual visit detections,
    and reference image as background.

    Parameters
    ----------
    field: :class:`~vasca.field.BaseField`
        VASCA field to be plotted.
    plot_sources : bool, optional
        Plot field sources, The default is True.
    plot_detections : bool, optional
        Plot visit detections. The default is True.
        Detections are only plotted if ``plot_sources`` is True.
    plot_map : bool, optional
        Plot reference image in the background. The default is True.
    kwargs_* : dict, optional
        Keyword arguments passed to respective plotting functions.

    Returns
    -------
    :py:class:`~matplotlib.figure.Figure`

    """

    logger.debug("Plotting sky map and/or sources'")

    # Figure settings
    kwargs_fig = {
        "figsize": (8, 7),
        "num": field.name,
        "constrained_layout": True,
        **(kwargs_fig if kwargs_fig is not None else dict()),
    }

    # Create new figure
    plt.close(kwargs_fig["num"])
    fig = plt.figure(**kwargs_fig)

    # Return empty figure if source nor map is selected
    if not plot_sources and not plot_map:
        return fig

    # Enforce Astropy WCS
    if isinstance(field.ref_wcs, WCS):
        ax = plt.subplot(projection=field.ref_wcs)
        ax.coords["ra"].set_major_formatter("d.dd")
        ax.coords["dec"].set_major_formatter("d.dd")
        ax.set_xlabel("Ra")
        ax.set_ylabel("Dec")
        ax.coords.grid(True, color="grey", ls="-", lw=0.5)
    else:
        raise ValueError(
            f"Expected WCS attribute for field, got {type(field.ref_wcs).__name__}"
        )

    if plot_map:
        graph = plot_field_sky_map(
            field, ax, **(kwargs_map if kwargs_map is not None else dict())
        )
        fig.colorbar(graph, label="Intensity [a.u.]", shrink=0.85)

    if plot_sources:
        plot_field_sky_sources(
            field,
            ax,
            plot_detections=plot_detections,
            **(kwargs_sources if kwargs_sources is not None else dict()),
            **(kwargs_detections if kwargs_detections is not None else dict()),
        )

    return fig


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
    ax : axes
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
    previosly created with region.add_coverage_hp.

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
    ax : axes
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
    vals: [float]
        Histogram bin values.
    bins : [float]
        Histogram bins, default is [selected, all] events.

    """
    logger.debug(f"Plotting histogram of variable '{var}'")

    if ax is None:
        ax = plt.gca()

    # Set marker properties for sources
    plot_kwargs = {
        "bins": "auto",
        "histtype": "barstacked",
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

    vals, bins, *_ = ax.hist(
        data,
        label=["unselected_" + str_nrnotsel, "selected_" + str_nrsel],
        **plot_kwargs,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Counts")

    return ax, vals, bins


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


# %% pipeline diagnostic plotting


# TODO: there are some GALEX specific variables
# will need to be adapted to other missions
def plot_pipe_diagnostic(tc, table_name, plot_type, fig_size=(12, 8)):
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

    Returns
    -------
    fig : matplotlib.figure
        Figure used for plotting
    var_plt : dict
        Dictionary with plot variables keys and plot settings values

    """

    var_plt = OrderedDict()
    if plot_type == "hist":
        # Detections diagnostic
        if table_name == "tt_detections":
            var_plt["s2n"] = {"logx": True}
            var_plt["mag"] = {"range": [14.5, 26.5]}
            var_plt["mag_err"] = {"range": [0.0, 3.0]}
            var_plt["r_fov"] = {"range": [0.0, 0.7]}
            var_plt["point_src_prob"] = {}
            var_plt["artifacts"] = {"histtype": "step"}
            var_plt["bright_match"] = {}
            var_plt["pos_err"] = {"range": [0.0, 0.0025]}
            fig, axs = plt.subplots(4, 2, figsize=fig_size, squeeze=False)
        elif table_name == "tt_sources":
            var_plt["nr_det"] = {}
            var_plt["flux_cpval"] = {}
            var_plt["mag_mean"] = {}
            var_plt["mag_dmax"] = {}
            var_plt["mag_dmax_sig"] = {"logx": True, "range": [-3, 2]}
            var_plt["nr_fd_srcs"] = {}
            fig, axs = plt.subplots(2, 3, figsize=fig_size, squeeze=False)
        else:
            logger.warning("Diegnostic for table '{table_name}' not defined")

    elif plot_type == "scatter":
        # Detections diagnostic
        if table_name == "tt_detections":
            var_plt[("s2n", "mag")] = {
                "invert_yaxis": True,
                "xlim": [1, 100],
                "ylim": [14.5, 26.5],
            }
            var_plt[("r_fov", "mag")] = {"invert_yaxis": True, "ylim": [15.5, 26.5]}
            var_plt[("point_src_prob", "mag")] = {
                "invert_yaxis": True,
                "ylim": [14.5, 26.5],
            }
            var_plt[("artifacts", "mag")] = {
                "invert_yaxis": True,
                "ylim": [14.5, 26.5],
            }
            var_plt[("r_fov", "artifacts")] = {}
            var_plt[("mag_err", "mag")] = {
                "xlim": [0.01, 3],
                "invert_yaxis": True,
                "ylim": [14.5, 26.5],
            }
            fig, axs = plt.subplots(3, 2, figsize=fig_size, squeeze=False)
        elif table_name == "tt_sources":
            var_plt[("flux_cpval", "mag_mean")] = {
                "xscale": "log",
                "xlim": [1e-23, 1.0],
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            var_plt[("flux_cpval", "mag_dmax_sig")] = {
                "xscale": "log",
                "xlim": [1e-23, 1.0],
                "ylim": [0, 15],
            }
            var_plt[("flux_cpval", "nr_det")] = {
                "xscale": "log",
                "xlim": [1e-23, 1.0],
            }
            var_plt[("flux_cpval", "flux_rchiq")] = {
                "xscale": "log",
                "yscale": "log",
                "xlim": [1e-23, 2],
            }
            var_plt[("nr_det", "mag_mean")] = {
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            var_plt[("nr_fd_srcs", "mag_mean")] = {
                "invert_yaxis": True,
                "ylim": [17.5, 24.5],
            }
            fig, axs = plt.subplots(2, 3, figsize=fig_size, squeeze=False)
        else:
            logger.warning("Diegnostic for table '{table_name}' not defined")
    else:
        logger.warning("Plot type '{plot_type}' unknown")

    axs = axs.flatten()
    ax_ctr = 0
    for var, plot_arg in var_plt.items():
        if plot_type == "hist":
            plot_table_hist(tc.__dict__[table_name], var, axs[ax_ctr], **plot_arg)
        elif plot_type == "scatter":
            plot_table_scatter(
                tc.__dict__[table_name], var[0], var[1], axs[ax_ctr], **plot_arg
            )
        ax_ctr += 1

    plt.tight_layout()
    plt.legend()
    return fig, var_plt


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
        Either field, region or source that contains the light curves.
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

    # Code snippets to be added later

    # from itertools import zip_longest
    # Draw selected  light curves
    # sel_srcs = field.tt_sources["sel"]
    # if sel_srcs.sum() > 0:
    #     all_srcs_ids = field.tt_sources[sel_srcs]["fd_src_id"].data

    #     # Loop over list in chuks of 14
    #     for fd_src_ids_chunk in zip_longest(
    #         *([np.nditer(all_srcs_ids)]) * 14, fillvalue=-1
    #     ):
    #         fig_lc = plt.figure(figsize=(10, 10))
    #         fd_src_ids_chunk = np.array(fd_src_ids_chunk, dtype=np.int64).flatten()
    #         fd_src_ids_chunk = np.delete(
    #             fd_src_ids_chunk, np.where(fd_src_ids_chunk == -1)
    #         )
    #         vvis.plot_light_curve(field, fd_src_ids=fd_src_ids_chunk, ylim=[24.5, 15.5])
    #         plt.tight_layout()
    #         srcs_name = "_".join([str(elem) for elem in fd_src_ids_chunk])

    #         fig_lc.savefig(
    #             field_dir + "lcs/lc_" + str(field.field_id) + "_" + srcs_name + ".png",
    #             dpi=150,
    #         )
    #         plt.close(fig_lc)
    #


def get_cutout_bounds(cutout, data_shape, out_frame="icrs", include_diag=False):
    """
    Computes the cutout boundaries as coordinates for upper right and lower left
    pixels in coordinate system of 'out_frame'.

    Parameters
    ----------
    cutout : astropy.nddata.Cutout2D
    data_shape : tuple
    out_frame : str, optional
    include_diag : bool

    Returns
    -------
    astropy.coordinates.SkyCoord or tuple(astropy.coordinates.SkyCoord, Angle)

    """
    # corner pixels
    # lower left
    x_ll = 0
    y_ll = 0
    # upper right
    x_ur = data_shape[1]
    y_ur = data_shape[0]
    # convert to world coordinates
    ll = cutout.wcs.pixel_to_world(x_ll, y_ll)
    ur = cutout.wcs.pixel_to_world(x_ur, y_ur)
    # separtion
    diag_sep = ll.separation(ur)
    # boundary coordinates
    bound_coords = SkyCoord(
        [ll.ra, ur.ra], [ll.dec, ur.dec], frame=ll.frame
    ).transform_to(out_frame)

    if include_diag:
        return (bound_coords, diag_sep)
    else:
        return bound_coords


def select_cutout(tt_cat, cutout_bounds, frame="icrs"):
    """
    Computes a bool array to select a dataset for given cutout bounds

    Parameters
    ----------
    tt_cat : astropy.table.Table
    cutout_bounds : astropy.coordinates.SkyCoords
    frame : str, optional

    Returns
    -------
    numpy.ndarray
    """
    x, y = None, None
    if frame == "icrs" or frame == "fk5":
        x = "ra"
        y = "dec"
        mask_cutout_ra = (tt_cat[x] <= cutout_bounds[0].ra) * (
            tt_cat[x] >= cutout_bounds[1].ra
        )
        mask_cutout_dec = (tt_cat[y] >= cutout_bounds[0].dec) * (
            tt_cat[y] <= cutout_bounds[1].dec
        )
    elif frame == "galctic":
        x = "lon"
        y = "lat"
        mask_cutout_ra = (tt_cat[x] <= cutout_bounds[0].l) * (
            tt_cat[x] >= cutout_bounds[1].l
        )
        mask_cutout_dec = (tt_cat[y] >= cutout_bounds[0].b) * (
            tt_cat[y] <= cutout_bounds[1].b
        )

    return mask_cutout_ra * mask_cutout_dec


def plot_source_tumbnail(
    coord_cat,
    rg,
    cutout_size=(30, 30),
    fig_title=None,
    marker_shows="field",
    color_shows="source",
):
    """
    Create detailed plot for a source closest to the target coordinate ``coord_cat``

    Parameters
    ----------
    coord_cat : astropy.coordinates.SkyCoord
    rg : vasca.Region
    cutout_size : tuple, optional
    fig_title : str, optional
    marker_shows : str, optional
    color_shows : str, optional
    """
    # Get Source instance by matching coordinates
    src, dist = rg.get_src_from_sky_pos(coord_cat.ra, coord_cat.dec)

    # Source coordinates
    coord_src = SkyCoord(
        src.tt_sources[0]["ra"] * u.deg, src.tt_sources[0]["dec"] * u.deg, frame="icrs"
    )

    # -> Todo: This should go to the logger
    print(f"Matched source at distance {dist.to('arcsec'):1.2f}")

    # Load reference map from field with longest total exposure time
    # -> Todo: Allow user to specify field
    fd_idx_max_xposure = np.argmax(src.tt_fields["time_bin_size_sum"])
    fd = rg.fields[src.tt_fields[fd_idx_max_xposure]["field_id"]]
    fd_src_id = src.tt_fields["fd_src_id"][fd_idx_max_xposure]

    # Create rectangular cutout
    fd_cutout = Cutout2D(
        fd.ref_img,
        position=coord_src,
        size=u.Quantity(cutout_size, u.arcsec).to(u.deg),
        wcs=fd.ref_wcs,
    )
    # Cutout boundary coordinates
    fd_cutout_bounds = get_cutout_bounds(fd_cutout, fd_cutout.data.shape)

    # IDs of fields and sources in cutout
    sel_cutout = select_cutout(rg.tt_detections, fd_cutout_bounds)
    field_ids = np.unique(rg.tt_detections[sel_cutout]["rg_fd_id"].data)
    source_ids = np.unique(rg.tt_detections[sel_cutout]["rg_src_id"].data)

    # Number of fields and sources in cutout
    n_fields = len(field_ids)
    n_sources = len(source_ids)

    # Detections plotting: Markers and colors according to field_id/src_id
    # (source colors/markers are fixed)
    markers = dict.fromkeys(["field", "source"])
    colors = dict.fromkeys(["field", "source"])
    for key, ids, n_ids in zip(
        ["field", "source"], [field_ids, source_ids], [n_fields, n_sources]
    ):
        markers[key] = {
            id: m for id, m in zip(ids, vutils.marker_set(n_ids, exclude=["*", "X"]))
        }
        colors[key] = {
            id: m for id, m in zip(ids, vutils.color_palette("turbo", n_ids))
        }

    def get_style(source_id, field_id):
        """Helper function returning a tuple with color and marker"""
        ids = {"field": field_id, "source": source_id}
        color = colors[color_shows][ids[color_shows]]
        marker = markers[marker_shows][ids[marker_shows]]
        return (color, marker)

    # Plot
    # -> Todo:  This is not optimal, method should accept axis as well
    plot_name = f"Source_{src.tt_sources['rg_src_id'][0]}"
    plt.close(plot_name)
    fig, ax = plt.subplots(
        num=plot_name,
        figsize=(6, 6),
        subplot_kw=dict(projection=fd_cutout.wcs),
        tight_layout=False,
    )
    if fig_title is not None:
        fig.suptitle(f"{fig_title}: {plot_name.replace('_', ' ')}")

    # Map
    ax.imshow(fd_cutout.data, cmap="gray", norm=LogNorm(), interpolation=None)
    ax.set_autoscale_on(False)

    # Style settings
    marker_size = 7**2
    marker_edge_lw = 1.5
    label_fontsize = 8

    # Source
    # Center
    ax.scatter(
        coord_src.ra,
        coord_src.dec,
        transform=ax.get_transform("world"),
        s=marker_size,
        edgecolor="red",
        facecolor="none",
        lw=marker_edge_lw,
        marker="*",
        label="Target source",
        zorder=5,
    )
    # Positional error
    s = SphericalCircle(
        coord_src,
        src.tt_sources["pos_err_mean"][0] * u.deg,
        edgecolor=(1, 0, 0, 0.5),
        facecolor=(1, 0, 0, 0.2),
        lw=1.5,
        transform=ax.get_transform("world"),
        zorder=5,
    )
    ax.add_patch(s)

    # Source detections
    for det in src.tt_detections:
        coord_det = SkyCoord(det["ra"] * u.deg, det["dec"] * u.deg, frame="icrs")

        field_id_det = det["rg_fd_id"]
        source_id_det = det["rg_src_id"]
        _, marker = get_style(source_id_det, field_id_det)

        # Center
        ax.scatter(
            coord_det.ra,
            coord_det.dec,
            transform=ax.get_transform("world"),
            s=marker_size,
            edgecolor="blue",
            facecolor="none",
            lw=1.5,
            marker=marker,
            label=f"src det. {field_id_det}",
            zorder=4,
        )
        # Positional error
        s_det = SphericalCircle(
            coord_det,
            det["pos_err"] * u.deg,
            edgecolor=(0, 0, 1, 0.5),
            facecolor=(0, 0, 1, 0.1),
            lw=0.5,
            transform=ax.get_transform("world"),
            zorder=4,
        )
        ax.add_patch(s_det)

    # Region sources
    sel_cutout = select_cutout(rg.tt_sources, fd_cutout_bounds)
    for rg_src in rg.tt_sources[sel_cutout]:
        # Skip targeted source
        rg_src_id = rg_src["rg_src_id"]
        if rg_src_id == src.tt_sources["rg_src_id"]:
            continue

        coord_rg_src = SkyCoord(
            rg_src["ra"] * u.deg, rg_src["dec"] * u.deg, frame="icrs"
        )
        # Fix source color
        color = vutils.color_palette("YlOrRd", 1)[0]
        # Center
        ax.scatter(
            coord_rg_src.ra,
            coord_rg_src.dec,
            transform=ax.get_transform("world"),
            s=marker_size,
            edgecolor=color,
            facecolor="none",
            lw=marker_edge_lw,
            marker="X",
            label=f"rg src {field_id_det}",
            zorder=5,
        )
        # Positional error
        s_det = SphericalCircle(
            coord_rg_src,
            rg_src["pos_err_mean"] * u.deg,
            edgecolor=(*color[:3], 0.5),
            facecolor=(*color[:3], 0.2),
            lw=1.5,
            transform=ax.get_transform("world"),
            zorder=5,
        )
        ax.add_patch(s_det)

    # Region detections
    sel_cutout = select_cutout(rg.tt_detections, fd_cutout_bounds)
    for det in rg.tt_detections[sel_cutout]:
        coord_det = SkyCoord(det["ra"] * u.deg, det["dec"] * u.deg, frame="icrs")

        field_id_det = det["rg_fd_id"]
        source_id_det = det["rg_src_id"]
        color, marker = get_style(source_id_det, field_id_det)

        # Center
        ax.scatter(
            coord_det.ra,
            coord_det.dec,
            transform=ax.get_transform("world"),
            s=marker_size,
            edgecolor=color,
            facecolor="none",
            lw=marker_edge_lw,
            marker=marker,
            label=f"rg det. {field_id_det}",
            alpha=0.5,
            zorder=3,
        )
        # Positional error
        s_det = SphericalCircle(
            coord_det,
            det["pos_err"] * u.deg,
            edgecolor=(*color[:3], 0.5),
            facecolor=(*color[:3], 0.1),
            lw=0.5,
            transform=ax.get_transform("world"),
            zorder=3,
        )
        ax.add_patch(s_det)

    # Modify grid
    ax.coords.grid(True, color="gray", ls="-", lw=0.75, alpha=0.3)

    # Set axis labels
    ra = ax.coords["ra"]
    dec = ax.coords["dec"]
    ra.set_major_formatter("d.ddd")
    dec.set_major_formatter("d.ddd")

    ax.set_xlabel("Ra", fontsize=label_fontsize)
    ax.set_ylabel("Dec", fontsize=label_fontsize)
    ax.xaxis.set_tick_params(labelsize=label_fontsize)
    ax.yaxis.set_tick_params(labelsize=label_fontsize)
    ax.tick_params(
        axis="x",
        labelsize=label_fontsize,
        bottom=True,
        top=True,
        direction="in",
        which="both",
        color="gray",
    )
    ax.tick_params(
        axis="y",
        labelsize=label_fontsize,
        left=True,
        right=True,
        direction="in",
        which="both",
        color="gray",
    )
