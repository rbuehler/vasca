import json
import os
from copy import deepcopy

import astropy.io.fits as fits
import astropy.units as u
import gPhoton as gp
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.wcs import WCS
from gPhoton.gphoton_utils import read_lc
from matplotlib.colors import LogNorm

from vasca.utils import tgalex_to_astrotime


def get_out_dir():
    # make sure output directory exists
    out_dir = "./out"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    return out_dir


def run_gFind(
    transient,
    ra,
    dec,
    band="NUV",
    detsize=1.1,
    maxgap=1500,
    refresh=False,
    overwrite=False,
    **gp_kwargs,
):
    # Directory to save results in
    out_dir = get_out_dir()

    # Default: reads cached gFind results from disc if it available
    if os.path.isfile(f"{out_dir}/{transient}_gFind_res.json") and not refresh:
        with open(f"{out_dir}/{transient}_gFind_res.json") as f:
            gfind_res = json.load(f)
    # Query gPhoton otherwise
    else:
        # load observations
        gfind_res = gp.gFind(
            band=band, skypos=[ra, dec], maxgap=maxgap, detsize=detsize, **gp_kwargs
        )

        # add time delta column
        dt = gfind_res[band]["t1"] - gfind_res[band]["t0"]
        # add gAperture compatible observation time ranges
        tranges = [
            [t0, t1] for t0, t1 in zip(gfind_res[band]["t0"], gfind_res[band]["t1"])
        ]
        # add query config
        config = {
            "transient": transient,
            "ra": ra,
            "dec": dec,
            "band": band,
            "detsize": detsize,
            "maxgap": maxgap,
            "refresh": refresh,
            "overwrite": overwrite,
        }
        config.update(gp_kwargs)

        gfind_res[band]["dt"] = dt
        gfind_res[band]["tranges"] = tranges
        gfind_res[band]["config"] = config

        # Save results to disc if not present before or overwrite existing
        if not os.path.isfile(f"{out_dir}/{transient}_gFind_res.json") or overwrite:
            with open(f"{out_dir}/{transient}_gFind_res.json", "w") as f:
                # convert to JSON-compatible types
                out_data = deepcopy(gfind_res)
                for key in out_data[band].keys():
                    if isinstance(out_data[band][key], np.ndarray):
                        out_data[band][key] = out_data[band][key].tolist()
                # write
                json.dump(out_data, f, indent=4)

    return gfind_res


def run_gMap(
    transient,
    ra,
    dec,
    trange,
    skyrange,
    band="NUV",
    detsize=1.1,
    counts=True,
    intensity=False,
    refresh=False,
    overwrite=True,
    return_data=True,
    **gp_kwargs,
):
    # Directory to save results in
    out_dir = get_out_dir()

    # create file names
    if counts:
        cntfile = f"{out_dir}/{transient}_counts.fits"
    else:
        cntfile = None
    if intensity:
        intfile = f"{out_dir}/{transient}_intensity.fits"
    else:
        intfile = None

    # load maps
    if refresh:
        gp.gMap(
            skypos=[ra, dec],
            trange=trange,
            skyrange=skyrange,
            band=band,
            detsize=detsize,
            cntfile=cntfile,
            intfile=intfile,
            overwrite=overwrite,
            **gp_kwargs,
        )

    # Return maps and wcs data in dictionary
    if return_data:
        gmap_res = dict.fromkeys(["counts", "intensity"])

        for k in gmap_res.keys():
            gmap_res[k] = dict.fromkeys(["data", "wcs"])

        if counts:
            with fits.open(f"{out_dir}/{transient}_counts.fits") as hdul:
                gmap_res["counts"]["data"] = hdul[0].data
                gmap_res["counts"]["wcs"] = WCS(hdul[0].header)
        if intensity:
            with fits.open(f"{out_dir}/{transient}_intensity.fits") as hdul:
                gmap_res["intensity"]["data"] = hdul[0].data
                gmap_res["intensity"]["wcs"] = WCS(hdul[0].header)

        return gmap_res
    else:
        return None


def plot_map(
    transient,
    map_type="counts",
    ra=None,
    dec=None,
    radius=None,
    annulus=None,
    save=False,
):
    # Loads map data
    out_dir = get_out_dir()
    map_file = f"{out_dir}/{transient}_{map_type}.fits"
    with fits.open(map_file) as hdul:
        map_wcs = WCS(hdul[0].header)
        map_data = hdul[0].data

    # Plot counts/intensity map

    # setup figure
    num = f"{transient}_{map_type}"
    plt.close(num)
    fig = plt.figure(num=num)
    ax = plt.subplot(projection=map_wcs)

    # plot map
    ax.imshow(map_data, norm=LogNorm())

    # plot photometry apertures
    phot_center = None
    if ra is not None and dec is not None:
        phot_center = SkyCoord(
            ra=ra * u.deg,
            dec=dec * u.deg,
            frame="icrs",
        )
    if radius is not None and phot_center is not None:
        s_phot = SphericalCircle(
            phot_center,
            radius * u.deg,
            edgecolor="white",
            facecolor="none",
            transform=ax.get_transform("world"),
        )
        ax.add_patch(s_phot)

    if annulus is not None and phot_center is not None:
        s_inner = SphericalCircle(
            phot_center,
            annulus[0] * u.deg,
            edgecolor="red",
            facecolor="none",
            transform=ax.get_transform("world"),
        )
        s_outer = SphericalCircle(
            phot_center,
            annulus[1] * u.deg,
            edgecolor="red",
            facecolor="none",
            transform=ax.get_transform("world"),
        )
        ax.add_patch(s_inner)
        ax.add_patch(s_outer)

    # style
    ax.coords["ra"].set_major_formatter("d.dddd")
    ax.coords["dec"].set_major_formatter("d.dddd")
    ax.set_xlabel("Ra")
    ax.set_ylabel("Dec")
    ax.coords.grid(True, color="grey", ls="-", lw=0.5)

    if save:
        out_dir = get_out_dir()
        plt.savefig(f"{out_dir}/{transient}_{map_type}.png", dpi=180)


def run_gAperture(
    transient,
    ra,
    dec,
    stepsize,
    radius,
    annulus,
    tranges,
    band="NUV",
    detsize=1.1,
    maxgap=1500.0,
    lc=True,
    photons=True,
    **gp_kwargs,
):
    # Directory to save results in
    out_dir = get_out_dir()

    # create file names
    if lc:
        csvfile = f"{out_dir}/{transient}_gAperture_lc.csv"
    if photons:
        photoncsvfile = f"{out_dir}/{transient}_gAperture_photons.csv"

    photon_events = gp.gAperture(
        skypos=[ra, dec],
        stepsz=stepsize,
        tranges=tranges,
        radius=radius,
        annulus=annulus,
        band=band,
        detsize=detsize,
        maxgap=maxgap,
        csvfile=csvfile,
        photoncsvfile=photoncsvfile,
        **gp_kwargs,
    )
    return photon_events


def plot_lc(transient, trange=None, xdates=True):
    # load light curve data
    out_dir = get_out_dir()
    lc_file = f"{out_dir}/{transient}_gAperture_lc.csv"
    lc_data = read_lc(lc_file)

    # filter bad time stamps
    if trange is None:
        sel = lc_data["t_mean"] > 1
    else:
        sel = (lc_data["t0"] >= trange[0]) & (lc_data["t1"] <= trange[1])

    # get human readable time stamps
    ts = np.asarray([tgalex_to_astrotime(t, "iso") for t in lc_data["t_mean"]])
    # setup figure
    num = f"{transient}_lc"
    plt.close(num)
    fig, ax = plt.subplots(num=num)

    # plot magnitude
    ax.errorbar(
        lc_data[sel]["t_mean"] if not xdates else ts[sel],
        lc_data[sel]["mag_bgsub"],
        yerr=[lc_data[sel]["mag_bgsub_err_1"], lc_data[sel]["mag_bgsub_err_2"]],
        ls="",
        color="k",
        marker="o",
        ms=3,
    )
    ax.set_title("calibrated light curve")
    ax.invert_yaxis()
    ax.set_xlabel("GALEX time")
    ax.set_ylabel("[mag]")
    # ax.set_xscale("log")
