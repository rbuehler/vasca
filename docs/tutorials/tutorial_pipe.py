# ---
# jupyter:
#   jupytext:
#     hide_notebook_metadata: true
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: vasca-github
#     language: python
#     name: vasca-github
# ---

# %% tags=["remove-cell"]
# ruff: noqa: T201

# %% tags=["remove-input"]
import pandas as pd
from IPython.display import HTML, display
from itables import init_notebook_mode, show

init_notebook_mode(all_interactive=True)

# Modify Table CSS so with colors that work ok in light and dark themes
class_specific_css = """
.dataTable th {
    font-weight: normal;
    background-color: #075;
    color: #fff;
}
.dataTable td {
        border-color: #f0f;
        background-color: #333;
        color: #fff;
}
.dt-container {
  font-size: small;
}

"""
display(HTML(f"<style>{class_specific_css}</style>"))

# %% [markdown]
# # Pipeline
#
# This is a tutorial showcasing VASCA's pipeline flow on a simple example. We will go
# through all the steps equivalent to what is done in [](#vasca_pipe.run_from_file).
# This is the same function that is called when starting the pipeline from the CLI using ``vasca-pipe``.
#
# The goal is to create a VASCA [](#Region) from multiple [](#GALEXField) for which we
# download the raw data online from [MAST](https://astroquery.readthedocs.io/en/latest/mast/mast.html).
# We apply quality cuts and do source clustering followed by variability analysis.
#
# For this tutorial we are interested in the near-ultraviolet (NUV) observations
# by GALEX. We are going to look at neighboring/overlapping fields all of which
# contain the location of famous Tidal Disruption Event [_PS1-10jh_](https://en.wikipedia.org/wiki/Pan-STARRS#Selected_discoveries)
# discovered by Pan-STARRS and observed by GALEX in 2010.
#
# :::{figure-md} galex-fields-ps1-10jh
# <img src="../images/GALEX_fields_ps1-10jh.jpg" alt="galex_fields_ps1-10jh" class="bg-primary mb-1" width="400px">
#
# GALEX sky map with field footprints of observations around the location of PS1-10jh (
# purple crosshair). Screenshot from [MAST Portal](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html)
# :::

# %% [markdown]
# ## General Configuration
#
# The standard pipeline processing starts by reading a yaml file. To keep this tutorial
# simple, we are going to introduce parts of the configuration step by step at the point
# where they are required in the pipeline flow.
#
# ```{note}
# An important premise of the configuration is that each parameter needs to be
# configured explicitly. This means even default values are specified all the time. This
# is a design decision purposefully made in order to ensure transparent and complete
# configurations. As a result, all possible parameters are always included when looking
# at configuration file.
# ```
#
# Let's begin with the ``general`` section. Here, basic information and functionality is
# configured. The ``name`` of the pipeline run specifies also the name of directory in
# which all results will be stored. The location of output directory is at ``out_dir_base``
# relative to the root directory of the package.
#
# VASCA uses the powerful logging system provided by [loguru](https://loguru.readthedocs.io/en/stable/index.html).
# The configuration specifies the [``log_level``](https://loguru.readthedocs.io/en/stable/api/logger.html#loguru._logger.Logger.level),
# which we are going to set to debugging mode here. By default VASCA is going to save
# all logging messages in a file stored in the output directory. ``log_file`` specifies
# the name of that file, while ``default`` tells the pipeline to use a default name.
#
# Parallel processing of the field-level analysis can be enabled when setting the number
# of CPU threads ``nr_cpus > 1``.
#
# VASCA can include field-averaged reference data, if such data is available additional
# to the visit-level data from the instruments mission pipeline. To save memory/storage
# and computation time it is configurable wether to include reference sources in the
# final [](#Region)-file (``save_ref_srcs``) and to repeat already processed fields that
# are included in the region (``run_fields``).

# %%
# Dictionary holding the configuration
config = {}

# General section of the configuration
config["general"] = {
    "name": "simple_pipe",
    "out_dir_base": "docs/tutorial_resources/vasca_pipeline",
    "log_level": "DEBUG",
    "log_file": "default",
    "nr_cpus": 3,
    "save_ref_srcs": True,
    "run_fields": True,
}

# %% [markdown]
# :::{tip}
# In case the output location is on a remote server and multiple users should be
# able to edit the data, i.e., reprocess data using an updated configruation, then
# one needs to manage user access priviliges. For convenience this can be done
# using [``umask``](https://en.wikipedia.org/wiki/Umask):
# ```Python
# import os
# os.umask("0o003", 0)
# ```
# This will grand user and group full permissions. The setting is only persistant
# throughout the Python session.
# :::

# %% [markdown]
# The pipeline begins with some prerequisites including enabling logging and creating
# the output directory

# %%
import sys
from loguru import logger
from pathlib import Path
from importlib.resources import files

# Setup output directory
out_dir_base = Path(files("vasca").parent / config["general"]["out_dir_base"])
out_dir_base.mkdir(parents=True, exist_ok=True)

pipe_dir = out_dir_base / config["general"]["name"]
pipe_dir.mkdir(parents=True, exist_ok=True)

# Path to log file
log_file = (
    pipe_dir / f'{config["general"]["name"]}.log'
)  # Pipeline name is used by default

# Logger configuration, both to stdout and .log file
log_cfg = {
    "handlers": [
        {
            "sink": sys.stdout,
            "format": "<green>{time:YYYY-MM-DD HH:mm:ss.SSSS}</green> "
            "<cyan>{name}</cyan>:<cyan>{line}</cyan> |"
            "<level>{level}:</level> {message}",
            "level": config["general"]["log_level"],
            "colorize": True,
            "backtrace": True,
            "diagnose": True,
        },
    ],
}
logger.configure(**log_cfg)  # Set config
logger.add(log_file)  # Add file logger
logger.enable("vasca")  # Enable logging

# Some initial logging messages
logger.info(f"Runing '{config['general']['name']}'")
logger.debug(f"Config. file: {log_file}")

# %% [markdown]
# ## Resources
# Next is the section about resource handling. This specifies the method used to load
# (``load_method``) field data and which data products should be included (tables,
# tables plus visit-level images, or only metadata ``load_products``). Additionally
# the ``coadd_exists`` flag tells the pipeline wether it can expect co-added (field-
# averaged) data. Finally, ``field_kwargs`` allows to pass parameters directly to
# the ``init`` function of a field class.
#
# Here we are going to initialize fields from local raw data, if present. Else the
# required data is downloaded from MAST. All data products will be included including
# co-add data. Using the [](#ResourceManager) we can tell the field class where to
# locate the data.

# %%
from vasca.resource_manager import ResourceManager

# Get the data locations using ResourceManager
rm = ResourceManager()
data_path = rm.get_path("docs_resources", "vasca")
visits_data_path = rm.get_path("gal_visits_list", "vasca")

# Resource section of the configuration
config["resources"] = {
    "load_method": "MAST_LOCAL",
    "load_products": "ALL",
    "coadd_exists": True,
    "field_kwargs": {
        "data_path": data_path,
        "visits_data_path": visits_data_path,
    },
}

# %% [markdown]
# ## Observations
# The observations section of the configuration is responsible for configuring
# which combination of instrument (``observatory``) and filter (``obs_filter``)to
# load data. Also it specifies the exact list of fields to load (``obs_field_ids``).
#
# In a later step we will also add here the selection parameters (``selection``) used
# for the quality cuts on the data and the field-level clustering settings
# (``cluster_det``).

# %%
config["observations"] = [
    {
        "observatory": "GALEX",
        "obs_filter": "NUV",
        "obs_field_ids": [
            3880803393752530944,  # MISGCSN2_10493_0117
            2529969795944677376,  # ELAISN1_09
            2597664528044916736,  # PS_ELAISN1_MOS15
        ],
        # "cluster_det": {},
        # "selection": {},
    },
    # More instruments/filters...
]

# %% [markdown]
# Find below the visit metadata about the fields under investigation.

# %% tags=["hide-input"]
from astropy.table import Table

df_gal_visits = (
    Table.read(visits_data_path)
    .to_pandas()
    .apply(lambda x: x.str.decode("utf-8") if x.dtype == "O" else x)
)
query = f"ParentImgRunID in {list(config['observations'][0]['obs_field_ids'])}"
df_select = df_gal_visits.query(query)
show(
    df_select,
    classes="display nowrap compact",
    scrollY="300px",
    scrollCollapse=True,
    paging=False,
    columnDefs=[{"className": "dt-body-left", "targets": "_all"}],
)

# %% [markdown]
# In the next step we will initialize a VASCA [](#Region) with all fields sequentially.
# [](#load_from_config) is a convenience function that acts as an interface between the
# region object and field-specific loading functions. This will downloads the data from
# MAST, it will detect if the data is already present on disc and loads the cashed
# files. To safe compute time later, a VASCA-field file is written to the download
# location so that one can use this file instead of creating a new field from raw data.
# This will be used during the field-level [processing](#field-level-analysis).

# %% tags=["hide-output"]
from vasca.region import Region

rg = Region.load_from_config(config)

# %% [markdown]
# This populates the region object with all specified fields, the relevant metadata is
# stored in {term}`tt_fields`.

# %%
rg.info()
# rg.tt_fields

# %% tags=["remove_input"]
df_tt_fields = rg.tt_fields.to_pandas().apply(
    lambda x: x.str.decode("utf-8") if x.dtype == "O" else x
)
show(
    df_tt_fields,
    classes="display nowrap compact",
    scrollY="300px",
    scrollCollapse=True,
    paging=False,
    columnDefs=[{"className": "dt-body-left", "targets": "_all"}],
)

# %% [markdown]
# ## Field-level analysis
#
# The field-level analysis incorporates, first, the data reduction and parameter mapping
# from raw data to VASCA field objects, second, the data quality selection and finally
# source clustering on the remaining high-quality detections.
#
# The first step is implicitly taken care of by the [](#GALEXField) class, where the raw
# data is loaded and only the column parameters are kept that are specified in the [](#tables_dict)
# module. A shortcut is provided through the [](#Region.get_field) method which is an
# interface to the ``load`` method of
# any field class.
#
# The configuration for the next two step requires the ``selection`` and ``cluster_det``
# entries under the observations section.
#
# ### Data selection
# ```{note}
# A crucial part of VASCA's flexibility to adapt to raw data of virtually any instrument
# comes from the fact that the parameter list used for data quality selection is not
# fixed and is allowed to vary for different instruments and filters. The only
# requirement is an existent entry in the [](#tables_dict) module for any parameter and
# a corresponding field class that includes these parameters in the {term}`tt_detections`
# table.
# ```
#
# The API for the data selection is provided by the [](#TableCollection.select_rows)
# method. Each entry under selection maps to this interface. The ``table`` parameters
# specifies which table to select on. Any selection operation modifies the ``sel``
# column of a given table. It contains boolean values so ``0`` means _unselected_ and
# ``1`` means _selected_.
#
# By specifying the ``presel_type``parameter, one controls the logic by which an
# existing selection is combined with a new one. The ``sel_type`` parameter specifies
# the logic by which the selection on a set of multiple column parameters is combined.
# Parameters ``range`` and ``bitmask`` provide the column parameter and artifact
# bitflag values that are used to make the selection. Using ``set_range`` on can choose
# to clip values of a certain column to minimum and maximum values.
#
# In combination with ``sel_type = "is_in"`` and ``var`` parameters, it is possible to
# select the rows of given column ``var``in the target table if a value is also present
# in the same column of a reference table (``ref_table``).

# %%
import numpy as np

# Updating the observations for GALEX-NUV observations
config["observations"][0].update(
    {
        "selection": {
            # Quality cuts on visit-level detections
            "det_quality": {
                "table": "tt_detections",
                "presel_type": "and",
                "sel_type": "and",
                "range": {
                    "s2n": [3.0, np.inf],
                    "r_fov": [0.0, 0.5],
                    "ellip_world": [0.0, 0.5],
                    "size_world": [0.0, 6.0],
                    "class_star": [0.15, 1.0],
                    "chkobj_type": [-0.5, 0.5],
                    "flux_app_ratio": [0.3, 1.05],
                },
                "bitmask": {
                    "artifacts": [2, 4, 8, 128, 256],
                },
                "set_range": {"pos_err": [0.5, 5]},
            },
            # Quality cuts on field-averaged detections
            "coadd_det_quality": {
                "table": "tt_detections",
                "presel_type": "and",
                "sel_type": "and",
                "range": {
                    "s2n": [5.0, np.inf],
                    "r_fov": [0.0, 0.5],
                    "ellip_world": [0.0, 0.5],
                    "size_world": [0.0, 6.0],
                    "class_star": [0.15, 1.0],
                    "chkobj_type": [-0.5, 0.5],
                    "flux_app_ratio": [0.3, 1.05],
                },
                "bitmask": {
                    "artifacts": [2, 4, 8, 128, 256],
                },
            },
            # Selection on only those detections wich are part of clusters
            "det_association": {
                "table": "tt_detections",
                "presel_type": "and",
                "sel_type": "is_in",
                "ref_table": "tt_sources",
                "var": "fd_src_id",
            },
        },
    }
)

# %% [markdown]
# ### Clustering
#
# Also the field-level clustering configuration showcases VASCA's modular
# approach. In the ``cluster_det`` section on specifies the clustering algorithm
# which, in principle, can be different for each instrument and filter. Although,
# at the moment only [mean-shift clustering](https://en.wikipedia.org/wiki/Mean_shift)
# is supported by VASCA.
#
# Again, the responsible API is provided by [](#TableCollection.cluster_meanshift).
# This method wraps a method provided by the [scikit-learn package](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MeanShift.html). The end result is that each field optains
# a new {term}`tt_visits` table that lists all identified sources as defined
# by their clustered detections. Sources have at the minimum one and as manny as ``n_vis``
# detections.
#
# Mean-shift is well suited for this use case due to several reasons. Most importantly
# it is that the algorithm doesn't require the total number of clusters as a parameter.
# In fact it is determining that number which would be otherwise very difficult to
# predict from the visit-level detections before having done the clustering.
#
# Another reason is its relatively simple algorithm where only one parameters is
# required. It is called the ``bandwidth`` which means, translated to the astronomy use
# case, the radial size of a typical source on the sky. It should be roughly chosen
# to match the instrument's PSF, which, for GALEX, is about 5 arcseconds. We set
# it slightly smaller to limit false associations also considering that the
# source center is usually much better constrained than the PSF might suggest.

# %%
# Updating the observations for GALEX-NUV observations continued...
config["observations"][0].update(
    {
        "cluster_det": {
            "meanshift": {
                "bandwidth": 4,
                "seeds": None,
                "bin_seeding": False,
                "min_bin_freq": 1,
                "cluster_all": True,
                "n_jobs": None,
                "max_iter": 300,
                "table_name": "tt_detections",
            },
        },
    },
)

# %% [markdown]
# ### Pipeline flow
# According to the configuration above, we can finally run the analysis. VASCA
# implements parallel processing ([](inv:#*.Pool.starmap)) for this part of the pipeline
# by applying the [](#run_field) method in parallel for each field.

# %% [markdown]
# First, the parameters for [](#run_field) are collected.

# %%
import vasca.utils as vutils

# Collect parameters from config
fd_pars: list[list[int, str, Region, dict]] = []
vobs: list[dict] = config["observations"]

obs_nr: int
field_nr: str
# Loop over observation list index
for obs_nr, _ in enumerate(vobs):
    # Loop over fields
    for field_nr in vobs[obs_nr]["obs_field_ids"]:
        # Construct VASCA field ID (prepending instrument/filter identifier)
        iprefix: str = vutils.dd_obs_id_add[
            vobs[obs_nr]["observatory"] + vobs[obs_nr]["obs_filter"]
        ]
        field_id: str = f"{iprefix}{field_nr}"
        # Save parameters outside loop
        fd_pars.append([obs_nr, field_id, rg, config])


# %% [markdown]
# Second, all fields are processed in parallel.

# %% tags=["hide-output"]
from multiprocessing.pool import Pool
import vasca.vasca_pipe as vpipe

# Run each field in a separate process in parallel
nr_cpus = config["general"]["nr_cpus"]
logger.info(f"Analyzing {len(fd_pars)} fields on {nr_cpus} parallel threads.")

with Pool(processes=nr_cpus) as pool:
    pool_return = pool.starmap(vpipe.run_field_docs, fd_pars)
pool.join()

logger.info("Done analyzing individual fields.")

# %% [markdown]
# Finally, the pool results are unpacked and the region object is updated with
# processed field information.
#
# A memory-saving procedure is used where first all fields are brought to the
# scope of the region object by filling the ``Region.field`` dictionary from
# which the field-level data is taken and stacked in respective region-owned
# tables using the [](#Region.add_table_from_fields) method. After this step
# all fields are discarded from the scope, to be deleted from the garbage collector.
#
# ```{hint}
# In VASCA a [](#Region) object keeps track of its fields in the ``Region.tt_fields``
# table. At any time one can load field data of a specific field using [](#Region.get_field).
# ```

# %% tags=["hide-output"]
from astropy.table import unique

# Loop over processed fields
for field in pool_return:
    # Populate region field dictionary with field data
    rg.fields[field.field_id] = field
    logger.info(f"Added field {field.field_id} from pool to region")

# Add field tables to region

# Visits metadata
rg.add_table_from_fields("tt_visits")
rg.tt_visits = unique(rg.tt_visits, keys="vis_id")

# Clustered sources
rg.add_table_from_fields("tt_sources")

# Visit-level detections
rg.add_table_from_fields("tt_detections", only_selected=False)

# Field-averaged detectsion
rg.add_table_from_fields("tt_coadd_detections")

# Discard fields. All that was needed has been transferred to region tables
del rg.fields

# %% [markdown]
# ## Region-level analysis
#
# In the final stage of the pipeline, all region-level analysis steps are performed.
# This stage encompasses three key tasks: first, managing sources located in
# overlapping sky regions; second, evaluating statistics for use in variability
# detection; and finally, preparing the pipeline results for writing to disk and
# generating the VASCA variable source catalog.

# %% [markdown]
# ### Overlapping fields
#
# VASCA merges field-level sources in overlapping sky regions in a second clustering
# step, where the same mean-shift algorithm is used but with a dedicated configuration.
#
# In case field-averaged (co-added) data exists, the field-level detections are merged
# in the same way, again, with a separate configuration.

# %%
config["cluster_src"] = {
    "meanshift": {
        "bandwidth": 4,
        "seeds": None,
        "bin_seeding": False,
        "min_bin_freq": 1,
        "cluster_all": True,
        "n_jobs": 1,
        "max_iter": 300,
        "table_name": "tt_sources",
    }
}

config["cluster_coadd_dets"] = {
    "meanshift": {
        "bandwidth": 4,
        "seeds": None,
        "bin_seeding": False,
        "min_bin_freq": 1,
        "cluster_all": True,
        "n_jobs": 1,
        "max_iter": 300,
        "table_name": "tt_coadd_detections",
    }
}

# %% tags=["hide-output"]
# Cluster field sources and co-adds in parallel
ll_cluster = [
    [
        config["cluster_src"]["meanshift"],
        rg.tt_sources,
        rg.tt_detections,
        False,
    ],
    [
        config["cluster_coadd_dets"]["meanshift"],
        rg.tt_coadd_detections,
        False,
        True,
    ],
]

with Pool(processes=2) as pool:
    pool_return = pool.starmap(vpipe.run_cluster_fields, ll_cluster)
pool.join()

# Copy parallelized results into region
for pool_rg in pool_return:
    if "tt_sources" in pool_rg._table_names:
        rg.tt_sources = pool_rg.tt_sources
        rg.tt_detections = pool_rg.tt_detections
    else:
        rg.add_table(pool_rg.tt_coadd_sources, "region:tt_coadd_sources")
        rg.tt_coadd_detections = pool_rg.tt_coadd_detections

# %% [markdown]
# ### Source statistics
#
# The primary statistic used by VASCA to detect variability is the probability of
# obtaining a flux with the observed fluctuations under the assumption that the null
# hypothesis is true, that is, constant flux (``flux_cpval``). The default selection is
# such that the chance for the observed variability being purely random can be ruled out
# at [5-sigma significance](https://en.wikipedia.org/wiki/Normal_distribution#Standard_deviation_and_coverage).
#
# Additionally the normalized excess variance (``flux_nxv``) and absolute flux limits
# is used to limit a contamination due to potential hidden systematic flux variations.
#
# Similarly a selection on variation of the spatial coordinates (``pos_cpval``. ``pos_xv``)
# is used to reduce the contamination du to false association of visit-levl detections
# in the clustering step.
#
# To ensure the statistical correctness only sources with more than three detections are
# considers (``n_det>3``).
#
# ```{note}
# The source selection is configured independently for each observational filter. This
# allows to adapt to potential systematics in a very flexible way.
# ```
#
# For a more concrete example on these calculations in VASCA, have a
# look at the tutorial on [Variability Statistics](tutorial_var_stat.md).

# %%
config["selection_src"] = {
    "src_variability_nuv": {
        "table": "tt_sources",
        "presel_type": "or",
        "sel_type": "and",
        "obs_filter": "NUV",
        "range": {
            "nr_det": [3, np.inf],
            "pos_cpval": [0.0001, np.inf],
            "pos_xv": [-np.inf, 2],
            "flux_cpval": [-0.5, 0.000000573303],
            "flux_nxv": [0.001, np.inf],
            "flux": [0.144543, 575.43],
        },
    },
    "src_coadd_diff_nuv": {
        "table": "tt_sources",
        "presel_type": "or",
        "sel_type": "and",
        "obs_filter": "NUV",
        "range": {
            "nr_det": [2, np.inf],
            "pos_cpval": [0.0001, np.inf],
            "pos_xv": [-np.inf, 2],
            "coadd_ffactor": [2.0, np.inf],
            "coadd_fdiff_s2n": [7, np.inf],
        },
    },
}

# %% [markdown]
#  An additional selection may be possible if co-added data is available. In this case
#  the association between VASCA sources and the field-averaged input data can be made.
#  This serves as cross-check since most static sources should be recovered in the
#  clustering step and match the field average.

# %%
config["assoc_src_coadd"] = {
    "dist_max": 1,  # Associate nearest source below this distance in arc_sec  OR
    "dist_s2n_max": 3,  # Associate nearest source with this distance in units of "squared summed position error"
}

# %% tags=["hide-output"]
import astropy.units as uu

# Calculate source statistics
rg.set_src_stats(src_id_name="rg_src_id")
rg.set_src_stats(src_id_name="coadd_src_id")

# Match sources to coadd sources
rg.cross_match(
    dist_max=config["assoc_src_coadd"]["dist_max"] * uu.arcsec,
    dist_s2n_max=config["assoc_src_coadd"]["dist_s2n_max"],
)

# Select variable sources, deselect all sources before and
# make sure all tables containting the region source ID
# are syncronized to this selection
rg.tt_sources["sel"] = False
rg.select_from_config(
    config["selection_src"]
)  # Loops over TableCollection.select_rows()
rg.synch_src_sel(remove_unselected=False)

# Set source ID mapping table
rg.set_src_id_info()

# %% tags=["hide-output"]
# View all table attributes that have been added to the region object
rg.info()

# %% # %% tags=["hide-input"]
# View all sources that passed the selection
df_sources = (
    vutils.select_obs_filter(rg.tt_sources, obs_filter_id=1)
    .to_pandas()
    .apply(lambda x: x.str.decode("utf-8") if x.dtype == "O" else x)
)
df_select = df_sources.query("sel")
show(
    df_select,
    classes="display nowrap compact",
    scrollY="300px",
    scrollCollapse=True,
    paging=False,
    columnDefs=[{"className": "dt-body-left", "targets": "_all"}],
)

# %% [markdown]
# ### Pipeline Output
#
# With a few simple export functions the full region file, the variable source catalog
# and its pipeline configuration are saved to the pipeline output directory.
#
# This concludes the tutorial. Readers are invited to look into the post-processing
# notebooks as listed [here](index.md).

# %% tags=["hide-output"]
import yaml

# Write region file
rg.write_to_fits(file_name=pipe_dir / f'region_{config["general"]["name"]}.fits')

# Export variable source catalog (only selected sources are included)
rc = rg.get_region_catalog()
rc.write_to_fits(file_name=pipe_dir / f'region_{config["general"]["name"]}_cat.fits')

# Write used config file
yaml_out_name = pipe_dir / f'cfg_ran_{config["general"]["name"]}.yaml'
with open(yaml_out_name, "w") as yaml_file:
    yaml.dump(config, yaml_file)

logger.info("Done running VASCA pipeline.")

# %%
# To be continued ...
