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
# ruff: noqa: T201, T203, E402

# %% [markdown]
# # Resource Manager
#
# With [](#ResourceManager) VASCA provides a utility that helps managing the raw input
# data. Data volumes processed by VASCA are generally pretty large and use cases as well
# as computation and storage resources can vary. [](#ResourceManager) adds an
# abstraction layer that is flexible enough to varying contexts while exposing a
# consistent API to the rest ofVASCA's pipeline functions.
#
# As an example, the processing of GALEX data for the proof-of-principle study was done
# by downloading raw data from [MAST](https://astroquery.readthedocs.io/en/latest/mast/mast.html)
# to a local directory, running the pipeline on a regular office laptop. This directory
# was then cloud-synced via DESY's NextCloud service to allow collaborative work with
# multiple users on the same dataset.
#
# Another use case is the unit testing for this package as well as this tutorial, wich
# should both work in a development environment and the GitHub continuous integration
# workflows.

# %% [markdown]
# ## Configuration
# Configuration files are used to specify file locations, environment variables and even
# specific data products that are relevant for the processing of a specific instrument's
# raw data. These can be freely edited by users to include data locations items as the
# use case requires. [](#ResourceManager) has the necessary consistency checks to warn
# if any miss-configuration has happened. So try it out!
#
# ### `.env`
# Text file located at the root directory of the package. This is read by the resource
# manager at initialization which uses [`dotenv`](https://saurabh-kumar.com/python-dotenv/)
# to set the environment variables temporarily during run time. See [`.env_template`](../../.env_template)
# when using VASCA for the first time.
#
# ### `resource_envs.yml`
# Configuration file specifying the required environment variables and associated
# attributes like a name, a project name and a short description to help other users to
# understand what variable is used for.
#
# ### `resource_catalag.yml`
# Configuration file that associates directory or file items to specific environment
# variables. Each item has a name, description, type, and path attribute.
#
# ```{note}
# The YAML configuration files are stored under the `vasca` module in a subdirectory
# named [`resource_metadata`](https://github.com/rbuehler/vasca/tree/main/vasca/resource_metadata)
# ```

# %% [markdown]
# ## Example

# %% [markdown]
# Initialize the ResourceManager and see what metadata it parsed from the config files.

# %%
from pprint import pprint
from vasca.resource_manager import ResourceManager

rm = ResourceManager()

# %% tags=["hide-output"]
# Resource item catalog
pprint(rm.metadata["catalog"], sort_dicts=False)

# %% tags=["hide-output"]
# Resource environment variables
pprint(rm.metadata["envs"], sort_dicts=False)

# %% [markdown]
# The main functionality: receiving paths from the ResourceManager to specific resource
# items:

# %%
rpath = rm.get_path(resource="gal_visits_list", storage="vasca")
print(rpath)

# %% [markdown]
# All paths returned by [](#get_path) are verified:

# %%
from pathlib import Path

Path(rpath).exists()

# %% [markdown]
# Otherwise actionable error messages are given.
#
# Resource name not found:

# %% tags=["raises-exception"]
rm.get_path(resource="foo", storage="vasca")

# %% [markdown]
# Storage system not recognized:

# %% tags=["raises-exception"]
rm.get_path(resource="gal_visits_list", storage="foo")
