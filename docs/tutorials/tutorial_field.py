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
# ruff: noqa: T201, T203, E402, B018

# %% [markdown]
# # GALEXField
#
# This is a short tutorial on the [](#GALEXField) class. This is an implementation of
# VASCA's [](#BaseField) tailored to the specific needs of GALEX archival data. It
# serves first and foremost as the interface between the instrument specifics, like raw
# data handling and nomenclature, and the rest of VASCA's functionality. Important are
# the loading functions that download raw data from [MAST](https://astroquery.readthedocs.io/en/latest/mast/mast.html)
# and load already processed field data into memory.

# %% [markdown]
# ## Example
#
# Let's load GALEX field data that is already downloaded (used for unit testing)

# %%
from pathlib import Path
from vasca.field import GALEXField
from vasca.resource_manager import ResourceManager

# %% [markdown]
# This shows the directories of raw data for two GALEX fields. The directories are named
# after their observation ID. This is used as a globally unique ID across all GALEX
# data.

# %%
rm = ResourceManager()
test_dir = rm.get_path("test_resources", "vasca")
[x.name for x in Path(test_dir).glob("*") if x.is_dir()]

# %% [markdown]
# We'll take the first one and insatiate a [](#GALEXField).

# %%
fd = GALEXField.load(
    gfield_id=6371091720297250816,
    obs_filter="NUV",
    method="MAST_LOCAL",
    load_products="TABLES",
    data_path=test_dir,
    visits_data_path=rm.get_path("gal_visits_list", "vasca"),
)

# %% [markdown]
# Using the [](#TableCollection.info) method one can see all data tables the field contains.

# %% tags=["hide-output"]
fd.info()

# %% [markdown]
# For convenience, some important attributes can be accessed directly:
#
# For example, the total exposure time over all visits or the center coordinate:

# %%
fd.time_bin_size_sum

# %%
fd.center

# %% [markdown]
# All available convenience attributes:

# %% tags=["hide-input"]
important_pars = {
    "field_id": "Unique field ID",
    "field_name": "Optional name if the field",
    "ra": "Right ascention coordinate (J2000)",
    "dec": "Declination coordinate (J2000)",
    "fov_diam": "Size of the field of view ",
    "observatory": "The instrument name",
    "obs_filter": "Observation band-pass filter",
    "center": "Center coordinate on the sky",
    "nr_vis": "Number of observational visits",
    "time_bin_size_sum": "Total exposure summed over all visits",
    "time_start": "Beggining of the first observation",
    "time_stop": "End of the last obseration",
}

for par, desc in important_pars.items():
    print(f"{par:<20} {str(getattr(fd, par)):>30} {desc:>60}")

# %%
