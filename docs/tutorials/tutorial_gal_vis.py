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
# # GALEX Visits
#
# In this tutorial we are going to inspect the _GALEX visits table_. This is the full
# list of observations by GALEX available on MAST. During initialization of a [](#GALEXField)
# with given field ID this list is searched for all visits belonging to this field.

# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import pandas as pd
from astropy.table import Table

from vasca.resource_manager import ResourceManager
import vasca.utils as vutils

# %%
# Initialize ResourceManager
rm = ResourceManager()
gal_visits = rm.get_path("gal_visits_list", "vasca")

# %%
# tt_gal_visits is a table containing info about all GALEX visits
tt_gal_visits = Table.read(gal_visits)

# %%
tt_gal_visits.info()

# %%
# Convert astropy Table to pandas DataFrame for better query/aggregation methods
# String columns are bytes in the Table representation and need to decoded in the DataFrame
df_gal_visits = tt_gal_visits.to_pandas().apply(
    lambda x: x.str.decode("utf-8") if x.dtype == "O" else x
)

# %% tags=["remove-input"]
show(
    df_gal_visits.describe(),
    classes="display nowrap compact",
    scrollY="300px",
    scrollCollapse=True,
    paging=False,
    columnDefs=[{"className": "dt-body-left", "targets": "_all"}],
)


# %% [markdown]
# ## Surveys
# Find more about GALEX surveys on the documentation [page](http://www.galex.caltech.edu/researcher/techdoc-ch2.html#2)
# **AIS**
# : All-Sky Imaging Survey. Exposure time of 100s over 26,000 square degrees of sky
# reaching a depth of mAB = 20-2 in both bands.
#
# **MIS**
# : Medium Sky Imaging Survey. Single orbit exposures (1,500s) of 1000 square degrees in
# positions that match the Sloan Digital Sky Survey (SDSS) spectroscopic footprint.
#
# **DIS**
# : Deep Sky Imaging Survey. Exposure goal of 30,000s over 80 square degrees of sky.
#
# **NGS**
# : Nearby Galaxy Survey. Nearby galaxies with a nominal exposure time of 1000s to 1500s.
#
# **GII**
# : Guest Investigator Imaging. 33% of observing time per year to peer reviewed
# proposals from the community.
#
# **CAI**
# : Calibration Imaging. White dwarf standards for the purposes of calibration.

# %% [markdown]
# ### Survey Stats
#
# Below are several statistics on the GALEX data set.
#
# ```{tip}
# One of the many useful utility functions for visualization is [](#color_palette).
# It returns a list of colors for given [color map](https://matplotlib.org/stable/users/explain/colors/colormaps.html)
# and specified number of elements
# ```

# %%
survey_names = ["AIS", "MIS", "DIS", "NGS", "GII", "CAI"]
n_surveys = len(survey_names)
survey_colors = vutils.color_palette("bright", n_surveys, show_in_notebook=True)

# %%
# Aggregate
df_survey_grpd = df_gal_visits.groupby("survey")
survey_stats = {}
for survey in survey_names:
    df_survey_visits = df_survey_grpd.get_group(survey)
    stats = {}
    stats["n_vis"] = len(df_survey_visits)
    stats["n_fd"] = len(df_survey_visits.ParentImgRunID.unique())
    stats["texp"] = df_survey_visits.nexptime.sum()
    survey_stats[survey] = stats

df_survey_stats = pd.DataFrame().from_dict(survey_stats).T
# df_survey_stats

# %% tags=["remove-input"]
show(
    df_survey_stats,
    classes="display nowrap compact",
    scrollY="300px",
    scrollCollapse=True,
    paging=False,
    columnDefs=[{"className": "dt-body-left", "targets": "_all"}],
)

# %% [markdown]
# **Number of visits**

# %% tags=["hide-input"]
plot_name = "n_vis"
fig, ax = plt.subplots(num=plot_name, figsize=(3, 3), layout="constrained")

total = df_survey_stats.n_vis.sum()
text = AnchoredText(f"Total: {total:1.0f}", loc="upper right")
ax.add_artist(text)

bars = ax.bar(
    survey_names,
    df_survey_stats.n_vis,
    color=survey_colors,
)
ax.bar_label(bars, labels=[f"{x:1.0%}" for x in df_survey_stats.n_vis / total])
ax.margins(y=0.2)

ax.set_ylabel("Number of visits")
ax.set_xlabel("Survey")
ax.set_yscale("log")
ax.grid(visible=True, linewidth=0.5, color="k", alpha=0.3, zorder=0)
ax.tick_params(axis="y", direction="in", left=True, right=True, which="both")
ax.tick_params(axis="x", direction="in", top=True, bottom=True, which="both")

# %% [markdown]
# **Number of fields**

# %% tags=["hide-input"]
plot_name = "n_fields"
fig, ax = plt.subplots(num=plot_name, figsize=(3, 3), layout="constrained")

total = df_survey_stats.n_fd.sum()
text = AnchoredText(f"Total: {total:1.0f}", loc="upper right")
ax.add_artist(text)

bars = ax.bar(
    survey_names,
    df_survey_stats.n_fd,
    color=survey_colors,
)
ax.bar_label(bars, labels=[f"{x:1.1%}" for x in df_survey_stats.n_fd / total])
ax.margins(y=0.2)
ax.margins(x=0.1)

ax.set_ylabel("Number of fields")
ax.set_xlabel("Survey")
ax.set_yscale("log")
ax.grid(visible=True, linewidth=0.5, color="k", alpha=0.3, zorder=0)
ax.tick_params(axis="y", direction="in", left=True, right=True, which="both")
ax.tick_params(axis="x", direction="in", top=True, bottom=True, which="both")

# %% [markdown]
# **Total exposure time**

# %% tags=["hide-input"]
plot_name = "texp"
fig, ax = plt.subplots(num=plot_name, figsize=(3, 3), layout="constrained")

total = df_survey_stats.texp.sum()
text = AnchoredText(f"Total: {total:1.0f}", loc="upper right")
ax.add_artist(text)

bars = ax.bar(
    survey_names,
    df_survey_stats.texp,
    color=survey_colors,
)
ax.bar_label(bars, labels=[f"{x:1.0%}" for x in df_survey_stats.texp / total])
ax.margins(y=0.2)

ax.set_ylabel("Total exposure time [s]")
ax.set_xlabel("Survey")
ax.set_yscale("log")
ax.grid(visible=True, linewidth=0.5, color="k", alpha=0.3, zorder=0)
ax.tick_params(axis="y", direction="in", left=True, right=True, which="both")
ax.tick_params(axis="x", direction="in", top=True, bottom=True, which="both")

# %% [markdown]
# **Number of visits per field distribution**
#
# ```{tip}
# Another useful utility for visualization is the [](#get_hist_bins) function which
# returns the bin edges for given bin size and a list of non-uniform 1-D data sets.
# ```

# %% tags=["hide-input"]
survey_nvis = {}
for survey in survey_names:
    df_survey_visits = df_gal_visits.query("survey==@survey")
    df_fd_grpd = df_survey_visits.groupby("ParentImgRunID")
    survey_nvis[survey] = [len(df_grp) for _, df_grp in df_fd_grpd]
bins = vutils.get_hist_bins(list(survey_nvis.values()), bin_size=1)

plot_name = "nvis_dist"
fig, axs = plt.subplot_mosaic(
    [[x] for x in survey_names],
    num=plot_name,
    figsize=(3, 5),
    layout="constrained",
    sharey=False,
    sharex=True,
)

for i, survey in enumerate(survey_names):
    ax = axs[survey]

    ax.hist(survey_nvis[survey], bins=bins, color=survey_colors[i])

    text = AnchoredText(
        f"{survey}, med. = {np.asarray(survey_nvis[survey]).mean():1.0f}",
        loc="upper right",
        prop={"fontsize": 7},
    )
    ax.add_artist(text)

    ax.set_ylim((0.7, 1e5))
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.grid(visible=True, linewidth=0.5, color="k", alpha=0.3, zorder=0)
    ax.tick_params(axis="y", direction="in", left=True, right=True, which="both")
    ax.tick_params(axis="x", direction="in", top=True, bottom=True, which="both")

fig.supylabel("Number of fields")
_ = fig.supxlabel("Number of visits")

# %% [markdown]
# **Exposure time distribution**

# %% tags=["hide-input"]
bins = vutils.get_hist_bins(
    [df_gal_visits.query("survey==@x").nexptime.tolist() for x in survey_names],
    bin_size=10,  # Seconds
)

plot_name = "texp_dist"
fig, axs = plt.subplot_mosaic(
    [[x] for x in survey_names],
    num=plot_name,
    figsize=(3, 5),
    layout="constrained",
    sharey=True,
    sharex=True,
)

for i, survey in enumerate(survey_names):
    ax = axs[survey]
    data = df_gal_visits.query("survey==@survey").nexptime
    ax.hist(data, color=survey_colors[i], bins=bins, label=survey)

    text = AnchoredText(
        f"{survey}, med. = {np.median(data):1.0f} s",
        loc="upper right",
        prop={"fontsize": 7},
    )
    ax.add_artist(text)

    ax.set_ylim((10, 5e4))
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.grid(visible=True, linewidth=0.5, color="k", alpha=0.3, zorder=0)
    ax.tick_params(axis="y", direction="in", left=True, right=True, which="both")
    ax.tick_params(axis="x", direction="in", top=True, bottom=True, which="both")

fig.supylabel("Number of visits")
_ = fig.supxlabel("Exposure time [s]")
