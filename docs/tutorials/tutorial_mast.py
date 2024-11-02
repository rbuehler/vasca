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
# # MAST Download
#
# In this short tutorial the MAST query functions of [](#GALEXField) are used to
# download fresh data.

# %%
from loguru import logger
from astropy.table import Table

from vasca.field import GALEXField
from vasca.resource_manager import ResourceManager

# %%
# Activate logging
logger.enable("vasca")

# %%
# Initialize ResourceManager
rm = ResourceManager()
docs_resources = rm.get_path("docs_resources", "vasca")
gal_visits = rm.get_path("gal_visits_list", "vasca")

# %%
# Let's look at a single field with two visits
field_name = "AIS_309_1_28"  # 2 visits, contains Crab pulsar
field_id = 6381787756527353856

# %%
# Show visits info about this field
# tt_gal_visits is a table containing info about all GALEX visits
tt_gal_visits = Table.read(gal_visits)
sel_fd = tt_gal_visits["ParentImgRunID"] == field_id
# tt_gal_visits[sel_fd]

# %% tags=["remove-input"]
show(
    tt_gal_visits[sel_fd].to_pandas(),
    classes="display nowrap compact",
    scrollY="300px",
    scrollCollapse=True,
    paging=False,
    columnDefs=[{"className": "dt-body-left", "targets": "_all"}],
)

# %% tags=["hide-output"]
# Initialize a new field and
# download data from MAST
fd = GALEXField.load(
    gfield_id=field_id,
    obs_filter="NUV",
    method="MAST_REMOTE",
    load_products="ALL",
    data_path=docs_resources,
    visits_data_path=gal_visits,
)
