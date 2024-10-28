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

# %% [markdown]
# # Table Display
# This test tutorial shows a Pandas DataFrame as an interactive table.

# %% tags=["remove-input"]
import pandas as pd
from IPython.display import HTML, display
from itables import init_notebook_mode, show

from vasca.tables_dict import dd_vasca_columns

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

# Get column definitions
df_cols = pd.DataFrame.from_dict(dd_vasca_columns)

show(
    df_cols.T,
    classes="display nowrap compact",
    scrollY="300px",
    scrollCollapse=True,
    paging=False,
    columnDefs=[{"className": "dt-body-left", "targets": "_all"}],
)
