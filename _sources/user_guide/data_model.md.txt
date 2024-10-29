---
file_format: mystnb
kernelspec:
    name: vasca-github
---
# Data Model

VASCA uses a hierarchical data model which wich defines cosmic _sources_, individual
_fields_ and whole _regions_ in the celestial sky.

:::{figure-md} data-model_
<img src="../images/VASCA_data_model_v2.jpg" alt="data_model" class="bg-primary mb-1" width="600px">

The VASCA data model.
::: 

This model aims to abstract the input data, astronomical photometric detection lists, in
in the most general way possible in order to allow integration of data from multiple
instruments and parallel processing of the data.

## Input Data
Typically, astronomical photometric surveys observe the sky by segmenting it into fields
which correspond to the field of view as defined by the instrument's telescope optics. A
field is then defined by a central coordinate and the diameter of the field of view.

VASCA relies as its input on the science data products that missions/organizations create
from their observational raw data. Specifically this means VASCA takes tables of
photometric detections that have a field and visit ID. For a reference of the full list
of required columns can be see here: [{class}`BaseField`](#base_field), [{class}`Region`](#vasca.tables_dict.region).

In the case of GALEX, the detection lists created by the mission pipeline ("mcat" files)
have a large number of different observables and parameters (columns) per detection (see
GALEX docs [here](http://www.galex.caltech.edu/wiki/Public:Documentation/Appendix_A.1)).
A subset is used by VASCA, where GALEX specific parameters are added to the ones used by
[](#BaseField). This is implemented in the instrument specific [](#TableCollection)
object [](#GALEXField) and its corresponding [column definitions](#vasca.tables_dict.galex_field).

## Data Structures
All data structures in VASCA are based on [](#TableCollection). As the name suggests
these objects describe a collection of astropy (inv:astropy:std:doc#*.Table) objects. A
general API is provided to [add](#add_table) and [remove](#remove_tables) tables to and
from these objects. This is used when, for instance, an individual VASCA [](#Source) is
extracted out of a [](#Region).

### VASCA Columns
See the full list of data columns defined by VASCA below:

```{code-cell}
:tags: [remove-input]

from itables import init_notebook_mode, show
import pandas as pd
from vasca.tables_dict import dd_vasca_columns
from IPython.display import HTML, display

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
    classes="display compact",
    scrollY="300px",
    scrollCollapse=True,
    paging=False,
    columnDefs=[{"className": "dt-body-left", "targets": "_all"}],
)
```

### VASCA Tables

Below are shown all tables that are used in VASCA in various contexts (e.g. input,
clustering, post-processing).

**`tt_visits`**
: Table that stores metadata about observational visits such as the visit ID (`vis_id`),
the associated field ID (`field_id`) and observation start time (`time_bin_start`) as
well as exposure time (`time_bin_size`).

```{note}
All parameters related to obersvational timing follow the naming scheme to be compatible
with astropy [](inv:astropy:std:doc#*.BinnedTimeSeries). This ist then utilized when
creating and analyzing light curves for specific sources.
```

**`tt_fields`**
: Table storing metadata about all fields such as the field ID (`field_id`), the center
coordinates (`ra`, `dec`), total exposure time (`time_bin_size_sum`) and the size of the
field of view (`r_fov`).

**`tt_coverage_hp`**
: Region observations properties in healpix binning. RING ordering and equatorial
coordinates

**`tt_coadd_detections`**
: Reference detections table. This data is provided at the input stage of the pipeline
and contains no intra-visit information.


**`tt_detections`**
: Visit-level detections table.

**`tt_sources`**
: Source information table.

**`tt_coadd_sources`**
: Source information table. 

**`tt_src_id_map`**
: Map between region and field source IDs (`rg_src_id`, `fd_src_id`).

**`tt_filters`**
: Filters, their IDs and index.

**`tt_sed`**
: Spectral Energy Distribution from VizieR database.

**`tt_gphoton_lc`**
: Light curve from [`gPhoton.gAperture`](https://github.com/cmillion/gPhoton/blob/master/docs/UserGuide.md#gaperturepy)

**`tt_spectrum`**
: Spectrum table

**`tt_lombscargle`**
: LombScargle results information

Additional tables (e.g. `tt_simbad` or `tt_gaiadr3`) are created during catalog
cross-matching in the [post-processing stage](https://github.com/rbuehler/vasca/blob/main/vasca/examples/vasca_pipe_post_process.ipynb)















