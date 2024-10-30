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
photometric detections that have a field and visit ID. A reference of the full list
of required columns can be found below [](#vasca-columns) or in the API reference:
[{class}`BaseField`](#base_field), [{class}`Region`](#vasca.tables_dict.region).

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

Below is a short description of the most important tables that are used in VASCA:

{.glossary}
**`tt_visits`**
: Table that stores metadata about observational visits such as the visit ID (`vis_id`),
the associated field ID (`field_id`) and observation start time (`time_bin_start`) as
well as exposure time (`time_bin_size`).

```{note}
All parameters related to obersvational timing follow the naming scheme to be compatible
with astropy [](inv:astropy:std:doc#*.BinnedTimeSeries). This ist then utilized when
creating and analyzing light curves for specific sources.
```

{.glossary}
**`tt_fields`**
: Table storing metadata about all fields such as the field ID (`field_id`), the center
coordinates (`ra`, `dec`), total exposure time (`time_bin_size_sum`) and the size of the
field of view (`r_fov`).

{.glossary}
**`tt_detections`**
: Visit-level detections table. All detections are listed that pass the quality cuts on
signal-to-noise (`s2n`), artifacts (`artifacts`) and more. Both [](#Region) and [](#BaseField)
hold this table. So all VASCA-generated ID parameters are limited to the respective
scope. Especially this also stores the region/field-level cluster ID for all detections,
i.e., the ID that associates all detections belonging to one cosmic source (`rg_src_id`,
`fd_src_id`).

{.glossary}
**`tt_sources`**
: Source information table where all cosmic sources are listed that pass the cuts applied
used for variability detection. Important variability parameters (`flux_cpval`,
`flux_nxv`) and cross-matching results (`gaiadr3_match_id`, `ogrp`) are given for each
source.

```{tip}
VASCA supports field-level reference data that might be provided by the
instrument's mission pipelines. This includes field-averaged (co-added) intensity sky
maps ([](#load_sky_map)) and detection lists. This data is stored in
`tt_coadd_detections` based on which `tt_coadd_sources` is created during the
region-level clustering stage. This information is useful when visualizing pipeline
results and for consistency checks on the source clustering.
```

**`tt_coadd_detections`**
: Reference detections table. This data is provided at the input stage of the pipeline
and contains no visit-to-visit information.

**`tt_coadd_sources`**
: Reference source information table. This table is created in the region-level
clustering stage where spatially overlapping detections in `tt_coadd_sources` are merged.

**`tt_src_id_map`**
: Map between region and field source IDs (`rg_src_id`, `fd_src_id`).

**`tt_filters`**
: Instrument filters, their IDs and index.

**`tt_coverage_hp`**
: Region observations properties in healpix binning. RING ordering and equatorial
coordinates

```{tip}
Additional tables are created during [catalog cross-matching](https://github.com/rbuehler/vasca/blob/main/vasca/examples/vasca_pipe_post_process.ipynb)
(`tt_simbad` or `tt_gaiadr3`) and [post-processing stages](https://github.com/rbuehler/vasca/blob/main/vasca/examples/inspect_matches.ipynb)
(`tt_lombscargle`, `tt_sed`).

VASCA, by design, handles visit-to-visit variability. In post-processing, users can
[investigate](https://github.com/rbuehler/vasca/blob/main/vasca/examples/run_gPhoton.py)
inter-visit variability with the help of [``gPhoton``](https://github.com/cmillion/gPhoton/tree/master)
```

**`tt_sed`**
: Spectral Energy Distribution from VizieR database.

**`tt_gphoton_lc`**
: Light curve from [`gPhoton.gAperture`](https://github.com/cmillion/gPhoton/blob/master/docs/UserGuide.md#gaperturepy)

**`tt_spectrum`**
: Spectrum table

**`tt_lombscargle`**
: LombScargle results information











