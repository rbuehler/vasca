# User Guide

This guide aims to introduce VASCA's basic functionality and central building blocks
of the package's internals. Good understanding of the latter will help expanding VASCA
to more instruments. For more detailed guides and descriptions see the respective [API reference](../api/index.rst),
[tutorials](../tutorials/tutorial_pipe.md) and [Jupyter examples](https://github.com/rbuehler/vasca/tree/main/vasca/examples)
on post-processing listed [here](../tutorials/index.md).

## Intro

VASCA is based on a very simple data model as you can see illustrated by the image below.
Photometric detections from repeated observations (visits) are taken as input. They must
be associated to a uniquely defined patch on the sky, a **_field_**. These detections are
then clustered on the field level to form uniquely defined **_sources_**. Multiple fields
are combined to form a **_region_**. Overlapping fields undergo a second clustering step
to assure the uniqueness of sources also on the region level. The final output of VASCA,
the variable source catalog, is based on variability detection, cross-matching and
classification of the region-level sources.

:::{figure-md} data-model
<img src="../images/VASCA_data_model_v2.jpg" alt="data_model" class="bg-primary mb-1" width="400px">

The VASCA data model.
:::

On the implementation side, VASCA provides three main data objects that inherit VASCA's
base data structure, the [](#TableCollection):

[](#Source)
: Unique cosmic source with an ID, sky coordinates, and a (multi-wavelength) [light curve](https://en.wikipedia.org/wiki/Light_curve)
from which variability parameters are computed in addition to possible IDs associating
the source to known objects from external catalogs.

[](#BaseField)
: Unique field defined on a sky area which holds all sources including their multi-visit
detections. Note that a field is uniquely specific in the instrument and the filter
(band-pass) it used to make the observation.

[](#Region)
: Region on the sky composed of multiple fields where sources are combined from
observations in different filters and instruments.

:::{admonition} Instrument independence
:class: tip
VASCA's instrument independence is largely owed to the fact that [](#BaseField) is indeed
specific to a given instrument and filter. This allows to treat field-level processing in
parallel and allows instrument- and filter-specific configuration of the pipeline.
:::

:::{admonition} Expanding VASCA to new instruments 
:class: tip
The only task required to expand VASCA to operate on data for a new instrument is to
create its own field class. This is necessary in order to map observational parameters
like flux, spatial coordinates and their uncertainties to the corresponding parameters in
VASCA. Fields my include co-added sky maps, or intensity images, as reference images and
also visit-level images. The field class must also provide a {meth}`load` method which
handles data I/O together with the [](#ResourceManager).
:::

All objects inheriting from [](#TableCollection) can be written to storage as [FITS](https://en.wikipedia.org/wiki/FITS)
files. These hold images and tables compatible with the FITS version 4.0 standard so that
users my use tools like [DS9](https://sites.google.com/cfa.harvard.edu/saoimageds9) and
[TOPCAT](https://www.star.bristol.ac.uk/mbt/topcat/) for data exploration and debugging.
To see what kind of Tables are stored in VASCA's data structures, see the glossary [here](data_model.md#vasca-tables).

## Using `VASCA`
More detailed information on using the package is provided on separate pages,
listed below.

```{toctree}
:maxdepth: 3
:titlesonly:
data_model
pipeline
post_processing
configuration
instruments
../api/index
```