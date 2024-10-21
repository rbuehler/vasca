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

VASCA relies as input on the science data products that missions/organizations create
from their observational raw data. Specifically this means VASCA takes tables of
photometric detections that have a field and visit ID. For a reference of the full list
of required columns can be see here: [{class}`BaseField`](#base_field), [{class}`Region`](#vasca.tables_dict.region)    