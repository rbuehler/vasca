---
jupytext:
  hide_notebook_metadata: true
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: vasca-github
  language: python
  name: vasca-github
---

```{code-cell}
:tags: [remove-cell]

# ruff: noqa: T201
```

# Pipeline

This is a tutorial showcasing VASCA's pipeline flow on a very simple example. The goal
is to create a VASCA [](#Region) from two [](#GALEXField) for which we download the
raw data online from [MAST](https://astroquery.readthedocs.io/en/latest/mast/mast.html).
We apply quality cuts and do source clustering followed by variability analysis and
finally source cross-matching.

+++

## General Configuration

The standard pipeline processing starts by reading a yaml file. To keep this tutorial
simple, we are going to introduce parts of the configuration step by step at the point
where they are required in the pipeline flow.

```{note}
An important premise of the configuration is that each parameter needs to be
configured explicitly. This means even default values are specified all the time. This
is a design decision purposefully made in order to ensure transparent and complete
configurations. As a result, all possible parameters are always included when looking
at configuration file.
```

Let's begin with the ``general`` section. Here, basic information and functionality is
configured. The ``name`` of the pipeline run specifies also the name of directory in
which all results will be stored. The location of output directory is at ``out_dir_base``
relative to the root directory of the package.

VASCA uses the powerful logging system provided by [logguru](https://loguru.readthedocs.io/en/stable/index.html).
The configuration specifies the [``log_level``](https://loguru.readthedocs.io/en/stable/api/logger.html#loguru._logger.Logger.level),
which we are going to set to debugging mode here. By default VASCA is going to save
all logging messages in a file stored in the output directory. ``log_file`` specifies
the name of that file, while ``default`` tells the pipeline to use a default name.

Parallel processing of the field-level analysis can be enabled when setting the number
of CPU threads ``nr_cpus > 1``.

VASCA can include field-averaged reference data, if such data is available additional
to the visit-level data from the instruments mission pipeline. To save memory/storage
and computation time it is configurable wether to include reference sources in the
final [](#Region)-file (``save_ref_srcs``) and to repeat already processed fields that
are included in the region (``run_fields``).

```{code-cell}
# Dictionary holding the configuration
config = {}

# General section of the configuration
config["general"] = {
    "name": "simple_pipe",
    "out_dir_base": "vasca_pipeline",
    "log_level": "DEBUG",
    "log_file": "default",
    "nr_cpus": 1,
    "save_ref_srcs": True,
    "run_fields": True,
}
```

```{code-cell}
# To be continued ...
```
