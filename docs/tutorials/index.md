---
file_format: mystnb
kernelspec:
    name: vasca-github
---
# Tutorials
In this section various tutorials are provided. This is a markdown file and this section
is still under development.

## Intro
All tutorials are jupyter-based. This documentation uses [`myst-nb`](https://myst-nb.readthedocs.io/en/latest/index.html)
and [`jupytext`](https://jupytext.readthedocs.io/en/latest/index.html) to execute and
render the content.

Additionally you can find various notebooks on post-processing VASCA's pipeline results
[here](https://github.com/rbuehler/vasca/blob/main/vasca/examples):

| Notebook                                                                                                                      | Description                                                      |
|-------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------|
| [vasca_pipe_post_process](https://github.com/rbuehler/vasca/blob/main/vasca/examples/vasca_pipe_post_process.ipynb)           | Post-process VASCA pipline results                               |
| [inspect_source_distributions](https://github.com/rbuehler/vasca/blob/main/vasca/examples/inspect_source_distributions.ipynb) | Show detailed distribution of the main VASCA selection variables |
| [select_sources](https://github.com/rbuehler/vasca/blob/main/vasca/examples/select_sources.ipynb)                             | Select VASCA sources based on chosen criteria                    |
| [inspect_sources](https://github.com/rbuehler/vasca/blob/main/vasca/examples/inspect_sources.ipynb)                           | Inspect source light curves and sky maps                         |
| [inspect_matches](https://github.com/rbuehler/vasca/blob/main/vasca/examples/inspect_matches.ipynb)                           | Analyse SIMBAD and Gaia counterparts                             |
| [inspect_match_distance](https://github.com/rbuehler/vasca/blob/main/vasca/examples/inspect_match_distance.ipynb)             | Compare match distance distribution to random                    |
| [create_pubcat](https://github.com/rbuehler/vasca/blob/main/vasca/examples/create_pubcat.ipynb)                               | Prepare a VASCA region catalog for upload to CDS                 |
| [inspect_artifacts](https://github.com/rbuehler/vasca/blob/main/vasca/examples/inspect_artifacts.ipynb)                       | Visualize artifacts and the detections on them                   |

## List of Tutorials
```{toctree}
:maxdepth: 2
:titlesonly:

tutorial_rm
tutorial_field
tutorial_pipe
simple_example
test_table
```