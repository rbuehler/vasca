# Getting started

## Installation

VASCA is tested to work with Python 3.10 and 3.11. Typically, you want to create a new
Python environment, e.g., with [pyenv-virtualenv](https://github.com/pyenv/pyenv-virtualenv)
or [mamba](https://github.com/mamba-org/mamba):

```bash
mamba create -n vasca python=3.11
mamba activate vasca
```

For standard usage, install VASCA from [PyPi](https://pypi.org/project/vasca/):

```bash
pip install vasca
```
Standard usage covers all pipeline functions for data from instruments that are already
implemented in VASCA. For more information consult the [list](user_guide/instruments.md#supported-instruments)
of supported instruments.

In order to extend VASCA to incorporate another instrument's data, install it directly
from the Github repository.

```bash
pip install -e git+https://github.com/rbuehler/vasca
```
This will ensure that all resources for testing and the jupyter examples are included in
the installation.

## Resource management

Management of the input observational data is handled in VASCA via the
[](#ResourceManager) class. Before first use, users need to edit environment variables
that specify the data storage locations in an `.env` file in VASCA's root directory.

:::{Tip}
The easiest way to set up the resource manager is to duplicate the `.env_template` and
rename it. Then edit your paths to the location of cloud-synced or local directories.
:::

The system is very flexible and can be tailored to your needs. New environment variables
specifying storage locations can be added in the [`resource_envs.yml`](https://github.com/rbuehler/vasca/blob/main/vasca/resource_metadata/resource_envs.yml)
file. New data items that might be required to support additional instruments can be
added in the [`resource_catalog.yml`](https://github.com/rbuehler/vasca/blob/main/vasca/resource_metadata/resource_catalog.yml)
file.

## Running the pipeline and post-processing

To start the pipeline processing go into the `vasca` directory and start the command line
script:
```shell
vasca-pipe <path-to-config.yml>
```
See the [user guide](user_guide/configuration.md#configuration) for more infos on how to
configure the pipeline in the yaml file.

We use [Jupyter Lab](https://github.com/jupyterlab/jupyterlab) for post-processing, with
functional examples provided in `vasca/examples`.

## Coding guidelines

We use the [PEP 8](https://realpython.com/python-pep8/) coding conventions. Before
contributing, please consider the use of automatic code formatting tools. We recommend [`ruff`](https://docs.astral.sh/ruff/).
It is a one-stop-shop that combines all style rules from tools like [`isort`](https://github.com/pycqa/isort),
[`flake8`](https://github.com/PyCQA/flake8), and [`black`](https://black.readthedocs.io/en/stable/#).
All ruff configurations can be found in the `pyproject.toml` file. We set [88 characters](https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html?highlight=88%20#line-length)
as the default line width. The recommended Python version to use is 3.11. For docstrings,
we use the [numpy format](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html).

## Build the documentation

For documentation, we use [Sphinx](https://www.sphinx-doc.org/en/master/). To build it
locally, run the following command from VASCA's `root` directory:
```shell
sphinx-build docs docs/_build
```

To create Unified Modeling Language diagrams, install [pyreverse](https://pylint.pycqa.org/en/latest/pyreverse.html)
and [graphviz](https://graphviz.org/), then run:

```shell
pyreverse vasca -o png -d ./docs/
```
