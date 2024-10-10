<!-- start installation -->
## Installation

Typically, you want to create a new Python 3.10 environment, e.g., with conda or mamba:

```bash
mamba create -n vasca python=3.10
mamba activate vasca
```

The installation steps are:

1. Clone this repository:

```bash
git clone git@github.com:rbuehler/vasca.git
```

2. Install the `VASCA` package:

```bash
cd ./vasca
pip install -e .
```

3. Set up the resource manager by including your cloud path location in the `.env_template` file and rename it to `.env`.
<!-- end installation -->

## Running the pipeline and post-processing

We use [Jupyter Lab](https://github.com/jupyterlab/jupyterlab) for post-processing, with functional examples provided in `vasca/examples`.

## Coding guidelines

We use the [PEP 8](https://realpython.com/python-pep8/) coding conventions. Before contributing, please consider the use of automatic code formatting tools like [isort](https://github.com/pycqa/isort), [flake8](https://github.com/PyCQA/flake8), and [black](https://black.readthedocs.io/en/stable/#). We set [88 characters](https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html?highlight=88%20#line-length) as the default line width. The recommended Python version to use is 3.10.x. For docstrings, we use the [numpy format](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html).

## Build the documentation

For documentation, we use [SPHINX](https://www.sphinx-doc.org/en/master/). To build it, run the following:

```bash
sphinx-apidoc -f -o docs vasca
cd docs/
make html
```

To create Unified Modeling Language diagrams, install [pyreverse](https://pylint.pycqa.org/en/latest/pyreverse.html) and [graphviz](https://graphviz.org/), then run:

```bash
pyreverse vasca -o png -d ./docs/
```
