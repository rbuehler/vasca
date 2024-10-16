# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from importlib.metadata import version as get_version

sys.path.insert(0, os.path.abspath(".."))


# -- Project information -----------------------------------------------------
#

project = "VASCA"
copyright = "BSD 3-Clause License"
author = "Rolf Buehler and Julian Schliwinski"

# The full version, including alpha/beta/rc/dev tags
release = get_version("vasca")
print(release)
print(type(release))
try:
    version = ".".join(release.split(".")[:3])
    if release.split(".")[3].startswith("dev"):
        version += "dev"
except:
    version = "1.0.13"


# -- General configuration ---------------------------------------------------
#

extensions = [
    # Sphinx
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    # External
    "myst_parser",
    "autodoc2",
    "sphinx_copybutton",
]

# -- Options for Autodoc --------------------------------------------------------------

autodoc2_packages = [
    {
        "path": "../vasca",
        "exclude_dirs": ["examples", "test"],
    }
]
autodoc2_output_dir = "api"

# -- Options for intersphinx -------------------------------------------------
#
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "loguru": ("https://loguru.readthedocs.io/en/stable", None),
    "astropy": ("http://docs.astropy.org/en/stable", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"
html_title = f"VASCA v{version}"

html_theme_options = {
    "source_repository": "https://github.com/rbuehler/vasca/",
    "source_branch": "main",
    "source_directory": "docs/",
}


# -- Options for Markdown files ----------------------------------------------
#

myst_heading_anchors = 3

# -- Options for coppybutton ----------------------------------------------
#

copybutton_exclude = ".linenos, .gp, .go"
