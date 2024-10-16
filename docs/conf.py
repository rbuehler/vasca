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

# -- Run Autodoc -------------------------------------------------------
#
# def run_apidoc(_) -> None:
#     ignore_paths = []
#
#     argv = [
#         "-f",
#         "-e",
#         "-M",
#         "-T",
#         # "--templatedir",
#         # "_templates",
#         "-o",
#         "api/",
#         "../vasca",
#     ] + ignore_paths
#
#     try:
#         import vasca
#     except ImportError:
#         raise ImportError(
#             "Package must first be installed before creating documentation"
#         )
#
#     from sphinx.ext import apidoc
#
#     apidoc.main(argv)
#
#
# def setup(app) -> None:
#     app.connect("builder-inited", run_apidoc)


# -- Project information -----------------------------------------------------
#

project = "VASCA"
copyright = "BSD 3-Clause License"
author = "Rolf Buehler and Julian Schliwinski"

# The full version, including alpha/beta/rc/dev tags
release = get_version("vasca")
version = ".".join(release.split(".")[:3])
if release.split(".")[3].startswith("dev"):
    version += "dev"


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
# autodoc_default_options = {
#     "members": True,
#     "member-order": "bysource",
#     "special-members": "__init__",
# }

# -- Options for Autodoc --------------------------------------------------------------

autosummary_generate = True

# autodoc_member_order = "bysource"
# autodoc_preserve_defaults = True
# autodoc_typehints = "description"

# numpydoc_show_class_members = False
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

# html_theme = "sphinx_rtd_theme"
# html_theme = "pydata_sphinx_theme"
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
