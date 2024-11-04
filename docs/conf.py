# Configuration file for the Sphinx documentation builder.
#

# -- Path setup --------------------------------------------------------------
#
import os
import sys
import jupytext
from importlib.metadata import version as get_version

sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
#
project = "VASCA"
copyright = "BSD 3-Clause License"
author = "Rolf Buehler and Julian Schliwinski"

# The full version, including alpha/beta/rc/dev tags
release = get_version("vasca")
try:
    version = ".".join(release.split(".")[:3])
    if release.split(".")[3].startswith("dev"):
        version += "dev"
except Exception as e:
    print(f"Failed to parse package version from `{release}`. Exception:\n {e}")  # noqa:T201
    version = "0.0.1"


# -- Convert py:percent files to markdown notebooks -----------------------------------
#
def convert_py_to_md() -> None:
    tutorials_dir = os.path.join(os.path.dirname(__file__), "tutorials")
    for filename in os.listdir(tutorials_dir):
        if filename.endswith(".py"):
            filepath = os.path.join(tutorials_dir, filename)
            md_filepath = filepath.replace(".py", ".md")
            jupytext.write(jupytext.read(filepath), md_filepath, fmt="md:myst")
            print(f"Converted {filepath} to {md_filepath}")


# Hook to run before the build process starts
def setup(app):
    app.connect("builder-inited", lambda app: convert_py_to_md())


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
    "myst_nb",
    "autodoc2",
    "sphinx_copybutton",
    "sphinx_tippy",
]

exclude_patterns = ["jupyter_execute", ".jupyter_cache"]

# -- Options for Autodoc --------------------------------------------------------------

autodoc2_packages = [
    {
        "path": "../vasca",
        "exclude_dirs": ["examples", "test"],
        "exclude_files": ["_version.py"],
    }
]
autodoc2_output_dir = "api"

# -- Options for intersphinx -------------------------------------------------
#
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "loguru": ("https://loguru.readthedocs.io/en/stable", None),
    "astropy": ("https://docs.astropy.org/en/stable", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

# -- Options for HTML output -------------------------------------------------
#
html_theme = "furo"
html_title = f"VASCA v{version}"
html_logo = "images/VASCA_icon.png"

html_theme_options = {
    "source_repository": "https://github.com/rbuehler/vasca/",
    "source_branch": "main",
    "source_directory": "docs/",
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/rbuehler/vasca",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
}

tippy_enable_mathjax = True
tippy_anchor_parent_selector = "div.content"

# -- Options for Markdown files ----------------------------------------------
#

myst_heading_anchors = 3
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "html_image",
    "attrs_block",
]
nb_execution_mode = "cache"
nb_execution_timeout = -1

# -- Options for coppybutton ----------------------------------------------
#

copybutton_exclude = ".linenos, .gp, .go"
