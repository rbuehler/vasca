
![VASCA icon](docs/images/VASCA_icon.png)
[![🧪 pytest](https://github.com/rbuehler/vasca/actions/workflows/ci.yml/badge.svg)](https://github.com/rbuehler/vasca/actions/workflows/ci.yml)
[![📚 docs](https://github.com/rbuehler/vasca/actions/workflows/docs.yml/badge.svg)](https://schliwiju.github.io/vasca-mirror/)
[![🚀 pypi](https://github.com/rbuehler/vasca/actions/workflows/pypi.yml/badge.svg)](https://pypi.org/project/vasca/)


# Variable Source Cluster Analysis (VASCA)

1. [Motivation](#motivation)
2. [Pipeline Overview](#pipeline-overview)
3. [Key Features](#key-features)
4. [Proof-of-Principle Study](#proof-of-principle-study)
5. [Documentation and Installation](#documentation-and-installation)
6. [Getting Started](docs/getting_started.md#getting-started)

## Motivation
VASCA is a high-performance software package developed to address the challenges of
time-domain astronomy, especially given the increasing volume of data from large-scale
surveys such as ZTF, LSST, and ULTRASAT. Designed to analyze time-variable cosmic sources
like active galactic nuclei, stars, and transient events, VASCA provides a modular,
scalable solution for integrating data from multiple instruments and conducting a
cohesive analysis.

## Pipeline Overview

The VASCA analysis pipeline consists of three primary steps:
1. **Spatial Clustering**: Associate detections from repeated observations to unique
cosmic sources using mean-shift clustering.
2. **Statistical Variability Detection**: Identify time-variable sources by testing flux
variations against a constant hypothesis at a 5-σ significance level.
3. **Source Classification**: Classify detected sources, including cross-matching with
external catalogs (e.g., SIMBAD, Gaia).

The main output of the pipeline is a catalog of time-variable cosmic
sources, including detailed classifications and cross-matches with existing astronomical
databases.

## Key Features

- **Simplicity and Modularity**: The software uses a hierarchical data model and modular
processing to ensure scalability and ease of use. It supports data from multiple
instruments seamlessly.
- **Proven Algorithms**: VASCA relies on established algorithms and statistical methods,
ensuring robustness and reducing the maintenance burden.
- **Focus on Specific Use Case**: Optimized for analyzing time-domain astronomical data,
VASCA keeps complexity low, simplifying auditing and debugging.
- **Standards Compliance**: Outputs are designed for publication readiness by adhering to
IAU and CDS standards, using widely-accepted, non-proprietary data formats. 
- **Customization and Extensibility**: VASCA allows flexible configuration, making it
adaptable to different datasets and instrument-specific requirements.

## Proof-of-Principle Study

VASCA was applied to a proof-of-principle study  using the Galaxy Evolution Explorer
(GALEX) archive (2003-2013). This study produced a catalog of over 4,000 UV-variable
sources, revealing UV variability across all classes of stars. Notably, a massive,
pulsating white dwarf exhibited unique long-term variability in the UV. The full article
including a description of VASCA's pipeline can be found here:
[The time-variable ultraviolet sky: Active galactic nuclei, stars, and white dwarfs](https://ui.adsabs.harvard.edu/abs/2024A%26A...687A.313B/abstract).

## Documentation and Installation

VASCA is distributed as an open-source package. Comprehensive documentation is available
[here](https://schliwiju.github.io/vasca-mirror/), including example notebooks and an API reference to help users get started.
For quick installation, VASCA can be installed via [PyPI](https://pypi.org/project/vasca/) using:
```shell
pip install vasca
```
For more info see the [installation guide](docs/installation_guide.md#installation).