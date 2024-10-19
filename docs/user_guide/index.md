# User Guide

## Intro
Detailed description of the code architecture and central data objects.

## Using `VASCA`
More detailed information on using the package is provided on separate pages,
listed below.

## Supported Instruments

| Instrument | Filters  | Data Description                                                                                                                                                                                                                                                                              |
|------------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| GALEX      | NUV, FUV | [GALEX](http://www.galex.caltech.edu) observations in the far and near ultra violet. Photometric data, processed by the mission pipeline, is accessed from MAST via the astroquery package. A manifest of all visits and all fields con be found here.                                                                        |
| GALEXDS    | NUV      | GALEX drift scan, was a special observation mode that was used at the end of the mission life time. Observations were funded by third parties and most data is not publicly accessible. One exception is data of the GCK survey. However, data must be downloaded locally to work with VASCA. |

## Reference/API
```{toctree}
:maxdepth: 3
:titlesonly:

../api/index
```