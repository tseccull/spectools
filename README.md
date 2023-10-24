# spectools

Developer:    Tom Seccull

Last Updated: 2023-10-24

A collection of scripts I that use for reduction, processing, and analysis of astronomical spectroscopic data.
This repo is not intended for distribution, so I cannot make any guarantees that this software will work,
or be useful to others. Those that choose to use these scripts, in whole or in part, do so at their own risk.
Support will not be provided.

If, against all odds, the software in this repo is helpful to your project in any way, please consider citing it
as shown in this repo's CITATION.cff file.


# scrap.py v1.0.0
This is essentially a Python wrapper for Curtis McCully's Astroscrappy, which is itself a Python implementation
of Pieter van Dokkum's LA Cosmic. This script is used for detecting, masking, and cleaning cosmic ray hits in
2D spectroscopic data. A modular design is intended to facilitate easy processing of data observed with a variety
of instruments. If Astroscrappy is used, both McCully et al., and van Dokkum should be cited:
[Astroscrappy Docs](https://astroscrappy.readthedocs.io/en/latest/index.html)
[McCully et al. 2018, Astropy/Astroscrappy: v1.0.5 Zenodo Release (v1.0.5). Zenodo](https://doi.org/10.5281/zenodo.1482019)
[van Dokkum 2001, PASP, 113, 1420](https://doi.org/10.1086/323894)

Requires: Astropy, Astroscrappy, NumPy, SciPy
