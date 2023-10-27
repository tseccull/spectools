# '''_-T-_'''

Developer:    Tom Seccull

This is a personal collection of scripts I that use for reduction, processing, and analysis of astronomical 
spectroscopic data. This repo is intended for storage and documentation, not distribution, so I cannot make any 
guarantees that this software will either work, or be useful to others. Those that choose to use these scripts, 
in whole or in part, do so at their own risk. Support will be very very limited. If, against all odds, any of the
scripts in this repo end up being useful in your project(s), please consider citing it.

The scripts in this repo rely on a lot of other Python packages that deserve recognition if they are used.
Please be sure to cite these as well where necessary:

[Astropy](https://www.astropy.org/acknowledging.html), [NumPy](https://numpy.org/citing-numpy/), [SciPy](https://scipy.org/citing-scipy/)


# scrap.py

v1.0.2

This is essentially a Python wrapper for Curtis McCully's Astroscrappy, which is itself a Python implementation
of Pieter van Dokkum's LA Cosmic. This script is used for detecting, masking, and cleaning cosmic ray hits in
2D spectroscopic data. A modular design is intended to facilitate easy processing of data observed with a variety
of instruments. If Astroscrappy is used, both McCully et al., and van Dokkum should be cited:

[Astroscrappy Docs](https://astroscrappy.readthedocs.io/en/latest/index.html)

[McCully et al. 2018, Astropy/Astroscrappy: v1.0.5 Zenodo Release (v1.0.5). Zenodo](https://doi.org/10.5281/zenodo.1482019)

[van Dokkum 2001, PASP, 113, 1420](https://doi.org/10.1086/323894)

Requires: [Astropy](https://www.astropy.org/), [Astroscrappy](https://doi.org/10.5281/zenodo.1482019), [NumPy](https://numpy.org/), [SciPy](https://scipy.org/)

Supported Instruments: [GMOS-N](https://www.gemini.edu/instrumentation/gmos), [GMOS-S](https://www.gemini.edu/instrumentation/gmos)

# License
All scripts in this repo are licensed under [GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) 

'''_-T-_''' looks like a pyramid to me. I couldn't think of a catchy name.

