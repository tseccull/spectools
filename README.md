# spectools

![GitHub Release](https://img.shields.io/github/v/release/tseccull/spectools)
[![DOI](https://zenodo.org/badge/709271938.svg)](https://zenodo.org/doi/10.5281/zenodo.12786056)

Author: Tom Seccull

This is a personal collection of scripts that I use for reduction, processing, 
and analysis of astronomical spectroscopic data. Although these scripts are 
still undergoing active development, this repo is intended more as a citable 
record than a public software distribution. By all means feel free to use any 
of these scripts, either in whole or in part, but note that you do so at your 
own risk. I cannot guarantee that this software will either work or be 
beneficial to you. If you download a copy of this repo, be sure to select the
most recent stable release at the link on the right.  

If, against all odds, any of the scripts in this repo end up helping you in 
your work, please consider citing the repo with the details provided via the 
Zenodo DOI above.

The scripts in this repo rely on multiple Python packages that deserve 
recognition if they are used. Please be sure to cite them where necessary:

[Astropy](https://www.astropy.org/acknowledging.html), [DRAGONS](https://www.gemini.edu/observing/phase-iii/reducing-data/dragons-data-reduction-software), [Matplotlib](https://matplotlib.org/stable/users/project/citing.html), [NumPy](https://numpy.org/citing-numpy/), [SciPy](https://scipy.org/citing-scipy/)


# dagrons.py 

v1.0.4

This script is a partial reduction pipeline for GMOS longslit spectroscopic 
data built with DRAGONS and based on its GMOS longslit spectroscopy reduction 
tutorial. Some aspects of the reduction have been automated to some degree, but
the bones are the same as what is presented in the DRAGONS documentation. The
recipe called by dagrons.py for reduction of science data and where it should
be pasted in DRAGONS is provided in `./dragons_reference_recipes/recipes_LS_SPECT.py`. 
The directory where this script is initiated is taken to be the working 
directory, with the directory containing data and calibrations for the 
reduction provided as an argument. The data and calibrators for only one on-sky
source are expected to be present in the data directory, such that all the data
and calibration frames have the same format/RoI. This means science targets and
their specphot standard stars may have to be stored and reduced separately, but
this script is less complex as a result. `dagrons.py` only calls a subset of 
the primitives provided by the full `reduceScience()` recipe that DRAGONS 
normally runs for GMOS longslit spectroscopic data including preparation, DQ 
and VAR frame addition, overscan correction, bias subtraction, ADU to e-
converion, flat-field correction, QE correction, and 2D spectrum distortion
correction (rectification). Cosmic ray flagging, fringe subtraction, sky 
subtraction, extraction, and stacking are all performed later by other scripts
in the spectools repo. DRAGONS should be cited if `dagrons.py` is used.

Requires: [DRAGONS](https://www.gemini.edu/observing/phase-iii/reducing-data/dragons-data-reduction-software)\
Supported Instruments: [GMOS-N](https://www.gemini.edu/instrumentation/gmos), [GMOS-S](https://www.gemini.edu/instrumentation/gmos)  


# scrap.py

v1.0.9

This is essentially a Python wrapper for Astropy's Astroscrappy, which is 
itself a Python implementation of Pieter van Dokkum's LA Cosmic. This script 
is used for detecting, masking, and cleaning cosmic ray hits in 2D 
spectroscopic data. A modular design is intended to facilitate easy processing 
of data observed with a variety of instruments. If Astroscrappy is used, both 
McCully et al., and van Dokkum should be cited:

[Astroscrappy Docs](https://astroscrappy.readthedocs.io/en/latest/index.html)

[McCully et al. 2018, Astropy/Astroscrappy: v1.0.5 Zenodo Release (v1.0.5). Zenodo](https://doi.org/10.5281/zenodo.1482019)

[van Dokkum 2001, PASP, 113, 1420](https://doi.org/10.1086/323894)

Requires: [Astropy](https://www.astropy.org/), [Astroscrappy](https://doi.org/10.5281/zenodo.1482019), [NumPy](https://numpy.org/), [SciPy](https://scipy.org/)\
Supported Instruments: [GMOS-N](https://www.gemini.edu/instrumentation/gmos), [GMOS-S](https://www.gemini.edu/instrumentation/gmos)


# fronge.py

v1.0.5

This script is designed to correct fringing in 2D spectroscopic data by 
creating a median fringe frame and subtracting it from science data. A modular 
design is intended to facilitate easy processing of data observed with a 
variety of instruments.

Requires: [Astropy](https://www.astropy.org/), [NumPy](https://numpy.org/)\
Supported Instruments: [GMOS-N](https://www.gemini.edu/instrumentation/gmos), [GMOS-S](https://www.gemini.edu/instrumentation/gmos)


# extinct.py

v1.0.4

This script applies an average extinction correction to input 1D astronomical 
spectroscopic data based on the measured atmospheric extinction curve of the 
observing site. The spectrum and its uncertainties are each multiplied by 
10 ^ (0.4 * airmass * k(lambda)), where airmass is the median airmass at which 
the spectrum was observed and k(lambda) is the wavelength dependent extinction 
curve at the observing site given in magnitudes per unit airmass. Note that
the extinction curves used here are averages that do not account
for variable atmospheric extinction due to variable concentrations of 
atmospheric water vapour or scattering particles. When comparing or calibrating
one corrected spectrum with another it is assumed that atmospheric conditions
were stable across the consecutive observations of both targets. This script
does not create new files, but instead just updates the input files. If this
script is used, be sure to cite the appropriate article or link for the
extinction curve. `extinct.py` is only readily compatible with spectra
extracted by [MOTES](https://github.com/tseccull/motes), and it should be run
on individual 1D spectra of a given target prior to stacking them.

Requires: [Astropy](https://www.astropy.org/), [Matplotlib](https://matplotlib.org/stable/users/project/citing.html), [NumPy](https://numpy.org/), [SciPy](https://scipy.org/)\
Supported Instruments: [GMOS-N](https://www.gemini.edu/instrumentation/gmos), [GMOS-S](https://www.gemini.edu/instrumentation/gmos)

Extinction Curve Sources:\
GMOS-N - [Buton et al. 2013, A&A, 549, A8](https://doi.org/10.1051/0004-6361/201219834)\
-------- see [Gemini website](https://www.gemini.edu/observing/telescopes-and-sites/sites#Extinction) for the 3100 angstrom point.\
GMOS-S - [Stritzinger et al. 2005, PASP, 117, 810](https://doi.org/10.1086/431468)
    

# stack.py

v1.0.4

This script takes multiple 1D spectra and combines them to produce a stacked 1D
spectrum with reduced noise. All input spectra are scaled to unity at either a
user-selected wavelength or the central wavelength of the spectra. The spectra 
are then all randomly resampled within their uncertainties to estimate the 
total combined distribution of possible values for each wavelength in the 
stacked spectrum. The uncertainty distribution of each data point in each
spectrum is assumed to be Gaussian with the mean defined by the value of 
the data point and the standard deviation defined by its uncertainty. For the 
resulting distribution of resampled points at each wavelength element, its 
mean and median are found to be within one standard error of the mean from each
other in almost all cases. Given its close proximity to the mean and its 
robustness against outliers, the median of the combined distribution in each 
wavelength element is taken as the value of the stacked spectrum at that 
wavelength. The standard error of the mean of the distribution is taken to be 
the uncertainty each stacked data point. `stack.py` is only readily compatible 
with spectra extracted by [MOTES](https://github.com/tseccull/motes).

Requires: [Astropy](https://www.astropy.org/), [Matplotlib](https://matplotlib.org/stable/users/project/citing.html), [NumPy](https://numpy.org/)\
Supported Instruments: [GMOS-N](https://www.gemini.edu/instrumentation/gmos), [GMOS-S](https://www.gemini.edu/instrumentation/gmos)


# stax.py

v1.1.1

`stax.py` works similarly to `stack.py` but only accepts longslit SpeX 
spectra observed in prism mode and reduced with [Spextool](https://irtfweb.ifa.hawaii.edu/~spex/Spextool.pdf).
The products of `stax.py` can be fed into scripts further down the reduction 
chain including `divide.py`, `bingrad.py`, `fits2dat.py`, and `qual.py`.

Unlike `stack.py`, `stax.py` can stack spectra in two different ways:

Default "Bootstrap Median" Operation:
See `stack.py` description; this mode works in an identical way to that.

"Summing" Operation":
For some bright targets it's necessary to take very short exposures to
avoid saturating the sensor. Exposures shorter than ~ 5 s are short
enough that frame-to-frame variation in the spectral continnum can be
caused by the movement of the PSF center across the slit that results
from tip/tilt atmospheric distortion. This especially needs to be
considered if the telescope used has no facility for tip/tilt 
correction. In summing mode the the spectra are stacked by simply
summing them all to average all spectra over a longer virtual
integration time that allows for the effects of tip/tilt atmospheric
turbulence to be averaged out. This should mitigate the effects of
slit losses on the stacked spectrum provided enough frames are stacked
to sum the total integration time up to > 5 s. Even longer total 
integration times should give better results. Uncertainties of each
stacked spectrum are summed in quadrature.   

Requires: [Astropy](https://www.astropy.org/), [Matplotlib](https://matplotlib.org/stable/users/project/citing.html), [NumPy](https://numpy.org/)\
Supported Instruments: [SpeX](https://irtfweb.ifa.hawaii.edu/~spex/index.html)


# divide.py

v1.0.6

This script divides one 1D spectrum by another. `divide.py` expects both spectra
to have wavelength axes of equal length and be the product of both extraction by MOTES
and stacking by `stack.py`. Typical use of this script is to calibrate a spectrum
of a minor planet with that of a solar twin or solar analog to derive the minor
planet's reflectance spectrum.

Requires: [Astropy](https://www.astropy.org/), [Matplotlib](https://matplotlib.org/stable/users/project/citing.html), [NumPy](https://numpy.org/)


# bingrad.py

v1.0.2

This script has two functions. Primarily it is used to bin spectroscopic data 
to boost its signal-to-noise ratio at the expense of spectral resolution. A 
binned spectrum can be plotted and/or saved to a new FITS file. The secondary
function is to allow the continuum gradient of the spectrum to be measured 
across a user-defined wavelength range via linear regression. The resulting 
linear fit can be plotted and its parameters will be printed in the terminal.
`bingrad.py` expects an input spectrum in the format produced by `divide.py`. 

Requires: [Astropy](https://www.astropy.org/), [Matplotlib](https://matplotlib.org/stable/users/project/citing.html), [NumPy](https://numpy.org/), [SciPy](https://scipy.org/)


# qual.py

v1.0.0

This script allows the user to manually set the quality flags for data points 
in a spectrum extracted with [MOTES](https://github.com/tseccull/motes).

Requires: [Astropy](https://www.astropy.org/), [Matplotlib](https://matplotlib.org/stable/users/project/citing.html), [NumPy](https://numpy.org/)


# fits2dat.py

v1.0.1

This script repackages 1D spectra extracted with [MOTES](https://github.com/tseccull/motes) 
and saved in .fits format into files in a csv format. Header metadata is stored
in at the top of the file with each line starting "#".

Requires: [Astropy](https://www.astropy.org/), [NumPy](https://numpy.org/)


# License
All scripts in this repo are licensed under [GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) 
