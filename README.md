# spectools

v24.194.0

Developer:    Tom Seccull

This is a personal collection of scripts I that use for reduction, processing, 
and analysis of astronomical spectroscopic data. I cannot make any guarantees 
that this software will either work, or be useful to others. Updates to these 
scripts may arrive in this repo randomly and without warning. Those that 
choose to use these scripts do so at their own risk. At very least, the 
shebang Python path at the top of each script should be updated to match the 
Python path in your system. Support from me may be limited. If, against all 
odds, any of the scripts in this repo end up being useful in your work, please 
consider citing the repo with the details provided in the CITATION.cff file. A
Zenodo DOI may be obtained for this repo at some point in the future.

The scripts in this repo rely on other Python packages that deserve recognition 
if they are used.
Please be sure to cite these as well where necessary:

[Astropy](https://www.astropy.org/acknowledging.html), [Matplotlib](https://matplotlib.org/stable/users/project/citing.html), [NumPy](https://numpy.org/citing-numpy/), [SciPy](https://scipy.org/citing-scipy/)


# scrap.py

v1.0.6

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

v1.0.3

This script is designed to correct fringing in 2D spectroscopic data by 
creating a median fringe frame and subtracting it from science data. A modular 
design is intended to facilitate easy processing of data observed with a 
variety of instruments.

Requires: [Astropy](https://www.astropy.org/), [NumPy](https://numpy.org/)\
Supported Instruments: [GMOS-N](https://www.gemini.edu/instrumentation/gmos), [GMOS-S](https://www.gemini.edu/instrumentation/gmos)


# extinct.py

v1.0.1

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

Requires: [Astropy](https://www.astropy.org/), [Matplotlib](https://matplotlib.org/stable/users/project/citing.html), [NumPy](https://numpy.org/)\
Supported Instruments: [GMOS-N](https://www.gemini.edu/instrumentation/gmos), [GMOS-S](https://www.gemini.edu/instrumentation/gmos)

Extinction Curve Sources:\
GMOS-N - [See the Gemini Observatory web pages.](https://www.gemini.edu/observing/telescopes-and-sites/sites#Extinction)\
GMOS-S - [Stone & Baldwin 1983, MNRAS, 204, 347](https://doi.org/10.1093/mnras/204.2.347)
    

# stack.py

v1.0.2

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
the uncertainty each stacked data point. stack.py is only readily compatible 
with spectra extracted by [MOTES](https://github.com/tseccull/motes).

Requires: [Astropy](https://www.astropy.org/), [Matplotlib](https://matplotlib.org/stable/users/project/citing.html), [NumPy](https://numpy.org/)\
Supported Instruments: [GMOS-N](https://www.gemini.edu/instrumentation/gmos), [GMOS-S](https://www.gemini.edu/instrumentation/gmos)


# divide.py

v1.0.2

This script divides one 1D spectrum by another. divide.py expects both spectra
to have wavelength axes of equal length and be the product of both extraction by MOTES
and stacking by stack.py. Typical use of this script is to calibrate a spectrum
of a minor planet with that of a solar twin or solar analog to derive the minor
planet's reflectance spectrum.

Requires: [Astropy](https://www.astropy.org/), [Matplotlib](https://matplotlib.org/stable/users/project/citing.html), [NumPy](https://numpy.org/)

# bingrad.py - under development

v0.0.6

This script has two functions. Primarily it is used to bin spectroscopic data 
to boost its signal-to-noise ratio at the expense of spectral resolution. A 
binned spectrum can be plotted and saved to a new FITS file. The secondary
function is to allow the continuum gradient of the spectrum to be measured 
across a user-defined wavelength range via linear regression. The resulting 
linear fit can be plotted and its parameters will be printed in the terminal.

# License
All scripts in this repo are licensed under [GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) 
