#! /home/tom/anaconda3/bin/python
"""
bingrad.py = written by Tom Seccull, 2024-07-05 - v0.0.0
	
	Last updated: 2024-07-05
	
	This script has two functions. Primarily it is used to bin spectroscopic 
	data to boost its signal-to-noise ratio at the expense of spectral 
	resolution. A binned spectrum can be plotted and saved to a new FITS file. 
	The secondary function is to allow the continuum gradient of the spectrum 
	to be measured across a user-defined wavelength range via linear 
	regression. The resulting linear fit can be plotted and its parameters 
	will be printed in the terminal.
"""

import argparse
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import linregress

parser = argparse.ArgumentParser(
	description="This script has two functions. Primarily it is used to bin\
	spectroscopic data to boost its signal-to-noise ratio at the expense of\
	spectral resolution. A binned spectrum can be plotted and saved to a new\
	FITS file. The secondary function is to allow the continuum\
	gradient of the spectrum to be measured across a user-defined wavelength\
	range via linear regression. The resulting linear fit can be plotted and\
	its parameters will be printed in the terminal."
)
parser.add_argument(
	"data_file", type=str,
	help="Input spectrum FITS file name."
)
parser.add_argument(
	"factor", type=int,
	help="Binning factor, a.k.a. bin width or number of points per bin."
)
parser.add_argument(
	"-r", "--rejection", type=float,
	help="This optional argument will be multiplied by the standard deviation\
	of the combined distribution of points present within each bin after each\
	original point is resampled within its Gaussian uncertainties. The\
	resulting number will be set as both the positive and negative distance\
	from the median of the distribution beyond which points will be rejected\
	from the final estimation of the median value of the bin. If this\
	argument is not set, the default behaviour of bingrad is to perform no\
	rejection."
)
parser.add_argument(
	"-g" ,"--gradient_wavelength_ranges",
	help="This is a comma-separated list of wavelengths that bound regions\
	where the spectral gradient will be measured from the binned spectrum.\
	Only one gradient measurement is done each time bingrad is run, but\
	multiple wavelength regions can be set to allow parts of the spectrum to\
	be skipped in the measurement if needed. The units of -gwr should match\
	those of the wavelength axis of the input spectrum. If np gradient\
	wavelength ranges are set, bingrad defaults to skipping the gradient\
	measurement."
)
parser.add_argument(
	"-p", "--plot", action="store_true",
	help="Plot the binned spectrum on screen. If a gradient has also been\
	measured, the linear fit to the chosen spectral regions will also be\
	presented."
)
parser.add_argument(
	"-s", "--save", action="store_true",
	help="Save the binned spectrum to a new FITS file."
)
args = parser.parse_args()


with fits.open(args.data_file) as data_file:
	headers = [x.header for x in data_file]
	frames = [x.data for x in data_file]
	
primary_head = headers[0]
primary_frame = frames[0]
qual = primary_frame[5]

# Cut off the qual==0 ends of the spectrum.
i = 0
while primary_frame[5][i] == 0: i += 1
j = np.shape(primary_frame)[1] - 1
while primary_frame[5][j] == 0: j -= 1
primary_frame = primary_frame[:, i:j+1]

# Make sure the number of unbinned points is even by removing the first or the
# last point in the spectrum based on whichever is less certain. Having an even
# number of data points in the spectrum simplifies the process of binning
# outward from the centre of the spectrum.
if np.shape(primary_frame)[1]%2 != 0:
	if primary_frame[2][0] > primary_frame[2][-1]:
		primary_frame = primary_frame[:, 1:]
	else:
		primary_frame = primary_frame[:, :-1]

wavelength_axis = primary_frame[0]
optimal_spectrum = primary_frame[1]
optimal_errors = primary_frame[2]
aperture_spectrum = primary_frame[3]
aperture_errors = primary_frame[4]
qual = primary_frame[5]

print(len(wavelength_axis)/args.factor, len(wavelength_axis)%args.factor)
