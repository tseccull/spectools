#! /home/tom/anaconda3/envs/work/bin/python
"""
stack.py - written by Tom Seccull, 2024-04-18 - v0.0.0

	Last updated: 2024-04-18
	
	This script takes multiple 1D spectra and combines them with a
	weighted bootstrapped median to produce a stacked 1D spectrum with
	reduced noise.
"""


import argparse
import astropy.io.fits as fits
import datetime
import glob
import matplotlib.pyplot as plt
import numpy as np
#import os
#import sys


def scale_spectrum(spectrum, scaling_index):
	scaling_region = spectrum[scaling_index-20:scaling_index+21]
	scaling_median = np.nanmedian(scaling_region)
	scaling_sigma = np.nanstd(scaling_region)
	upper_scaling_limit = scaling_median+(3*scaling_sigma)
	lower_scaling_limit = scaling_median-(3*scaling_sigma)
	scaling_region = scaling_region[
		np.where(
			np.logical_and(
				scaling_region < upper_scaling_limit,
				scaling_region > lower_scaling_limit
			)
		)
	]
	scale = np.median(scaling_region)
	scaled_spectrum = spectrum/scale
	
	return scaled_spectrum, scale


parser = argparse.ArgumentParser(
	description="This script takes multiple 1D spectra and combines\
	them with a weighted median to produce a stacked 1D spectrum with\
	reduced noise."
)
parser.add_argument(
	"-p", "--plot", action="store_true", 
	help="[boolean] Plot the stacked spectrum on screen."
)
parser.add_argument(
	"-sa", "--save", action="store_true",
	help="[boolean] Save the stacked spectrum to a new .fits file."
)
parser.add_argument(
	"-sc", "--scaling_wavelength", default=0.0,
	help="[float] Wavelength at which all spectra are scaled to unity\
	before stacking. Must have the same units as the wavelength axis of\
	the spectrum. The median of the twenty data points both sides of\
	the chosen wavelength (40 points total) will be used as the scaling\
	factor. Defaults to '0.0' which will result in the spectra being\
	scaled to unity at their central wavelength."
)
args = parser.parse_args()

files = sorted(glob.glob("*.fits"))

headers = {}
scaled_optimal_spectra = {}
scaled_aperture_spectra = {}
optimal_scales = {}
aperture_scales = {}

for f in files:
	with fits.open(f) as file_hdu_list:
		headers[f[:-5]] = file_hdu_list[0].header
		optimal_frame = file_hdu_list[0].data
		aperture_frame = file_hdu_list[1].data

	wavelength_axis = optimal_frame[0]
	optimal_spectrum = optimal_frame[1]
	aperture_spectrum = aperture_frame[1]
	
	if args.scaling_wavelength == 0.0:
		args.scaling_wavelength += (wavelength_axis[0]+wavelength_axis[-1])*0.5
	
	scaling_index = np.argmin(np.abs(wavelength_axis-args.scaling_wavelength))
	
	scaled_optimal_spectra[f[:-5]], optimal_scales[f[:-5]] = scale_spectrum(
		optimal_spectrum, scaling_index
	)
	scaled_aperture_spectra[f[:-5]], aperture_scales[f[:-5]] = scale_spectrum(
		aperture_spectrum, scaling_index
	)
