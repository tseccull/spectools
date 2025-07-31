#!/usr/bin/env python3

"""
	stax.py

	Copyright (C) 2025-07-29 Tom Seccull
	
	This script is part of the spectools repo hosted at 
	https://github.com/tseccull/spectools
	https://doi.org/10.5281/zenodo.12786056
	
	If used, please cite the spectools DOI above.
	
	This script is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.

	Last updated - 2025-07-31

	Description --------------------------------------------------------
	This script is similar to stack.py, but handles spectra observed 
	with SpeX at the NASA Infrared Telescope Facility and then reduced 
	with Spextool. Like stack.py, it takes multiple 1D spectra and 
	combines them to produce a stacked 1D spectrum with reduced noise. 
	
	- Default "Bootstrap Median" Operation -  
	All input spectra are scaled to unity at either a user-selected 
	wavelength or the central wavelength of the spectra. The spectra are
	then all randomly resampled within their uncertainties to estimate
	the total combined distribution of possible values for each 
	wavelength in the stacked spectrum. The uncertainty distribution of
	each data point in each spectrum is assumed to be Gaussian with the
	mean defined by the value of the data point and the standard
	deviation defined by its uncertainty. For the resulting distribution
	of resampled points at each wavelength element, its mean and median
	are found to be within one standard error of the mean from each
	other in almost all cases. Given its close proximity to the mean and
	its robustness against outliers, the median of the combined
	distribution in each wavelength element is taken as the value of the
	stacked spectrum at that wavelength. The standard error of the mean
	of the distribution is taken to be the uncertainty each stacked
	data point.
	
	- "Summing" Operation" -
	For some brigh targets its necessary to take very short exposures to
	avoid saturating the sensor. Exposures shorter than ~ 5 s are short
	enough that frame-to-frame variation in the spectral continnum can
	be caused by the movement of the PSF center across the slit that
	results from tip/tilt atmospheric distortion. This especially needs
	to be considered if the telescope used has no facility for tip/tilt
	correction. In summing mode the the spectra are stacked by simply
	summing them all to average all spectra over a longer virtual
	integration time that allows for the effects of tip/tilt
	atmospheric turbulence to be averaged out. This should mitigate the
	effects of slit losses on the stacked spectrum provided enough frames
	are stacked to sum the total integration time up to > 5 s. Even 
	longer total integration times should give better results.
	Uncertainties of each stacked spectrum are summed in quadrature.   
"""

__author__ = "Tom Seccull"
__version__ = "1.1.0"

import argparse
import astropy.io.fits as fits
import copy
import datetime
import glob
import stax.spex_stax_head as spexsh
import matplotlib.pyplot as plt
import numpy as np


def bootstrap_median(files, scaling_wavelength):
	"""
	This function is used to stack the spectrum in the default
	bootstrapped median mode of operation.
	
	Args:
	 -- files (list)
			List of filenames for files containing spectra to be
			stacked.
	 -- scaling_wavelength (float)
			Wavelength at which the stacked spectrum will be scaled to
			unity.
	
	Returns
	 -- stacking_results (list)
			List of objects resulting from the stacking procedure.
	"""
	
	headers = {}
	frames = {}
	scaled_spectra = []
	scaled_errors = []
	scales = {}

	for f in files:
		with fits.open(f) as file_hdu_list:
			headers[f[:-5]] = file_hdu_list[0].header
			frames[f[:-5]] = file_hdu_list[0].data
			frame = file_hdu_list[0].data
	
		wavelength_axis = frame[0]
		spectrum = frame[1]
		error = frame[2]
		
		scaling_wavelength = float(scaling_wavelength)
		
		if scaling_wavelength == 0.0:
			scaling_wavelength += (wavelength_axis[0]+wavelength_axis[-1])*0.5
		
		scaling_index = np.argmin(np.abs(wavelength_axis-scaling_wavelength))
		
		scaled_data = scale_spectrum(frame, scaling_index)
		scaled_spectra.append(scaled_data[0])
		scaled_errors.append(scaled_data[1])
		scales[f[:-5]] = scaled_data[2]
	
	scaled_spectra = np.array(scaled_spectra).T
	scaled_errors = np.array(scaled_errors).T
	
	stacked_data = stacking(
		scaled_spectra,
		scaled_errors,
		wavelength_axis,
		len(files)
	)
	
	stacking_results = [wavelength_axis, stacked_data, scales, headers, frames]

	return stacking_results


def scale_spectrum(spectrum_data, scaling_index):
	"""
	This function scales a spectrum to unity at the wavelength
	element indicated by scaling_index.
	
	Args:
	 -- spectrum_data (numpy.array)
			Input spectrum
	 -- scaling index (int)
			Element of the spectrum where input spectrum wil be scaled
			to unity.
	
	Returns:
	 -- scaled_data (list)
			List containing the rescaled spectrum, the associated
			rescaled errors, and the scaling factor used during scaling.
	"""
	scaling_region = spectrum_data[1][scaling_index-20:scaling_index+21]
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
	scaled_spectrum = spectrum_data[1]/scale
	scaled_error = spectrum_data[2]/scale
	
	scaled_data = [scaled_spectrum, scaled_error, scale]

	return scaled_data


def stacking(spectra, errors, wavelength_axis, number_of_spectra):
	"""
	When supplied with spectra, their errors, a common wavelength axis,
	and the number of spectra to be stacked, this function median
	resamples the spectra within their uncertainties, median stacks
	them, and estimates the standard error of the resulting data points
	of the stacked spectrum.
	
	Args:
	 -- spectra (numpy.array)
			Array of spectra to be stacked.
	 -- errors (numpy.arrray)
			Array of errors for the spectra to be stacked.
	 -- wavelength_axis (numpy.array)
			Common wavelength axis for the spectra to be stacked.
	 -- number_of_spectra (int)
			The integer number of spectra being stacked together.
	
	Returns:
	 -- stacked_frame (numpy.array)
			Array containing the wavelengths axis of the stacked
			spectrum, the stacked spectral data, and the estimated
			errors of the stacked spectrum.
	"""
	
	stacked_spectrum = []
	stacked_errors = []
	for i in range(len(wavelength_axis)):
		uncertain_samples = []
		for j in range(number_of_spectra):
			uncertain_samples.append(
				np.random.normal(
					loc=spectra[i][j],
					scale=errors[i][j],
					size=100
				)
			)
			
		uncertain_samples = np.array(uncertain_samples)

		stacked_spectrum.append(np.median(uncertain_samples))
		stacked_errors.append(np.std(uncertain_samples)/((len(files)-1)**0.5))
	
	stacked_spectrum = np.array(stacked_spectrum)
	stacked_errors = np.array(stacked_errors)
	stacked_frame = np.array(
		[wavelength_axis, stacked_spectrum, stacked_errors]
	)
	
	return stacked_frame	


def sum_std(files, scaling_wavelength):
	"""
	This function is used to stack the spectrum in the summing mode of
	operation.
	
	Args:
	 -- files (list)
			List of filenames for files containing spectra to be
			stacked.
	 -- scaling_wavelength (float)
			Wavelength at which the stacked spectrum will be scaled to
			unity.
	
	Returns
	 -- stacking_results (list)
			List of objects resulting from the stacking procedure.
	"""
	
	headers = {}
	frames = {}
	scales = {}

	for f in files:
		with fits.open(f) as file_hdu_list:
			headers[f[:-5]] = file_hdu_list[0].header
			frames[f[:-5]] = file_hdu_list[0].data
			frame = file_hdu_list[0].data
		
		if f == files[0]:
			wavelength_axis = frame[0]
			stacked_frame = copy.deepcopy(frame)
			stacked_frame[2] *= stacked_frame[2]
		else:
			stacked_frame[1] += frame[1]
			stacked_frame[2] += (frame[2]*frame[2])
		
	stacked_frame[2] = stacked_frame[2]**0.5
	
	scaling_wavelength = float(scaling_wavelength)
	
	if scaling_wavelength == 0.0:
		scaling_wavelength += (wavelength_axis[0]+wavelength_axis[-1])*0.5
	
	scaling_index = np.argmin(np.abs(wavelength_axis-scaling_wavelength))
	
	scaled_data = scale_spectrum(stacked_frame, scaling_index)
	scaled_stack = [wavelength_axis, scaled_data[0], scaled_data[1]]
	scales["sum"] = scaled_data[2]
	
	stacking_results = [wavelength_axis, scaled_stack, scales, headers, frames]
	
	return stacking_results


###############################################################################
#### SCRIPT STARTS HERE  # # # # # # # # # # # # # # # # # # # # # # # # # #### 
###############################################################################	

parser = argparse.ArgumentParser(
	description="This script takes multiple 1D SpeX spectra and\
	combines them with a weighted median to produce a stacked 1D\
	spectrum with reduced noise."
)
parser.add_argument(
	"-p", "--plot", action="store_true", 
	help="[boolean] Plot the stacked spectrum on screen."
)
parser.add_argument(
	"-s", "--save", action="store_true",
	help="[boolean] Save the stacked spectrum to a new .fits file."
)
parser.add_argument(
	"-w", "--scaling_wavelength", default=0.0,
	help="[float] Wavelength at which all spectra are scaled to unity\
	before stacking. Must have the same units as the wavelength axis of\
	the spectrum. The median of the twenty data points both sides of\
	the chosen wavelength (40 points total) will be used as the scaling\
	factor. Defaults to '0.0' which will result in the spectra being\
	scaled to unity at their central wavelength."
)
parser.add_argument(
	"-u", "--summing", action="store_true",
	help="[boolean] Stack the spectra by summing them and then scaling.\
	This kind of stacking is appropriate to use when exposure times are\
	shorter than ~5 s."
)
args = parser.parse_args()

files = sorted(glob.glob("*.fits"))

if args.summing:
	stack_results = sum_std(files, args.scaling_wavelength)
else:
	stack_results = bootstrap_median(files, args.scaling_wavelength)

wavelength_axis = stack_results[0]
stacked_data = stack_results[1]
scales = stack_results[2]
headers = stack_results[3]
frames = stack_results[4]

if args.plot:
	fig = plt.figure(figsize=(10,6))
	ax = plt.subplot()
	ax.errorbar(
		wavelength_axis,
		stacked_data[1],
		yerr=stacked_data[2],
		marker='.',
		linestyle='',
		color="k",
		label="Stack"
	)
	ax.grid(linestyle="dashed")
	ax.set_xlabel("Wavelength, " + headers[files[0][:-5]]["XUNITS"])
	ax.set_ylabel("Scaled DN/s")
	ax.legend()
	plt.show()

if args.save:
	save_dict = {"stacked" : stacked_data}
	
	for i in headers:
		save_dict[i] = frames[i]
		
	instrument = headers[files[0][:-5]]["INSTRUME"]
	
	stack_header_dict = {
		"SpeX Spectrograph" : spexsh.stack_header_spex
	}
		
	new_hdu = fits.PrimaryHDU()
	new_header = new_hdu.header
	new_header["DATE"] = (
		datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%S"),
		"UT File creation date."
	)
	new_header["FITSDOI"] = (
		"10.1051/0004-6361:20010923", "FITS format definition paper DOI"
	)
	new_header["ORIGIN"] = ("stax.py v" + __version__, "Script that created this file")
	new_header["STAXDOI"] = ("10.5281/zenodo.12786056", "Script repository DOI")
	if args.summing:
		new_header["STAXMODE"] = ("Sum", "stax.py stacking method")
		new_header["STXSCALE"] = (scales["sum"], "Spectrum stacking scaling factor")
	else:
		new_header["STAXMODE"] = ("Bootstrap Median", "stax.py stacking method")	

	stack_hdu = stack_header_dict[instrument](
		new_hdu, headers, files, args.scaling_wavelength
	)
	
	stack_hdu.data = save_dict["stacked"]
	hdu_list = [stack_hdu]
	hdu_list_keys = list(save_dict.keys())
	for key in hdu_list_keys[1:]:
		if len(scales) > 1:
			headers[key]["SCALFACT"] = (
				scales[key], "Spectrum stacking scaling factor"
			)
		headers[key]["EXTNAME"] = key
		hdu_list.append(fits.ImageHDU(save_dict[key], header=headers[key]))
	
	hdu_list = fits.HDUList(hdu_list)
	hdu_list.writeto(new_header["OBJECT"].replace(" ","_")+"_STACK.fits")
	hdu_list.close()
