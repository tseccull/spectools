#!/usr/bin/env python3

"""
	divide.py

	Copyright (C) 2024-05-08 Tom Seccull
	
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

	Last updated - 2025-03-19

	Description --------------------------------------------------------
	This script divides one 1D spectrum by another. divide.py expects
	both spectra to have a common wavelength axis, be the product
	of extraction by MOTES, and stacking by stack.py. Typical use of
	this script is to calibrate a spectrum of a minor planet with that
	of a solar twin or solar analog to derive the minor planet's
	reflectance spectrum.
"""

__author__ = "Tom Seccull"
__version__ = "1.0.5"

import argparse
import astropy.io.fits as fits
import datetime
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.widgets import Slider


def division(wave_axis, shifted_wave_axis, spectra, uncertainties):
	"""
	This function takes two spectra that are offset by a uniform
	wavelength shift, interpolates the shifted second spectrum to the
	wavelength axis of the first and divides the first spectrum by the 
	second spectrum. Uncertainties in the spectral data are also 
	propagated.
	
	Args:
	 -- wave_axis (numpy.array)
			Unshifted wavelength axis of the first spectrum.
	 -- shifted_wave_axis (numpy.array)
			Shifted wavelength axis of the second spectrum.
	 -- spectra (numpy.array)
			2D array containing the spectroscopic data of both the first
			and second spectra.
	 -- uncertainties (numpy.array)
			2D array containing the uncertainties of the spectroscopic
			data of the first and second spectra.
			
	Returns:
	 -- output_array (numpy.array)
			2D array containing the wavelength axis of the ratio
			spectrum, the ratio spectrum, the uncertainties of the ratio
			spectrum, and a row containing quality flags for each point
			in the ratio spectrum where 1=GOOD and 0=BAD.
	"""
	
	shifted_spectra = [spectra[0]]
	shifted_uncertainties = [uncertainties[0]]

	shifted_spectra.append(
		np.interp(wave_axis, shifted_wave_axis, spectra[1], left=0., right=0.)
	)

	shifted_uncertainties.append(
		(
			np.interp(
				wave_axis,
				shifted_wave_axis,
				spectra[1] + uncertainties[1],
				left=0., 
				right=0.
			)
			- shifted_spectra[1]
		)
	)

	mult_spec = shifted_spectra[0] * shifted_spectra[1]
	qual_1d = [0 if x!=0. else 1 for x in mult_spec]

	if offset > 0:
		zero_edges = [
			x+1 for x in range(len(qual_1d)) if qual_1d[x:x+2] == [1,0]
		]
		zero_edges = np.array(zero_edges)
		qual_1d = np.array(qual_1d)	
		qual_1d[zero_edges] = 1
	elif offset < 0:
		zero_edges = [
			x for x in range(len(qual_1d)) if qual_1d[x:x+2] == [0,1]
		]
		zero_edges = np.array(zero_edges)
		qual_1d = np.array(qual_1d)	
		qual_1d[zero_edges] = 1
	else:
		qual_1d = np.array(qual_1d)

	shifted_spectra[0][qual_1d==1] = 0
	shifted_spectra[1][qual_1d==1] = 0
	shifted_uncertainties[0][qual_1d==1] = 0
	shifted_uncertainties[1][qual_1d==1] = 0

	shifted_spectra[0][shifted_spectra[1]==0] = 0
	shifted_spectra[1][shifted_spectra[0]==0] = 1
	
	ratio_spectrum = shifted_spectra[0]/shifted_spectra[1]
	
	shifted_uncertainties[0][shifted_spectra[0]==0] = 0
	shifted_spectra[0][shifted_uncertainties[0]==0] = 1
	
	shifted_uncertainties[1][shifted_spectra[1]==0] = 0
	shifted_spectra[1][shifted_uncertainties[1]==0] = 1
	
	shifted_variance_0 = (shifted_uncertainties[0]/shifted_spectra[0])**2
	shifted_variance_1 = (shifted_uncertainties[1]/shifted_spectra[1])**2
	
	ratio_uncertainties = (
		np.abs(ratio_spectrum) 
	  * ((shifted_variance_0 + shifted_variance_1)**0.5)
	)
	
	output_array = np.array(
		[wavelength_axis, ratio_spectrum, ratio_uncertainties, qual_1d]
	)
	
	return output_array


def update_wave(val):
	"""
	This function controls updates to the wavelength shift plot based on
	the input from the plot's slider.
	
	Args:
	 -- val (numpy.float)
			The value of the wavelength shift slider.
	
	Returns: None
	"""
	
	global offset
	offset = wave_slider.val
	two_spec.set_data(wavelength_axis+offset, optimal_spectra[1])
	fig.canvas.draw()


###############################################################################
#### SCRIPT STARTS HERE  # # # # # # # # # # # # # # # # # # # # # # # # # ####
###############################################################################

parser = argparse.ArgumentParser(
	description="This script divides one 1D spectrum by another.\
	divide.py expects both spectra to have a common wavelength axis, be\
	the product of extraction by MOTES, and stacking by stack.py."
)
parser.add_argument("spec_file_one", type=str,
	help="[str] File name of the first spectrum, which takes the 'numerator'\
	role in the division." 
)
parser.add_argument("spec_file_two", type=str,
	help="[str] File name of the second spectrum, which takes the\
	'denominator' role in the division." 
)
parser.add_argument("-p", "--plot", action="store_true",
	help="[boolean] Plot the ratio spectrum on screen."
)
parser.add_argument("-s", "--save", action="store_true",
	help="[boolean] Save the ratio spectrum to a new .fits file."
)
args = parser.parse_args()

with fits.open(args.spec_file_one) as one:
	one_headers = [x.header for x in one]
	one_frames = [x.data for x in one]

with fits.open(args.spec_file_two) as two:
	two_headers = [x.header for x in two]
	two_frames = [x.data for x in two]

primary_head_one = one_headers[0]
primary_head_two = two_headers[0]

if one_headers[0]["NAXIS1"] < two_headers[0]["NAXIS1"]:
	diff = two_headers[0]["NAXIS1"] - one_headers[0]["NAXIS1"]
	two_frames = [x[:,:-diff] for x in two_frames]
elif one_headers[0]["NAXIS1"] > two_headers[0]["NAXIS1"]:
	diff = one_headers[0]["NAXIS1"] - two_headers[0]["NAXIS1"]
	one_frames = [x[:,:-diff] for x in one_frames]

wavelength_axis = one_frames[0][0]

optimal_spectra = np.array([one_frames[0][1], two_frames[0][1]])
optimal_uncertainties = np.array([one_frames[0][2], two_frames[0][2]])
aperture_spectra = np.array([one_frames[0][3], two_frames[0][3]])
aperture_uncertainties = np.array([one_frames[0][4], two_frames[0][4]])

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.19)

one_spec = ax.errorbar(
	wavelength_axis,
	optimal_spectra[0],
	yerr = optimal_uncertainties[0],
	marker='.',
	label=primary_head_one["OBJECT"]
)

two_spec, = ax.plot(
	wavelength_axis,
	optimal_spectra[1],
	marker='.',
	zorder=2,
	label=primary_head_two["OBJECT"]
)

ax.set_xlabel(r"$\lambda$, " + primary_head_one["WAVU"])
ax.set_ylabel("Relative Counts")
ax.legend()

ax_wave_slider = plt.axes([0.25, 0.05, 0.5, 0.03])
wave_slider = Slider(
	ax_wave_slider,
	r"$\lambda$ Offset, " + primary_head_two["WAVU"],
	-5,
	5,
	valinit=0,
	valfmt="%.2f",
	valstep=.005
)
wave_slider.on_changed(update_wave)
plt.show()

try:
	offset = round(offset,1)
except NameError:
	offset = 0.0

shifted_wavelength_axis = wavelength_axis + offset

optimal_ratio_frame = division(
	wavelength_axis,
	shifted_wavelength_axis,
	optimal_spectra,
	optimal_uncertainties
)

aperture_ratio_frame = division(
	wavelength_axis,
	shifted_wavelength_axis,
	aperture_spectra,
	aperture_uncertainties
)

aperture_ratio_frame[1][optimal_ratio_frame[3]==1] = 0
aperture_ratio_frame[2][optimal_ratio_frame[3]==1] = 0
aperture_ratio_frame[3][optimal_ratio_frame[3]==1] = 0

if args.plot:
	plt.errorbar(
		optimal_ratio_frame[0],
		optimal_ratio_frame[1],
		yerr=optimal_ratio_frame[2],
		marker='.',
		color='k',
		label=(
			primary_head_one["OBJECT"] 
		  + "/" 
		  + primary_head_two["OBJECT"] 
		  + " - Optimal"
		), 
		zorder=1
	)
	plt.errorbar(
		aperture_ratio_frame[0],
		aperture_ratio_frame[1],
		yerr=aperture_ratio_frame[2],
		marker='.',
		color='r',
		label=(
			primary_head_one["OBJECT"] 
		  + "/" 
		  + primary_head_two["OBJECT"] 
		  + " - Aperture"
		),
		zorder=0
	)
	
	ylim_bottom, ylim_top = plt.ylim()
	qual_shade_bottom = np.ones(len(aperture_ratio_frame[1])) * ylim_bottom
	qual_shade_top = np.ones(len(aperture_ratio_frame[1])) * ylim_top
	plt.fill_between(
		optimal_ratio_frame[0],
		qual_shade_bottom,
		qual_shade_top,
		where=optimal_ratio_frame[3]==1,
		facecolor="red",
		alpha=0.5,
		label="QUAL=1"
	)
	
	plt.ylim(ylim_bottom, ylim_top)
	plt.xlabel(r"$\lambda$, " + primary_head_one["WAVU"])
	plt.ylabel("Relative Reflectance")
	plt.grid(linestyle="dashed", alpha=0.5)
	plt.legend()
	plt.show()

if args.save:
	
	combi_ratio_frame = np.array(
		[
			optimal_ratio_frame[0],
			optimal_ratio_frame[1],
			optimal_ratio_frame[2],
			aperture_ratio_frame[1],
			aperture_ratio_frame[2],
			optimal_ratio_frame[3]
		]
	)
	
	ratio_hdu = fits.PrimaryHDU(combi_ratio_frame)
	ratio_header = ratio_hdu.header
	ratio_header["DATE"] = (datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), "UT file creation date.")
	ratio_header["FITSDOI"] = ("10.1051/0004-6361:20010923", "FITS format definition paper DOI")
	ratio_header["ORIGIN"] = ("divide.py v" + __version__, "Script that created this file.")
	ratio_header["DIVDOI"] = ("10.5281/zenodo.12786056", "Script repository DOI")
	ratio_header["INPUT1"] = (args.spec_file_one, "First input spectrum file")
	ratio_header["INPUT2"] = (args.spec_file_two, "Second input spectrum file")
	ratio_header["OBJECT1"] = (primary_head_one["OBJECT"], "Name of object in first spectrum")
	ratio_header["OBJECT2"] = (primary_head_two["OBJECT"], "Name of object in second spectrum")
	ratio_header["EXTENS1"] = (args.spec_file_one[:-5], "Extension of first spectrum data frame")
	ratio_header["EXTENS2"] = (args.spec_file_two[:-5], "Extension of second spectrum data frame")
	ratio_header["WAVEOFFS"] = (offset, "Wavelength offset applied to spectrum two")
	ratio_header["WAVOFFSU"] = (primary_head_two["WAVU"], "Unit of wavelength offset")
	ratio_header["WAVU"] = (primary_head_one["WAVU"], "Wavelength axis units")
	ratio_header["HDUROW0"] = "Wavelength axis, " + primary_head_one["WAVU"]
	ratio_header["HDUROW1"] = "Optimal ratio spectrum"
	ratio_header["HDUROW2"] = "Optimal ratio spectrum uncertainty"
	ratio_header["HDUROW3"] = "Aperture ratio spectrum"
	ratio_header["HDUROW4"] = "Aperture ratio spectrum uncertainty"
	ratio_header["HDUROW5"] = "Quality flags: 0=GOOD, 1=BAD"
	ratio_header["EXTNAME"] = (
		primary_head_one["OBJECT"].replace(" ","_")
	  + "/" 
	  + primary_head_two["OBJECT"].replace(" ","_")
	)
	
	listed_hdus = [ratio_hdu]
	
	for i in range(len(one_headers)):
		if "STACK" in one_headers[i]["EXTNAME"]:
			one_headers[i]["EXTNAME"] = (
				one_headers[i]["OBJECT"].replace(" ","_") + "_STACK"
			)
		listed_hdus.append(fits.ImageHDU(one_frames[i], header=one_headers[i]))
	for i in range(len(two_headers)):
		if "STACK" in two_headers[i]["EXTNAME"]:
			two_headers[i]["EXTNAME"] = (
				two_headers[i]["OBJECT"].replace(" ","_") + "_STACK"
			)
		listed_hdus.append(fits.ImageHDU(two_frames[i], header=two_headers[i]))
	
	new_file_name = (
		primary_head_one["OBJECT"].replace(" ","_") 
	  + "%" 
	  + primary_head_two["OBJECT"].replace(" ","_") 
	  + ".fits"
	)
	
	hdu_list = fits.HDUList(listed_hdus)
	hdu_list.writeto(new_file_name)
	hdu_list.close()
