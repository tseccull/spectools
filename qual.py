#!/usr/bin/env python3

"""
	qual.py

	Copyright (C) 2025-04-09 Tom Seccull
	
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

	Last updated - 2025-04-09

	Description --------------------------------------------------------	
	This script allows the user to manually set the quality flags for
	data points in a spectrum.
"""

__version__ = "1.0.0"
__author__ = "Tom Seccull"

import argparse
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np

###############################################################################
#### SCRIPT STARTS HERE  # # # # # # # # # # # # # # # # # # # # # # # # # ####
###############################################################################

parser = argparse.ArgumentParser(
	description="This script allows the user to manually set the\
	quality flags for data points in a spectrum."
)
parser.add_argument(
	"data_file", type=str,
	help="Input spectrum FITS file name."
)
parser.add_argument(
	"-g", "--good_points", type=str,
	help="Wavelength range(s) where quality flags should be set to 0\
	(GOOD). Start and end of each wavelength range is separated by\
	'-' and ranges are separated by commas. e.g. 500-560,570,600."
)
parser.add_argument(
	"-b", "--bad_points", type=str,
	help="Wavelength range(s) where quality flags should be set to 1\
	(BAD). Start and end of each wavelength range is separated by '-'\
	and ranges are separated by commas. e.g. 500-560,570,600."
)
parser.add_argument(
	"-p", "--plot", action="store_true",
	help="Plot the new spectrum on screen and display points with\
	updated qual flags. If a gradient has also been measured, the\
	linear fit to the chosen spectral regions will also be presented."
)
parser.add_argument(
	"-s", "--save", action="store_true",
	help="Update the original spectrum file with new qual flags."
)
args = parser.parse_args()

with fits.open(args.data_file) as hdu_list:
	data_frame = hdu_list[0].data

wave_axis = data_frame[0]
spectrum = data_frame[1]
uncertainty = data_frame[2]
qual = data_frame[5]

if args.plot:
	waves = [wave_axis, wave_axis[qual==0]]
	specs = [spectrum, spectrum[qual==0]]
	uncs = [uncertainty, uncertainty[qual==0]]
	colors = ["k", "k"]
	alphas = [0.3, 1.]

if args.good_points:
	good_points = args.good_points.split(",")
	good_points = [[float(x.split("-")[0]),float(x.split("-")[1])] for x in good_points]
	good_inds = [
		np.where(np.logical_and(wave_axis>x[0],wave_axis<x[1])) for x in good_points
	]
	qual[good_inds] = 0.0

if args.bad_points:
	bad_points = args.bad_points.split(",")
	bad_points = [[float(x.split("-")[0]),float(x.split("-")[1])] for x in bad_points]
	bad_inds = [
		np.where(np.logical_and(wave_axis>x[0],wave_axis<x[1])) for x in bad_points
	]
	qual[bad_inds] = 1.0

if args.plot:
	if args.good_points:
		waves.append(np.squeeze(wave_axis[good_inds]))
		specs.append(np.squeeze(spectrum[good_inds]))
		uncs.append(np.squeeze(uncertainty[good_inds]))
		colors.append("#1f77b4")
		alphas.append(1.)
	if args.bad_points:
		waves.append(np.squeeze(wave_axis[bad_inds]))
		specs.append(np.squeeze(spectrum[bad_inds]))
		uncs.append(np.squeeze(uncertainty[bad_inds]))
		colors.append("#ff7f0e")
		alphas.append(1.)
	
	for i in range(len(waves)):
		plt.errorbar(
			waves[i], specs[i], uncs[i], color=colors[i], linestyle="", marker=".", alpha=alphas[i]
		)
	plt.show()

if args.save:
	with fits.open(args.data_file, mode="update") as hdu_list:
		hdu_list[0].data[5] = qual
		hdu_list.flush()
