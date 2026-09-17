#!/usr/bin/env python3

"""
iooio.py

	Copyright (C) 2024-05-06 Tom Seccull
	
	This script is part scrap.py, in the spectools repo hosted at 
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

	Last updated - 2026-09-17

	Description --------------------------------------------------------
	This module contains the preparation and save functions called by
	scrap.py for IO:O images.
"""

__author__ = "Tom Seccull"

import astropy.io.fits as fits
import astroscrappy as asc
import copy
import datetime
import numpy as np


def prep_ioo(in_file, fine_structure_mode):
	"""
	Combines all necessary data for detect_cosmics() into a dictionary
	for an image observed with IO:O.
	
	Args:
	 -- in_file (.fits HDU list)
			The object produced by using fits.open on the fits file
			currently being processed. It contains all the dataframes
			and headers for the current spectrum.
	 -- fine_structure_mode (str)
			A string keyword to tell prep_gmos how detect_cosmics() will
		    generate the fine structure image.
		
	Returns:
	 -- detect_cosmics_input (dict)
			A dictionary of dataframes and parameters collected from the
			input file or calculated from it. Items in this dictionary
			are all inputs for detect_cosmics().
	 -- detect_cosmics_input["in_data_frame"] (numpy.ndarray)
			The 2D science dataframe. The value of each pixel is in 
		    Analog-to-Digital Units (ADU)
	 -- detect_cosmics_input["in_quality_frame"] (numpy.ndarray)
			2D frame flagging bad pixels in the science and variance
			frames
	 -- detect_cosmics_input["in_background_frame"] (numpy.ndarray)
			estimate of the 2D background in the science frame. Units
			are ADUs
	 -- detect_cosmics_input["in_variance_frame"] (numpy.ndarray)
			2D frame containing the variance of each pixel in the
			science frame; units are ADU^2
	 -- detect_cosmics_input["detector_gain"] (float)
			average CCD detector gain for this data in e-/ADU
	 -- detect_cosmics_input["read_noise"] (float)
			average detector readout noise e- rms
	 -- detect_cosmics_input["psf_model"] (str)
			notes the Point Spread Function model adopted by
			detect_cosmics() when building the fine structure image.
			"gaussy" is used here because there is no option for a
			directional Moffat profile. 
	 -- detect_cosmics_input["fwhm"] (float)
			Full Width at Half Maximum measured for the median spatial
			profile of the spectrum. Units are pixels.
	 -- detect_cosmics_input["psf_size"] (int)
			size of the PSF model in pixels that will be convolved with
			the data by detect_cosmics() to produce the fine structure
			model. This value must be odd.
	-- detect_cosmics_input["saturation_level"] (int)
			number of counts at which the detector saturates.
	"""
	
	in_header = in_file["PRIMARY"].header
	in_data = in_file["PRIMARY"].data
	
	# I don't know why NAXIS2 and NAXIS1 are this way round.
	frame_shape = (in_header["NAXIS2"], in_header["NAXIS1"])
	
	background_frame = np.tile(in_header["BACKGRD"], frame_shape)
	
	count_read_noise = (in_header["READNOIS"] / in_header["GAIN"])
	count_dark_signal = 0.002 * in_header["EXPTIME"]
	variance_frame = np.abs(in_data) + (count_read_noise**2) + count_dark_signal
	
	# Create the output dictionary and fill it with relevant dataframes
	# and values from the input file.
	detect_cosmics_input = {
		      "in_data_frame" : in_file["PRIMARY"].data,
		   "in_quality_frame" : np.zeros(frame_shape),
		"in_background_frame" : background_frame,
		  "in_variance_frame" : variance_frame,
		      "detector_gain" : in_header["GAIN"],
		         "read_noise" : in_header["READNOIS"],
		          "psf_model" : "gauss",
		   "saturation_level" : 65536*in_header["GAIN"]	
	}

	# If the fine stucture image is to be generated with a convolution
	# of a model PSF, estimate the PSF's FWHM and size in pixels.
	# It is not always trivial to estimate FWHM from IO:O data, so
	# assume 1 arcsec seeing
	if fine_structure_mode == "convolve":
		fwhm = 1/in_header["CCDSCALE"]
		detect_cosmics_input["fwhm"] = fwhm
		psf_scale = np.ceil(2.8*fwhm)
		if psf_scale % 2 == 0:
			psf_scale += 1
		
		detect_cosmics_input["psf_size"] = psf_scale
		
	# If a median filter is being used to generate the fine structure
	# image, fwhm and psf_size aren't needed. In this case we set
	# their values to the defaults for detect_cosmics() with the
	# knowledge that they won't be used.
	else:
		detect_cosmics_input["fwhm"] = 2.5
		detect_cosmics_input["psf_size"] = 7
	
	return detect_cosmics_input


def save_ioo(
	file_name,
	in_file,
	primary_header,
	cosmic_ray_mask,
	clean_science_frame,
	detect_cosmics_parameters,
	command_line_arguments,
	scrap_version
):
	'''
	Constructs and saves a new .fits file combining the original input
	dataframes and headers with the cleaned science data and updated
	quality mask.
	
	Args:
	 -- file_name (str)
			The name of the input file.
	 -- in_file (.fits HDU list)
			The object produced by using fits.open on the fits file
			currently being processed. It contains all the dataframes
			and headers for the current spectrum.
	 -- primary_header (.fits header)
			The header of the primary header data unit in in_file.
	 -- cosmic_ray_mask (numpy.ndarray)
			A 2D array flagging the location of cosmic ray detections.
	 -- clean_science_frame (numpy.ndarray)
			The 2D science data array after cosmic rays have been
			cleaned.
	 -- detect_cosmics_parameters (dict)
			A dictionary of data and parameters fed to Astroscrappy 
			detect_cosmics().
	 -- command_line_arguments (class)
			The scrap.py command line argument namespace.
	 -- scrap_version (str)
	        String noting the current version of scrap.py.
	Returns:
	 -- None
	'''
	
	# Update Primary header of the fits file.
	primary_header["CRSCRIPT"] = (
		"scrap.py v" + scrap_version, "Cosmic ray masking/cleaning script")
	primary_header["SCRAPDOI"] = ("10.5281/zenodo.12786056", "Script repository DOI")
	primary_header["CRMETHOD"] = (
		"Astroscrappy v" + asc.__version__,
		"Cosmic ray masking/cleaning method"
	)
	primary_header["ASCDOI"] = (
		"10.5281/zenodo.1482019", "Astroscrappy Zenodo DOI"
	)
	primary_header["VDOKDOI"] = (
		"10.1086/323894", "van Dokkum 2001 PASP paper DOI"
	)
	primary_header["CRDATE"] = (
		datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
		"UT timestamp for Astroscrappy"
	)
	primary_header["CRMETHOD"] = (
		"Astroscrappy v" + asc.__version__,
		"Cosmic ray masking/cleaning method"
	)
	primary_header["CRDATE"]   = (
		datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
		"UT timestamp for Astroscrappy"
	)
	primary_header["CRBKGD"]   = (
		"BACKGRD", 
		"Background estimation method for Astroscrappy"
	)
	primary_header["CRSIGCLP"] = (
		command_line_arguments.sigma_clip,
		"Astroscrappy sigclip value"
	)
	primary_header["CRSIGFRC"] = (
		command_line_arguments.sigma_frac,
		"Astroscrappy sigfrac value"
	)
	primary_header["CROBJLIM"] = (
		command_line_arguments.obj_limit,
		"Astroscrappy objlim value"
	)
	primary_header["CRDETSAT"] = (
		detect_cosmics_parameters["saturation_level"],
		"Astroscrappy satlevel value (e-)"
	)
	primary_header["CRNITER"]  = (
		command_line_arguments.iteration_number,
		"Astroscrappy niter value"
	)
	primary_header["CRSEPMED"] = (
		command_line_arguments.separable_median,
		"Astroscrappy sepmed value"
	)
	primary_header["CRDCTYPE"] = (
		command_line_arguments.data_clean_type,
		"Astroscrappy cleantype value"
	)
	primary_header["CRFSMODE"] = (
		command_line_arguments.fine_structure_mode,
		"Astroscrappy fsmode value"
	)
	
	quality_header = copy.deepcopy(primary_header)
	quality_header["EXTNAME"] = "COSMIC_MASK"
	
	# If fsmode is "median", then no psf parameters.
	if command_line_arguments.fine_structure_mode == "convolve":
		primary_header["CRPSFMOD"] = (
			detect_cosmics_parameters["psf_model"],
			"Astroscrappy psfmodel value"
		)
		primary_header["CRPSFWHM"] = (
			detect_cosmics_parameters["fwhm"],
			"Astroscrappy psffwhm value (pix)"
		)
		primary_header["CRPSFSIZ"] = (
			detect_cosmics_parameters["psf_size"],
			"Astroscrappy psfsize value (pix)"
		)
	
	# Construct the output .fits file.
	primary_hdu = fits.PrimaryHDU(clean_science_frame, header=primary_header)
	cr_hdu = fits.ImageHDU(cosmic_ray_mask, header=quality_header)
	original_science_hdu = fits.ImageHDU(in_file["PRIMARY"].data)
	original_science_hdu.header["EXTNAME"] = "OG_SCI"
	hdu_list=fits.HDUList(
		[
			primary_hdu,
			cr_hdu,
			original_science_hdu
		]
	)
	hdu_list.writeto("c" + file_name)
	hdu_list.close()
