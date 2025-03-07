#!/usr/bin/env python3

"""
gmosio.py - written by T. Seccull, 2024-05-06
	
	Called by: scrap.py
	Last updated: - 2025-03-07

	This module contains the preparation and save functions called by
	scrap.py for GMOS spectra.
"""

__author__ = "Tom Seccull"

import astropy.io.fits as fits
import astroscrappy as asc
import copy
import datetime
import numpy as np

from .moffat import moffat_least_squares


def prep_gmos(in_file, fine_structure_mode):
	"""
	Combines all necessary data for detect_cosmics() into a dictionary
	for a spectrum observed with GMOS-N or GMOS-S.
	
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
	
	# Take a median of the 2D science frame along the spatial axis to
	# get a median sky spectrum, then tile that to be the same shape as
	# the original data frame. This will be used as the estimated
	# background.
	background_frame = np.tile(
		np.median(in_file["SCI"].data, axis=0), 
		(np.shape(in_file["SCI"])[0], 1)
	)
	
	# DQ=64 marks unilluminated regions of the frame, but does not
	# indicate bad data.
	in_file["DQ"].data[in_file["DQ"].data==64] = 0
	
	# Create the output dictionary and fill it with relevant dataframes
	# and values from the input file.
	detect_cosmics_input = {
		      "in_data_frame" : in_file["SCI"].data,
		   "in_quality_frame" : in_file["DQ"].data,
		"in_background_frame" : background_frame,
		  "in_variance_frame" : in_file["VAR"].data,
		      "detector_gain" : in_file["SCI"].header["GAIN"],
		         "read_noise" : in_file["SCI"].header["RDNOISE"],
		          "psf_model" : "gaussy"	
	}
	
	# The detectors of GMOS-N and GMOS-S are slightly different, become
	# non-linear at different points, and therefore have different full
	# well depths.
	if in_file[0].header["INSTRUME"] == "GMOS-N":
		detect_cosmics_input["saturation_level"] = 106822
	else:
		detect_cosmics_input["saturation_level"] = 117963

	# If the fine stucture image is to be generated with a convolution
	# of a model PSF, estimate the PSF's FWHM and size in pixels.
	if fine_structure_mode == "convolve":
		# All this is to get an initial estimate of the IQ. Tables below
		# are based on the condition constraints used by Gemini. See web
		# page below.
		# https://www.gemini.edu/observing/telescopes-and-sites/sites#ImageQuality
		iq_dict = {
			"20-percentile": 0,
			"70-percentile": 1,
			"85-percentile": 2,
			"100-percentile": 3,
			"Any": 3,
			"UNKNOWN": 3,
		}
	
		wav_tab = np.array(
			[
				[0000.0, 4000.0, 0],
				[4000.0, 5500.0, 1],
				[5500.0, 7000.0, 2],
				[7000.0, 8500.0, 3],
				[8500.0, 9750.0, 4],
				[9750.0, 11000.0, 5],
			]
		)
	
		iq_tab = np.array(
			[
				[0.6, 0.90, 1.20, 2.00],
				[0.6, 0.85, 1.10, 1.90],
				[0.5, 0.75, 1.05, 1.80],
				[0.5, 0.75, 1.05, 1.70],
				[0.5, 0.70, 0.95, 1.70],
				[0.4, 0.70, 0.95, 1.65],
			]
		)
	
		iq = in_file[0].header["RAWIQ"]
	
		short_wav = in_file["SCI"].header["CRVAL1"]
	
		for i in wav_tab:
			if short_wav > i[0] and short_wav < i[1]:
				seeing = float(iq_tab[int(i[2])][int(iq_dict[iq])])
				break
	
		pixel_resolution = in_file["SCI"].header["PIXSCALE"]
		
		# Measure the FWHM of the spectrum's spatial profile. A Moffat
		# profile is fitted to the median spatial profile to get an
		# accurate measure of the Full Width at Half Maximum, even
		# though the psf model used by detect_cosmics() is a gaussian.
		
		# Calculate median spatial profile of the spectrum.
		median_profile = np.nanmedian(in_file["SCI"].data, axis=1)
		
		# Scipy least squares doesn't like really tiny numbers like
		# fluxes in erg/s/cm^2/Angstrom, so it's necessary to scale the
		# data to a size that least squares can handle. The shape of the
		# profile fitted to the scaled spatial profile is the same as
		# the unscaled, so FWHM is unaffected.
		absolute_median_counts = np.abs(np.nanmedian(median_profile))
		data_scale = 10 ** np.abs(np.floor(np.log10(absolute_median_counts)))
		
		# Fit the median spatial profile with a Moffat function.
		moffat_parameters = moffat_least_squares(
			range(np.shape(in_file["SCI"].data)[0]),
			median_profile * data_scale,
			seeing,
			pixel_resolution,
			50
		)
		
		# Get an improved estimate of the FWHM of the spectrum from the
		# best fit Moffat profile.
		sqrt_beta_factor = np.sqrt((2**(1/moffat_parameters[3])) - 1)
		fwhm = 2 * moffat_parameters[2] * sqrt_beta_factor
		
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


def save_gmos(
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
	primary_header.comments["GEM-TLM"] = "Time of last modification by GEMINI"
	primary_header.comments["DATE"] = "Date FITS file was generated by IRAF"
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
	
	# Update the SCI frame header and copy it to create the OG_SCI header. 
	science_header = in_file["SCI"].header
	clean_science_header = copy.deepcopy(science_header)
	science_header["EXTNAME"] = "OG_SCI"
	
	quality_header = in_file["DQ"].header
	quality_frame = in_file["DQ"].data + cosmic_ray_mask
	
	update_headers = [primary_header, clean_science_header, quality_header]
	
	for head in update_headers:
		# Update the header of the new SCI header to add details on the
		# cosmic ray detection, masking, and cleaning process.
		head["CRMETHOD"] = (
			"Astroscrappy v" + asc.__version__,
			"Cosmic ray masking/cleaning method"
		)
		head["CRDATE"]   = (
			datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
			"UT timestamp for Astroscrappy"
		)
		head["CRBKGD"]   = (
			"Median", 
			"Background estimation method for Astroscrappy"
		)
		head["CRSIGCLP"] = (
			command_line_arguments.sigma_clip,
			"Astroscrappy sigclip value"
		)
		head["CRSIGFRC"] = (
			command_line_arguments.sigma_frac,
			"Astroscrappy sigfrac value"
		)
		head["CROBJLIM"] = (
			command_line_arguments.obj_limit,
			"Astroscrappy objlim value"
		)
		head["CRDETSAT"] = (
			detect_cosmics_parameters["saturation_level"],
			"Astroscrappy satlevel value (e-)"
		)
		head["CRNITER"]  = (
			command_line_arguments.iteration_number,
			"Astroscrappy niter value"
		)
		head["CRSEPMED"] = (
			command_line_arguments.separable_median,
			"Astroscrappy sepmed value"
		)
		head["CRDCTYPE"] = (
			command_line_arguments.data_clean_type,
			"Astroscrappy cleantype value"
		)
		head["CRFSMODE"] = (
			command_line_arguments.fine_structure_mode,
			"Astroscrappy fsmode value"
		)
		
		# If fsmode is "median", then no psf parameters.
		if command_line_arguments.fine_structure_mode == "convolve":
			head["CRPSFMOD"] = (
				detect_cosmics_parameters["psf_model"],
				"Astroscrappy psfmodel value"
			)
			head["CRPSFWHM"] = (
				detect_cosmics_parameters["fwhm"],
				"Astroscrappy psffwhm value (pix)"
			)
			head["CRPSFSIZ"] = (
				detect_cosmics_parameters["psf_size"],
				"Astroscrappy psfsize value (pix)"
			)
	
	# Construct the output .fits file.
	primary_hdu = fits.PrimaryHDU(header=primary_header)
	mdf_hdu = in_file["MDF"]
	clean_science_hdu = fits.ImageHDU(
		clean_science_frame, 
		header=clean_science_header
	)
	variance_hdu = in_file["VAR"]
	quality_hdu = fits.ImageHDU(quality_frame, header=quality_header)
	original_science_hdu = in_file["OG_SCI"]
	hdu_list=fits.HDUList(
		[
			primary_hdu,
			mdf_hdu,
			clean_science_hdu,
			variance_hdu,
			quality_hdu,
			original_science_hdu
		]
	)
	hdu_list.writeto("c" + file_name)
	hdu_list.close()
