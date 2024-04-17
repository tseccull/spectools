#! /home/tom/anaconda3/envs/work/bin/python

"""
scrap.py - written by T. Seccull, 2023-10-24 - v1.0.4

	Last updated - 2024-04-16

	Locate, mask, and clean cosmic ray hits in 2D spectroscopic data
	with Astroscrappy/LA Cosmic. This script runs the detect_cosmics()
	function from Astroscrappy on supplied astronomical spectroscopic
	data. Astroscrappy is a Python implentation of Pieter van Dokkum's
	LA Cosmic. Cite both Astroscrappy and LA Cosmic if used. scrap.py
	will assume all files in the current directory are .fits formatted
	2D spectra that need their cosmic rays masked and will try to apply
	detect_cosmics() to each in turn. scrap.py will replace the primary
	input data frame with the cleaned data array in the output .fits
	file; a copy of the original input data will be stored in a new
	Header Data Unit (HDU) of the output .fits file. The boolean cosmic
	ray mask that scrap.py creates will be either added to an existing
	quality frame (e.g. flagging bad pixels), or if a quality frame
	doesn't exist in the original file the crmask will be added to the
	output file as the new quality frame in a new HDU. A gaussian
	psfmodel is assumed, so in detect_cosmics() psfk=None by default,
	and psfbeta is ignored. Other detect_cosmics() parameters are either
	taken directly from the fits data and headers, or they can be set 
	optionally when scrap.py is run.
    
    Original paper describing LA Cosmic 
		van Dokkum 2001, PASP, 113, 1420
		https://doi.org/10.1086/323894
    Astroscrappy Zenodo DOI
		McCully et al. 2018, Astropy/Astroscrappy: v1.0.5 Zenodo Release
		https://doi.org/10.5281/zenodo.1482019
    Astroscrappy docs - https://astroscrappy.readthedocs.io/en/latest/
"""


import argparse
import astropy.io.fits as fits
import astroscrappy as asc
import copy
import datetime
import glob
import numpy as np

from scipy.optimize import least_squares


###############################################################################
def moffat_least_squares(
	spatial_axis, column, seeing, pixel_resolution, end_clip
):
    """
    Takes a data column, spatial axis and seeing of the observation and
    fits a Moffat function to the column using a least squares method.
    Returns the best fit parameters of the Moffat function.

    Args:
     -- spatial_axis (numpy.ndarray)
			The spatial axis of the data being fit.
     -- column (numpy.ndarray)
			The data being fitted.
     -- seeing (float)
			The estimated FWHM of the spatial profile.
     -- pixel_resolution (float)
			The spatial resolution of each pixel in arcsec/pixel.
     -- end_clip (int)
			The number of pixels at each end of the spatial profile
			array to ignore when fitting the Moffat profile.

    Returns:
     -- parameter_list (list)
			The list of best fit output parameters returned by the least
			squares routine.
    """

	# Clip the median spatial profile to be fitted based on the value of
	# end_clip
    column[:end_clip] = np.median(column)
    column[-end_clip:] = np.median(column)

    # Set up initial conditions for the least squares fit.
    # x0 = [
    #	 amplitude,
    # 	 centre,
    #	 alpha,
    #	 beta,
    #	 background gradient,
    #	 background level
    # ]
    # Initial beta estimate comes from optimal value from atmospheric
    # turbulence theory as described in 
    # Trujillo, I. et al. (2001), MNRAS, 328, 977-985
    # See https://ui.adsabs.harvard.edu/abs/2001MNRAS.328..977T/abstract
    x0 = [
        np.nanmedian(np.sort(column)[-3:]),
        np.argmax(column),
        seeing / pixel_resolution,
        4.765,
        0.0,
        np.median(np.concatenate((column[:5], column[-5:]))),
    ]

    # Run the least squares fit.
    res_lsq = least_squares(
        moffat_resid,
        x0,
        bounds=(
            [
                0.0, 
                np.argmax(column) - 1.0,
                0.0, 
                0.0, 
                0.0, 
                -np.inf
            ],
            [
                np.inf,
                np.argmax(column) + 1,
                (5 * seeing / pixel_resolution),
                5.0,
                np.inf,
                np.inf,
            ],
        ),
        args=(spatial_axis, column),
        method="trf",
        ftol=1e-12,
    )
    
    parameter_list = [
        res_lsq.x[0],
        res_lsq.x[1],
        res_lsq.x[2],
        res_lsq.x[3],
        res_lsq.x[4],
        res_lsq.x[5]
    ]
    return parameter_list


###############################################################################
def moffat_resid(x, spatial_axis, data):
    """
    Calculates residuals of fitted moffat profile and the data for the
    least squares fitting.

    Description:
        A = x[0]
        c = x[1]
        alpha = x[2]
        beta = x[3]
        B = x[4]
        m = x[5]

    Args:
     -- x (numpy.ndarray)
			An array of parameters defining the shape of the model
			moffat profile.
     -- spatial_axis (numpy.ndarray)
			The spatial axis of the data.
     -- data (numpy.ndarray)
			The data.

    Returns:
     -- residual (numpy.ndarray)
			The residual array between the model moffat profile and the
			data.
    """
    
    moffat_profile = x[0]*(
		(1+((spatial_axis-x[1])*(spatial_axis-x[1]))/(x[2]*x[2]))**-x[3]
	)
    
    residual = moffat_profile + x[4] + (spatial_axis*x[5]) - data
    
    return residual

 
###############################################################################
def prep_gmos(in_file, primary_header, fine_structure_mode):
	"""
	Combines all necessary data for detect_cosmics() into a dictionary
	for a spectrum observed with GMOS-N or GMOS-S.
	
	Args:
	 -- in_file (.fits HDU list)
			The object produced by using fits.open on the fits file
			currently being processed. It contains all the dataframes
			and headers for the current spectrum.
	 -- primary_header (.fits header)
			The header of the primary header data unit in in_file
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
	
	# Create the output dictionary and fill it with relevant dataframes
	# and values from the input file.
	detect_cosmics_input = {
		      "in_data_frame" : in_file["SCI"].data,
		   "in_quality_frame" : in_file["DQ"].data,
		"in_background_frame" : background_frame,
		  "in_variance_frame" : in_file["VAR"].data,
		      "detector_gain" : primary_header["GAINMULT"],
		         "read_noise" : primary_header["RDNOISE"],
		          "psf_model" : "gaussy"	
	}
	
	# The detectors of GMOS-N and GMOS-S are slightly different, become
	# non-linear at different points, and therefore have different full
	# well depths.
	if primary_header["INSTRUME"] == "GMOS-N":
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
	
		iq = primary_header["RAWIQ"]
	
		short_wav = in_file["SCI"].header["CRVAL1"]
	
		for i in wav_tab:
			if short_wav > i[0] and short_wav < i[1]:
				seeing = float(iq_tab[int(i[2])][int(iq_dict[iq])])
				break
	
		pixel_resolution = primary_header["PIXSCALE"]
		
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


###############################################################################
def save_gmos(
	file_name,
	in_file,
	primary_header,
	cosmic_ray_mask,
	clean_science_frame,
	detect_cosmics_parameters,
	command_line_arguments
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
	Returns:
	 -- None
	'''
	
	# Update Primary header of the fits file.
	primary_header.comments["IRAF-TLM"] = "Time of last modification by IRAF"
	primary_header.comments["DATE"] = "Date FITS file was generated by IRAF"
	primary_header["CRMETHOD"] = (
		"Astroscrappy v" + asc.__version__,
		"Cosmic ray masking/cleaning method"
	)
	primary_header["CRDATE"] = (
		datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
		"UT timestamp for Astroscrappy"
	)
	
	# Update the SCI frame header and copy it to create the OG_SCI header. 
	science_header = in_file["SCI"].header
	science_header.comments["IRAF-TLM"] = "Time of last modification by IRAF"
	science_header.comments["DATE"] = "Date FITS file was generated by IRAF"
	clean_science_header = copy.deepcopy(science_header)
	science_header["EXTNAME"] = "OG_SCI"
	
	# Update the header of the new SCI header to add details on the
	# cosmic ray detection, masking, and cleaning process.
	clean_science_header["CRMETHOD"] = (
		"Astroscrappy v" + asc.__version__,
		"Cosmic ray masking/cleaning method"
	)
	clean_science_header["CRDATE"]   = (
		datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
		"UT timestamp for Astroscrappy"
	)
	clean_science_header["CRBKGD"]   = (
		"Median", 
		"Background estimation method for Astroscrappy"
	)
	clean_science_header["CRSIGCLP"] = (
		command_line_arguments.sigma_clip,
		"Astroscrappy sigclip value"
	)
	clean_science_header["CRSIGFRC"] = (
		command_line_arguments.sigma_frac,
		"Astroscrappy sigfrac value"
	)
	clean_science_header["CROBJLIM"] = (
		command_line_arguments.obj_limit,
		"Astroscrappy objlim value"
	)
	clean_science_header["CRDETSAT"] = (
		detect_cosmics_parameters["saturation_level"],
		"Astroscrappy satlevel value (e-)"
	)
	clean_science_header["CRNITER"]  = (
		command_line_arguments.iteration_number,
		"Astroscrappy niter value"
	)
	clean_science_header["CRSEPMED"] = (
		command_line_arguments.seperable_median,
		"Astroscrappy sepmed value"
	)
	clean_science_header["CRDCTYPE"] = (
		command_line_arguments.data_clean_type,
		"Astroscrappy cleantype value"
	)
	clean_science_header["CRFSMODE"] = (
		command_line_arguments.fine_structure_mode,
		"Astroscrappy fsmode value"
	)
	
	# If fsmode is "median", then no psf parameters.
	if command_line_arguments.fine_structure_mode == "convolve":
		clean_science_header["CRPSFMOD"] = (
			detect_cosmics_parameters["psf_model"],
			"Astroscrappy psfmodel value"
		)
		clean_science_header["CRPSFWHM"] = (
			detect_cosmics_parameters["fwhm"],
			"Astroscrappy psffwhm value (pix)"
		)
		clean_science_header["CRPSFSIZ"] = (
			detect_cosmics_parameters["psf_size"],
			"Astroscrappy psfsize value (pix)"
		)
	
	# Update the header of the new DQ frame to add details on the
	# cosmic ray detection, masking, and cleaning process. Add the
	# cosmic ray mask to the DQ frame. 
	quality_header = in_file["DQ"].header
	quality_frame = in_file["DQ"].data + cosmic_ray_mask
	
	quality_header.comments["IRAF-TLM"] = "Time of last modification by IRAF"
	quality_header.comments["DATE"] = "Date FITS file was generated by IRAF"
	quality_header["CRMETHOD"] = (
		"Astroscrappy v" + asc.__version__,
		"Cosmic ray masking/cleaning method"
	)
	quality_header["CRDATE"]   = (
		datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
		"UT timestamp for Astroscrappy"
	)
	quality_header["CRBKGD"]   = (
		"Median",
		"Background estimation method for Astroscrappy"
	)
	quality_header["CRSIGCLP"] = (
		command_line_arguments.sigma_clip,
		"Astroscrappy sigclip value"
	)
	quality_header["CRSIGFRC"] = (
		command_line_arguments.sigma_frac,
		"Astroscrappy sigfrac value"
	)
	quality_header["CROBJLIM"] = (
		command_line_arguments.obj_limit,
		"Astroscrappy objlim value"
	)
	quality_header["CRDETSAT"] = (
		detect_cosmics_parameters["saturation_level"],
		"Astroscrappy satlevel value (e-)"
	)
	quality_header["CRNITER"]  = (
		command_line_arguments.iteration_number,
		"Astroscrappy niter value"
	)
	quality_header["CRSEPMED"] = (
		command_line_arguments.seperable_median,
		"Astroscrappy sepmed value"
	)
	quality_header["CRDCTYPE"] = (
		command_line_arguments.data_clean_type,
		"Astroscrappy cleantype value"
	)
	quality_header["CRFSMODE"] = (
		command_line_arguments.fine_structure_mode,
		"Astroscrappy fsmode value"
	)
	
	# If fsmode is "median", then no psf parameters.
	if command_line_arguments.fine_structure_mode == "convolve":
		quality_header["CRPSFMOD"] = (
			detect_cosmics_parameters["psf_model"],
			"Astroscrappy psfmodel value"
		)
		quality_header["CRPSFWHM"] = (
			detect_cosmics_parameters["fwhm"],
			"Astroscrappy psffwhm value (pix)"
		)
		quality_header["CRPSFSIZ"] = (
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


###############################################################################
#### SCRIPT STARTS HERE  # # # # # # # # # # # # # # # # # # # # # # # # # #### 
###############################################################################

# Parse scrap.py arguments
parser = argparse.ArgumentParser(
	description="Locate, mask, and clean cosmic ray hits in 2D\
	spectroscopic data with Astroscrappy/LACosmic. This script runs the\
	detect_cosmics() function from Astroscrappy on supplied\
	astronomical spectroscopic data. Astroscrappy is a Python\
	implentation of Pieter van Dokkum's LACosmic. Cite both\
	Astroscrappy and LACosmic if used. scrap.py will assume all files\
	in the current directory are .fits formatted 2D spectra that need\
	their cosmic rays masked and will try to apply detect_cosmics() to\
	each in turn. scrap.py will replace the primary input data frame\
	with the cleaned data array in the output .fits file; a copy of the\
	original input data will be stored in a new Header Data Unit (HDU)\
	of the output .fits file. The boolean cosmic ray mask that scrap.py\
	creates will be either added to an existing quality frame (e.g.\
	flagging bad pixels), or if a quality frame doesn't exist in the\
	original file the crmask will be added to the output file as the\
	new quality frame in a new HDU. A gaussian psfmodel is assumed, so\
	in detect_cosmics() psfk=None by default, and psfbeta is ignored.\
	Other detect_cosmics() parameters are either taken directly from\
	the fits data and headers, or they can be set optionally when\
	scrap.py is run. Astroscrappy docs and links to citables ->\
	https://astroscrappy.readthedocs.io/en/latest/"
)
parser.add_argument("-sc", "--sigma_clip", default=4.5, type=float,
	help="[float] Laplacian-to-noise limit for cosmic ray detection.\
	Lower values will flag more pixels as cosmic rays. Default: 4.5"
)
parser.add_argument("-sf", "--sigma_frac", default=0.3, type=float,
	help="[float] Fractional detection limit for neighboring pixels.\
	For cosmic ray neighbor pixels, a lapacian-to-noise detection\
	limit of sigfrac * sigclip will be used. Default: 0.3"
)
parser.add_argument("-ol", "--obj_limit", default=5.0, type=float,
	help="[float] Minimum contrast between Laplacian image and the fine\
	structure image. Increase this value if cores of bright stars are\
	flagged as cosmic rays. Default: 5.0"
)
parser.add_argument("-in", "--iteration_number", default=4, type=int,
	help="[int] Number of iterations of the LA Cosmic algorithm to\
	perform. Default: 4"
)
parser.add_argument("-sm", "--seperable_median", default=True, type=bool,
	help="[boolean] Use the separable median filter instead of the full\
	median filter. The separable median is not identical to the full\
	median filter, but they are approximately the same and the\
	separable median filter is significantly faster and still detects\
	cosmic rays well. Default: True"
)
parser.add_argument("-ct", "--data_clean_type", default="meanmask", type=str, 
	help="[str] {'median', 'medmask', 'meanmask', 'idw'} Set which\
	clean algorithm is used: ('median': An umasked 5x5 median filter),\
	('medmask': A masked 5x5 median filter), ('meanmask': A masked 5x5\
	mean filter), ('idw': A masked 5x5 inverse distance weighted\
	interpolation). Default: 'meanmask'"
)
parser.add_argument("-fsm", "--fine_structure_mode", default="convolve", 
	type=str, help="[str] {'median', 'convolve'} Method to build the\
	fine structure image: ('median': Use the median filter in the\
	standard LA Cosmic algorithm), ('convolve': Convolve the image with\
	the psf kernel to calculate the fine structure image).\
	Default: 'convolve'"
)
parser.add_argument("-v", "--verbose", default=False, type=bool,
	help="[boolean] Set whether or not detect_cosmics() will print to\
	screen. Default: False"
)
args = parser.parse_args()


# List .fits files in current directory.
files = sorted(glob.glob('*.fits'))

# Create dictionary that scrap.py will use to call instrument specific
# data preparation functions.
instrument_prep = {
	"GMOS-N": prep_gmos,
	"GMOS-S": prep_gmos
}

# Create dictionary that scrap.py will use to call instrument specific
# data saving functions.
instrument_save = {
	"GMOS-N": save_gmos,
	"GMOS-S": save_gmos
}

# Run detect_cosmics on every fits file listed in files.
for f in files:
	with fits.open(f) as spectrum_file:
		primary_header = spectrum_file[0].header
		instrument = primary_header["INSTRUME"]
		
		# Based on the value of instrument, this calls a prep_instrument
		# function
		detect_cosmics_parameters = instrument_prep[instrument](
			spectrum_file, primary_header, args.fine_structure_mode
		)
		
		# Use Astroscrappy to detect, mask, and clean cosmic rays.
		cosmic_ray_mask, clean_science_data = asc.detect_cosmics(
			detect_cosmics_parameters["in_data_frame"],
			inmask    = detect_cosmics_parameters["in_quality_frame"],
			inbkg     = detect_cosmics_parameters["inbkgd"],
			invar     = detect_cosmics_parameters["invari"],
			sigclip   = args.sigmaClip,
			sigfrac   = args.sigmaFrac,
			objlim    = args.objLimit,
			gain      = detect_cosmics_parameters["adgain"],
			readnoise = detect_cosmics_parameters["readns"],
			satlevel  = detect_cosmics_parameters["satlvl"],
			niter     = args.numIter,
			sepmed    = args.separatedMed,
			cleantype = args.dataCleanType,
			fsmode    = args.fine_structure_mode,
			psfmodel  = detect_cosmics_parameters["pmodel"],
			psffwhm   = detect_cosmics_parameters["fwhm"],
			psfsize   = detect_cosmics_parameters["psfsiz"],
			verbose   = args.verbose
		)
		
		# Save cleaned science frame, and save cosmic_ray_mask as part
		# of the quality frame.
		instrument_save[instrument](
			f, spectrum_file, primary_header, cosmic_ray_mask*1, clean_science_data, detect_cosmics_parameters, args
		)
