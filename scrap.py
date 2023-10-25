#! /home/tom/anaconda3/envs/work/bin/python

"""
scrap.py - written by T. Seccull, 2023-10-24 - v1.0.1 

	Locate, mask, and clean cosmic ray hits in 2D spectroscopic data with Astroscrappy/LA Cosmic. A 
	script that runs the detect_cosmics() function from Astroscrappy on supplied astronomical 
	spectroscopic data. Astroscrappy is a Python implentation of Pieter van Dokkum's LA Cosmic. 
	Cite both Astroscrappy and LA Cosmic if used. scrap.py will assume all files in the current 
	directory are .fits formatted 2D spectra that need their cosmic rays masked and will try to 
	apply detect_cosmics() to each in turn. scrap.py will replace the primary input data frame with 
	the cleaned data array in the output .fits file; a copy of the original input data will be 
	stored in a new Header Data Unit (HDU) of the output .fits file. The boolean cosmic ray mask 
	that scrap.py creates will be either added to an existing quality frame (e.g. flagging bad 
	pixels), or if a quality frame doesn't exist in the original file the crmask will be added to 
	the output file as the new quality frame in a new HDU. A gaussian psfmodel is assumed, so in 
	detect_cosmics() psfk=None by default, and psfbeta is ignored. Other detect_cosmics() 
	parameters are either taken directly from the fits data and headers, or they can be set 
	optionally when scrap.py is run.
    
    Original paper describing LA Cosmic 
		van Dokkum 2001, PASP, 113, 1420 - https://doi.org/10.1086/323894
    Astroscrappy Zenodo DOI
		McCully et al. 2018, Astropy/Astroscrappy: v1.0.5 Zenodo Release (v1.0.5). Zenodo
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


####################################################################################################
def moffat_least_squares(r, col, seeing, pixres, eClip):
    """
    Takes a data column, spatial axis and seeing of the observation and fits a Moffat function to 
    the column using a least squares method. Returns the best fit parameters of the Moffat 
    function.

    Args:
        r (numpy.ndarray)  : spatial axis of the data being fit
        col (numpy.ndarray): data being fitted
        seeing (float)     : estimated FWHM of the spatial profile
        pixres (float)     : spatial resolution of each pixel in arcsec/pixel
        eClip (int)        : number of pixels at each end of the spatial profile array to ignore 
                             when fitting the Moffat profile.

    Returns:
        param_list (list)  : list of best fit output parameters returned by the 
                             least squares routine.
    """

	# Clip the median spatial profile to be fitted based on the value of eClip
    col[:eClip] = np.median(col)
    col[-eClip:] = np.median(col)

    # Set up initial conditions for the least squares fit.
    # x0 = [amplitude, centre, alpha, beta, background gradient, background level]
    # Initial beta estimate comes from optimal value from atmospheric turbulence theory as 
    # described in Trujillo, I. et al. (2001), MNRAS, 328, 977-985
    # See https://ui.adsabs.harvard.edu/abs/2001MNRAS.328..977T/abstract
    x0 = [
        np.nanmedian(np.sort(col)[-3:]),
        np.argmax(col),
        seeing / pixres,
        4.765,
        0.0,
        np.median(np.concatenate((col[:5], col[-5:]))),
    ]

    # Run the least squares fit.
    res_lsq = least_squares(
        moffat_resid,
        x0,
        bounds=(
            [
                0.0, 
                np.argmax(col) - 1.0,
                0.0, 
                0.0, 
                0.0, 
                -np.inf
            ],
            [
                np.inf,
                np.argmax(col) + 1,
                (5 * seeing / pixres),
                5.0,
                np.inf,
                np.inf,
            ],
        ),
        args=(r, col),
        method="trf",
        ftol=1e-12,
    )
    
    param_list = [
        res_lsq.x[0],
        res_lsq.x[1],
        res_lsq.x[2],
        res_lsq.x[3],
        res_lsq.x[4],
        res_lsq.x[5]
    ]
    return param_list


####################################################################################################
def moffat_resid(x, datarange, data):
    """
    Calculates residuals of fitted moffat profile and the data for the least squares fitting.

    Description:
        A = x[0]
        c = x[1]
        alpha = x[2]
        beta = x[3]
        B = x[4]
        m = x[5]

    Args:
        x (numpy.ndarray)        : an array of parameters defining the shape of the model moffat 
                                   profile
        datarange (numpy.ndarray): spatial axis of the data
        data (numpy.ndarray)     : the data

    Returns:
        residual (numpy.ndarray) : the residual array between the model moffat profile and the data
    """

    moff = x[0] * ((1 + ((datarange - x[1]) * (datarange - x[1])) / (x[2] * x[2])) ** -x[3])
    
    residual = moff + x[4] + (datarange * x[5]) - data
    
    return residual

 
####################################################################################################
def prep_gmos(iFile, iHead, fSMode):
	"""
	Combines all necessary data for detect_cosmics() into a dictionary for a spectrum observed with 
	GMOS-N or GMOS-S.
	
	Args:
		iFile (.fits HDU list): the object produced by using fits.open on the fits file currently 
		                        being processed. It contains all the dataframes and headers for the 
		                        current spectrum.
		iHead (.fits header)  : the header of the primary header data unit in iFile
		fSMode (str)          : string keyword to tell prep_gmos how detect_cosmics() will
		                        generate the fine structure image.
		
	Returns:
		detCosmicsInput (dict): dictionary of dataframes and parameters collected from the input 
		                        file or calculated from it. Items in this dictionary are all inputs
		                        for detect_cosmics():
		indata (numpy.ndarray): the 2D science dataframe. The value of each pixel is in 
		                        Analog-to-Digital Units (ADU)
		inqual (numpy.ndarray): 2D frame flagging bad pixels in the science and variance frames
		inbkgd (numpy.ndarray): estimate of the 2D background in the science frame. Units are ADUs
		invari (numpy.ndarray): 2D frame containing the variance of each pixel in the science 
		                        frame; units are ADU^2
		adgain (float)        : average CCD detector gain for this data in e-/ADU
		readns (float)        : average detector readout noise e- rms
		pmodel (str)          : notes the Point Spread Function model adopted by detect_cosmics() 
		                        when building the fine structure image. "gaussy" is used here 
		                        because there is no option for a directional Moffat profile. 
		fwhm (float)          : Full Width at Half Maximum measured for the median spatial profile 
		                        of the spectrum. Units are pixels.
		psfsiz (int)          : size of the PSF model in pixels that will be convolved with the 
		                        data by detect_cosmics() to produce the fine structure model. This 
		                        value must be odd.
	"""
	
	# Take a median of the 2D science frame along the spatial axis to get a median sky spectrum, 
	# then tile that to be the same shape as the original data frame. This will be used as the 
	# estimated background.
	bgFrame = np.tile(np.median(iFile["SCI"].data, axis=0), (np.shape(iFile["SCI"])[0], 1))
	
	# Create the output dictionary and fill it with relevant dataframes and values from the input 
	# file.
	detCosmicsInput = {
		"indata": iFile["SCI"].data,
		"inqual": iFile["DQ"].data,
		"inbkgd": bgFrame,
		"invari": iFile["VAR"].data,
		"adgain": iHead["GAINMULT"],
		"readns": iHead["RDNOISE"],
		"pmodel": "gaussy"	
	}
	
	# The detectors of GMOS-N and GMOS-S are slightly different, become non-linear at different 
	# points, and therefore have different full well depths.
	if iHead["INSTRUME"] == "GMOS-N":
		detCosmicsInput["satlvl"] = 106822
	else:
		detCosmicsInput["satlvl"] = 117963

	# If the fine stucture image is to be generated with a convolution of a model PSF, estimate
	# the PSF's FWHM and size in pixels.
	if fSMode == "convolve":
		# All this is to get an initial estimate of the IQ. Tables below are based on the condition 
		# constraints used by Gemini. See web page below.
		# https://www.gemini.edu/observing/telescopes-and-sites/sites#ImageQuality
		IQ_dict = {
			"20-percentile": 0,
			"70-percentile": 1,
			"85-percentile": 2,
			"100-percentile": 3,
			"Any": 3,
			"UNKNOWN": 3,
		}
	
		WavTab = np.array(
			[
				[0000.0, 4000.0, 0],
				[4000.0, 5500.0, 1],
				[5500.0, 7000.0, 2],
				[7000.0, 8500.0, 3],
				[8500.0, 9750.0, 4],
				[9750.0, 11000.0, 5],
			]
		)
	
		IQTab = np.array(
			[
				[0.6, 0.90, 1.20, 2.00],
				[0.6, 0.85, 1.10, 1.90],
				[0.5, 0.75, 1.05, 1.80],
				[0.5, 0.75, 1.05, 1.70],
				[0.5, 0.70, 0.95, 1.70],
				[0.4, 0.70, 0.95, 1.65],
			]
		)
	
		iq = iHead["RAWIQ"]
	
		shortWav = iFile["SCI"].header["CRVAL1"]
	
		for i in WavTab:
			if shortWav > i[0] and shortWav < i[1]:
				seeing = float(IQTab[int(i[2])][int(IQ_dict[iq])])
				break
	
		pixres = iHead["PIXSCALE"]
		
		# Measure the FWHM of the spectrum's spatial profile. A Moffat profile is fitted to the 
		# median spatial profile to get an accurate measure of the Full Width at Half Maximum, even 
		# though the psf model used by detect_cosmics() is a gaussian.
		
		# Calculate median spatial profile of the spectrum.
		medProfile = np.nanmedian(iFile["SCI"].data, axis=1)
		
		# Scipy least squares doesn't like really tiny numbers like fluxes in erg/s/cm^2/Angstrom, 
		# so it's necessary to scale the data to a size that least squares can handle. The shape of 
		# the profile fitted to the scaled spatial profile is the same as the unscaled, so FWHM is 
		# unaffected.
		datascale = 10 ** np.abs(np.floor(np.log10(np.abs(np.nanmedian(medProfile)))))
		
		# Fit the median spatial profile with a Moffat function.
		moffparams = moffat_least_squares(
			range(np.shape(iFile["SCI"].data)[0]),
			medProfile * datascale,
			seeing,
			pixres,
			50
		)
		
		# Get an improved estimate of the FWHM of the spectrum from the best fit Moffat profile.
		fwhm = 2 * moffparams[2] * np.sqrt((2 ** (1 / moffparams[3])) - 1)
		
		detCosmicsInput["fwhm"] = fwhm
		
		psfSize = np.ceil(2.8*fwhm)
		if psfSize % 2 == 0:
			psfSize += 1
		
		detCosmicsInput["psfsiz"] = psfSize
		
	# If a median filter is being used to generate the fine structure image, psffwhm and psfsize
	# aren't needed. In this case we set their values to the defaults for detect_cosmics() with
	# the knowledge that they won't be used.
	else:
		detCosmicsInput["fwhm"] = 2.5
		detCosmicsInput["psfsiz"] = 7
	
	return detCosmicsInput


####################################################################################################
def save_gmos(fileName, iFile, iHead, cMask, cleanSci, dCParams, commandArgs):
	'''
	Constructs and saves a new .fits file combining the original input dataframes and headers with 
	the cleaned science data and updated quality mask.
	
	Args:
		fileName (str)           : Name of the input file
		iFile (.fits HDU list)   : the object produced by using fits.open on the fits file 
								   currently being processed. It contains all the dataframes and 
								   headers for the current spectrum.
		iHead (.fits header)     : the header of the primary header data unit in iFile
		cMask (numpy.ndarray)    : 2D array flagging the location of cosmic ray detections
		cleanSci (numpy.ndarray) : 2D science data array after cosmic rays have been cleaned
		dCParams (dict)          : dictionary of data and parameters fed to Astroscrappy 
		                           detect_cosmics()
		commandArgs (class)      : scrap.py command line argument namespace
	Returns:
		None
	'''
	
	# Update Primary header of the fits file.
	iHead.comments["IRAF-TLM"] = "Time of last modification by IRAF"
	iHead.comments["DATE"] = "Date FITS file was generated by IRAF"
	iHead["CRMETHOD"] = ("Astroscrappy v" + asc.__version__, "Cosmic ray masking/cleaning method")
	iHead["CRDATE"]   = (
		datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
		"UT timestamp for Astroscrappy"
	)
	
	# Update the SCI frame header and copy it to create the OG_SCI header. 
	sciHead = iFile["SCI"].header
	sciHead.comments["IRAF-TLM"] = "Time of last modification by IRAF"
	sciHead.comments["DATE"] = "Date FITS file was generated by IRAF"
	cleanSciHead = copy.deepcopy(sciHead)
	sciHead["EXTNAME"] = "OG_SCI"
	
	# Update the header of the new SCI header to add details on the cosmic ray detection, masking, 
	# and cleaning process.
	cleanSciHead["CRMETHOD"] = (
		"Astroscrappy v" + asc.__version__, "Cosmic ray masking/cleaning method"
	)
	cleanSciHead["CRDATE"]   = (
		datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
		"UT timestamp for Astroscrappy"
	)
	cleanSciHead["CRBKGD"]   = ("Median", "Background estimation method for Astroscrappy")
	cleanSciHead["CRSIGCLP"] = (commandArgs.sigmaClip, "Astroscrappy sigclip value")
	cleanSciHead["CRSIGFRC"] = (commandArgs.sigmaFrac, "Astroscrappy sigfrac value")
	cleanSciHead["CROBJLIM"] = (commandArgs.objLimit, "Astroscrappy objlim value")
	cleanSciHead["CRDETSAT"] = (dCParams["satlvl"], "Astroscrappy satlevel value (e-)")
	cleanSciHead["CRNITER"]  = (commandArgs.numIter, "Astroscrappy niter value")
	cleanSciHead["CRSEPMED"] = (commandArgs.separatedMed, "Astroscrappy sepmed value")
	cleanSciHead["CRDCTYPE"] = (commandArgs.dataCleanType, "Astroscrappy cleantype value")
	cleanSciHead["CRFSMODE"] = (commandArgs.finStrucMode, "Astroscrappy fsmode value")
	if commandArgs.finStrucMode == "convolve": # If fsmode is "median", then no psf parameters.
		cleanSciHead["CRPSFMOD"] = (dCParams["pmodel"], "Astroscrappy psfmodel value")
		cleanSciHead["CRPSFWHM"] = (dCParams["fwhm"], "Astroscrappy psffwhm value (pix)")
		cleanSciHead["CRPSFSIZ"] = (dCParams["psfsiz"], "Astroscrappy psfsize value (pix)")
	
	# Update the header of the new DQ header to add details on the cosmic ray detection, masking, 
	# and cleaning process. Add the cosmic ray mask to the DQ frame. 
	qualHead = iFile["DQ"].header
	qualData = iFile["DQ"].data + cMask
	
	qualHead.comments["IRAF-TLM"] = "Time of last modification by IRAF"
	qualHead.comments["DATE"] = "Date FITS file was generated by IRAF"
	qualHead["CRMETHOD"] = (
		"Astroscrappy v" + asc.__version__, "Cosmic ray masking/cleaning method"
	)
	qualHead["CRDATE"]   = (
		datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
		"UT timestamp for Astroscrappy"
	)
	qualHead["CRBKGD"]   = ("Median", "Background estimation method for Astroscrappy")
	qualHead["CRSIGCLP"] = (commandArgs.sigmaClip, "Astroscrappy sigclip value")
	qualHead["CRSIGFRC"] = (commandArgs.sigmaFrac, "Astroscrappy sigfrac value")
	qualHead["CROBJLIM"] = (commandArgs.objLimit, "Astroscrappy objlim value")
	qualHead["CRDETSAT"] = (dCParams["satlvl"], "Astroscrappy satlevel value (e-)")
	qualHead["CRNITER"]  = (commandArgs.numIter, "Astroscrappy niter value")
	qualHead["CRSEPMED"] = (commandArgs.separatedMed, "Astroscrappy sepmed value")
	qualHead["CRDCTYPE"] = (commandArgs.dataCleanType, "Astroscrappy cleantype value")
	qualHead["CRFSMODE"] = (commandArgs.finStrucMode, "Astroscrappy fsmode value")
	if commandArgs.finStrucMode == "convolve": # If fsmode is "median", then no psf parameters.
		qualHead["CRPSFMOD"] = (dCParams["pmodel"], "Astroscrappy psfmodel value")
		qualHead["CRPSFWHM"] = (dCParams["fwhm"], "Astroscrappy psffwhm value (pix)")
		qualHead["CRPSFSIZ"] = (dCParams["psfsiz"], "Astroscrappy psfsize value (pix)")
	
	# Construct the output .fits file.
	iHDU = fits.PrimaryHDU(header=iHead)
	mdfHDU = iFile["MDF"]
	cleanSciHDU = fits.ImageHDU(cleanSci, header=cleanSciHead)
	varHDU = iFile["VAR"]
	qualHDU = fits.ImageHDU(qualData, header=qualHead)
	sciHDU = iFile["OG_SCI"]
	hduList=fits.HDUList([iHDU, mdfHDU, cleanSciHDU, varHDU, qualHDU, sciHDU])
	hduList.writeto("c" + fileName)
	hduList.close()


####################################################################################################
#### SCRIPT STARTS HERE # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #### 
####################################################################################################

# Parse scrap.py arguments
parser = argparse.ArgumentParser(
	description="Locate, mask, and clean cosmic ray hits in 2D spectroscopic data with\
	             Astroscrappy/LACosmic. A script that runs the detect_cosmics() function from\
	             Astroscrappy on supplied astronomical spectroscopic data. Astroscrappy is a\
	             Python implentation of Pieter van Dokkum's LACosmic. Cite both Astroscrappy and\
	             LACosmic if used. scrap.py will assume all files in the current directory\
	             are .fits formatted 2D spectra that need their cosmic rays masked and will try to\
	             apply detect_cosmics() to each in turn. scrap.py will replace the primary input\
	             data frame with the cleaned data array in the output .fits file; a copy of the\
	             original input data will be stored in a new Header Data Unit (HDU) of the output\
	             .fits file. The boolean cosmic ray mask that scrap.py creates will be either\
	             added to an existing quality frame (e.g. flagging bad pixels), or if a quality\
	             frame doesn't exist in the original file the crmask will be added to the output\
	             file as the new quality frame in a new HDU. A gaussian psfmodel is assumed, so in\
	             detect_cosmics() psfk=None by default, and psfbeta is ignored. Other\
	             detect_cosmics() parameters are either taken directly from the fits data and\
	             headers, or they can be set optionally when scrap.py is run. Astroscrappy docs\
	             and links to citables -> https://astroscrappy.readthedocs.io/en/latest/"
)
parser.add_argument("-sc", "--sigmaClip", default=4.5, type=float, help="[float]\
    Laplacian-to-noise limit for cosmic ray detection. Lower values will flag more pixels as\
    cosmic rays. Default: 4.5"
)
parser.add_argument("-sf", "--sigmaFrac", default=0.3, type=float, help="[float] Fractional\
    detection limit for neighboring pixels. For cosmic ray neighbor pixels, a lapacian-to-noise\
    detection limit of sigfrac * sigclip will be used. Default: 0.3"
)
parser.add_argument("-ol", "--objLimit", default=5.0, type=float, help="[float] Minimum contrast\
    between Laplacian image and the fine structure image. Increase this value if cores of bright\
    stars are flagged as cosmic rays. Default: 5.0"
)
parser.add_argument("-ni", "--numIter", default=4, type=int, help="[int] Number of iterations of\
    the LA Cosmic algorithm to perform. Default: 4"
)
parser.add_argument("-sm", "--separatedMed", default=True, type=bool, help="[boolean] Use the\
    separable median filter instead of the full median filter. The separable median is not\
    identical to the full median filter, but they are approximately the same and the separable\
    median filter is significantly faster and still detects cosmic rays well. Default: True"
)
parser.add_argument("-ct", "--dataCleanType", default="meanmask", type=str, help="\
    [str] {'median', 'medmask', 'meanmask', 'idw'} Set which clean algorithm is used:\
    ('median': An umasked 5x5 median filter), ('medmask': A masked 5x5 median filter),\
    ('meanmask': A masked 5x5 mean filter),\
    ('idw': A masked 5x5 inverse distance weighted interpolation.)\
    Default: 'meanmask'"
)
parser.add_argument("-fsm", "--finStrucMode", default="convolve", type=str, help="\
    [str] {'median', 'convolve'} Method to build the fine structure image:\
    ('median': Use the median filter in the standard LA Cosmic algorithm),\
    ('convolve': Convolve the image with the psf kernel to calculate the fine structure image),\
    Default: 'convolve'"
)
parser.add_argument("-v", "--verbose", default=False, type=bool, help="[boolean] Set whether or\
	not detect_cosmics() will print to screen. Default: False"
)
args = parser.parse_args()


# List .fits files in current directory.
#files = sorted(glob.glob('*.fits'))
files = glob.glob('*.fits')

# Create dictionary that scrap.py will use to call instrument specific data preparation functions.
instrument_prep = {
	"GMOS-N": prep_gmos,
	"GMOS-S": prep_gmos
}

# Create dictionary that scrap.py will use to call instrument specific data saving functions.
instrument_save = {
	"GMOS-N": save_gmos,
	"GMOS-S": save_gmos
}

# Run detect_cosmics on every fits file listed in files.
for f in files:
	with fits.open(f) as imgFile:
		imgHead = imgFile[0].header
		inst = imgHead["INSTRUME"]
		
		# Based on the value of inst, this calls a prep_instrument function
		detCosmicParams = instrument_prep[inst](imgFile, imgHead, args.finStrucMode)
		
		# Use Astroscrappy to detect, mask, and clean cosmic rays.
		crMask, cleanData = asc.detect_cosmics(
			detCosmicParams["indata"],
			inmask    = detCosmicParams["inqual"],
			inbkg     = detCosmicParams["inbkgd"],
			invar     = detCosmicParams["invari"],
			sigclip   = args.sigmaClip,
			sigfrac   = args.sigmaFrac,
			objlim    = args.objLimit,
			gain      = detCosmicParams["adgain"],
			readnoise = detCosmicParams["readns"],
			satlevel  = detCosmicParams["satlvl"],
			niter     = args.numIter,
			sepmed    = args.separatedMed,
			cleantype = args.dataCleanType,
			fsmode    = args.finStrucMode,
			psfmodel  = detCosmicParams["pmodel"],
			psffwhm   = detCosmicParams["fwhm"],
			psfsize   = detCosmicParams["psfsiz"],
			verbose   = args.verbose
		)
		
		# Save cleaned science frame, and save crMask as part of the quality frame
		instrument_save[inst](f, imgFile, imgHead, crMask*1, cleanData, detCosmicParams, args)
	
