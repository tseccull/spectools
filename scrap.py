#! /home/tom/miniforge3/envs/work/bin/python

"""
scrap.py - written by T. Seccull, 2024-05-06 - v1.0.7

	Last updated - 2024-08-10

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
import glob
import scrap.gmosio as gmosio


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
parser.add_argument("-sm", "--separable_median", default=True, type=bool,
	help="[boolean] Use the separable median filter instead of the full\
	median filter. The separable median is not identical to the full\
	median filter, but they are approximately the same and the\
	separable median filter is significantly faster and still detects\
	cosmic rays well. Default: True"
)
parser.add_argument("-ct", "--data_clean_type", default="meanmask", type=str, 
	help="[str] {'median', 'medmask', 'meanmask', 'idw'} Set which\
	clean algorithm is used: ('median': An unmasked 5x5 median filter),\
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
	"GMOS-N": gmosio.prep_gmos,
	"GMOS-S": gmosio.prep_gmos
}

# Create dictionary that scrap.py will use to call instrument specific
# data saving functions.
instrument_save = {
	"GMOS-N": gmosio.save_gmos,
	"GMOS-S": gmosio.save_gmos
}

# Run detect_cosmics on every fits file listed in files.
for f in files:
	with fits.open(f) as spectrum_file:
		primary_header = spectrum_file[0].header
		instrument = primary_header["INSTRUME"]
		
		# Based on the value of instrument, this calls a prep_instrument
		# function
		detect_cosmics_parameters = instrument_prep[instrument](
			spectrum_file, args.fine_structure_mode
		)
		
		# Use Astroscrappy to detect, mask, and clean cosmic rays.
		cosmic_ray_mask, clean_science_data = asc.detect_cosmics(
			detect_cosmics_parameters["in_data_frame"],
			inmask    = detect_cosmics_parameters["in_quality_frame"],
			inbkg     = detect_cosmics_parameters["in_background_frame"],
			invar     = detect_cosmics_parameters["in_variance_frame"],
			sigclip   = args.sigma_clip,
			sigfrac   = args.sigma_frac,
			objlim    = args.obj_limit,
			gain      = detect_cosmics_parameters["detector_gain"],
			readnoise = detect_cosmics_parameters["read_noise"],
			satlevel  = detect_cosmics_parameters["saturation_level"],
			niter     = args.iteration_number,
			sepmed    = args.separable_median,
			cleantype = args.data_clean_type,
			fsmode    = args.fine_structure_mode,
			psfmodel  = detect_cosmics_parameters["psf_model"],
			psffwhm   = detect_cosmics_parameters["fwhm"],
			psfsize   = detect_cosmics_parameters["psf_size"],
			verbose   = args.verbose
		)
		
		# Save cleaned science frame, and save cosmic_ray_mask as part
		# of the quality frame.
		instrument_save[instrument](
			f,
			spectrum_file,
			primary_header,
			cosmic_ray_mask*1,
			clean_science_data,
			detect_cosmics_parameters,
			args
		)
