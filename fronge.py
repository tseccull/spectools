#! /home/tom/anaconda3/envs/work/bin/python
"""
fronge.py - written by T. Seccull, 2023-10-27 - v1.0.1

    Last updated: 2024-04-13

	This script is designed to act on 2D spectroscopic data where the
	target has been observed in a sequence of offset dither patterns
	along the spectrographic slit. For the science frame in each file,
	fronge.py will search for all science frames that have a different
	dither position. The median of these science frames is taken to
	construct a fringe frame that is then subtracted from the science
	frame currently being processed. The uncertainties of the values in
	the fringe frame are estimated by taking the median absolute
	deviation of the set of values median-combined into each pixel of
	the fringe frame. fringe.py will save a new file for each that it
	processes that will contain the fringe corrected science frame,
	updated uncertainty frame, the fringe frame, and the fringe frame's
	uncertainty frame. fringe.py will assume that all .fits files in the
	current directory are intended for processing, and will attempt to
	fringe correct each one. This script has not optional arguments
	except -h, which will show this description.
"""


import argparse
import astropy.io.fits as fits
import copy
import datetime
import glob
import numpy as np


###############################################################################
def save_gmos(file_string, fringe_frame, mad_frame, other_dither_data):
	"""
	Constructs and saves a new .fits file combining the original input
	dataframes and headers with the defringed science frame, updated
	variance frame, the new fringe frame, and the fringe frame's
	uncertainty frame.
	
	Args:
	 --	file_string (str) - A name for the current frame (i.e. the file
		    name with the ".fits" cut off the end).	
	 --	fringe_frame (numpy.ndarray) - A 2D fringe frame array.
	 --	mad_frame (numpy.ndarray) - A 2D array of median absolute 
		    deviation (i.e. estimated uncertainty) values for the fringe
		    frame array.
	 --	other_dither_data (list) - A list of file keywords for all the
		    files that were median combined to make the fringe frame
	
	Returns:
	 --	None
	"""
	
	with fits.open(file_string+".fits") as in_file_hdu_list:
		new_file_hdu_list = copy.deepcopy(in_file_hdu_list)
	
	extensions = []
	for hdu in new_file_hdu_list:
		if "XTENSION" in hdu.header:
			extensions.append(hdu.header["EXTNAME"])
			
	if "OG_SCI" not in extensions:
		new_file_hdu_list.append(copy.deepcopy(new_file_hdu_list["SCI"]))
		new_file_hdu_list[-1].header["EXTNAME"] = "OG_SCI"
	
	if "OG_VAR" not in extensions:
		new_file_hdu_list.append(copy.deepcopy(new_file_hdu_list["VAR"]))
		new_file_hdu_list[-1].header["EXTNAME"] = "OG_VAR"
	
	new_file_hdu_list["SCI"].data -= fringe_frame
	new_file_hdu_list.append(fits.ImageHDU(fringe_frame))
	new_file_hdu_list[-1].header["EXTNAME"] = "FRINGE_FRAME"
	
	new_file_hdu_list["VAR"].data += (mad_frame*mad_frame)
	new_file_hdu_list.append(fits.ImageHDU(mad_frame))
	new_file_hdu_list[-1].header["EXTNAME"] = "MAD_FRINGE_FRAME"
	
	new_file_hdu_list[0].header["FRNGDATE"] = (
		datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
		"UT timestamp for fringe frame subtraction"
	)
	
	frame_list = ["SCI", "VAR", "FRINGE_FRAME", "MAD_FRINGE_FRAME"]
	for frame in frame_list:
		new_file_hdu_list[frame].header["FRNGDATE"] = (
			datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
			"UT timestamp for Fringeframe subtraction"
		)
		for i, file_name in enumerate(other_dither_data):
			new_file_hdu_list[frame].header["FRNGIN"+str(i+1)] = (
				file_name+".fits", "SCI medianed into fringe frame"
			)
	
	new_file_hdu_list.writeto("f" + file_string + ".fits")
	new_file_hdu_list.close()


###############################################################################
#### SCRIPT STARTS HERE  # # # # # # # # # # # # # # # # # # # # # # # # # #### 
###############################################################################

# Parse fronge.py help message.
parser = argparse.ArgumentParser(
	description="This script is designed to act on 2D spectroscopic\
	data where the target has been observed in a sequence of offset\
	dither patterns along the spectrographic slit. For the science\
	frame in each file, fronge.py will search for all science frames\
	that have a different dither position. The median of these science\
	frames is taken to construct a fringe frame that is then subtracted\
	from the science frame currently being processed. The uncertainties\
	of the values in the fringe frame are estimated by taking the\
	median absolute deviation of the set of values median-combined into\
	each pixel of the fringe frame. fringe.py will save a new file for\
	each that it processes that will contain the fringe corrected\
	science frame, updated uncertainty frame, the fringe frame, and the\
	fringe frame's uncertainty frame. fringe.py will assume that all\
	.fits files in the current directory are intended for processing,\
	and will attempt to fringe correct each one. This script has no\
	optional arguments except -h, which will show this description."
)
args = parser.parse_args()

# Grab all fits files in the current directory.
files = sorted(glob.glob("*.fits"))

# Set the science extension name in the headers of fits files from
# different instruments.
instrument_data_hdu = {
	"GMOS-N": "SCI",
	"GMOS-S": "SCI"
}

# Set the keyword that notes the offset in arcseconds along the slit
# relative to the zero position for data observed with different
# instruments.
instrument_offset_keyword = {
	"GMOS-N": "YOFFSET",
	"GMOS-S": "YOFFSET"
}

# Create dictionary that scrap.py will use to call instrument specific
# data saving functions.
instrument_save = {
	"GMOS-N": save_gmos,
	"GMOS-S": save_gmos
}

# All these dictionaries use the same keywords, but point to different
# variables. In all cases the keywords are just the name of the .fits
# file with the ".fits" removed from the end.
data_frames = {}        # Dictionary of 2D science frames.
dither_points  = {}     # Dictionary of offset values along the slit.

# Dictionary containing a list of all other offset values relative to
# that of the current file.
other_dither_points = {}

# Loop through input files, collecting the relevant data frames and
# metadata needed for creation of a fringe frame for each science frame.
for file_name in files:
	with fits.open(file_name) as image_file:
		# Get science frame and header metadata.
		primary_header = image_file[0].header
		instrument = primary_header["INSTRUME"]
		science_frame = image_file[instrument_data_hdu[instrument]].data
		
	# Subtract median spatial background from the science frame to
	# remove sky emissionlines.
	background_1d = np.nanmedian(science_frame, axis=0)
	background_2d = np.tile(background_1d, (np.shape(science_frame)[0], 1))
	data_frames[file_name[:-5]] = science_frame - background_2d
	
	# Save the offset of the current science frame along the slit.
	dither_points[file_name[:-5]] = round(
		primary_header[instrument_offset_keyword[instrument]], 1
	)

# Make a list containing all dither positions along the slit represented
# by the current file list.
all_dithers = []
for dither_point_value in list(dither_points.values()):
	if dither_point_value not in all_dithers:
		all_dithers.append(dither_point_value)

# For each science frame:
for k in data_frames:
	# Get the dither points in all dither points that aren't the dither
	# value of the current frame.
	other_dither_points[k] = [x for x in all_dithers if x!=dither_points[k]]
	
	# Collect all science frames with dithers other than that of the
	# current frame.
	other_dither_data = [
		x for x in data_frames if dither_points[x] in other_dither_points[k]
	]
	other_dither_frames = np.array([data_frames[x] for x in other_dither_data])
	
	# Create the fringe frame by median combining the frames in the
	# other_dither_frames
	fringe_frame = np.nanmedian(other_dither_frames, axis=0)
	# Create the uncertainty frame of the fringe frame by estimating its
	# median absolute deviation
	mad_dither_frames = np.abs(other_dither_frames - fringe_frame)
	mad_frame = np.nanmedian(mad_dither_frames, axis=0)
	
	# Subtract the fringe frame from the science frame and save the
	# fringe-corrected data.	
	instrument_save[instrument](k, fringe_frame, mad_frame, other_dither_data)
