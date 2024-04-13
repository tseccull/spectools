#! /home/tom/anaconda3/envs/work/bin/python
"""
fronge.py - written by T. Seccull, 2023-10-27 - v1.0.0

	This script is designed to act on 2D spectroscopic data where the target has been observed in
	a sequence of offset dither patterns along the spectrographic slit. For the science frame in
	each file, fronge.py will search for all science frames that have a different dither position.
	The median of these science frames is taken to construct a fringe frame that is then subtracted
	from the science frame currently being processed. The uncertainties of the values in the fringe
	frame are estimated by taking the median absolute deviation of the set of values
	median-combined into each pixel of the fringe frame. fringe.py will save a new file for each
	that it processes that will contain the fringe corrected science frame, updated uncertainty
	frame, the fringe frame, and the fringe frame's uncertainty frame. fringe.py will assume that
	all .fits files in the current directory are intended for processing, and will attempt to 
	fringe correct each one. This script has not optional arguments except -h, which will show this
	description.
"""


import argparse
import astropy.io.fits as fits
import copy
import datetime
import glob
import numpy as np


####################################################################################################
def save_gmos(file_string, fringe_frame, mad_frame, other_dither_data):
	"""
	Constructs and saves a new .fits file combining the original input dataframes and headers with 
	the defringed science frame, updated variance frame, the new fringe frame, and the fringe
	frame's uncertainty frame.
	
	Args:
		file_string (str)                : name for the current frame (i.e. the file name with the 
		                           ".fits" cut off the end)
		fringe_frame (numpy.ndarray)  : 2D fringe frame array
		mad_frame (numpy.ndarray)  : 2D array of median absolute deviation (i.e. estimated 
		                           uncertainty) values for the fringe frame array 
		other_dither_data (list)            : list of file keywords for all the files that were median 
		                           combined to make the fringe frame
	Returns:
		None
	"""
	
	# Open the current fits file, extract relevant data frames and metadata, combine it with the
	# new data frames and metadata and construct a new fits HDUList object that will be written to
	# a new .fits file.
	with fits.open(file_string+".fits") as iFile:
		varHead = iFile["VAR"].header
		ogVarHead = copy.deepcopy(varHead)
		ogVarHead["EXTNAME"] = "OG_VAR"
		iFile.append(fits.ImageHDU(copy.deepcopy(iFile["VAR"].data), header=ogVarHead))
		iFile["SCI"].data -= fringe_frame
		iFile["VAR"].data += (mad_frame*mad_frame)
		iFile.append(fits.ImageHDU(fringe_frame))
		iFile[-1].header = copy.deepcopy(iFile["OG_SCI"].header)
		iFile[-1].header["EXTNAME"] = "FRINGE_FRAME"
		iFile.append(fits.ImageHDU(mad_frame))
		iFile[-1].header = copy.deepcopy(iFile["OG_SCI"].header)
		iFile[-1].header["EXTNAME"] = "MAD_FRINGE_FRAME"
		newFile = copy.deepcopy(iFile)
	
		
	newFile[0].header["FRNGDATE"] = (datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
		"UT timestamp for fringe frame subtraction"
	)
	fFrames = ["SCI", "VAR", "FRINGE_FRAME", "MAD_FRINGE_FRAME"]
	for frame in fFrames:
		newFile[frame].header["FRNGDATE"] = (
			datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
			"UT timestamp for Fringeframe subtraction"
		)
		for i, fName in enumerate(other_dither_data):
			newFile[frame].header["FRNGIN"+str(i+1)] = (
				fName+".fits", "SCI medianed into fringe frame"
			)
	
	newFile.writeto("f" + file_string + ".fits")
	newFile.close()


####################################################################################################
#### SCRIPT STARTS HERE # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #### 
####################################################################################################

# Parse fronge.py help message.
parser = argparse.ArgumentParser(
	description="This script is designed to act on 2D spectroscopic data where the target has been\
	observed in a sequence of offset dither patterns along the spectrographic slit. For the\
	science frame in each file, fronge.py will search for all science frames that have a different\
	dither position. The median of these science frames is taken to construct a fringe frame that\
	is then subtracted from the science frame currently being processed. The uncertainties of the\
	values in the fringe frame are estimated by taking the median absolute deviation of the set of\
	values median-combined into each pixel of the fringe frame. fringe.py will save a new file for\
	each that it processes that will contain the fringe corrected science frame, updated\
	uncertainty frame, the fringe frame, and the fringe frame's uncertainty frame. fringe.py will\
	assume that all .fits files in the current directory are intended for processing, and will\
	attempt to fringe correct each one. This script has no optional arguments except -h, which\
	will show this description."
)
args = parser.parse_args()

# Grab all fits files in the current directory.
files = sorted(glob.glob("*.fits"))

# Set the science extension name in the headers of fits files from different instruments.
instDataHDU = {
	"GMOS-N": "SCI",
	"GMOS-S": "SCI"
}

# Set the keyword that notes the offset in arcseconds along the slit relative to the zero position
# for data observed with different instruments.
instOffsetKeyword = {
	"GMOS-N": "YOFFSET",
	"GMOS-S": "YOFFSET"
}

# Create dictionary that scrap.py will use to call instrument specific data saving functions.
instrument_save = {
	"GMOS-N": save_gmos,
	"GMOS-S": save_gmos
}

# All these dictionaries use the same keywords, but point to different variables. In all cases the
# keywords are just the name of the .fits file with the ".fits" removed from the end.
dataFrames = {}        # Dictionary of 2D science frames
ditherPoints  = {}     # Dictionary of offset values along the slit.
otherDitherPoints = {} # Dictionary containing a list of all other offset values relative to that
                       # of the current file.

# Loop through input files, collecting the relevant data frames and metadata needed for creation of
# a fringe frame for each science frame.
for f in files:
	with fits.open(f) as imgFile:
		# Get science frame and header metadata.
		imgHead = imgFile[0].header
		inst = imgHead["INSTRUME"]
		sciFrame = imgFile[instDataHDU[inst]].data
		
		# Subtract median spatial background from the science frame to remove sky emissionlines.
		medBackground = np.tile(np.nanmedian(sciFrame, axis=0), (np.shape(sciFrame)[0], 1))
		dataFrames[f[:-5]] = sciFrame - medBackground
		
		# Save the offset of the current science frame along the slit.
		ditherPoints[f[:-5]] = round(imgHead[instOffsetKeyword[inst]], 1)

# Make a list containing all dither positions along the slit represented by the current file list.
allDithers = []
[allDithers.append(x) for x in list(ditherPoints.values()) if x not in allDithers]

# For each science frame:
for k in dataFrames:
	# Get the dither points in all dither points that aren't the dither value of the current frame.
	otherDitherPoints[k] = [x for x in allDithers if x!=ditherPoints[k]]
	
	# Collect all science frames with dithers other than that of the current frame.
	otherDitherData = [x for x in dataFrames if ditherPoints[x] in otherDitherPoints[k]]
	otherDitherFrames = np.array([dataFrames[x] for x in otherDitherData])
	
	# Create the fringe frame by median combining the frames in the otherDitherFrames
	fringeFrame = np.nanmedian(otherDitherFrames, axis=0)
	# Create the uncertainty frame of the fringe frame by estimating its median absolute deviation
	madDitherFrames = np.abs(otherDitherFrames - fringeFrame)
	madFrame = np.nanmedian(madDitherFrames, axis=0)
	
	# Subtract the fringe frame from the science frame and save the fringe-corrected data.	
	instrument_save[inst](k, fringeFrame, madFrame, otherDitherData)
	
