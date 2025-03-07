#!/usr/bin/env python3

"""
	fronge.py

	Copyright (C) 2024-05-06 Tom Seccull
	
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

	Last updated - 2025-03-07

	Description --------------------------------------------------------

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
	fringe correct each one. This script has no optional arguments
	except -h, which will show this description.
"""

__version__ = "1.0.5"
__author__ = "Tom Seccull"

import argparse
import astropy.io.fits as fits
import copy
import datetime
import fronge.gmosio as gmosio
import glob
import numpy as np


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
	"GMOS-N": gmosio.save_gmos,
	"GMOS-S": gmosio.save_gmos
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
	instrument_save[instrument](
		k, fringe_frame, mad_frame, other_dither_data, __version__
	)
