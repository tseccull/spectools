#! /home/tom/anaconda3/envs/work/bin/python
"""
gmosio.py - written by T. Seccull, 2024-05-06

	Called by fronge.py
    Last updated: 2024-05-09

	This module contains the save_gmos() function called by fronge.py.
"""

import astropy.io.fits as fits
import copy
import datetime


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
	
	new_file_hdu_list[0].header["FRNGSCPT"] = (
		"fronge.py v1.0.3", "Script used to perform fringe correction"
	)
	new_file_hdu_list[0].header["FRNGDOI"] = (
		"UNKNOWN", "DOI of fronge.py repository"
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
