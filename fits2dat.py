#!/usr/bin/env python3

"""
	fits2dat.py

	Copyright (C) 2025-04-02 Tom Seccull
	
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

	Last updated - 2025-04-02

	Description --------------------------------------------------------	
	This script repackages 1D spectra extracted with MOTES and saved in 
	.fits format into files in a csv format. Header metadata is stored 
	in at the top of the file with each line starting "#".
"""

__version__ = "1.0.0"
__author__ = "Tom Seccull"

import astropy.io.fits as fits
import glob
import numpy as np

def write_cards(txt_file, header):
	"""
	This function converts each FITS header card to a string that is
	printed to the header of the csv formatted .dat file.
	
	Args:
	 -- txt_file (_io.TextIOWrapper)
			Open file object for the .dat file currently being written.
	 -- header (astropy.io.fits.header.Header)
			Astropy header object containing header cards to be written.
	
	Returns: None
	"""
	
	for card in header:
		card_value = str(header[card])
		card_comment = header.comments[card]
		card_string = "# " + card + " = " + card_value + " / " + card_comment + "\n"
		txt_file.write(card_string)
	
	return None

###############################################################################
#### SCRIPT STARTS HERE  # # # # # # # # # # # # # # # # # # # # # # # # # #### 
###############################################################################	

files = sorted(glob.glob("*.fits"))
hdu_extnames = ["/", "STACK"]
dashes = "# -------- "
end = " Header\n"
cols = "# Data columns correspond to HDUROW values in the top header of this file.\n"

# This bit is fugly, but it works.
for fits_file in files:
	with fits.open(fits_file) as fits_hdu_list:
		for hdu in fits_hdu_list:
			if any([x in hdu.header["EXTNAME"] for x in hdu_extnames]):
				extname = hdu.header["EXTNAME"]
				with open(fits_file[:-5] + ".dat", "a") as txt_file:
					txt_file.write(dashes + fits_file + " - " + extname + end)
					write_cards(txt_file, hdu.header)
					txt_file.write("#\n")

		with open(fits_file[:-5] + ".dat", "a") as txt_file:
			txt_file.write("# " + fits_file + " - Data\n")
			txt_file.write(cols)
			for row in fits_hdu_list[0].data.T:
				row_string = []
				[row_string.append(str(x) + ", ") for x in row]
				row_string = "".join(row_string)[:-2] + "\n"
				txt_file.write(row_string)
