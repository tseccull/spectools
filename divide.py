#! /home/tom/anaconda3/envs/work/bin/python
"""
divide.py - written by Tom Seccull, 2024-05-08 - v0.0.1

	Last updated: 2024-05-09
	
	This script divides one 1D spectrum by another. divide.py expects
	both spectra to have a common wavelength axis, be the product
	of extraction by MOTES, and stacking by stack.py. 
"""

import argparse
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(
	description="This script divides one 1D spectrum by another.\
	divide.py expects both spectra to have a common wavelength axis, be\
	the product of extraction by MOTES, and stacking by stack.py."
)
parser.add_argument("spec_file_one", type=str,
	help="[str] File name of the first spectrum, which takes the 'numerator'\
	role in the division." 
)
parser.add_argument("spec_file_two", type=str,
	help="[str] File name of the second spectrum, which takes the\
	'denominator' role in the division." 
)
parser.add_argument("-p", "--plot", action="store_true",
	help="[boolean] Plot the ratio spectrum on screen."
)
parser.add_argument("-s", "--save", action="store_true",
	help="[boolean] Save the ratio spectrum to a new .fits file."
)
args = parser.parse_args()

with fits.open(args.spec_file_one) as one:
	one_hdu_list = [x for x in one]

with fits.open(args.spec_file_two) as two:
	two_hdu_list = [x for x in two]
	
[print(x.header["EXTNAME"]) for x in one_hdu_list]
[print(x.header["EXTNAME"]) for x in two_hdu_list]
