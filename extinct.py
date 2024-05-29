#! /home/tom/anaconda3/envs/work/bin/python
"""
extinct.py - written by Tom Seccull, 2024-05-29 - v0.0.1

	Last updated: 2024-05-29
	
	This script applies a correction to input 1D astronomical 
	spectroscopic data based on the measured atmospheric extinction
	curve of the observing site. The spectrum and its uncertainties are
	each multiplied by 10 ^ (0.4 * airmass * k(lambda)), where airmass
	is the median airmass at which the spectrum was observed and
	k(lambda) is the wavelength dependent extinction curve of the
	observing site.
"""

import argparse
import astropy.io.fits as fits
import glob
import matplotlib.pyplot as plt
import numpy as np
import os

parser = argparse.ArgumentParser(
	description="This script applies a correction to input 1D\
	astronomical spectroscopic data based on the measured atmospheric\
	extinction curve of the observing site. The spectrum and its\
	uncertainties are each multiplied by\
	10 ^ (0.4 * airmass * k(lambda)), where airmass is the median\
	airmass at which the spectrum was observed and k(lambda) is the\
	wavelength dependent extinction curve of the observing site."
)
parser.add_argument("-p", "--plot", action="store_true",
	help="[boolean] Plot each spectrum before and after extinction\
	correction on the same plot."
)
parser.add_argument("-s", "--save", action="store_true",
	help="[boolean] Save the extinction corrected spectrum to a new\
	file."
)
args = parser.parse_args()

script_directory = os.path.abspath(os.path.dirname(__file__))
extinction_curve_directory = script_directory + "/extinct/"
extinction_functions = {
	"GMOS-N" : "mkoextinct.dat",
	"GMOS-S" : "ctioextinct.dat"
}

files = sorted(glob.glob("*.fits"))

instrument = None

for i, file_name in enumerate(files):
	with fits.open(file_name) as file_hdu:
		
