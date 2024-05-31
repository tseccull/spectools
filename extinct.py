#! /home/tom/anaconda3/envs/work/bin/python
"""
extinct.py - written by Tom Seccull, 2024-05-29 - v1.0.0

	Last updated: 2024-05-31
	
	This script applies a correction to input 1D astronomical 
	spectroscopic data based on the measured atmospheric extinction
	curve of the observing site. The 1D spectrum and its uncertainties
	are each multiplied by 10 ^ (0.4 * airmass * k(lambda)), where
	airmass is the mean airmass at which the spectrum was observed and
	k(lambda) is the wavelength dependent extinction curve of the
	observing site given in magnitudes per unit airmass. Note that the
	extinction curves used here are averages that do not account for
	variable atmospheric extinction due to variable concentrations of
	atmospheric water vapour or scattering particles. When comparing or
	calibrating one corrected spectrum with another it is assumed that 
	atmospheric conditions were stable across the consecutive 
	observations of both targets.This script does not create new files,
	but instead just updates the input files.
"""


import argparse
import astropy.io.fits as fits
import datetime
import glob
import matplotlib.pyplot as plt
import numpy as np
import os


def get_extinction(extinct_curve_file):
	with open(extinct_curve_file) as ext:
		lines = ext.readlines()
	lines = [x.strip() for x in lines]
	lines = [x.split(" ") for x in lines]
	
	extinct_wave = []
	extinct_curve = []
	for i in lines:
		extinct_wave.append(float(i[0]))
		extinct_curve.append(float(i[1]))
		
	return np.array(extinct_wave), np.array(extinct_curve)
	

parser = argparse.ArgumentParser(
	description="This script applies a correction to input 1D\
	astronomical spectroscopic data based on the measured atmospheric\
	extinction curve of the observing site. The spectrum and its\
	uncertainties are each multiplied by\
	10 ^ (0.4 * airmass * k(lambda)), where airmass is the median\
	airmass at which the spectrum was observed and k(lambda) is the\
	wavelength dependent extinction curve of the observing site given\
	in magnitudes per unit airmass. Note that the extinction curves\
	used here are averages that do not account for variable atmospheric\
	extinction due to variable concentrations of atmospheric water\
	vapour or scattering particles. When comparing or calibrating one\
	corrected spectrum with another it is assumed that atmospheric\
	conditions were stable across the consecutive observations of both\
	targets. This script does not create new files, but instead just\
	updates the input files."
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

airmass_keys = {
	"GMOS-N" : "AIRMASS",
	"GMOS-S" : "AIRMASS"
}

files = sorted(glob.glob("*.fits"))

instrument = None

for i, file_name in enumerate(files):
	with fits.open(file_name, mode='update') as file_hdu:
		primary_header = file_hdu[0].header
		
		if primary_header["INSTRUME"] != instrument:
			instrument = primary_header["INSTRUME"]
			extinct_wavelength, extinct_curve = get_extinction(
				extinction_curve_directory + extinction_functions[instrument]
			)
		
		airmass = primary_header[airmass_keys[instrument]]
		wavelength_axis = file_hdu[0].data[0]
		spectrum = file_hdu[0].data[1]
		uncertainties = file_hdu[0].data[2]
		interp_extinction = np.interp(
			wavelength_axis, extinct_wavelength, extinct_curve
		)
		corrected_spectrum = spectrum * (
			10 ** (0.4 * airmass * interp_extinction)
		)
		corrected_uncertainties = uncertainties * (
			10 ** (0.4 * airmass * interp_extinction)
		)
		
		if args.plot:
			plt.errorbar(
				wavelength_axis,
				spectrum,
				yerr=uncertainties,
				label="Airmass: " + str(round(airmass, 3))
			)
			plt.errorbar(
				wavelength_axis,
				corrected_spectrum,
				yerr=corrected_uncertainties,
				label="Airmass: 1.000"
			)
			plt.grid(linestyle="dashed", alpha=0.5)
			plt.xlabel(r"$\lambda$, " + primary_header["WAVU"])
			plt.ylabel(primary_header["HDUROW1"])
			plt.legend()
			plt.title(file_name)
			plt.show() 
		
		if args.save:
			file_hdu[0].data[1] = corrected_spectrum
			file_hdu[0].data[2] = corrected_uncertainties
			file_hdu[0].header["EXTIDATE"] = (
				datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 
				"UT timestamp for extinction correction"
			)
			file_hdu[0].header["EXTISCPT"] = (
				"extinct.py v1.0.0",
				"Script used to perform extinction correction"
			)
			file_hdu[0].header["EXTIDOI"] = (
				"UNKNOWN", "DOI of extinct.py repository"
			)
			file_hdu[0].header["EXTICURV"] = (
				extinction_functions[instrument],
				"Extinction curve used in extinction correction."
			)
			file_hdu.flush()
			os.rename(file_name, "e" + file_name)		
