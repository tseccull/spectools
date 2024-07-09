#! /home/tom/anaconda3/bin/python
"""
bingrad.py = written by Tom Seccull, 2024-07-05 - v0.0.3
	
	Last updated: 2024-07-09
	
	This script has two functions. Primarily it is used to bin spectroscopic 
	data to boost its signal-to-noise ratio at the expense of spectral 
	resolution. A binned spectrum can be plotted and saved to a new FITS file. 
	The secondary function is to allow the continuum gradient of the spectrum 
	to be measured across a user-defined wavelength range via linear 
	regression. The resulting linear fit can be plotted and its parameters 
	will be printed in the terminal.
"""

import argparse
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import linregress


def binning(spectrum_data, error_data):
	
	uncertain_samples = []
	
	# If all data points in bin fall in a chip gap mark the binned point and 
	# its uncertainty as zero
	if all(x==0 for x in spectrum_data):
		return 0., 0.
	
	error_data = error_data[spectrum_data!=0]
	spectrum_data = spectrum_data[spectrum_data!=0]
	for j in range(len(spectrum_data)):
		uncertain_samples.append(
			np.random.normal(
				loc=spectrum_data[j],
				scale=error_data[j],
				size=100
			)
		)
	
	# Median is within the standard error of the mean from the mean in all
	# cases. Taking the median as the binned value and the standard error
	# of the mean as the uncertainty.
	uncertain_samples = np.array(uncertain_samples).flatten()
	binned_point_data = np.median(uncertain_samples)
	binned_point_error = (
		  np.std(uncertain_samples)
		/ ((len(spectrum_data)-1)**0.5)
	)
	#plt.hist(uncertain_samples, bins=100)
	#plt.axvline(binned_point_data, color='k')
	#plt.axvline(binned_point_data-binned_point_error, color='r')
	#plt.axvline(binned_point_data+binned_point_error, color='r')
	#plt.show()
	
	return binned_point_data, binned_point_error


parser = argparse.ArgumentParser(
	description="This script has two functions. Primarily it is used to bin\
	spectroscopic data to boost its signal-to-noise ratio at the expense of\
	spectral resolution. A binned spectrum can be plotted and saved to a new\
	FITS file. The secondary function is to allow the continuum\
	gradient of the spectrum to be measured across a user-defined wavelength\
	range via linear regression. The resulting linear fit can be plotted and\
	its parameters will be printed in the terminal."
)
parser.add_argument(
	"data_file", type=str,
	help="Input spectrum FITS file name."
)
parser.add_argument(
	"factor", type=int,
	help="Binning factor, a.k.a. bin width or number of points per bin."
)
parser.add_argument(
	"-r", "--rejection", type=float,
	help="This optional argument will be multiplied by the standard deviation\
	of the combined distribution of points present within each bin after each\
	original point is resampled within its Gaussian uncertainties. The\
	resulting number will be set as both the positive and negative distance\
	from the median of the distribution beyond which points will be rejected\
	from the final estimation of the median value of the bin. If this\
	argument is not set, the default behaviour of bingrad is to perform no\
	rejection."
)
parser.add_argument(
	"-g" ,"--gradient_wavelength_ranges",
	help="This is a comma-separated list of wavelengths that bound regions\
	where the spectral gradient will be measured from the binned spectrum.\
	Only one gradient measurement is done each time bingrad is run, but\
	multiple wavelength regions can be set to allow parts of the spectrum to\
	be skipped in the measurement if needed. The units of -gwr should match\
	those of the wavelength axis of the input spectrum. If np gradient\
	wavelength ranges are set, bingrad defaults to skipping the gradient\
	measurement."
)
parser.add_argument(
	"-p", "--plot", action="store_true",
	help="Plot the binned spectrum on screen. If a gradient has also been\
	measured, the linear fit to the chosen spectral regions will also be\
	presented."
)
parser.add_argument(
	"-s", "--save", action="store_true",
	help="Save the binned spectrum to a new FITS file."
)
args = parser.parse_args()


with fits.open(args.data_file) as data_file:
	headers = [x.header for x in data_file]
	frames = [x.data for x in data_file]
	
primary_head = headers[0]
primary_frame = frames[0]
qual = primary_frame[5]

# Cut off the qual==0 ends of the spectrum.
i = 0
while primary_frame[5][i] == 1: i += 1
j = np.shape(primary_frame)[1] - 1
while primary_frame[5][j] == 1: j -= 1
primary_frame = primary_frame[:, i:j+1]

# Make sure the number of unbinned points is even by removing the first or the
# last point in the spectrum based on whichever has a larger uncertainty. 
# Having an even number of data points in the spectrum simplifies the process
# of binning outward from the centre of the spectrum.
if np.shape(primary_frame)[1]%2 != 0:
	if primary_frame[2][0] > primary_frame[2][-1]:
		primary_frame = primary_frame[:, 1:]
	else:
		primary_frame = primary_frame[:, :-1]

wavelength_axis = primary_frame[0]
optimal_spectrum = primary_frame[1]
optimal_errors = primary_frame[2]
aperture_spectrum = primary_frame[3]
aperture_errors = primary_frame[4]
qual = primary_frame[5]

upper_bins = np.arange(
	int(len(wavelength_axis)*0.5) + args.factor,
	len(wavelength_axis), 
	args.factor
)
upper_bins = np.append(upper_bins, [len(wavelength_axis)])

lower_bins = np.arange(
	int(len(wavelength_axis)*0.5), 0, -args.factor
)
lower_bins = np.append(lower_bins, [0])

bins = np.concatenate([np.flip(lower_bins), upper_bins])

binned_wavelength_axis = []
binned_optimal_spectrum = []
binned_optimal_errors = []
binned_aperture_spectrum = []
binned_aperture_errors = []
binned_qual = []

i = 0
wavelength_step = wavelength_axis[1]-wavelength_axis[0]
end_points_needed = args.factor - bins[1]
while bins[i] < bins[-1]:
	if bins[i+1]-bins[i] < args.factor and bins[i] == 0:
		wavelength_extension = np.arange(
			wavelength_axis[0]-(end_points_needed * wavelength_step),
			wavelength_axis[0], 
			wavelength_step
		)
		bin_wavelengths = np.concatenate(
			[wavelength_extension, wavelength_axis[:bins[1]]]
		)
	elif bins[i+1]-bins[i] < args.factor and bins[i+1] == len(wavelength_axis):
		wavelength_extension = np.arange(
			wavelength_axis[-1] + wavelength_step,
			wavelength_axis[-1] + ((end_points_needed + 1) * wavelength_step),
			wavelength_step
		)
		bin_wavelengths = np.concatenate(
			[wavelength_axis[bins[i]:], wavelength_extension]
		)
	else:
		bin_wavelengths = wavelength_axis[bins[i]:bins[i+1]]
	binned_wavelength_axis.append(np.median(bin_wavelengths))
	
	bin_optimal_data = optimal_spectrum[bins[i]:bins[i+1]]
	bin_optimal_errors = optimal_errors[bins[i]:bins[i+1]]
	bin_aperture_data = aperture_spectrum[bins[i]:bins[i+1]]
	bin_aperture_errors = aperture_errors[bins[i]:bins[i+1]]
	
	binned_optimal_datapoint, binned_optimal_error = binning(
		bin_optimal_data, bin_optimal_errors
	)
	binned_aperture_datapoint, binned_aperture_error = binning(
		bin_aperture_data, bin_aperture_errors
	)
	
	binned_optimal_spectrum.append(binned_optimal_datapoint)
	binned_optimal_errors.append(binned_optimal_error)
	binned_aperture_spectrum.append(binned_aperture_datapoint)
	binned_aperture_errors.append(binned_aperture_error)
	
	binned_points = [
		binned_optimal_datapoint, 
		binned_optimal_error, 
		binned_aperture_datapoint,
		binned_aperture_error
	]
	if all(x==0 for x in binned_points):
		binned_qual.append(1)
	else:
		binned_qual.append(0)
	
	i += 1

binned_wavelength_axis = np.array(binned_wavelength_axis)
binned_optimal_spectrum = np.array(binned_optimal_spectrum)
binned_optimal_errors = np.array(binned_optimal_errors)
binned_aperture_spectrum = np.array(binned_aperture_spectrum)
binned_aperture_errors = np.array(binned_aperture_errors)
binned_qual = np.array(binned_qual)

#if args.gradient_wavelength_ranges:

if args.plot:
	fig = plt.figure(figsize=(10,6))
	ax = plt.subplot()
	ax.errorbar(
		binned_wavelength_axis,
		binned_optimal_spectrum,
		yerr=binned_optimal_errors,
		marker='.',
		linestyle='',
		color="k",
		label="Optimal Spectrum, Binned x" + str(args.factor)
	)
	ax.errorbar(
		binned_wavelength_axis,
		binned_aperture_spectrum + 1,
		yerr=binned_aperture_errors,
		marker='.',
		linestyle='',
		color="dimgray",
		label="Aperture Spectrum +1, Binned x" + str(args.factor)
	)
	ylim_bottom, ylim_top = plt.ylim()
	qual_shade_bottom = np.ones(len(binned_wavelength_axis)) * ylim_bottom
	qual_shade_top = np.ones(len(binned_wavelength_axis)) * ylim_top
	plt.fill_between(
		binned_wavelength_axis,
		qual_shade_bottom,
		qual_shade_top,
		where=binned_qual==1,
		facecolor="red",
		alpha=0.5,
		label="QUAL=1"
	)
	ax.set_ylim(ylim_bottom, ylim_top)
	ax.grid(linestyle="dashed")
	ax.set_xlabel("Wavelength, " + headers[1]["WAVU"])
	ax.set_ylabel("Relative Reflectance")
	ax.legend()
	plt.show()

#if args.save:
