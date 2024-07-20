#! /home/tom/anaconda3/envs/work/bin/python
"""
bingrad.py = written by Tom Seccull, 2024-07-05 - v1.0.0
	
	Last updated: 2024-07-16
	
	This script has two functions. Primarily it is used to bin 
	spectroscopic data to boost its signal-to-noise ratio at the expense
	of spectral resolution. A binned spectrum can be plotted and saved
	to a new FITS file. The secondary function is to allow the continuum
	gradient of the spectrum to be measured across a user-defined
	wavelength range via linear regression. The resulting linear fit can
	be plotted and its parameters will be printed in the terminal.
"""

import argparse
import astropy.io.fits as fits
import bingrad.gradient as gradient
import datetime
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.gridspec import GridSpec


def binning(spectrum_data, error_data):
	"""
	This function takes the spectroscopic data and uncertainties
	contained within a section of a larger spectrum that is to be
	binned. Each of the supplied datapoints is resampled from a gaussian
	distribution defined with the mean as the point's value and sigma as
	the point's uncertainty. The resulting distribution of resampled
	points from each of the datapoints in the bin are combined into a
	single common distribution. The median of this new distribution is
	returned as the value of the binned point while the standard error
	of the distribution's mean is returned as the binned point's
	uncertainty.
	
	Args:
	 -- spectrum_data (numpy.array)
			Array containing a subsection of the full spectrum that is
			to be binned into a single point.
	 -- error_data (numpy.array)
			Array containing the uncertainties associated with the
			spectrum values supplied in spectrum_data
	
	Returns:
	 -- binned_point_data (float)
			The data point derived from binning spectrum_data.
	 -- binned_point_error (float)
			An estimate of the uncertainty of binned_point_data. 
	"""
	
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
				size=5000
			)
		)
	
	# Median is within the standard error of the mean from the mean in all
	# cases. Taking the median as the binned value and the standard error
	# of the mean as the uncertainty.
	uncertain_samples = np.array(uncertain_samples).flatten()
	binned_point_data = np.median(uncertain_samples)
	if len(spectrum_data) > 1:
		binned_point_error = (
			np.std(uncertain_samples)
			/ ((len(spectrum_data)-1)**0.5)
		)
	else: binned_point_error = np.std(uncertain_samples)
	
	return binned_point_data, binned_point_error


###############################################################################
#### SCRIPT STARTS HERE  # # # # # # # # # # # # # # # # # # # # # # # # # ####
###############################################################################

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

combi_binned_frame = np.array(
		[
			binned_wavelength_axis,
			binned_optimal_spectrum,
			binned_optimal_errors,
			binned_aperture_spectrum,
			binned_aperture_errors,
			binned_qual
		]
	)

if args.gradient_wavelength_ranges:
	opt_params, ape_params, wavelength_floats = gradient.grad_measurement(
		combi_binned_frame,
		args.factor,
		args.gradient_wavelength_ranges,
		headers
	)

if args.plot:
	fig = plt.figure(figsize=(6,10))
	gs = GridSpec(2,1, figure=fig)
	ax1 = plt.subplot(gs[0,0])
	ax2 = plt.subplot(gs[1,0], sharex=ax1, sharey=ax1)
	
	axes = [ax1, ax2]
	specs = [binned_optimal_spectrum, binned_aperture_spectrum]
	errs = [binned_optimal_errors, binned_aperture_errors]
	cols = ["k","dimgray"]
	labels = [
		"Optimal Spectrum, Binned x" + str(args.factor),
		"Aperture Spectrum, Binned x" + str(args.factor),
	]
	
	for i, ax in enumerate(axes):
		ax.errorbar(
			binned_wavelength_axis,
			specs[i],
			yerr=errs[i],
			marker='.',
			linestyle='',
			color=cols[i],
			label=labels[i]
		)
	
		ylim_bottom, ylim_top = ax.get_ylim()
		qual_shade_bottom = np.ones(len(binned_wavelength_axis)) * ylim_bottom
		qual_shade_top = np.ones(len(binned_wavelength_axis)) * ylim_top
		ax.fill_between(
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
		if i==1:
			ax.set_xlabel("Wavelength, " + headers[1]["WAVU"])

	plt.figtext(0.03, 0.435, "Relative Reflectance", rotation=90)
	plt.setp(ax1.get_xticklabels(), visible=False)
	
	if args.gradient_wavelength_ranges:
		gradient.plot_gradient(
			opt_params,
			ape_params,
			combi_binned_frame,
			wavelength_floats,
			axes
		)
	
	else:
		ax1.legend(loc=4)
		ax2.legend(loc=4)

	plt.subplots_adjust(hspace=0)
	plt.show()

if args.save:
	binned_hdu = fits.PrimaryHDU(combi_binned_frame)
	binned_header = binned_hdu.header
	binned_header["DATE"] = (datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), "UT file creation date.")
	binned_header["FITSDOI"] = ("10.1051/0004-6361:20010923", "FITS format definition paper DOI")
	binned_header["ORIGIN"] = ("bingrad.py v0.0.7", "Script that created this file.")
	binned_header["DIVDOI"] = ("10.5281/zenodo.12786056", "Script repository DOI")
	binned_header["INPUT1"] = (args.data_file, "Input spectrum file")
	binned_header["OBJECT1"] = (primary_head["OBJECT1"], "Name of first object in ratio spectrum")
	binned_header["OBJECT2"] = (primary_head["OBJECT2"], "Name of second object in ratio spectrum")
	binned_header["EXTENS"] = (primary_head["EXTNAME"], "Extension of unbinned data frame")
	binned_header["BINFACTR"] = (args.factor, "Binning factor used to bin spectrum")
	binned_header["HDUROW0"] = "Wavelength axis, " + primary_head["WAVU"]
	binned_header["HDUROW1"] = "Optimal binned spectrum"
	binned_header["HDUROW2"] = "Optimal binned spectrum uncertainty"
	binned_header["HDUROW3"] = "Aperture binned spectrum"
	binned_header["HDUROW4"] = "Aperture binned spectrum uncertainty"
	binned_header["HDUROW5"] = "Quality flags: 0=GOOD, 1=BAD"
	binned_header["EXTNAME"] = primary_head["EXTNAME"].split("STACK")[0] + "BINNED"
	
	listed_hdus = [binned_hdu]
	
	for i in range(len(headers)):
		listed_hdus.append(fits.ImageHDU(frames[i], header=headers[i]))
		
	hdu_list = fits.HDUList(listed_hdus)
	hdu_list.writeto("b" + args.data_file)
	hdu_list.close()
