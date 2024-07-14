#! /home/tom/anaconda3/bin/python
"""
bingrad.py = written by Tom Seccull, 2024-07-05 - v0.0.7
	
	Last updated: 2024-07-14
	
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

from matplotlib.gridspec import GridSpec
from scipy.stats import linregress


def binning(spectrum_data, error_data):
	"""
	This function takes the spectroscopic data and uncertainties contained
	within a section of a larger spectrum that is to be binned. Each of the
	supplied datapoints is resampled from a gaussian distribution defined with
	the mean as the point's value and sigma as the point's uncertainty. The
	resulting distribution of resampled points from each of the datapoints in
	the bin are combined into a single common distribution. The median of this
	new distribution is returned as the value of the binned point while the 
	standard error of the distributions mean is returned as the binned point's
	uncertainty.
	
	Args:
	 -- spectrum_data (numpy.array)
			Array containing a subsection of the full spectrum that is to be
			binned into a single point.
	 -- error_data (numpy.array)
			Array containing the uncertainties associated with the spectrum
			values supplied in spectrum_data
	
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


def grad_measurement(point_dist, grad_indices, grad_wavelengths, size, p100nm):
	"""
	This script takes in a large 2D array of data. Each row in this ray is 
	spectroscopic data produced by resampling a common origin spectrum within
	its uncertainties. To measure the gradient of the origin spectrum, a linear
	regression is performed on each resampled spectrum to produce a 
	distribution of gradients from which a median and standard error of the
	mean can be calculated as the respective value of the gradient of the
	origin spectrum and the uncertainty of that gradient.
	
	Args:
	 -- point_dist (numpy.array)
		2D array of spectroscopic data where each row is version of a common 
		original spectrum that has been sampled from within the original
		spectrums uncertainties.
	 -- grad_indices (numpy.array)
		Array of indices marking the region of the spectroscopic data to be 
		measured. This array will skip gaps in the data caused by CCD chip gaps
		or those marked by the user.
	 -- grad_wavelengths (numpy.array)
		The section of the full spectrum wavelength axis that is marked by
		grad_indices.
	 -- size (int)
		The number of times the spectrum has been resampled to produce 
		point_dist.
	 -- p100nm (float)
		Conversion factor that accounts for wavelength units of the data when
		converting the raw slope measurement to units of %/100 nm.
	
	Returns:
	 -- median_grad (float)
		Converted spectral gradient in units of %/100 nm
	 -- grad_error (float)
		Uncertainty of spectral gradient in units of %/100 nm
	 -- median_slope (float)
		Raw median of the gradients returned by linregress. This is used to
		plot the linear fit if args.plot is called.
	 -- median_intercept (float)
		Raw median of the intercepts returned by linregress. This is used to
		plot the linear fit if args.plot is called.
	"""
	
	gradients = []
	fit_slopes = []
	fit_intercepts = []
	for i in range(size):
		grad_spectrum = (point_dist[i][grad_indices])
		linear_fit = linregress(grad_wavelengths, grad_spectrum)
		x1 = gradient_wavelengths[0]
		x2 = gradient_wavelengths[-1]
		y1 = (x1*linear_fit.slope) + linear_fit.intercept
		y2 = (x2*linear_fit.slope) + linear_fit.intercept
		gradients.append(100*((y2-y1)/((x2-x1)/p100nm))*(2/(y1+y2)))
		fit_slopes.append(linear_fit.slope)
		fit_intercepts.append(linear_fit.intercept)
	
	gradients = np.array(gradients)
	fit_slopes = np.array(fit_slopes)
	fit_intercepts = np.array(fit_intercepts)
	
	median_grad = np.median(gradients)
	grad_error = np.std(gradients)/((size-1)**0.5)
	median_slope = np.median(fit_slopes)
	median_intercept = np.median(fit_intercepts)

	return median_grad, grad_error, median_slope, median_intercept


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

if args.gradient_wavelength_ranges:
	sample_size = 5000
	
	optimal_point_distributions = []
	aperture_point_distributions = []
	for i in range(len(binned_wavelength_axis)):
		optimal_point_distributions.append(
			np.random.normal(
				loc=binned_optimal_spectrum[i],
				scale=binned_optimal_errors[i]*((args.factor-1)**0.5),
				size=(sample_size)
			)
		)
		aperture_point_distributions.append(
			np.random.normal(
				loc=binned_aperture_spectrum[i],
				scale=binned_aperture_errors[i]*((args.factor-1)**0.5),
				size=(sample_size)
			)
		)
	optimal_point_distributions = np.array(optimal_point_distributions).T
	aperture_point_distributions = np.array(aperture_point_distributions).T
	
	wavelength_strings = args.gradient_wavelength_ranges.split(",")
	wavelength_floats = [float(x) for x in wavelength_strings]
	wavelength_floats = np.array(wavelength_floats)
	
	hi_wavelength = wavelength_floats[-1]
	lo_wavelength = wavelength_floats[0]
	
	gradient_indices = []
	for i in range(int(len(wavelength_floats)*0.5)):
		grad_region_indices = np.where(
			np.logical_and(
				binned_wavelength_axis > wavelength_floats[2*i],
				binned_wavelength_axis < wavelength_floats[(2*i)+1]
			)
		)
		valid_grad_region_indices = (
			[x for x in grad_region_indices[0] if binned_qual[x]==0]
		)
		gradient_indices.append(valid_grad_region_indices)
	
	gradient_indices = np.concatenate(gradient_indices)
	gradient_wavelengths = binned_wavelength_axis[gradient_indices]
	
	per_hundred_nm_conversions = {
		"angstroms" : 1000.
	}
	
	per_100_nm_factor = per_hundred_nm_conversions[primary_head["WAVU"]]
	
	opt_gradient, opt_gradient_error, opt_slope, opt_inter = grad_measurement(
		optimal_point_distributions,
		gradient_indices,
		gradient_wavelengths,
		sample_size,
		per_100_nm_factor
	)
	ape_gradient, ape_gradient_error, ape_slope, ape_inter = grad_measurement(
		aperture_point_distributions,
		gradient_indices,
		gradient_wavelengths,
		sample_size,
		per_100_nm_factor
	)
	stack_headers = [x for x in headers if "STACK" in x["EXTNAME"]]
	airmass_difference = np.abs(
		stack_headers[0]["MDAIRMSS"] - stack_headers[1]["MDAIRMSS"]
	)
	opt_grad_uncertainty = np.sqrt(
		opt_gradient_error**2
		+ np.abs(opt_gradient - ape_gradient)**2
		+ airmass_difference**2
	)
	ape_grad_uncertainty = np.sqrt(
		ape_gradient_error**2
		+ np.abs(opt_gradient - ape_gradient)**2
		+ airmass_difference**2
	)
	
	opt_grad = "%.2f" % round(opt_gradient,2)
	u_opt_grad = "%.2f" % round(opt_grad_uncertainty,2)
	ape_grad = "%.2f" % round(ape_gradient,2)
	u_ape_grad = "%.2f" % round(ape_grad_uncertainty,2)
	print(
		"Optimal Spectrum Gradient: " 
		+ opt_grad 
		+ " +/- " 
		+ u_opt_grad
		+ " %/100 nm"
	)
	print(
		"Aperture Spectrum Gradient: " 
		+ ape_grad 
		+ " +/- " 
		+ u_ape_grad
		+ " %/100 nm"
	)
	
if args.plot:
	fig = plt.figure(figsize=(6,10))
	gs = GridSpec(2,1, figure=fig)
	ax1 = plt.subplot(gs[0,0])
	ax2 = plt.subplot(gs[1,0])
	
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
		if i==0:
			ax.set_xticklabels([])
		if i==1:
			ax.set_xlabel("Wavelength, " + headers[1]["WAVU"])

	plt.figtext(0.03, 0.435, "Relative Reflectance", rotation=90)
	
	if args.gradient_wavelength_ranges:
		opt_line = (binned_wavelength_axis * opt_slope) + opt_inter
		ape_line = (binned_wavelength_axis * ape_slope) + ape_inter
		opt_input_points = (wavelength_floats * opt_slope) + opt_inter
		ape_input_points = (wavelength_floats * ape_slope) + ape_inter
		lines = [opt_line, ape_line]
		cols = ["goldenrod", "magenta"]
		grads = [opt_grad, ape_grad]
		u_grads = [u_opt_grad, u_ape_grad]
		input_points = [opt_input_points, ape_input_points]
		for i , ax in enumerate(axes):
			ax.plot(
				binned_wavelength_axis,
				lines[i],
				linestyle="dashed",
				color=cols[i],
				linewidth=2.5,
				zorder=9,
				label=(
					"Linear Fit: Gradient = "
					+ grads[i]
					+ " +/- "
					+ u_grads[i]
					+ " %/100 nm"
				)
			)
			ax.plot(
				wavelength_floats,
				input_points[i],
				linestyle="",
				marker='.',
				markersize=12,
				mec=cols[i],
				mfc="white",
				zorder=9,
				label="Gradient Input Wavelengths"
			)
			ax.legend(loc=4)
	
	else:
		ax1.legend(loc=4)
		ax2.legend(loc=4)

	print(args.data_file)
	plt.subplots_adjust(hspace=0)
	plt.show()

if args.save:
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
	
	binned_hdu = fits.PrimaryHDU(combi_binned_frame)
	binned_header = binned_hdu.header
	binned_header["DATE"] = (datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), "UT file creation date.")
	binned_header["FITSDOI"] = ("10.1051/0004-6361:20010923", "FITS format definition paper DOI")
	binned_header["ORIGIN"] = ("bingrad.py v0.0.7", "Script that created this file.")
	binned_header["DIVDOI"] = ("UNKNOWN", "Script repository DOI")
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
	binned_header["EXTNAME"] = (
		primary_head["EXTNAME"].split("STACK")[0] + "BINNED"
	)
	
	listed_hdus = [binned_hdu]
	
	for i in range(len(headers)):
		listed_hdus.append(fits.ImageDU(frames[i], header=headers[i]))
		
	hdu_list = fits.HDUList(listed_hdus)
	hdu_list.writeto("b" + args.data_file)
	hdu_list.close()
	
	#### Test save function and add wavelength output to printout from gradient measurement.
