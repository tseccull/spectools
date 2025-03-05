#!/usr/bin/env python3

"""
gradient.py - written by Tom Seccull, 2024-07-16

	Called by: bingrad.py
	Last updated: 2025-03-05
	
	This file contains all the functions needed by bingrad.py to perform
	a linear regression of the supplied spectral data and measure the
	value and uncertainty of its gradient within a user-defined
	wavelength region. 
"""


import numpy as np

from scipy.stats import linregress


def grad_measurement(binned_frame, factor, grad_wavelengths, headers):
	"""
	Prepares the binned spectrum for gradient measurement and supplies
	the prepared spectrum to the line_fitting function. Takes the
	output of the line fitting function to determine the final value
	and uncertainty of the spectral gradient in %/100 nm and prints
	the results to the terminal.
	
	Args:
	 -- binned_frame (numpy.array)
			2D array containing all the binned spectroscopic data. Array
			rows are as follows:
			binned_frame[0]: binned wavelength axis
			binned_frame[1]: binned optimally extracted spectrum
			binned_frame[2]: uncertainties of binned_frame[1]
			binned_frame[3]: binned aperture extracted spectrum
			binned_frame[4]: uncertainties of binned_frame[3]
			binned_frame[5]: binned_quality frame - 0 = GOOD, 1 = BAD
	 -- factor (int)
			User-defined binning factor (number of points per bin) used
			to bin the original spectrum.
	 -- grad_wavelengths (str)
			User-defined wavelength regions where the spectrum's
			gradient will be measured.
	 -- headers (list)
			List of the headers of all the extensions found in the 
			file used as input to bingrad.py
			
	Returns:
	 -- optimal_line_parameters (list)
			list of parameters determined from gradient measurement of
			the optimally extracted spectrum:
			opt_slope: median slope returned by linregress
			opt_inter: median intercept returned by linregress
			opt_grad: median spectral gradient in %/100 nm
			u_opt_grad: uncertainty of opt_grad in %/100 nm
	 -- aperture_line_parameters (list)
			list of parameters determined from gradient measurement of
			the aperture extracted spectrum:
			ape_slope: median slope returned by linregress
			ape_inter: median intercept returned by linregress
			ape_grad: median spectral gradient in %/100 nm
			u_ape_grad: uncertainty of ape_grad in %/100 nm
	 -- wavelength_floats (numpy.array)
			1D array containing the user-defined wavelengths used to
			define the region of the spectrum where the spectral
			gradient is measured.
	"""
	sample_size = 5000
	
	optimal_point_distributions = []
	aperture_point_distributions = []
	for i in range(len(binned_frame[0])):
		optimal_point_distributions.append(
			np.random.normal(
				loc=binned_frame[1][i],
				scale=binned_frame[2][i]*((factor-1)**0.5),
				size=(sample_size)
			)
		)
		aperture_point_distributions.append(
			np.random.normal(
				loc=binned_frame[3][i],
				scale=binned_frame[4][i]*((factor-1)**0.5),
				size=(sample_size)
			)
		)
	optimal_point_distributions = np.array(optimal_point_distributions).T
	aperture_point_distributions = np.array(aperture_point_distributions).T
	
	wavelength_strings = grad_wavelengths.split(",")
	wavelength_floats = [float(x) for x in wavelength_strings]
	wavelength_floats = np.array(wavelength_floats)
	
	hi_wavelength = wavelength_floats[-1]
	lo_wavelength = wavelength_floats[0]
	
	gradient_indices = []
	for i in range(int(len(wavelength_floats)*0.5)):
		grad_region_indices = np.where(
			np.logical_and(
				binned_frame[0] > wavelength_floats[2*i],
				binned_frame[0] < wavelength_floats[(2*i)+1]
			)
		)
		valid_grad_region_indices = (
			[x for x in grad_region_indices[0] if binned_frame[5][x]==0]
		)
		gradient_indices.append(valid_grad_region_indices)
	
	gradient_indices = np.concatenate(gradient_indices)
	gradient_wavelengths = binned_frame[0][gradient_indices]
	
	per_hundred_nm_conversions = {
		"angstroms" : 1000.
	}
	
	per_100_nm_factor = per_hundred_nm_conversions[headers[0]["WAVU"]]
	
	opt_gradient, opt_gradient_error, opt_slope, opt_inter = line_fitting(
		optimal_point_distributions,
		gradient_indices,
		gradient_wavelengths,
		sample_size,
		per_100_nm_factor
	)
	ape_gradient, ape_gradient_error, ape_slope, ape_inter = line_fitting(
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
	print(
		"Reference wavelength: "
		+ str((wavelength_floats[0]
		+ wavelength_floats[-1])*0.5)
		+ " "
		+ headers[0]["WAVU"]
	)
	print("Measured Wavelength Range(s):")
	for i in range(int(len(wavelength_floats)*0.5)):
		print(
			"    "
			+ str(wavelength_floats[i*2])
			+ "--"
			+ str(wavelength_floats[(i*2)+1])
			+ " "
			+ headers[0]["WAVU"]
		)
	
	optimal_line_parameters = [opt_slope, opt_inter, opt_grad, u_opt_grad]
	aperture_line_parameters = [ape_slope, ape_inter, ape_grad, u_ape_grad]
	
	return optimal_line_parameters, aperture_line_parameters, wavelength_floats


def line_fitting(point_dist, grad_indices, grad_wavelengths, size, p100nm):
	"""
	This script takes in a large 2D array of data. Each row in this ray
	is spectroscopic data produced by resampling a common origin
	spectrum within its uncertainties. To measure the gradient of the
	origin spectrum, a linear regression is performed on each resampled
	spectrum to produce a distribution of gradients from which a median
	and standard error of the mean can be calculated as the respective
	value of the gradient of the origin spectrum and the uncertainty of
	that gradient.
	
	Args:
	 -- point_dist (numpy.array)
		2D array of spectroscopic data where each row is version of a
		common original spectrum that has been sampled from within the
		original spectrums uncertainties.
	 -- grad_indices (numpy.array)
		Array of indices marking the region of the spectroscopic data to
		be measured. This array will skip gaps in the data caused by CCD
		chip gaps or those marked by the user.
	 -- grad_wavelengths (numpy.array)
		The section of the full spectrum wavelength axis that is marked
		by grad_indices.
	 -- size (int)
		The number of times the spectrum has been resampled to produce 
		point_dist.
	 -- p100nm (float)
		Conversion factor that accounts for wavelength units of the data
		when converting the raw slope measurement to units of %/100 nm.
	
	Returns:
	 -- median_grad (float)
		Converted spectral gradient in units of %/100 nm
	 -- grad_error (float)
		Uncertainty of spectral gradient in units of %/100 nm
	 -- median_slope (float)
		Raw median of the gradients returned by linregress. This is used
		to plot the linear fit if args.plot is called.
	 -- median_intercept (float)
		Raw median of the intercepts returned by linregress. This is used
		to plot the linear fit if args.plot is called.
	"""
	
	gradients = []
	fit_slopes = []
	fit_intercepts = []
	for i in range(size):
		grad_spectrum = (point_dist[i][grad_indices])
		linear_fit = linregress(grad_wavelengths, grad_spectrum)
		x1 = grad_wavelengths[0]
		x2 = grad_wavelengths[-1]
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


def plot_gradient(
	opt_params, ape_params, binned_frame, wavelength_floats, axes
):
	"""
	Takes the output of the gradient measurement process and plots it
	over the binned spectrum when plotting is called for.
	
	Args:
	 -- optimal_line_parameters (list)
			list of parameters determined from gradient measurement of
			the optimally extracted spectrum:
			opt_slope: median slope returned by linregress
			opt_inter: median intercept returned by linregress
			opt_grad: median spectral gradient in %/100 nm
			u_opt_grad: uncertainty of opt_grad in %/100 nm
	 -- aperture_line_parameters (list)
			list of parameters determined from gradient measurement of
			the aperture extracted spectrum:
			ape_slope: median slope returned by linregress
			ape_inter: median intercept returned by linregress
			ape_grad: median spectral gradient in %/100 nm
			u_ape_grad: uncertainty of ape_grad in %/100 nm
	 -- binned_frame (numpy.array)
			2D array containing all the binned spectroscopic data. Array
			rows are as follows:
			binned_frame[0]: binned wavelength axis
			binned_frame[1]: binned optimally extracted spectrum
			binned_frame[2]: uncertainties of binned_frame[1]
			binned_frame[3]: binned aperture extracted spectrum
			binned_frame[4]: uncertainties of binned_frame[3]
			binned_frame[5]: binned_quality frame - 0 = GOOD, 1 = BAD
	 -- wavelength_floats (numpy.array)
			1D array containing the user-defined wavelengths used to
			define the region of the spectrum where the spectral
			gradient is measured.
	 -- axes (list)
			A list containing the two axes present in the figure plotted
			by bingrad.py
	
	Returns:
	None
	"""
	
	opt_line = (binned_frame[0] * opt_params[0]) + opt_params[1]
	ape_line = (binned_frame[0] * ape_params[0]) + ape_params[1]
	opt_input_points = (wavelength_floats * opt_params[0]) + opt_params[1]
	ape_input_points = (wavelength_floats * ape_params[0]) + ape_params[1]
	lines = [opt_line, ape_line]
	cols = ["goldenrod", "magenta"]
	grads = [opt_params[2], ape_params[2]]
	u_grads = [opt_params[3], ape_params[3]]
	input_points = [opt_input_points, ape_input_points]
	for i , ax in enumerate(axes):
		ax.plot(
			binned_frame[0],
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

