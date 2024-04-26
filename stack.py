#! /home/tom/anaconda3/envs/work/bin/python
"""
stack.py - written by Tom Seccull, 2024-04-18 - v0.0.1

	Last updated: 2024-04-26
	
	This script takes multiple 1D spectra and combines them with a
	weighted bootstrapped median to produce a stacked 1D spectrum with
	reduced noise.
"""


import argparse
import astropy.io.fits as fits
import datetime
import glob
import matplotlib.pyplot as plt
import numpy as np


def stack_header_gmos(new_hdu, heads, files):
	new_head = new_hdu.header
	key_0 = files[0][:-5]
	
	# Observation metadata
	new_head["SEQEXVER"] = (heads[key_0]["SEQEXVER"], "Seqexec version")
	new_head["GEMPRGID"] = (heads[key_0]["GEMPRGID"], "Gemini program ID")
	new_head["OBSID"] = (heads[key_0]["OBSID"], "Gemini observation ID")
	new_head["OBSERVER"] = (heads[key_0]["OBSERVER"], "Observer")
	new_head["OPERATOR"] = (heads[key_0]["SSA"], "Telescope operator")
	new_head["OBSERVAT"] = ("Gemini Observatory")
	new_head["LT"] = (heads[key_0]["LT"], "Local time at sequence start")
	new_head["UT"] = (heads[key_0]["UT"], "UT at sequence start")
	new_head["ST"] = (heads[key_0]["ST"], "Sidereal time at sequence start")
	new_head["RAWIQ"] = (heads[key_0]["RAWIQ"], "Raw image quality")
	new_head["RAWCC"] = (heads[key_0]["RAWCC"], "Raw cloud cover")
	new_head["RAWWV"] = (heads[key_0]["RAWWV"], "Raw water vapor")
	new_head["RAWBG"] = (heads[key_0]["RAWBG"], "Raw sky brightness")
	new_head["RAWPIREQ"] = (heads[key_0]["RAWPIREQ"], "PI's condition requirements met?")
	
	# Weather metadata
	humidities = np.array([heads[x[:-5]]["HUMIDITY"] for x in files])
	min_humidity = np.min(humidities)
	max_humidity = np.max(humidities)
	med_humidity = np.median(humidities)
	new_head["LOHUMID"] = (min_humidity, "Minimum relative humidity, %")
	new_head["HIHUMID"] = (max_humidity, "Maximum relative humidity, %")
	new_head["MEDHUMID"] = (med_humidity, "Median relative humidity, %")
	
	temperatures = np.array([heads[x[:-5]]["TAMBIENT"] for x in files])
	min_temp = np.min(temperatures)
	max_temp = np.max(temperatures)
	med_temp = np.median(temperatures)
	new_head["LOTAMBI"] = (min_temp, "Minimum ambient temperature, C")
	new_head["HITAMBI"] = (max_temp, "Maximum ambient temperature, C")
	new_head["MEDTAMBI"] = (med_temp, "Median ambient temperature, C")
	
	pressures = np.array([heads[x[:-5]]["PRESSUR2"] for x in files])
	min_pressure = np.min(pressures)
	max_pressure = np.max(pressures)
	med_pressure = np.median(pressures)
	new_head["LOPRESS"] = (min_pressure, "Minimum atmospheric pressure, Pa")
	new_head["HIPRESS"] = (max_pressure, "Maximum atmospheric pressure, Pa")
	new_head["MEDPRESS"] = (med_pressure, "Median atmospheric pressure, Pa")
	
	dew_points = np.array([heads[x[:-5]]["DEWPOINT"] for x in files])
	min_dew_points = round(np.min(dew_points), 2)
	max_dew_points = round(np.max(dew_points), 2)
	med_dew_points = round(np.median(dew_points), 2)
	new_head["LODEWPNT"] = (min_dew_points, "Minimum dew point, C")
	new_head["HIDEWPNT"] = (max_dew_points, "Maximum dew point, C")
	new_head["MEDDWPNT"] = (med_dew_points, "Median dew point, C")
	
	wind_speeds_ms = np.array([heads[x[:-5]]["WINDSPEE"] for x in files])
	min_wind_speed_ms = round(np.min(wind_speeds_ms), 2)
	max_wind_speed_ms = round(np.max(wind_speeds_ms), 2)
	med_wind_speed_ms = round(np.median(wind_speeds_ms), 2)
	new_head["LWNDVLMS"] = (min_wind_speed_ms, "Minimum wind speed, m/s")
	new_head["HWNDVLMS"] = (max_wind_speed_ms, "Maximum wind speed, m/s")
	new_head["MWNDVLMS"] = (med_wind_speed_ms, "Median wind speed, m/s")
	
	wind_speeds_mph = np.array([heads[x[:-5]]["WINDSPE2"] for x in files])
	min_wind_speed_mph = round(np.min(wind_speeds_ms), 2)
	max_wind_speed_mph = round(np.max(wind_speeds_ms), 2)
	med_wind_speed_mph = round(np.median(wind_speeds_ms), 2)
	new_head["LWNDVLMH"] = (min_wind_speed_mph, "Minimum wind speed, mph")
	new_head["HWNDVLMH"] = (max_wind_speed_mph, "Maximum wind speed, mph")
	new_head["MWNDVLMH"] = (med_wind_speed_mph, "Median wind speed, mph")
	
	wind_directions = np.array([heads[x[:-5]]["WINDDIRE"] for x in files])
	radian_wind_directions = wind_directions*(np.pi/180)
	sum_sin_wind_directions = np.sum(np.sin(radian_wind_directions))
	sum_cos_wind_directions = np.sum(np.cos(radian_wind_directions))
	radian_mean_wind_direction = np.arctan2(
		sum_sin_wind_directions, sum_cos_wind_directions
	)
	mean_wind_direction = round(radian_mean_wind_direction/(np.pi/180), 2)
	new_head["MEDWDDIR"] = (mean_wind_direction, "Circular mean wind direction, deg")
	
	# Target and guiding metadata
	new_head["OBJECT"] = (heads[key_0]["OBJECT"], "Target name")
	new_head["PARALLAX"] = (heads[key_0]["PARALLAX"], "Parallax of target")
	new_head["RADVEL"] = (heads[key_0]["RADVEL"], "Heliocentric radial velocity")
	new_head["EQUINOX"] = (heads[key_0]["EQUINOX"], "Equinox of coordinate system")
	new_head["TRKEQUIN"] = (heads[key_0]["TRKEQUIN"], "Tracking equinox")
	new_head["FRAME"] = (heads[key_0]["FRAME"], "Target coordinate system")
	
	if all(heads[key_0][x]==0 for x in ["RATRACK", "DECTRACK"]):
		new_head["SIDEREAL"] = (True, "Is target sidereal? (True/False)")
		new_head["EPOCH"] = (heads[key_0]["EPOCH"], "Epoch for target coordinates")
		
		round_ra = round(heads[key_0]["RA"], 5)
		new_head["RA"] = (round_ra, "Target RA, deg")
		
		round_dec = round(heads[key_0]["DEC"], 5)
		new_head["DEC"] = (round_dec, "Target Dec, deg")
		
		ra_rate_deg = heads[key_0]["DECTRACK"] * 15
		new_head["RATRACK"] = (ra_rate_deg, "RA differential tracking rate, deg/hr")
		
		dec_rate_deg = heads[key_0]["DECTRACK"]
		new_head["DECTRACK"] = (dec_rate_deg, "Dec differential tracking rate, deg/hr")
	
	else:
		new_head["SIDEREAL"] = (False, "Is target sidereal? (True/False)")
		
		tracking_epochs = np.array([heads[x[:-5]]["TRKEPOCH"] for x in files])
		tracking_intervals = tracking_epochs - np.roll(tracking_epochs, 1)
		tracking_intervals = np.roll(tracking_intervals, -1)
		tracking_intervals[-1] = np.median(tracking_intervals)
		tracking_intervals *= 24
		ras = [heads[key_0]["RA"]]
		decs = [heads[key_0]["DEC"]]
		initial_epoch = heads[key_0]["EPOCH"]
		low_middle = int(np.floor(len(files)/2))
		high_middle = int(np.ceil(len(files)/2))
		middle = int(len(files)/2)
		for i in range(1, high_middle+1):
			ra_addition = (
				heads[files[i-1][:-5]]["RATRACK"] * tracking_intervals[i] * 15
			)
			ras.append(ras[-1] + ra_addition)
			
			dec_addition = (
				heads[files[i-1][:-5]]["DECTRACK"] * tracking_intervals[i]
			)
			decs.append(decs[-1] + dec_addition)
		
		if len(files)%2 == 0:
			half_duration = np.sum(tracking_intervals[:middle]/24)
			median_epoch = round(initial_epoch + half_duration, 7)
			median_ra = round(ras[middle], 5)
			median_dec = round(decs[middle], 5)
			new_head["EPOCHMED"] = (median_epoch, "Epoch for median target coordinates")
			new_head["RAMED"] = (median_ra, "Target median RA, deg")
			new_head["DECMED"] = (median_dec, "Target median Dec, deg")
		
		else:
			half_duration_1 = np.sum(tracking_intervals[:low_middle]/24)
			half_duration_2 = np.sum(tracking_intervals[:high_middle]/24)
			half_duration = (half_duration_1 + half_duration_2) * 0.5
			median_epoch = round(initial_epoch + half_duration, 7)
			new_head["EPOCHMED"] = (median_epoch, "Epoch for median target coordinates")
			
			median_ra = round((ras[low_middle] + ras[high_middle]) * 0.5, 5)
			new_head["RAMED"] = (median_ra, "Target median RA, deg")
			
			median_dec = round((decs[low_middle] + decs[high_middle]) * 0.5, 5)
			new_head["DECMED"] = (median_dec,"Target median Dec, deg")
		
		ra_rates = np.array([heads[x[:-5]]["RATRACK"] for x in files])
		median_ra_rates = round(np.median(ra_rates)*15, 10)
		new_head["MEDRATRA"] = (median_ra_rates, "Median RA differential tracking rate, deg/hr")
		
		dec_rates = np.array([heads[x[:-5]]["DECTRACK"] for x in files])
		median_dec_rates = round(np.median(dec_rates), 10)
		new_head["MEDRATDE"] = (median_dec_rates, "Median Dec differential tracking rate, deg/hr")
	
	proper_motion_ra = round(heads[key_0]["PMRA"]*15, 6)
	proper_motion_dec = round(heads[key_0]["PMDEC"], 6)	
	new_head["PMRA"] = (proper_motion_ra, "Target proper motion, ''/yr")
	new_head["PMDEC"] = (proper_motion_dec, "Target proper motion, ''/yr")
	
	new_head["TRACKWAV"] = (heads[key_0]["WAVELENG"], "Effective target tracking wavelength, angstrom")
	new_head["PWFS1_ST"] = (heads[key_0]["PWFS1_ST"], "PWFS1 probe state (frozen/guiding/parked)")
	new_head["PWFS2_ST"] = (heads[key_0]["PWFS2_ST"], "PWFS2 probe state (frozen/guiding/parked)")
	new_head["OIWFS_ST"] = (heads[key_0]["OIWFS_ST"], "OIWFS probe state (frozen/guiding/parked)")
	new_head["AOWFS_ST"] = (heads[key_0]["AOWFS_ST"], "AOWFS probe state (frozen/guiding/parked)")
	
	guide_probes = {"PWFS1_ST":"P1", "PWFS2_ST":"P2", "OIWFS_ST":"OI", "AOWFS_ST":"AO"}
	guide_probe_list = [x for x, y in new_head.items() if y=="guiding"]
	guide_probe_key = guide_probe_list[0]
	guide_probe = guide_probe_key[:-3]
	guider = guide_probes[guide_probe_key]
	
	guide_star_name = heads[key_0][guider+"AOBJEC"]
	comment_string = guide_probe + " guide star name"
	new_head[guider+"AOBJEC"] = (guide_star_name, comment_string)
	
	guide_star_ra = round(heads[key_0][guider+"ARA"], 5)
	comment_string = guide_probe + " guide star RA"
	new_head[guider+"ARA"] = (guide_star_ra, comment_string)
	
	guide_star_dec = round(heads[key_0][guider+"ADEC"], 5)
	comment_string = guide_probe + " guide star DEC"
	new_head[guider+"ADEC"] = (guide_star_dec, comment_string)
	
	guide_star_rv = heads[key_0][guider+"ARV"]
	comment_string = guide_probe + " guide star heliocentric radial velocity"
	new_head[guider+"ARV"] = (guide_star_rv, comment_string)
	
	guiding_wavelength = heads[key_0][guider+"AWAVEL"]
	comment_string = guide_probe + " effective target tracking wavelength, ang"
	new_head[guider+"AWAVEL"] = (guiding_wavelength, comment_string)
	
	guide_star_epoch = heads[key_0][guider+"AEPOCH"]
	comment_string = guide_probe + " guide star coordinate epoch"
	new_head[guider+"AEPOCH"] = (guide_star_epoch, comment_string)
	
	guide_star_equinox = heads[key_0][guider+"AEQUIN"]
	comment_string = guide_probe + " guide star coordinate equinox"
	new_head[guider+"AEQUIN"] = (guide_star_equinox, comment_string)
	
	guide_star_frame = heads[key_0][guider+"AFRAME"]
	comment_string = guide_probe + " guide star coordinate system"
	new_head[guider+"AFRAME"] = (guide_star_frame, comment_string)
	
	guide_star_pmra = round(heads[key_0][guider+"APMRA"]*15, 6)
	comment_string = guide_probe + " guide star proper motion in RA, ''/yr"
	new_head[guider+"APMRA"] = (guide_star_pmra, comment_string)
	
	guide_star_pmdec = round(heads[key_0][guider+"APMDEC"], 6)
	comment_string = guide_probe + " guide star proper motion in Dec, ''/yr"
	new_head[guider+"APMDEC"] = (guide_star_pmdec, comment_string)
	
	guide_star_parallax = heads[key_0][guider+"APARAL"]
	comment_string = guide_probe + " guide star parallax"
	new_head[guider+"APARAL"] = (guide_star_parallax, comment_string)
	
	guider_focus_offset = heads[key_0][guider+"FOCUS"]
	comment_string = guide_probe + " focus offset, mm"
	new_head[guider+"FOCUS"] = (guider_focus_offset, comment_string)
	
	guider_frequency = heads[key_0][guider+"FREQ"]
	comment_string = guide_probe + " sampling  frequency, Hz"
	new_head[guider+"FREQ"] = (guider_frequency, comment_string)
	
	# Telescope metadata
	new_head["SCOPE"] = (heads[key_0]["TELESCOP"], "Telescope")
	new_head["CGUIDMOD"] = (heads[key_0]["CGUIDMOD"], "Carousel guide mode")
	new_head["M2BAFFLE"] = (heads[key_0]["M2BAFFLE"], "Position of M2 baffle")
	new_head["M2CENBAF"] = (heads[key_0]["M2CENBAF"], "Position of M2 central hole baffle")
	dithers = []
	for file_name in files:
		y_offset = round(heads[file_name[:-5]]["YOFFSET"])
		if y_offset in dithers:
			continue
		i = len(dithers)
		dithers.append(y_offset)
		dither_key_string = "DITHER" + ("0"*(8-len(str(i+1))-6)) + str(i+1)
		comment_string = "Spatial dither pattern offset point #" + str(i+1) + ", arcsec"
		new_head[dither_key_string] = (y_offset, comment_string)
	
	new_head["PA"] = (heads[key_0]["PA"], "Sky position angle at sequence start")
	new_head["IAA"] = (heads[key_0]["IAA"], "Instrument Alignment Angle")
	
	science_fold_rots = np.array([heads[x[:-5]]["SFRT2"] for x in files])
	median_sf_angle = round(np.median(science_fold_rots), 12)
	new_head["MEDSFRT"] = (median_sf_angle, "Median science fold rotation angle, deg")
	
	science_fold_tilts = np.array([heads[x[:-5]]["SFTILT"] for x in files])
	median_sf_tilt = round(np.median(science_fold_tilts), 4)
	new_head["MEDSFTIL"] = (median_sf_tilt, "Median science fold tilt angle, deg")
	
	science_fold_linears = np.array([heads[x[:-5]]["SFLINEAR"] for x in files])
	median_sf_linear = np.median(science_fold_linears)
	new_head["MEDSFLIN"] = (median_sf_linear, "Median science fold linear position, mm")
	
	new_head["AOFOLD"] = (heads[key_0]["AOFOLD"], "AO pick-off mirror position")
		
	# Instrument metadata
	new_head["INSTRUME"] = (heads[key_0]["INSTRUME"], "Instrument used in observation")
	
	# Stacking metadata
	new_head["STACKED"] = (len(files), "Number of stacked frames")
	for i ,file_name in enumerate(files):
		stack_key_string = "STACK" + ("0"*(8-len(str(i+1))-5)) + str(i+1)
		comment_string = "Stacked file #" + str(i+1)
		new_head[stack_key_string] = (file_name, comment_string)
	
	for i in files:
		print(headers[i[:-5]]["HUMIDITY"])
	return new_hdu


def scale_spectrum(spectrum_data, scaling_index):
	scaling_region = spectrum_data[1][scaling_index-20:scaling_index+21]
	scaling_median = np.nanmedian(scaling_region)
	scaling_sigma = np.nanstd(scaling_region)
	upper_scaling_limit = scaling_median+(3*scaling_sigma)
	lower_scaling_limit = scaling_median-(3*scaling_sigma)
	scaling_region = scaling_region[
		np.where(
			np.logical_and(
				scaling_region < upper_scaling_limit,
				scaling_region > lower_scaling_limit
			)
		)
	]
	scale = np.median(scaling_region)
	scaled_spectrum = spectrum_data[1]/scale
	scaled_error = spectrum_data[2]/scale
	
	return [scaled_spectrum, scaled_error, scale]


def stacking(spectra, errors, wavelength_axis, number_of_spectra):
	stacked_spectrum = []
	stacked_errors = []
	for i in range(len(wavelength_axis)):
		uncertain_samples = []
		for j in range(number_of_spectra):
			uncertain_samples.append(
				np.random.normal(
					loc=spectra[i][j],
					scale=errors[i][j],
					size=100
				)
			)
			
		uncertain_samples = np.array(uncertain_samples)

		stacked_spectrum.append(np.median(uncertain_samples))
		stacked_errors.append(np.std(uncertain_samples)/((len(files)-1)**0.5))
	
	stacked_spectrum = np.array(stacked_spectrum)
	stacked_errors = np.array(stacked_errors)
	
	return np.array([wavelength_axis, stacked_spectrum, stacked_errors])	
	

parser = argparse.ArgumentParser(
	description="This script takes multiple 1D spectra and combines\
	them with a weighted median to produce a stacked 1D spectrum with\
	reduced noise."
)
parser.add_argument(
	"-p", "--plot", action="store_true", 
	help="[boolean] Plot the stacked spectrum on screen."
)
parser.add_argument(
	"-s", "--save", action="store_true",
	help="[boolean] Save the stacked spectrum to a new .fits file."
)
parser.add_argument(
	"-w", "--scaling_wavelength", default=0.0,
	help="[float] Wavelength at which all spectra are scaled to unity\
	before stacking. Must have the same units as the wavelength axis of\
	the spectrum. The median of the twenty data points both sides of\
	the chosen wavelength (40 points total) will be used as the scaling\
	factor. Defaults to '0.0' which will result in the spectra being\
	scaled to unity at their central wavelength."
)

args = parser.parse_args()

files = sorted(glob.glob("*.fits"))

headers = {}
optimal_frames = {}
aperture_frames = {}
scaled_optimal_spectra = []
scaled_optimal_errors = []
scaled_aperture_spectra = []
scaled_aperture_errors = []
optimal_scales = {}
aperture_scales = {}

for f in files:
	with fits.open(f) as file_hdu_list:
		headers[f[:-5]] = file_hdu_list[0].header
		optimal_frames[f[:-5]] = file_hdu_list[0].data
		aperture_frames[f[:-5]] = file_hdu_list[1].data
		optimal_frame = file_hdu_list[0].data
		aperture_frame = file_hdu_list[1].data

	wavelength_axis = optimal_frame[0]
	optimal_spectrum = optimal_frame[1]
	optimal_error = optimal_frame[2]
	aperture_spectrum = aperture_frame[1]
	aperture_error = aperture_frame[2]
	
	if args.scaling_wavelength == 0.0:
		args.scaling_wavelength += (wavelength_axis[0]+wavelength_axis[-1])*0.5
	
	scaling_index = np.argmin(np.abs(wavelength_axis-args.scaling_wavelength))
	
	scaled_optimal_data = scale_spectrum(optimal_frame, scaling_index)
	scaled_optimal_spectra.append(scaled_optimal_data[0])
	scaled_optimal_errors.append(scaled_optimal_data[1])
	optimal_scales[f[:-5]] = scaled_optimal_data[2]
	
	scaled_aperture_data = scale_spectrum(aperture_frame, scaling_index)
	scaled_aperture_spectra.append(scaled_aperture_data[0])
	scaled_aperture_errors.append(scaled_aperture_data[1])
	aperture_scales[f[:-5]] = scaled_aperture_data[2]

scaled_optimal_spectra = np.array(scaled_optimal_spectra).T
scaled_optimal_errors = np.array(scaled_optimal_errors).T
scaled_aperture_spectra = np.array(scaled_aperture_spectra).T
scaled_aperture_errors = np.array(scaled_aperture_errors).T

stacked_optimal_data = stacking(
	scaled_optimal_spectra,
	scaled_optimal_errors,
	wavelength_axis,
	len(files)
)

stacked_aperture_data = stacking(
	scaled_aperture_spectra,
	scaled_aperture_errors,
	wavelength_axis,
	len(files)
)

if args.plot:
	fig = plt.figure(figsize=(10,6))
	ax = plt.subplot()
	ax.errorbar(
		wavelength_axis,
		stacked_optimal_data[1],
		yerr=stacked_optimal_data[2],
		marker='.',
		linestyle='',
		color="k",
		label="Optimal Stack"
	)
	ax.errorbar(
		wavelength_axis,
		stacked_aperture_data[1] + 1,
		yerr=stacked_aperture_data[2],
		marker='.',
		linestyle='',
		color="dimgray",
		label="Aperture Stack +1"
	)
	ax.grid(linestyle="dashed")
	ax.set_xlabel("Wavelength, " + headers[f[:-5]]["HIERARCH WAVU"])
	ax.set_ylabel("Scaled Counts")
	ax.legend()
	plt.show()

if args.save:
	save_dict = {"stacked" : np.concatenate(
			[stacked_optimal_data, stacked_aperture_data[1:]]
		)}
	
	for i in headers:
		save_dict[i] = np.concatenate(
			[optimal_frames[i], aperture_frames[i][1:]]
		)
		
	instrument = headers[files[0][:-5]]["INSTRUME"]
	
	stack_header_dict = {
		"GMOS-N" : stack_header_gmos,
		"GMOS-S" : stack_header_gmos
	}
	
	new_hdu = fits.PrimaryHDU()
	new_header = new_hdu.header
	new_header["DATE"] = (
		datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"),
		"File creation date."
	)
	new_header["FITSDOI"] = (
		"10.1051/0004-6361:20010923", "FITS format definition paper DOI"
	)
	new_header["ORIGIN"] = ("stack.py", "Script that created this file")
	new_header["VERSION"] = ("v0.0.0", "stack.py version")
	new_header["DOI"] = ("UNKNOWN", "Script repository DOI")
	stack_hdu = stack_header_dict[instrument](new_hdu, headers, files)
	stack_hdu.data = save_dict["stacked"]
	print(stack_hdu.header)
