#! /home/tom/anaconda3/envs/work/bin/python
"""
gmos_stack_head.py - written by Tom Seccull, 2024-05-06

	Called by: stack.py
	Last updated: 2024-05-06
	
	This file contains the stack_header_gmos() function, which creates
	the astropy header for the HDU containing the stacked GMOS spectra
	produced by stack.py.
"""

import copy
import datetime
import numpy as np


def stack_header_gmos(new_hdu, heads, files, scale_wavelength):
	"""
	Creates the header for the fits Header Data Unit containing the
	stacked spectra. Returns new_hdu with that new header.
	
	Args:
	 -- new_hdu (astropy header data unit)
			The new HDU containing to contain the stacked spectra and
			their metadata.
	 -- heads (dict)
			Dictionary containing the headers of all the spectra
			combined to make the stacked spectrum.
	 -- files (list)
			List of the filenames of the spectra being stacked.
	 -- scale_wavelength (float)
			The wavelength at which the spectra were scaled prior to
			being stacked.
	
	Returns:
	 -- new_hdu (astropy header data unit)
			Function returns the HDU that was entered with an updated
			and expanded header.
	"""
	
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
	
	# Weather/condition metadata
	mean_airmasses = np.array([heads[x[:-5]]["AIRMASS"] for x in files])
	start_airmasses = np.array([heads[x[:-5]]["AMSTART"] for x in files])
	end_airmasses = np.array([heads[x[:-5]]["AMEND"] for x in files])
	airmasses = np.concatenate(
		[start_airmasses, mean_airmasses, end_airmasses]
	)
	min_airmass = round(np.min(airmasses), 3)
	max_airmass = round(np.max(airmasses), 3)
	new_head["LOAIRMSS"] = (min_airmass, "Minimum sequence airmass")
	new_head["HIAIRMSS"] = (max_airmass, "Maximum sequence airmass")
	
	condition_keys = [
		"HUMIDITY", "TAMBIENT", "PRESSUR2", "DEWPOINT", "WINDSPEE", "WINDSPE2"
	]
	new_condition_keys = ["HUMID", "TAMBI", "PRESS", "DWPNT", "WNDMH", "WNDMS"]
	condition_comments = [
		"relative humidity, %",
		"ambient temperature, C",
		"atmospheric pressure, Pa", 
		"dew point, C", "wind speed, m/s", 
		"wind speed, mph"
	] 
	
	for i in range(len(condition_keys)):
		conditions = np.array(
			[heads[x[:-5]][condition_keys[i]] for x in files]
		)
		min_condition = np.min(conditions)
		max_condition = np.max(conditions)
		med_condition = np.median(conditions)
		
		new_head["LO"+new_condition_keys[i]] = (
			min_condition, "Minimum " + condition_comments[i]
		)
		new_head["HI"+new_condition_keys[i]] = (
			max_condition, "Maximum " + condition_comments[i]
		)
		new_head["MED"+new_condition_keys[i]] = (
			med_condition, "Median " + condition_comments[i]
		)
	
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
	
	guide_keys = [
		"AOBJEC", "ARA", "ADEC", "ARV", "AWAVEL", "AEPOCH", "AEQUIN",
		"AFRAME", "APMRA", "APMDEC", "APARAL", "FOCUS", "FREQ"
	]
	guide_comments = [
		" guide star name",
		" guide star RA",
		" guide star DEC",
		" guide star heliocentric radial velocity",
		" effective target tracking wavelength, ang",
		" guide star coordinate epoch",
		" guide star coordinate equinox",
		" guide star coordinate system",
		" guide star proper motion in RA, ''/yr",
		" guide star proper motion in DEC, ''/yr",
		" guide star parallax",
		" focus offset, mm",
		" sampling frequency, Hz"
	]
	
	for i in range(len(guide_keys)):
		guide_star_parameter = heads[key_0][guider+guide_keys[i]]
		comment_string = guide_probe + guide_comments[i]
		new_head[guider+guide_keys[i]] = (guide_star_parameter, comment_string)
	
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
		comment_string = (
			"Spatial dither pattern offset point #" + str(i+1) + ", arcsec"
		)
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
	new_head["INPORT"] = (heads[key_0]["INPORT"], "ISS port where GMOS was attached")
	new_head["GMOSCC"] = (heads[key_0]["GMOSCC"], "GMOS components controller software")
	new_head["CONID"] = (heads[key_0]["CONID"], "Detector controller ID")
	new_head["DETECTOR"] = (heads[key_0]["DETECTOR"], "Detector name")
	new_head["DEWAR"] = (heads[key_0]["DEWAR"], "Dewar name")
	new_head["ARRYTSET"] = (heads[key_0]["ARRYTSET"], "Array temperature setpoint, C")
	
	temps = ["A", "B", "C", "D"]
	for t in temps:
		array_temps = np.array([heads[x[:-5]]["ARRYTMP"+t] for x in files])
		median_temp = round(np.median(array_temps), 3)
		comment_string = "Median array temperature " + t + ", C"
		new_head["MEDARTP"+t] = (median_temp, comment_string)
	
	new_head["DETSIZE"] = (heads[key_0]["DETSIZE"], "Detector size, pixels")
	new_head["NCCDS"] = (heads[key_0]["NCCDS"], "Number of CCD chips")
	new_head["NAMPS"] = (heads[key_0]["NAMPS"], "Number of amplifiers")
	new_head["AMPINTEG"] = (heads[key_0]["AMPINTEG"], "Amplifier integration time")
	new_head["NSUBEXP"] = (heads[key_0]["NSUBEXP"], "Number of sub exposures")
	new_head["EXPTIME"] = (heads[key_0]["EXPTIME"], "Exposure time, s")
	new_head["DARKTIME"] = (round(heads[key_0]["DARKTIME"], 2), "Dark current integration, s")
	new_head["MASKID"] = (heads[key_0]["MASKID"], "Mask/IFU barcode")
	new_head["MASKNAME"] = (heads[key_0]["MASKNAME"], "Mask name")
	new_head["MASKTYPE"] = (heads[key_0]["MASKTYP"], "Mask/IFU type (-1=IFU/0=none/1=mask)")
	new_head["MASKLOC"] = (heads[key_0]["MASKLOC"], "Mask/IFU location (-1=unknown/0=FP/1=cassette)")
	new_head["FILTER1"] = (heads[key_0]["FILTER1"], "Filter 1 name")
	new_head["FILTID1"] = (heads[key_0]["FILTID1"], "Filter 1 barcode")
	new_head["FILTER2"] = (heads[key_0]["FILTER2"], "Filter 2 name")
	new_head["FILTID2"] = (heads[key_0]["FILTID2"], "Filter 2 barcode")
	new_head["GRATING"] = (heads[key_0]["GRATING"], "Grating name")
	new_head["GRATID"] = (heads[key_0]["GRATID"], "Grating barcode")
	new_head["CENTWAVE"] = (heads[key_0]["CENTWAVE"], "Central wavelength, nm")
	new_head["GRWLEN"] = (heads[key_0]["GRWLEN"], "Grating wavelength at slit, nm")
	new_head["GRORDER"] = (heads[key_0]["GRORDER"], "Grating order")
	new_head["GRTILT"] = (heads[key_0]["GRTILT"], "Grating tilt angle, deg")
	new_head["DTMODE"] = (heads[key_0]["DTMODE"], "Detector translation stage mode")
	new_head["GMOSDC"] = (heads[key_0]["GMOSDC"], "GMOS detector controller s/w")
	new_head["DETTYPE"] = (heads[key_0]["DETTYPE"], "Detector array type")
	new_head["DETID"] = (heads[key_0]["DETID"], "Chip IDs")
	new_head["DETRO1X"] = (heads[key_0]["DETRO1X"], "ROI 1 X start")
	new_head["DETRO1XS"] = (heads[key_0]["DETRO1XS"], "ROI 1 X size")
	new_head["DETRO1Y"] = (heads[key_0]["DETRO1Y"], "ROI 1 Y start")
	new_head["DETRO1YS"] = (heads[key_0]["DETRO1YS"], "ROI 1 Y size")
	new_head["PIXSCALE"] = (heads[key_0]["PIXSCALE"], "Pixel scale in Y, ''/pixel")
	new_head["CCDSUM"] = (heads[key_0]["CCDSUM"], "Detector binning, pixels")
	new_head["GAIN"] = (heads[key_0]["GAIN"], "Amplifier gain")
	new_head["GAINMULT"] = (heads[key_0]["GAINMULT"], "Gain multiplication")
	new_head["RDNOISE"] = (heads[key_0]["RDNOISE"], "Readout noise")
	
	# Data reduction metadata
	new_head["OBSMODE"] = (heads[key_0]["OBSMODE"], "Observing mode (Image/IFU/MOS/LONGSLIT)")
	
	reduction_datetimes = [heads[x[:-5]]["GEM-TLM"] for x in files]
	reduction_datetimes = [
		datetime.datetime.fromisoformat(x) for x in reduction_datetimes
	]
	last_reduction_timestamp = max(reduction_datetimes).strftime("%Y-%m-%dT%H:%M:%S")
	new_head["GEM-TLM"] = (last_reduction_timestamp, "Latest UT data reduction timestamp")
	new_head["CRMETHOD"] = (heads[key_0]["CRMETHOD"], "Cosmic ray masking/cleaning method in scrap.py")
	
	reduction_keys = [
		"GPREPARE", "GGAIN", "GIREDUCE", "GMOSAIC", "GSAPPWAV",
		"GSREDUCE", "WMEF", "GSTRANSF", "CRDATE", "FRNGDATE", "MOTESDAT"
	]
	reduction_comments = copy.deepcopy(reduction_keys)
	reduction_comments[-3] = "scrap.py"
	reduction_comments[-2] = "fronge.py"
	reduction_comments[-1] = "MOTES"
	for i, k in enumerate(reduction_keys):
		date_times = [heads[x[:-5]][k] for x in files]
		date_times = [datetime.datetime.fromisoformat(x) for x in date_times]
		last_time_stamp = max(date_times).strftime("%Y-%m-%dT%H:%M:%S")
		timestamp_comment = "Latest UT timestamp for " + reduction_comments[i]
		new_head[k] = (last_time_stamp, timestamp_comment)
	
	new_head["WAVU"] = (heads[key_0]["WAVU"], "Wavelength unit")
	low_wave_comment = "Lower limit of wavelength range, " + new_head["WAVU"]
	new_head["WAVL"] = (heads[key_0]["WAVL"], low_wave_comment)
	high_wave_comment = "Upper limit of wavelength range, " + new_head["WAVU"]
	new_head["WAVH"] = (heads[key_0]["WAVH"], high_wave_comment)
	
	iqs = np.array([heads[x[:-5]]["IQ"] for x in files])
	min_iq = np.min(iqs)
	max_iq = np.max(iqs)
	med_iq = np.median(iqs)
	new_head["LOIQ"] = (min_iq, "Best image quality - lowest FWHM, arcsec ")
	new_head["HIIQ"] = (max_iq, "Worst image quality - highest FWHM, arcsec")
	new_head["MEDIQ"] = (med_iq, "Median image quality - median FWHM, arcsec")
	new_head["SNBINLIM"] = (heads[key_0]["SNBINLIM"], "maximum SNR per bin for extraction")
	new_head["COLBNLIM"] = (heads[key_0]["COLBNLIM"], "minimum number of columns per bin")
	new_head["FWHMMULT"] = (heads[key_0]["FWHMMULT"], "FWHM multiplier used to bound extraction region")
	new_head["SKYBINLM"] = (heads[key_0]["SKYBINLM"], "max SNR per bin for sky subtraction")
	new_head["SKYFWHM"] = (heads[key_0]["SKYFWHM"], "FWHM multiplier used to define the backround")
	
	# Stacking metadata
	new_head["STACKED"] = (len(files), "Number of stacked frames")
	for i ,file_name in enumerate(files):
		stack_key_string = "STACK" + ("0"*(8-len(str(i+1))-5)) + str(i+1)
		comment_string = "Stacked file #" + str(i+1)
		new_head[stack_key_string] = (file_name, comment_string)
	
	new_head["SCALEWAV"] = (round(scale_wavelength, 2), "Wavelength where spectra are scaled to 1.")
	new_head["HDUROW0"] = "Wavelength Axis, angstroms"
	new_head["HDUROW1"] = "Stacked Scaled Optimal Spectrum"
	new_head["HDUROW2"] = "Stacked Scaled Optimal Spectrum Uncertainty"
	new_head["HDUROW3"] = "Stacked Scaled Aperture Spectrum"
	new_head["HDUROW4"] = "Stacked Scaled Aperture Spectrum Uncertainty"
	new_head["EXTNAME"] = "STACKED_1D_SPEC"
		
	return new_hdu
