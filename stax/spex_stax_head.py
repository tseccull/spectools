#!/usr/bin/env python3

"""
	spex_stack_head.py

	Copyright (C) 2025-07-29 Tom Seccull
	
	This script is part stax.py, in the spectools repo hosted at 
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

	Last updated - 2025-08-06

	Description --------------------------------------------------------
	This file contains the stack_header_spex() function, which creates
	the astropy header for the HDU containing the stacked SpeX spectra
	produced by stax.py.
"""

import copy
import datetime
import numpy as np

from astropy import units as u
from astropy.coordinates import Angle


def last_reduction_modifications(new_head, heads, files):
	"""
	Takes header timestamps for last interaction with IDL/Spextool from 
	a set of spectra to be stacked, finds the latest time for each, and 
	assigns that value to the "IDLSTAMP" timestamp in the new stacked
	header.
	
	Args:
	 -- new_head (astropy fits header)
	        The header being constructed for the new stacked spectrum.
	 -- heads (list)
	        A list of the headers being combined to make new_head.
	 -- files (list)
	        List of the filenames of the spectra being stacked.
	
	Returns:
	 -- new_head (astropy fits header)
	        The header being constructed for the new stacked spectrum.
	"""
	
	date_times = []
	month_dict = {
		"Jan" : "01", "Feb" : "02", "Mar" : "03", "Apr" : "04",
		"May" : "05", "Jun" : "06", "Jul" : "07", "Aug" : "08",
		"Sep" : "09", "Oct" : "10", "Nov" : "11", "Dec" : "12"
	}
	
	for file_name in files:
		idl_timestamp_list = heads[file_name[:-5]].comments["SIMPLE"].split(" ")
		idl_timestamp_list = [x for x in idl_timestamp_list if x != '']
		idl_timestamp_string = (
			idl_timestamp_list[7] + "-" +
			month_dict[idl_timestamp_list[4]] + "-" +
			idl_timestamp_list[5].zfill(2) + "T" +
			idl_timestamp_list[6] 
		)
		date_times.append(idl_timestamp_string)
	
	date_times = [datetime.datetime.fromisoformat(x) for x in date_times]
	last_time_stamp = max(date_times).strftime("%Y-%m-%dT%H:%M:%S")
	
	new_head["IDLSTAMP"] = (last_time_stamp, "Timestamp for last change by IDL/Spextool")
	
	return new_head


def median_angle(heads, files, angle_flag=None):
	"""
	Calculates median Right Ascension, Declination, or Hour Angle 
	(depending on angle flag) from header keywords in heads.
	
	Args:
	 -- heads (list)
	        A list of the headers being combined during stacking.
	 -- files (list)
	        List of the filenames of the spectra being stacked.
	
	Returns:
	 -- median_tcs_angle (numpy.float64)
	        Median value for the angles being queried.
	"""
	
	unit_dict = {
		"ra" : u.hourangle,
		"dec" : u.deg,
		"ha" : u.deg
	}
	
	angle_keys = {
		"ra" : "TCS_RA",
		"dec" : "TCS_DEC",
		"ha" : "TCS_HA"
	}
	
	angs = []
	for file_name in files:
		angle = Angle(
			heads[file_name[:-5]][angle_keys[angle_flag]], 
			unit_dict[angle_flag]
		)
		angs.append(angle.deg)
	tcs_angles = np.array(angs)
	if np.max(tcs_angles)-np.min(tcs_angles) > 180 and angle_flag == "ra":
		[x-360 for x in tcs_angles if x > 180]
	median_tcs_angle = np.median(tcs_angles)
	if median_tcs_angle < 0 and angle_flag == "ra":
		median_tcs_angle += 360
	
	return median_tcs_angle
	

def median_tcs_value(heads, files, head_key):
	"""
	Returns the median value for a particular card in the headers of the
	stacked files.
	
	Args:
	-- heads (dict)
			Dictionary containing the headers of all the spectra
			combined to make the stacked spectrum.
	-- files (list)
			List of the filenames of the spectra being stacked.
	-- head_key (string)
			Keyword for the header card being processed.
			
	Returns:
	 -- median value (np.float64, string)
			If the card values are successfully medianed, then this
			is the median value. If not, an "Insufficient Data" string
			will be returned.
	"""
	try:
		tcs_values = np.array([heads[x[:-5]][head_key] for x in files])
		median_value = np.median(tcs_values)
	except:
		median_value = "Insufficient Data"
	
	return median_value


def stack_header_spex(new_hdu, heads, files, scale_wavelength):
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
	
	last_reduction_modifications(new_head, heads, files)
	new_head["OPERINST"] = ("Uni. of Hawaii Inst. for Astronomy", "Observatory operating institute")
	new_head["SCOPE"] = (heads[key_0]["TELESCOP"], "Telescope")
	new_head["INSTRUME"] = (heads[key_0]["INSTRUME"], "Instrument used in observation")
	new_head["OBSERVER"] = (heads[key_0]["OBSERVER"], "Observer")
	new_head["OBJECT"] = (heads[key_0]["OBJECT"], "Target name")
	observation_datetime = heads[key_0]["DATE_OBS"] + "T" + heads[key_0]["TIME_OBS"]
	new_head["UT"] = (observation_datetime, "UT at sequence start")
	new_head["MJD"] = (heads[key_0]["MJD_OBS"], "MJD at sequence start")
	new_head["INTTIME"] = (heads[key_0]["ITIME"], "Integration time, s")
	new_head["COADDS"] = (heads[key_0]["CO_ADDS"], "Number of coadds")
	new_head["BEAMPAT"] = (heads[key_0]["BEAMPAT"], "Beam pattern")
	new_head["NDR"] = (heads[key_0]["NDR"], "Number of non-destructive reads")
	new_head["TABLE_MS"] = (heads[key_0]["TABLE_MS"])
	new_head["QTHLAMP"] = (heads[key_0]["QTH_LAMP"], "Quart 0.8-2.5 10W lamp")
	new_head["INCLAMP"] = (heads[key_0]["INC_LAMP"], "Inc 0.8-2.5 0.1W lamp")
	new_head["IR_SRC"] = (heads[key_0]["IR_SRC"])
	new_head["ARG_SRC"] = (heads[key_0]["ARG_SRC"])
	new_head["CALMIR"] = (heads[key_0]["CALMIR"], heads[key_0].comments["CALMIR"])
	new_head["DIT"] = (heads[key_0]["DIT"], heads[key_0].comments["DIT"])
	new_head["OSF"] = (heads[key_0]["OSF"], heads[key_0].comments["OSF"])
	new_head["ROTATOR"] = (heads[key_0]["ROTATOR"], heads[key_0].comments["ROTATOR"])
	new_head["POSANGLE"] = (heads[key_0]["POSANGLE"], heads[key_0].comments["POSANGLE"])
	new_head["SLIT"] = (heads[key_0]["SLIT"], heads[key_0].comments["SLIT"])
	new_head["GRAT"] = (heads[key_0]["GRAT"], heads[key_0].comments["GRAT"])
	new_head["GFLT"] = (heads[key_0]["GFLT"], heads[key_0].comments["GFLT"])
	new_head["PLATE_SC"] = (heads[key_0]["PLATE_SC"], heads[key_0].comments["PLATE_SC"])
	new_head["DATATYPE"] = (heads[key_0]["DATATYPE"], heads[key_0].comments["DATATYPE"])
	
	median_tcs_ut = median_tcs_value(heads, files, "TCS_UTC")
	new_head["MEDTCSUT"] = (median_tcs_ut, "Median TCS UT MJD time")
	
	median_tcs_ra = median_angle(heads, files, angle_flag="ra")
	median_tcs_dec = median_angle(heads, files, angle_flag="dec")
	median_tcs_ha = median_angle(heads, files, angle_flag="ha")
	new_head["MEDTCSRA"] = (median_tcs_ra, "Median TCS Right Ascension, deg FK5 J2000")
	new_head["MEDTCSDC"] = (median_tcs_dec, "Median TCS Declination, deg FK5 J2000")
	new_head["MEDTCSHA"] = (median_tcs_ha, "Median TCS Hour Angle, deg")
	
	median_tcs_az = median_tcs_value(heads, files, "TCS_AZ")
	new_head["MEDTCSAZ"] = (median_tcs_az, "Median TCS Azimuth, deg")
	
	median_tcs_el = median_tcs_value(heads, files, "TCS_EL")
	new_head["MEDTCSEL"] = (median_tcs_el, "Median TCS Elevation, deg")
	
	airmasses = np.array([heads[x[:-5]]["TCS_AM"] for x in files])
	median_tcs_am = np.median(airmasses)
	min_tcs_am = np.min(airmasses)
	max_tcs_am = np.max(airmasses)
	new_head["MINTCSAM"] = (min_tcs_am, "Minimum TCS Airmass")
	new_head["MEDTCSAM"] = (median_tcs_am, "Median TCS Airmass")
	new_head["MAXTCSAM"] = (max_tcs_am, "Maximum TCS Airmass")
	
	median_tcs_focus = median_tcs_value(heads, files, "TCS_FOC")
	new_head["MEDTCSFC"] = (median_tcs_focus, "Median TCS actual focus")

	median_tcs_air_temp = median_tcs_value(heads, files, "TCS_AIRT")
	new_head["MEDTCSAT"] = (median_tcs_air_temp, "Median TCS outside air temp, deg C")
	
	median_tcs_humidity = median_tcs_value(heads, files, "TCS_HUM")
	new_head["MEDTCSHM"] = (median_tcs_humidity, "Median TCS humiditiy, 0-100")
	
	median_tcs_wind_speed = median_tcs_value(heads, files, "TCS_WSP")
	new_head["MEDTCSWS"] = (median_tcs_wind_speed, "Median TCS wind speed")
	
	median_tcs_wind_direction = median_tcs_value(heads, files, "TCS_WDIR")
	new_head["MEDTCSWD"] = (median_tcs_wind_direction, "TCS wind direction, deg")
	
	new_head["INSTR"] = (heads[key_0]["INSTR"], heads[key_0].comments["INSTR"])
	new_head["MODENAME"] = (heads[key_0]["MODENAME"], heads[key_0].comments["MODENAME"])
	new_head["ITOT"] = (heads[key_0]["ITOT"], heads[key_0].comments["ITOT"])
	new_head["LINCRMAX"] = (heads[key_0]["LINCRMAX"], heads[key_0].comments["LINCRMAX"])
	new_head["FLAT"] = (heads[key_0]["FLAT"], heads[key_0].comments["FLAT"])
	new_head["WAVECAL"] = (heads[key_0]["WAVECAL"], heads[key_0].comments["WAVECAL"])
	new_head["WAVETYPE"] = (heads[key_0]["WAVETYPE"], heads[key_0].comments["WAVETYPE"])
	new_head["NAPS"] = (heads[key_0]["NAPS"], heads[key_0].comments["NAPS"])
	new_head["NORDERS"] = (heads[key_0]["NORDERS"], heads[key_0].comments["NORDERS"])
	new_head["ORDERS"] = (heads[key_0]["ORDERS"], heads[key_0].comments["ORDERS"])
	new_head["SLTH_PIX"] = (heads[key_0]["SLTH_PIX"], heads[key_0].comments["SLTH_PIX"])
	new_head["SLTH_ARC"] = (heads[key_0]["SLTH_ARC"], heads[key_0].comments["SLTH_ARC"])
	new_head["SLTW_PIX"] = (heads[key_0]["SLTW_PIX"], heads[key_0].comments["SLTW_PIX"])
	new_head["SLTW_ARC"] = (heads[key_0]["SLTW_ARC"], heads[key_0].comments["SLTW_ARC"])
	new_head["RP"] = (heads[key_0]["RP"], heads[key_0].comments["RP"])
	new_head["SLTH_PIX"] = (heads[key_0]["SLTH_PIX"], heads[key_0].comments["SLTH_PIX"])
	new_head["AP01RAD"] = (heads[key_0]["AP01RAD"], heads[key_0].comments["AP01RAD"])
	new_head["DISPO01"] = (heads[key_0]["DISPO01"], heads[key_0].comments["DISPO01"])
	new_head["BGSTART"] = (heads[key_0]["BGSTART"], heads[key_0].comments["BGSTART"])
	new_head["BGWIDTH"] = (heads[key_0]["BGWIDTH"], heads[key_0].comments["BGWIDTH"])
	new_head["BGORDER"] = (heads[key_0]["BGORDER"], heads[key_0].comments["BGORDER"])
	new_head["WAVU"] = (heads[key_0]["XUNITS"], heads[key_0].comments["XUNITS"])
	new_head["YUNITS"] = (heads[key_0]["YUNITS"], heads[key_0].comments["YUNITS"])
	new_head["SPXTLVER"] = (heads[key_0]["VERSION"], heads[key_0].comments["VERSION"])
	
	# Stacking metadata
	new_head["STACKED"] = (len(files), "Number of stacked frames")
	for i ,file_name in enumerate(files):
		stack_key_string = "STACK" + ("0"*(8-len(str(i+1))-5)) + str(i+1)
		comment_string = "Stack file"
		new_head[stack_key_string] = (file_name, comment_string)
	
	new_head["SCALEWAV"] = (round(scale_wavelength, 2), "Wavelength where spectra are scaled to 1.")
	new_head["HDUROW0"] = "Wavelength Axis, " + new_head["WAVU"]
	new_head["HDUROW1"] = "Stacked Scaled Spectrum"
	new_head["HDUROW2"] = "Stacked Scaled Spectrum Uncertainty"
	new_head["EXTNAME"] = "STACKED_1D_SPEC"
	
	return new_hdu
