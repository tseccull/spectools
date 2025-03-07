#!/usr/bin/env python3

"""
moffat.py

	Copyright (C) 2024-05-06 Tom Seccull
	
	This script is part scrap.py, in the spectools repo hosted at 
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

	Last updated - 2025-03-07

	Description --------------------------------------------------------
	This file contains the moffat_least_squares() and moffat_resid()
	functions which are needed to determine the FWHM of spectra when
	scrap.py is run. moffat_least_squares is called by prep_gmos()
	functions in the instrument io modules.
"""

__author__ = "Tom Seccull"

import numpy as np

from scipy.optimize import least_squares


def moffat_least_squares(
	spatial_axis, column, seeing, pixel_resolution, end_clip
):
    """
    Takes a data column, spatial axis and seeing of the observation and
    fits a Moffat function to the column using a least squares method.
    Returns the best fit parameters of the Moffat function.

    Args:
     -- spatial_axis (numpy.ndarray)
			The spatial axis of the data being fit.
     -- column (numpy.ndarray)
			The data being fitted.
     -- seeing (float)
			The estimated FWHM of the spatial profile.
     -- pixel_resolution (float)
			The spatial resolution of each pixel in arcsec/pixel.
     -- end_clip (int)
			The number of pixels at each end of the spatial profile
			array to ignore when fitting the Moffat profile.

    Returns:
     -- parameter_list (list)
			The list of best fit output parameters returned by the least
			squares routine.
    """

	# Clip the median spatial profile to be fitted based on the value of
	# end_clip
    column[:end_clip] = np.median(column)
    column[-end_clip:] = np.median(column)

    # Set up initial conditions for the least squares fit.
    # x0 = [
    #	 amplitude,
    # 	 centre,
    #	 alpha,
    #	 beta,
    #	 background gradient,
    #	 background level
    # ]
    # Initial beta estimate comes from optimal value from atmospheric
    # turbulence theory as described in 
    # Trujillo, I. et al. (2001), MNRAS, 328, 977-985
    # See https://ui.adsabs.harvard.edu/abs/2001MNRAS.328..977T/abstract
    x0 = [
        np.nanmedian(np.sort(column)[-3:]),
        np.argmax(column),
        seeing / pixel_resolution,
        4.765,
        0.0,
        np.median(np.concatenate((column[:5], column[-5:]))),
    ]

    # Run the least squares fit.
    res_lsq = least_squares(
        moffat_resid,
        x0,
        bounds=(
            [
                0.0, 
                np.argmax(column) - 1.0,
                0.0, 
                0.0, 
                0.0, 
                -np.inf
            ],
            [
                np.inf,
                np.argmax(column) + 1,
                (5 * seeing / pixel_resolution),
                5.0,
                np.inf,
                np.inf,
            ],
        ),
        args=(spatial_axis, column),
        method="trf",
        ftol=1e-12,
    )
    
    parameter_list = [
        res_lsq.x[0],
        res_lsq.x[1],
        res_lsq.x[2],
        res_lsq.x[3],
        res_lsq.x[4],
        res_lsq.x[5]
    ]
    return parameter_list


def moffat_resid(x, spatial_axis, data):
    """
    Calculates residuals of fitted moffat profile and the data for the
    least squares fitting.

    Description:
        A = x[0]
        c = x[1]
        alpha = x[2]
        beta = x[3]
        B = x[4]
        m = x[5]

    Args:
     -- x (numpy.ndarray)
			An array of parameters defining the shape of the model
			moffat profile.
     -- spatial_axis (numpy.ndarray)
			The spatial axis of the data.
     -- data (numpy.ndarray)
			The data.

    Returns:
     -- residual (numpy.ndarray)
			The residual array between the model moffat profile and the
			data.
    """
    
    moffat_profile = x[0]*(
		(1+((spatial_axis-x[1])*(spatial_axis-x[1]))/(x[2]*x[2]))**-x[3]
	)
    
    residual = moffat_profile + x[4] + (spatial_axis*x[5]) - data
    
    return residual
