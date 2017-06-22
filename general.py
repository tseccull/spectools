#! /Users/tseccull/anaconda2/bin/python

import numpy as np


# Take .fits header parameters and use them to create a wavelength axis.
def wavaxis(header):
    numpix = header['NAXIS1']
    startwavlen = header['CRVAL1']
    increment = header['CDELT1']

    return startwavlen + (np.arange(numpix) * increment)
