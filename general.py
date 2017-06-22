#! /Users/tseccull/anaconda2/bin/python

import numpy as np


# Take .fits header parameters and use them to create a wavelength axis.
def wavaxis(header, increment):
    numpix = header['NAXIS1']
    startwavlen = header['CRVAL1']

    return startwavlen + (np.arange(numpix) * increment)
