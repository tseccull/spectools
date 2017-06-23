#! /Users/tseccull/anaconda2/bin/python

# USES A 1 X 7 BOXCAR TO CONVOLVE A SPECTRUM OBSERVED WITH A 0.7" SLIT SOP THAT IT MATCHES THE RESOLTUION OF A SPECTRUM
# OBSERVED WITH A 1.0" SLIT.
import astropy.io.fits as fits
import glob
import numpy as np
import scipy.signal as signal

files = glob.glob('*.fits')
box = np.repeat(1./7., 10)

for i in files:
    with fits.open(i) as han:
        head = han[0].header
        wav = han[0].data[0]
        data = han[0].data[1]

    conv_data = signal.convolve(data, box, mode='same')

    hdu = fits.PrimaryHDU([wav, conv_data], header=head)
    hdulist = fits.HDUList(hdu)
    hdulist.writeto('_'.join(i.split('_')[0:-1]) + '_CON_' + i.split('_')[-1])
    hdulist.close()
