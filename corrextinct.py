#! /Users/tseccull/anaconda2/bin/python


# SCRIPT THAT USES THE EXTINCTION TABLE PROVIDED FOR THE FORS@ INSTRUMENT TO CORRECT ANY EXTICNTION IN THE SPECTRUM
# CAUSED BY AIRMASS. THIS SCRIPT REQUIRES THE PRESENCE OF THE 'fors2_extinct_table.fits'FILE IN  THE DIRECTORY
# CONTAINING THE SPECTRA TO BE CORRECTED.
import astropy.io.fits as fits
import glob
import matplotlib.pyplot as plt
import scipy.interpolate as interp

with fits.open('fors2_extinct_table.fits') as han:
    print han.info()
    extab = han[1].data

files = glob.glob('*LSS*.fits')

for i in files:
    with fits.open(i) as han:
        head = han[0].header
        wav = han[0].data[0]
        flux = han[0].data[1]
        if head['INSTRUME'] == 'FORS2' and head['HIERARCH GME UTC EXTRACTION DATE']:
            quals = han[0].data[3]

    interptab = interp.interp1d(extab['WAVE'], extab['EXTINCTION'])
    interpextab = interptab(wav)

    corrspec = flux * (10 ** (0.4 * head['HIERARCH ESO TEL AIRM START'] * interpextab))

    plt.plot(wav, flux, wav, corrspec)
    plt.show()

    hdu = fits.PrimaryHDU([wav, corrspec, quals], header=head)
    hdulist = fits.HDUList(hdu)
    hdulist.writeto('_'.join(i.split('_')[0:-1]) + '_EXT_' + i.split('_')[-1])
    hdulist.close()
