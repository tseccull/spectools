#! /Users/tseccull/anaconda2/bin/python

# A SCRIPT FOR MEDIAN OR MEAN STACKING OF MULTIPLE SPECTRA INTO A SINGLE SPECTRUM. THIS SCRIPT WILL ONLY WORK IN A
# DIRECTORY CONTAINING ONLY THE SPECTRA TO BE STACKED, WHICH MUST ALL SHARE A COMMON WAVELENGTH AXIS (I.E. FOR XSHOOTER
# SPECTRA THEY MUST ALL BE FROM THE SAME ARM). IT'S ALSO ASSUMED ALL THE SPECTRA HAVE BEEN OBSERVED USING THE SAME
# INSTRUMENT.

import argparse
import astropy.io.fits as fits
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import general

# PARSE COMMAND LINE ARGUMENTS
parser = argparse.ArgumentParser(description='A script for stacking multiple spectra into a single spectrum. Must '
                                             'be performed in a directory containing only the spectra to be stacked.')
parser.add_argument('-m', '--method', default='median', help='stacking method, takes mean or median, default is median')
parser.add_argument('-p', '--plot', help='plot the stacked spectrum and show on screen', action='store_true')
parser.add_argument('-s', '--save', help='save the stacked spectrum to a new .fits file', action='store_true')
args = parser.parse_args()

# GLOB PATHS TO SPECTRA TO BE STACKED
specpaths = glob.glob('*.fits')

# SET UP BLANK LIST WHICH WILL HAVE THE FLUX DATA APPENDED TO IT FOR STACKING
datalist = []

# EXTRACT INFORMATION FROM THE HEADERS OF THE .FITS FILES AS WELL AS DATA WHICH WILL BE STORED IN DATALIST
for i, path in enumerate(specpaths):
    with fits.open(path) as han:
        head = han[0].header
        datalist.append(han[0].data[1])

        # OBTAIN/MODIFY THE COSMIC RAY MASK SO IT MATCHES THE STACKED SPECTRUM AND OBTAIN THE WAVELENGTH AXIS
        if head['INSTRUME'] == 'FORS2':
            if 'HIERARCH GME UTC EXTRACTION DATE' in head:
                if i == 0:
                    wav = han[0].data[0]
                    masterquals = han[0].data[2]
                else:
                    masterquals *= han[0].data[2]
            else:
                increment = head['CD1_1']
                wav = general.wavaxis(head, increment)

        elif head['INSTRUME'] == 'XSHOOTER':
            wav = han[0].data[0]

    # EXTRACT HEADER INFORMATION
    if i == 0:
        if head['INSTRUME'] == 'XSHOOTER':
            arm = head['HIERARCH ESO SEQ ARM']

        obj = head['OBJECT'].replace(' ', '_')

# NORMALISE THE SPECTRA
normdatalist = []
for k in datalist:
    if head['INSTRUME'] == 'FORS2':
        normwav = 5000

    elif head['INSTRUME'] == 'XSHOOTER':
        if arm == 'UVB':
            normwav = 540
        elif arm == 'VIS':
            normwav = 654
        elif arm == 'NIR':
            normwav = 1050

    normidx = np.argmin(np.abs(wav - normwav))
    k /= np.mean(k[normidx - 50:normidx + 49])

    normdatalist.append(k)

# STACK THE SPECTRA USING A MEAN
if args.method == 'mean':
    print '[STACKSPEC] Stack method = mean'
    sumspec = np.nansum(normdatalist, axis=0)
    stackspec = np.divide(sumspec, len(normdatalist))

# STACK THE SPECTRA BY USING A MEDIAN
if args.method == 'median':
    print '[STACKSPEC] Stack method = median'

    stackspec = []
    median_array = np.zeros((len(normdatalist), len(normdatalist[0])))

    for j in range(len(normdatalist)):
        median_array[j] = normdatalist[j]

    columnmax = np.amax(median_array, axis=0)
    columnmin = np.amin(median_array, axis=0)
    columnmaxmin = np.array([columnmax, columnmin])
    mediansum = np.nansum(columnmaxmin, axis=0)
    stackspec = np.divide(mediansum, 2)

# PLOT THE SPECTRUM IF SELECTED
if args.plot:
    plt.plot(wav, stackspec)
    plt.title('STACK_' + str.upper(args.method) + '_' + obj)
    plt.show()

# SAVE THE SPECTRUM TO A NEW FITS FILE IN THE PARENT DIRECTORY IF SELECTED
if args.save:
    if head['INSTRUME'] == 'FORS2' and head['HIERARCH GME UTC EXTRACTION DATE']:
        hdu = fits.PrimaryHDU([wav, stackspec, masterquals.astype(bool)], header=head)
    else:
        hdu = fits.PrimaryHDU([wav, stackspec], header=head)

    hdulist = fits.HDUList(hdu)
    cwd = os.getcwd()

    if os.path.isfile(os.path.dirname(cwd) + '/' + obj + '_' + str.upper(args.method) + '_STACK_' + '.fits'):
        print '[STACKSPEC] ' + os.path.dirname(cwd) + '/' + obj + '_' + str.upper(args.method) + '_STACK_' + '.fits ' \
              'already exists. New stack will not be saved. Terminating.'
        sys.exit()
    else:
        hdulist.writeto(os.path.dirname(cwd) + '/' + obj + '_' + str.upper(args.method) + '_STACK_' + '.fits')
        hdulist.close()
        print '[STACKSPEC] Stacked spectrum saved at ' + os.path.dirname(cwd) + '/' + obj + '_' + \
              str.upper(args.method) + '_STACK_' + '.fits'
