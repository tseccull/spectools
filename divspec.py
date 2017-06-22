#! /Users/tseccull/anaconda2/bin/python

# A SCRIPT FOR DIVIDING ONE SPECTRUM BY ANOTHER. THEY MUST SHARE A COMMON WAVELENGTH. AXIS


import argparse
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as num
import os
import sys

# PARSE ARGUMENTS
parser = argparse.ArgumentParser(description='A script for dividing one spectrum by another. They must share a common'
                                             'wavelength axis.')
parser.add_argument('targFile', type=str, help='target spectrum to be divided')
parser.add_argument('starFile', type=str, help='star spectrum to be divided')
parser.add_argument('-p', '--plot', help='plot the stacked spectrum and show on screen', action='store_true')
parser.add_argument('-s', '--save', help='save the stacked spectrum to a new .fits file', action='store_true')
args = parser.parse_args()

# READ IN HEADER DATA AND SPECTRUM DATA FOR THE FIRST SPECTRUM
with fits.open(args.targFile) as han:
    hdata = han[0].data[1]
    wav = han[0].data[0]
    hhead = han[0].header
    if 'HIERARCH GME UTC EXTRACTION DATE' in hhead:
        hquals = han[0].data[2]
    if hhead['INSTRUME'] == 'XSHOOTER':
        arm = hhead['HIERARCH ESO SEQ ARM']

# READ IN HEADER DATA AND SPECTRUM DATA FOR THE SECOND SPECTRUM
with fits.open(args.starFile) as solo:
    shead = solo[0].header
    sdata = solo[0].data[1]
    if 'HIERARCH GME UTC EXTRACTION DATE' in shead:
        squals = solo[0].data[2]

targobj = hhead['OBJECT'].replace(' ', '_')
starobj = shead['OBJECT'].replace(' ', '_')

# DIVIDE THE FIRST SPECTRUM BY THE SECOND SPECTRUM
divSpec = num.divide(hdata, sdata)
divquals = num.multiply(hquals, squals)

# PLOT THE DIVIDED SPECTRUM IF SELECTED
if args.plot:
    if 'HIERARCH GME UTC EXTRACTION DATE' in hhead:
        plt.plot(wav, divSpec, wav, divquals)
    else:
        plt.plot(wav, divSpec)

    plt.title(targobj + '_' + '_DIV_' + starobj)
    plt.ylabel('Relative Reflectance')

    if hhead['INSTRUME'] == 'FORS2':
        plt.xlabel('Wavelength, Angstroms')
    if hhead['INSTRUME'] == 'XSHOOTER':
        plt.xlabel('Wavelength, nm')

    plt.show()

# SAVE THE DIVIDED SPECTRUM TO A NEW FITS FILE IF SELECTED
if args.save:
    if 'HIERARCH GME UTC EXTRACTION DATE' in hhead:
        hdu = fits.PrimaryHDU([wav, divSpec, divquals], header=hhead)
    else:
        hdu = fits.PrimaryHDU([wav, divSpec], header=hhead)

    hdulist = fits.HDUList(hdu)

    if hhead['INSTRUME'] == 'FORS2':
        if os.path.isfile(targobj + '_DIV_' + starobj + '.fits'):
            print '[DIVSPEC] ' + targobj + '_DIV_' + starobj + '.fits ' \
                  'already exists. New spectrum will not be saved. Terminating.'
            sys.exit()
        else:
            hdulist.writeto(targobj +'_DIV_' + starobj + '.fits')
            hdulist.close()
            cwd = os.getcwd()
            print '[DIVSPEC] Divided spectrum saved at ' + cwd + '/' + targobj + '_DIV_' + starobj + '.fits'

    if hhead['INSTRUME'] == 'XSHOOTER':
        if os.path.isfile(targobj + '_' + targexptime + '_DIV_' + starobj + '_' + starexptime + '_' + arm + '.fits'):
            print '[DIVSPEC] ' + targobj + '_' + targexptime + '_DIV_' + starobj + '_' + starexptime + '_' + arm + \
                  '.fits already exists. New spectrum will not be saved. Terminating.'
            sys.exit()
        else:
            hdulist.writeto(targobj + '_' + targexptime + '_DIV_' + starobj + '_' + starexptime + '_' + arm + '.fits')
            hdulist.close()
            cwd = os.getcwd()
            print '[DIVSPEC] Divided spectrum saved at ' + cwd + '/' + targobj + '_' + targexptime + '_DIV_' + \
                  starobj + '_' + starexptime + '_' + arm + '.fits'
