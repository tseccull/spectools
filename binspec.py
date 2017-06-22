#! /Users/tseccull/anaconda2/bin/python
# Author T. Seccull 2016
# A script for binning a spectrum.

from spectools import wavaxis
import argparse
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as num
import os
import sys

from matplotlib.backends.backend_pdf import PdfPages

# Parse arguments
parser = argparse.ArgumentParser(description='A script for binning a spectrum.')
parser.add_argument('datafile', type=str, help='spectrum file to be binned')
parser.add_argument('width', type=int, help='bin width, number of points per bin')
parser.add_argument('-r', '--rejection', type=float, help='Number that when multiplied by the standard '
                                                                     'deviation defines the boundaries for rejection of points from a bin.')
parser.add_argument('-p', '--plot', help='plot the binned spectrum and show on screen as a continuous line', action='store_true')
parser.add_argument('-e', '--errs', help='plot the binned spectrum and show on screen as a scatter plot showing errorbars', action='store_true')
parser.add_argument('-s', '--save', help='save the stacked spectrum to a new .fits file', action='store_true')
parser.add_argument('-sfp', '--save_plot', help='save a plot of the normalised binned spectrum, only compatible with X-Shooter spectra which are all arms stitched together', action='store_true')
parser.add_argument('-n', '--norm', help='normalise the spectrum', action='store_true')
args = parser.parse_args()

# Load in the spectrum to be binned.
with fits.open(args.datafile) as han:
    handata = han[0].data
    hhead = han[0].header

# Create wavelength axis
try:
    wav = handata[0]
    hdata = handata[1]
except:
    wav = wavaxis(hhead)

# If the spectrum is from FORS observations reduced with GME, mask out pixels with cosmic rays based on the mask provided.
if hhead['INSTRUME'] == 'FORS2' and 'HIERARCH GME UTC EXTRACTION DATE' in hhead:
    quals = handata[2]
    hdata[num.where(quals < 1)] = num.nan

# Bin the spectrum
outw = []
outs = []
uncFlux = []

for i in range(0, len(wav), args.width):
    # Quit the loop at the end of the wavelength range being binned.
    if i + args.width / 2 >= len(wav):
        break

    # Bin the wavelength axis
    outw.append(wav[i + args.width / 2])

    # If rejection is required remove points from the bin which are beyond (rejection_factor * sigma) from the mean.
    if args.rejection:
        x = hdata[(1 * max(i - args.width, 0)):(1 * min(i + args.width, len(hdata)))]
        m = num.nanmean(x)
        std = num.nanstd(x)
        pointskept = num.where(num.logical_and((x - m) < args.rejection * std, (x - m) > (args.rejection / -1) * std))
        xkept = x[pointskept]
        outs.append(num.median(xkept))
        if len(xkept) > 0:
            uncFlux.append((num.nanmin(xkept) - num.nanmax(xkept))/2)
        else:
            uncFlux.append(num.nan)
    # Else bin the spectrum with a median.
    else:
        outs.append(num.median(hdata[i:min(i + args.width, len(hdata))]))

outw = num.array(outw)
uncWav = num.array([args.width/2] * len(outw)) * 0.02
outs = num.array(outs)

# Plot the spectrum if required.
if args.plot:
    if args.errs:
        plt.errorbar(outw, outs, xerr=uncWav, yerr=uncFlux, fmt='.')
    elif args.norm:
        normIdx = num.argmin(num.abs(outw - 5630))
        outs /= num.mean(outs[normIdx - 50:normIdx + 49])
        plt.plot(outw, outs)
    else:
        plt.plot(outw, outs)
    #plt.ylim(-0.001, 0.001)
    plt.title(args.datafile[:-5] + '_BIN_B' + str(args.width) + '_R' + str(args.rejection))
    plt.show()

# Save the spectrum to a new fits file if required.
if args.save:
    if args.norm:
        normIdx = num.argmin(num.abs(outw - 5630))
        outs /= num.mean(outs[normIdx - 50:normIdx + 49])

    arrayList = [outw, outs]
    hdu = fits.PrimaryHDU(arrayList, header=hhead)
    hdulist = fits.HDUList(hdu)

    if os.path.isfile(args.datafile[:-5] + '_BIN_B' + str(args.width) + '_R' + str(args.rejection) + '.fits'):
        print '[BINSPEC] ' + args.datafile[:-5] + '_BIN_B' + str(args.width) + '_R' + str(args.rejection) + '.fits ' \
              'already exists. New stack will not be saved. Terminating.'
        sys.exit()
    else:
        hdulist.writeto(args.datafile[:-5] + '_BIN_B' + str(args.width) + '_R' + str(args.rejection) + '.fits')
        hdulist.close()
        cwd = os.getcwd()
        print '[BINSPEC] Binned spectrum saved at ' + cwd + '/' + args.datafile[:-5] + '_BIN_B' + str(args.width) \
              + '_R' + str(args.rejection) + '.fits. Errs not saved.'

if args.save_plot:
    pp = PdfPages(args.datafile[:-5] + '.pdf')
    normIdx = num.argmin(num.abs(outw - 5630))
    outs /= outs[normIdx]
    plt.plot(outw, outs)
    plt.xlabel('Wavelength, nm')
    plt.ylabel('Reflectance normalised at 658 nm')
    plt.title(hhead['OBJECT'])
    pp.savefig()
    pp.close()
