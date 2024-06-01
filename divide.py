#! /home/tom/anaconda3/envs/work/bin/python
"""
divide.py - written by Tom Seccull, 2024-05-08 - v0.0.3

	Last updated: 2024-06-01
	
	This script divides one 1D spectrum by another. divide.py expects
	both spectra to have a common wavelength axis, be the product
	of extraction by MOTES, and stacking by stack.py. 
"""

import argparse
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.widgets import Slider


def division(wave_axis, shifted_wave_axis, spectra, uncertainties):
	shifted_spectra = [spectra[0]]
	shifted_uncertainties = [uncertainties[0]]

	shifted_spectra.append(
		np.interp(wave_axis, shifted_wave_axis, spectra[1], left=0., right=0.)
	)

	shifted_uncertainties.append(
		(
			np.interp(
				wave_axis,
				shifted_wave_axis,
				spectra[1] + uncertainties[1],
				left=0., 
				right=0.
			)
			- shifted_spectra[1]
		)
	)

	mult_spec = shifted_spectra[0] * shifted_spectra[1]
	qual_1d = [1 if x!=0. else 0 for x in mult_spec]

	if offset > 0:
		zero_edges = [x+1 for x in range(len(qual_1d)) if qual_1d[x:x+2] == [0,1]]
		zero_edges = np.array(zero_edges)
		qual_1d = np.array(qual_1d)	
		qual_1d[zero_edges] = 0
	elif offset < 0:
		zero_edges = [x for x in range(len(qual_1d)) if qual_1d[x:x+2] == [1,0]]
		zero_edges = np.array(zero_edges)
		qual_1d = np.array(qual_1d)	
		qual_1d[zero_edges] = 0
	else:
		qual_1d = np.array(qual_1d)

	shifted_spectra[0][qual_1d==0] = 0
	shifted_spectra[1][qual_1d==0] = 0
	shifted_uncertainties[0][qual_1d==0] = 0
	shifted_uncertainties[1][qual_1d==0] = 0

	shifted_spectra[0][shifted_spectra[1]==0] = 0
	shifted_spectra[1][shifted_spectra[0]==0] = 1
	
	ratio_spectrum = shifted_spectra[0]/shifted_spectra[1]
	
	shifted_uncertainties[0][shifted_spectra[0]==0] = 0
	shifted_spectra[0][shifted_uncertainties[0]==0] = 1
	
	shifted_uncertainties[1][shifted_spectra[1]==0] = 0
	shifted_spectra[1][shifted_uncertainties[1]==0] = 1
	
	ratio_uncertainties = np.abs(ratio_spectrum)*((((shifted_uncertainties[0]/shifted_spectra[0])**2)+((shifted_uncertainties[1]/shifted_spectra[1])**2))**0.5)
	
	return np.array([wavelength_axis, ratio_spectrum, ratio_uncertainties, qual_1d])

def update_wave(val):
	global offset
	offset = wave_slider.val
	two_spec.set_data(wavelength_axis+offset, optimal_spectra[1])
	fig.canvas.draw()


parser = argparse.ArgumentParser(
	description="This script divides one 1D spectrum by another.\
	divide.py expects both spectra to have a common wavelength axis, be\
	the product of extraction by MOTES, and stacking by stack.py."
)
parser.add_argument("spec_file_one", type=str,
	help="[str] File name of the first spectrum, which takes the 'numerator'\
	role in the division." 
)
parser.add_argument("spec_file_two", type=str,
	help="[str] File name of the second spectrum, which takes the\
	'denominator' role in the division." 
)
parser.add_argument("-p", "--plot", action="store_true",
	help="[boolean] Plot the ratio spectrum on screen."
)
parser.add_argument("-s", "--save", action="store_true",
	help="[boolean] Save the ratio spectrum to a new .fits file."
)
args = parser.parse_args()

with fits.open(args.spec_file_one) as one:
	one_headers = [x.header for x in one]
	one_frames = [x.data for x in one]

with fits.open(args.spec_file_two) as two:
	two_headers = [x.header for x in two]
	two_frames = [x.data for x in two]

primary_head_one = one_headers[0]
primary_head_two = two_headers[0]
wavelength_axis = one_frames[0][0]

optimal_spectra = np.array([one_frames[0][1], two_frames[0][1]])
optimal_uncertainties = np.array([one_frames[0][2], two_frames[0][2]])
aperture_spectra = np.array([one_frames[0][3], two_frames[0][3]])
aperture_uncertainties = np.array([one_frames[0][4], two_frames[0][4]])

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.19)
one_spec = ax.errorbar(wavelength_axis, optimal_spectra[0], yerr = optimal_uncertainties[0], marker='.')
two_spec, = ax.plot(wavelength_axis, optimal_spectra[1], marker='.', zorder=2)
ax.set_xlabel(r"$\lambda$, " + primary_head_one["WAVU"])
ax.set_ylabel("Relative Counts")

ax_wave_slider = plt.axes([0.25, 0.05, 0.5, 0.03])
wave_slider = Slider(
	ax_wave_slider,
	r"$\lambda$ Offset, " + primary_head_two["WAVU"],
	-20,
	20,
	valinit=0,
	valfmt="%.1f",
	valstep=.1
)
wave_slider.on_changed(update_wave)

plt.show()

try:
	offset = round(offset,1)
except NameError:
	offset = 0.0

shifted_wavelength_axis = wavelength_axis + offset

optimal_ratio_frame = division(
	wavelength_axis,
	shifted_wavelength_axis,
	optimal_spectra,
	optimal_uncertainties
)

plt.errorbar(optimal_ratio_frame[0], optimal_ratio_frame[1], yerr=optimal_ratio_frame[2])
plt.plot(optimal_ratio_frame[0], optimal_ratio_frame[3])
plt.show()
exit()

optimal_shifted_spectra = [optimal_spectra[0]]
optimal_shifted_uncertainties = [optimal_uncertainties[0]]

optimal_shifted_spectra.append(
	np.interp(
		wavelength_axis,
		shifted_wavelength_axis,
		optimal_spectra[1],
		left=0., 
		right=0.
	)
)

optimal_shifted_uncertainties.append(
	(
		np.interp(
			wavelength_axis,
			shifted_wavelength_axis,
			optimal_spectra[1] + optimal_uncertainties[1],
			left=0., 
			right=0.
		)
		- optimal_shifted_spectra[1]
	)
)

mult_spec = optimal_shifted_spectra[0] * optimal_shifted_spectra[1]
qual_1d = [1 if x!=0. else 0 for x in mult_spec]

if offset > 0:
	zero_edges = [x+1 for x in range(len(qual_1d)) if qual_1d[x:x+2] == [0,1]]
	zero_edges = np.array(zero_edges)
	qual_1d = np.array(qual_1d)	
	qual_1d[zero_edges] = 0
elif offset < 0:
	zero_edges = [x for x in range(len(qual_1d)) if qual_1d[x:x+2] == [1,0]]
	zero_edges = np.array(zero_edges)
	qual_1d = np.array(qual_1d)	
	qual_1d[zero_edges] = 0
else:
	qual_1d = np.array(qual_1d)

optimal_shifted_spectra[0][qual_1d==0] = 0
optimal_shifted_spectra[1][qual_1d==0] = 0
optimal_shifted_uncertainties[0][qual_1d==0] = 0
optimal_shifted_uncertainties[1][qual_1d==0] = 0

optimal_shifted_spectra[0][optimal_shifted_spectra[1]==0] = 0
optimal_shifted_spectra[1][optimal_shifted_spectra[0]==0] = 1
qual_1d[optimal_shifted_spectra[0]==0] = 0

divided_optimal_spectrum = optimal_shifted_spectra[0]/optimal_shifted_spectra[1]
plt.plot(wavelength_axis, divided_optimal_spectrum)
plt.plot(wavelength_axis, qual_1d)
plt.plot(wavelength_axis, optimal_shifted_spectra[0])
plt.plot(wavelength_axis, optimal_shifted_spectra[1])
plt.show()

print(offset)
