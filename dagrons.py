#!/home/tom/miniforge3/envs/dragons/bin/python
"""
dagrons.py = written by Tom Seccull, 2024-08-11 - v1.0.3
	
	Last updated: 2025-02-25
	
	This script is a partial reduction pipeline for GMOS longslit
	spectroscopic data built with DRAGONS and based on its GMOS longslit
	spectroscopy reduction tutorial. Some aspects of the reduction have
	been automated to some degree, but the bones are the same as what is
	presented in the DRAGONS documentation. The recipe called by
	dagrons.py for reduction of science data and where it should be
	pasted in DRAGONS is provided in
	`./dragons_reference_recipes/recipes_LS_SPECT.py`. The directory
	where this script is initiated is taken to be the working directory,
	with the directory containing data and calibrations for the
	reduction provided as an argument. The data and calibrators for only
	one on-sky source are expected to be present in the data directory,
	such that all the data and calibration frames have the same
	format/RoI. This means science targets and their specphot standard
	stars may have to be stored and reduced separately, but this script
	is less complex as a result. dagrons.py only calls a subset of the
	primitives provided by the full reduceScience() recipe that DRAGONS
	normally runs for GMOS longslit spectroscopic data including
	preparation, DQ and VAR frame addition, overscan correction, bias 
	subtraction, ADU to e- converion, flat-field correction, QE
	correction, and 2D spectrum distortion correction (rectification).
	Cosmic ray flagging, fringe subtraction, sky subtraction,
	extraction, and stacking are all performed later by other scripts in
	the spectools repo.DRAGONS should be cited if dagrons.py is used.
"""


# Check DRAGONS miniforge env has been activated. This script assumes that 
# DRAGONS has been properly installed and configured before use.
import os
if os.environ["CONDA_PREFIX"] != "/home/tom/miniforge3/envs/dragons":
	print("\n...Run 'conda activate dragons' and try again...\n")
	exit()


# Import libraries
import argparse
import astrodata
import gemini_instruments
import glob
import shutil

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system import cal_service
from recipe_system.reduction.coreReduce import Reduce
from recipe_system.utils.reduce_utils import normalize_ucals


# Parse command line arguments.
parser = argparse.ArgumentParser(
	description="This script is a partial reduction pipeline for GMOS\
	longslit spectroscopic data built around DRAGONS and its GMOS\
	longslit spectroscopy reduction tutorial. Some aspects of the\
	reduction have been automated to a greater degree, but the bones\
	are the same as what is presented in the DRAGONS documentation. The\
	directory where this script is initiated is taken to be the working\
	directory, with the directory containing data and calibrations for\
	the reduction provided as an argument. The data and calibrators for\
	a only one on-sky source are expected to be present in the data\
	directory, such that all the data and calibration frames have the\
	same format/RoI. This means science targets and their specphot\
	standard stars may have to be stored and reduced separately, but\
	this script is less complex as a result. dagrons.py only calls a\
	subset of the primitives provided by the full reduceScience()\
	recipe that DRAGONS normally runs for GMOS longslit spectroscopic\
	data including preparation, DQ and VAR frame addition, overscan\
	correction, bias subtraction, ADU to e- converion, flat-field\
	correction, QE correction, and 2D spectrum distortion correction\
	(rectification). Cosmic ray flagging, fringe subtraction, sky\
	subtraction, extraction, and stacking are all performed later by\
	other scripts in the spectools repo."
)
parser.add_argument(
	"data_directory", type=str,
	help="The path to the directory containing the science and\
	calibration data to be reduced."
)
parser.add_argument(
	"-i", "--interactive", action="store_true",
	help="When activated this argument will execute a run of this\
	script where all interactive elements are switched on."
)
args = parser.parse_args()


# Set up the calibration service
# https://dragons.readthedocs.io/projects/gmosls-drtutorial/en/release-3.2.x/cal_service.html#cal-service
dragonsrc_path = "/home/tom/.dragons/dragonsrc"
working_directory = os.getcwd()
with open(dragonsrc_path, "r+") as dragonsrc_file:
	dragonsrc_lines = dragonsrc_file.readlines()
	dragonsrc_lines[-1] = (
		"databases= " + working_directory + "/cal_manager.db get store\n"
	)
	dragonsrc_file.seek(0)
	dragonsrc_file.writelines(dragonsrc_lines)
	dragonsrc_file.truncate()

# Delete the log file from a previous run and set up the DRAGONS logger.
if os.path.isfile(working_directory + "/dagrons.log"):
	os.remove(working_directory + "/dagrons.log")
logutils.config(file_name='dagrons.log')

# If calibration files are still present from a previous run, delete them.
if os.path.isdir(working_directory + "/calibrations/"):
	shutil.rmtree(working_directory + "/calibrations/")

old_file_subpaths = [
	"/*_bias.fits",
	"/*_flat.fits",
	"/*_arc.fits",
	"/*_mosaic.pdf",
	"/*_prepared.fits",
	"/*_dqAdded.fits",
	"/*_varAdded.fits",
	"/*_overscanCorrected.fits",
	"/*_biasCorrected.fits",
	"/*_ADUToElectrons.fits",
	"/*_varAdded_crash.fits",
	"/*_wavelengthSolutionAttached.fits",
	"/*_flatCorrected.fits",
	"/*_QECorrected.fits",
	"/*_distortionCorrected.fits",
	"/*_2D.fits"
]
for subpath in old_file_subpaths:
	old_files = glob.glob(working_directory + subpath)
	if len(old_files) > 0:
		[os.remove(x) for x in old_files]

# (Re)Initialize the calibration database. If a calibration database already
# exists from a previous run of this script it will be deleted and recreated.
if os.path.isfile(working_directory + "/cal_manager.db"):
	os.remove(working_directory + "/cal_manager.db")
caldb  = cal_service.set_local_database()
caldb.init()

# Create a list of all the files in the data directory.
if args.data_directory[-1] == "/":
	all_files = glob.glob(args.data_directory + "*.fits")
else:
	all_files = glob.glob(args.data_directory + "/*.fits")
all_files.sort()

# Make file lists for all file types
biases = dataselect.select_data(all_files, ["BIAS"])
flats = dataselect.select_data(all_files, ["FLAT"])
arcs = dataselect.select_data(all_files, ["ARC"])
target_frames = dataselect.select_data(all_files, [], ["CAL"])
target_frames = [[x] for x in target_frames]

# Store the bad pixel mask in the calibration database
[caldb.add_cal(bpm) for bpm in dataselect.select_data(all_files, ["BPM"])]

# Make Master Bias
reduce_bias = Reduce()
reduce_bias.files.extend(biases)
reduce_bias.runr()
master_bias = reduce_bias.output_filenames[0]
master_bias_path = (
	working_directory + "/calibrations/processed_bias/" + master_bias
)

# Master Flat Field
reduce_flats = Reduce()
reduce_flats.files.extend(flats)
if args.interactive:
	reduce_flats.uparms = [("interactive", True)]
reduce_flats.recipename = "makeProcessedFlatStack"
reduce_flats.runr()
master_flat = reduce_flats.output_filenames[0]
master_flat_path = (
	working_directory + "/calibrations/processed_flat/" + master_flat
)

# Calculate Wavelength Solution
reduce_arc = Reduce()
reduce_arc.files.extend(arcs)
if args.interactive:
	reduce_arc.uparms = [("interactive", True)]
reduce_arc.runr()
master_arc = reduce_arc.output_filenames[0]
master_arc_path = (
	working_directory + "/calibrations/processed_arc/" + master_arc
)

# Collect processed calibs for target frame reduction.
processed_calibs = [
	"processed_bias:" + master_bias_path,
	"processed_flat:" + master_flat_path,
	"processed_arc:" + master_arc_path
]

# Reduce target frames. A reference copy of the reduceMinorPlanetScience() 
# recipe is stored in ./dragons_reference_recipes/recipes_LS_SPECT.py
for frame in target_frames:
	reduce_target = Reduce()
	reduce_target.files.extend(frame)
	reduce_target.ucals = normalize_ucals(processed_calibs)
	reduce_target.recipename = "reduceMinorPlanetScience"
	reduce_target.runr()
