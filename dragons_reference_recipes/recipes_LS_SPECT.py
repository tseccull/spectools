# Reference file showing the reduceMinorPlanetScience() recipe I've added to
# my copy of recipes_LS_SPECT.py in DRAGONS. This function is added to 
# ~/miniforge3/envs/dragons/lib/python3.10/site-packages/geminidr/gmos/recipes/sq/recipes_LS_SPECT.py


def reduceMinorPlanetScience(p):
    """
    This is a reduced version of reduceScience() that calls a subset of
    the primitives provided by that recipe. This recipe only runs
    preparation, DQ and VAR frame addition, overscan correction, bias
    subtraction, ADU to e- converion, flat-field correction, QE
    correction, and 2D spectrum distortion correction (rectification).
    Cosmic ray flagging, fringe subtraction, sky subtraction,
    extraction, and stacking are all performed later.
    
    This recipe function is called when DRAGONS is run as part of the
    dagrons.py script.

    Parameters
    ----------
    p : :class:`geminidr.gmos.primitives_gmos_longslit.GMOSLongslit`

    """
    p.prepare()
    p.addDQ()
    p.addVAR(read_noise=True)
    p.overscanCorrect()
    p.biasCorrect()
    p.ADUToElectrons()
    p.addVAR(poisson_noise=True)
    p.attachWavelengthSolution()
    p.flatCorrect()
    p.QECorrect()
    p.distortionCorrect()
    p.storeProcessedScience(suffix="_2D")
