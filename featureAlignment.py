#!/usr/bin/python

import sys, re, numpy as np
import statsmodels as stat
import utils
from pyteomics import mzxml

def calibrateFeatures(ref, comp, parameters):
    initMzTol = float(params['tol_initial'])
    sdWidth = float(params['sd_width'])
    ref = np.sort(ref, order = "Intensity")[::-1]   # Sort features in descending order of intensity
    comp = np.sort(comp, order = "Intensity")[::-1]

    # Calibration of RT and m/z globally
    print ("Global calibration of feature is being performed")
    rtShifts, mzShifts = globalCalibration(ref, comp, initMzTol)


def globalCalibration(ref, comp, mzTol = 20):
    nPeaks = round(0.05 * ref.shape[0])    # Number of peaks to be considered for global calibration
    rtShifts = []   # Array for RT-shifts between reference and compared runs
    mzShifts = []   # Array for mz-shifts (ppm)
    i, j = 0, 1
    while (j <= nPeaks):
        z = ref['z'][i] # From the 1st feature of reference run (i.e. strongest feature)
        mz = ref['mz'][i]
        lL = mz - mz * mzTol / 1e6
        uL = mz + mz * mzTol / 1e6
        rt = ref['RT'][i]
        intensity = ref['Intensity'][i]
        if z == 0:
            # For a reference feature with undetermined charge, consider all possible charges in compared feature
            rowInd = np.where((comp['mz'] >= lL) & (comp['mz'] < uL))[0]
        else:
            rowInd = np.where((comp['mz'] >= lL) & (comp['mz'] < uL) & (comp['z'] == z))[0]

        if len(rowInd) > 0:
            rowInd = rowInd[0]
            rtShifts.append(comp['RT'][rowInd] - rt)
            mzShifts.append((comp['mz'][rowInd] - mz) / mz * 1e6)
            comp = np.delete(comp, rowInd, 0)
            j += 1

        i += 1

    # For more robust calculation, top and bottom 10% values are trimmed
    rtShifts = np.array(rtShifts)
    rtShifts = rtShifts[(rtShifts > np.percentile(rtShifts, 10)) &
                     (rtShifts < np.percentile(rtShifts, 90))]
    mzShifts = np.array(mzShifts)
    mzShifts = mzShifts[(mzShifts > np.percentile(mzShifts, 10)) &
                     (mzShifts < np.percentile(mzShifts, 90))]

    return rtShifts, mzShifts

###########################################
################ Main part ################
###########################################

paramFile = r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\jumpm_negative.params"
featureFiles = [r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.1.feature",
                r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.1.feature",
                r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.1.feature"]

nFiles = len(featureFiles)

################################
# Load parameters and features #
################################
params = utils.readParams(paramFile)
# Features from .feature files are stored in fArray. For example,
# featureFiles = [file1, file2, file3]
# fArray[0] = features from file1 (which has column names like 'index', 'mz', etc.)
# fArray[1] = features from file2
# ...
# The array of m/z values from the first feature file can be accessed by fArray[0]['mz']
# The header of .feature file is used as column names of the array
# Note that "/" (slash) is ignored when the file is loaded through "genfromtxt"
fArray = []
for i in range(0, len(featureFiles)):
    data = np.genfromtxt(featureFiles[i], delimiter = "\t", dtype = None, names = True)
    fArray.append(data)

if nFiles > 1: # Multiple feature files -> alignment is required
    print ("Feature calibration")

    ###################################
    # Selection of a reference sample #
    ###################################
    if params['reference_feature'] == "0":
        # A sample with the largest median of top 100 intensities is set to a reference run
        refNo = 0
        refIntensity = 0
        for i in range(0, len(fArray)):
            tmpIntensity = np.median(sorted(fArray[i]['Intensity'], reverse = True))
            if tmpIntensity >= refIntensity:
                refNo = i
                refIntensity = tmpIntensity
    else:
        try:
            refNo = featureFiles.index(params['reference_feature'])
        except:
            sys.exit("'reference_feature' parameter should be correctly specified")
    print ("%s is chosen as the reference run" % featureFiles[refNo])

    ############################################################
    # Calibration of features against those in a reference run #
    ############################################################
    for i in range(0, len(featureFiles)):
        if i != refNo:
            fArray[i] = calibrateFeatures(fArray[refNo], fArray[i], params)






