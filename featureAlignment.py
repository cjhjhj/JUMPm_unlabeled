#!/usr/bin/python

import sys, os, subprocess, numpy as np
import utils
from pyteomics import mzxml

def calibrateFeatures(ref, comp, parameters):
    # ref and comp are "ndarray"s with the following column name
    # 'index' = nominal feature index
    # 'mz' = m/z value of a feature
    # 'z' = charge of the feature (0 for undetermined)
    # 'MS1ScanNumber' = representative MS1 scan number of the feature
    # 'minMS1ScanNumber' = minimum MS1 scan number of the feature (spanned)
    # 'maxMS1ScanNumber' = maximum MS1 scan number of the feature (spanned)
    # 'RT' = representative RT (should correspond to MS1ScanNumber)
    # 'minRT' = minimum RT
    # 'maxRT' = maximum RT
    # 'intensity' = intensity of the feature
    # 'SN' = signal-to-noise ratio
    # 'PercentageTF' = percentage of true feature (= feature width over RT/scan)

    initMzTol = int(params['tol_initial'])
    sdWidth = int(params['sd_width'])
    ref = np.sort(ref, order = "Intensity")[::-1]   # Sort features in descending order of intensity
    comp = np.sort(comp, order = "Intensity")[::-1]

    # Calibration of RT and m/z globally
    print ("Global calibration of features is being performed")
    rtShifts, mzShifts = globalCalibration(ref, comp, initMzTol)
    print ("Based on the matched features within %d ppm" % initMzTol)
    print ("The global RT-shift is %.4f second" % np.median(rtShifts))
    print ("The global m/z-shift is %.4f ppm" % np.median(mzShifts))
    print ("RT and m/z of the compared features are calibrated according to the above information")
    comp['RT'] = comp['RT'] - np.median(rtShifts)
    comp['mz'] = comp['mz'] / (1 + np.median(mzShifts) / 1e6)

    # Calibration of RT and m/z locally using LOWESS
    print ()
    print ("Local calibration of features is being performed (through LOESS modeling")
    print ("RT- and m/z-tolerances will be dynamically estimated over RT- and m/z-range as follows,")
    print ("  RT- and m/z-tolerance = %d x dynamically estimated SD of RT- and m/z-shifts" % sdWidth)
    print ("LOESS modeling may take some time. Please be patient ...")
    rtSd = np.maximum(1e-3, np.std(rtShifts))
    mzSd = np.maximum(1e-3, np.std(mzShifts))
    rtShifts, mzShifts = localCalibration(ref, comp, rtSd, mzSd, params, "RT")


def localCalibration(ref, comp,rtSd, mzSd, params, type):
    n = ref.shape[0]    # Number of rows in "ref"
    if ~isinstance(rtSd, (list, np.ndarray)):
        rtSd = np.repeat(rtSd, n)
    if ~isinstance(mzSd, (list, np.ndarray)):
        mzSd = np.repeat(mzSd, n)
    # if len(rtSd) == 1:
    #     rtSd = np.repeat(rtSd, n)
    # if len(mzSd) == 1:
    #     mzSd = np.repeat(mzSd, n)
    sdWidth = float(params['sd_width'])
    rtTol = rtSd * sdWidth
    mzTol = mzSd * sdWidth

    # Look for comparable features between "ref" and "comp"
    # and prepare LOWESS modeling
    refRt, compRt, rtShifts = [], [], []
    refMz, compMz, mzShifts = [], [], []
    for i in range(0, n):   # For each feature in "ref", look for a matching one in "comp"
        z = ref['z'][i]
        mz = ref['mz'][i]
        rt = ref['RT'][i]
        intensity = ref['Intensity'][i]
        rtDev = comp['RT'] - rt
        mzDev = (comp['mz'] - mz) / comp['mz'] * 1e6    # Unit of PPM
        if z == 0:  # Undetermined charge
            # For the feature with undetermined charge,
            # look for a matching one without considering charge state
            rowInd = np.where((abs(rtDev) <= rtTol[i]) &
                              (abs(mzDev) <= mzTol[i]))[0]
        else:
            # For the feature with a charge state,
            # look for a matching one with considering charge state
            rowInd = np.where((abs(rtDev) <= rtTol[i]) &
                              (abs(mzDev) <= mzTol[i]) &
                              (comp['z'] == z))[0]
        if len(rowInd) > 0:
            # When multiple features in "comp" are matched to a feature in "ref",
            # choose the one with the highest intensity
            # Since "comp" is sorted by descending order of intensity,
            # the first one has the highest intensity
            rowInd = rowInd[0]
            refRt.append(rt)
            refMz.append(mz)
            compRt.append(comp['RT'][rowInd])
            compMz.append(comp['mz'][rowInd])
            rtShifts.append(rtDev[rowInd])
            mzShifts.append(mzDev[rowInd])

    # Perform LOWESS regression to calibrate RT and m/z
    # Note that "rpy2" does not support Windows system officially
    np.savetxt('refRt.txt', refRt)
    np.savetxt('compRt.txt', compRt)
    np.savetxt('compRt_new.txt', comp['RT'])
    rPath = "C:\\Program Files\\R\\R-3.6.2\\bin\\Rscript.exe"
    script = "lowess.R"
    cmd = [rPath, script, 'refRt.txt', 'compRt.txt', 'compRt_new.txt']
    subprocess.call(cmd, shell = False)
    print ("Check")

    return




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






