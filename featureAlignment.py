#!/usr/bin/python

import sys, os, numpy as np
import rpy2.robjects as ro
import matplotlib.pyplot as plt
from rpy2.robjects.vectors import IntVector, FloatVector
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
    # "intensity" = intensity of the feature
    # 'SN' = signal-to-noise ratio
    # 'PercentageTF' = percentage of true feature (= feature width over RT/scan)

    initMzTol = int(params['tol_initial'])
    sdWidth = int(params['sd_width'])
    ref = np.sort(ref, order = "intensity")[::-1]   # Sort features in descending order of intensity
    comp = np.sort(comp, order = "intensity")[::-1]

    # Calibration of RT and m/z globally
    print ("  Global calibration of features is being performed")
    rtShifts, mzShifts = globalCalibration(ref, comp, initMzTol)
    print ("    Based on the matched features within %d ppm" % initMzTol)
    print ("    The global RT-shift is %.4f second" % np.median(rtShifts))
    print ("    The global m/z-shift is %.4f ppm" % np.median(mzShifts))
    print ("    RT and m/z of the compared features are calibrated according to the above information")
    comp["RT"] = comp["RT"] - np.median(rtShifts)
    comp["mz"] = comp["mz"] / (1 + np.median(mzShifts) / 1e6)

    # Calibration of RT and m/z locally using LOESS (stepwise)
    print ("  Local calibration of features is being performed (through LOESS modeling")
    print ("    RT- and m/z-tolerances will be dynamically estimated over RT- and m/z-range as follows,")
    print ("    RT- and m/z-tolerance = %d x dynamically estimated SD of RT- and m/z-shifts" % sdWidth)
    print ("    LOESS modeling may take some time. Please be patient ...")
    rtSd = np.maximum(1e-3, np.std(rtShifts))
    mzSd = np.maximum(1e-3, np.std(mzShifts))
    ref, comp, rtSd, mzSd = localCalibration(ref, comp, rtSd, mzSd, params, "RT")
    print ("    The 1st round of RT-calibration is done")
    print ("      min SD of RT-shifts = %.4f second" % np.amin(rtSd))
    print ("      max SD of RT-shifts = %.4f second" % np.amax(rtSd))
    ref, comp, rtSd, mzSd = localCalibration(ref, comp, rtSd, mzSd, params, "RT")
    print ("    The 2nd round of RT-calibration is done")
    print ("      min SD of RT-shifts = %.4f second" % np.amin(rtSd))
    print ("      max SD of RT-shifts = %.4f second" % np.amax(rtSd))
    ref, comp, rtSd, mzSd = localCalibration(ref, comp, rtSd, mzSd, params, "mz")
    print ("    The 1st round of m/z-calibration is done")
    print ("      min SD of m/z-shifts = %.4f second" % np.amin(mzSd))
    print ("      max SD of m/z-shifts = %.4f second" % np.amax(mzSd))
    ref, comp, rtSd, mzSd = localCalibration(ref, comp, rtSd, mzSd, params, "mz")
    print ("    The 2nd round of m/z-calibration is done")
    print ("      min SD of m/z-shifts = %.4f second" % np.amin(mzSd))
    print ("      max SD of m/z-shifts = %.4f second" % np.amax(mzSd))
    print ()
    return comp, rtSd, mzSd # "comp" is the set of calibrated features


def globalCalibration(ref, comp, mzTol = 20):
    nPeaks = round(0.05 * ref.shape[0])    # Number of peaks to be considered for global calibration
    rtShifts = []   # Array for RT-shifts between reference and compared runs
    mzShifts = []   # Array for mz-shifts (ppm)
    i, j = 0, 1
    while j <= nPeaks:
        z = ref["z"][i] # From the 1st feature of reference run (i.e. strongest feature)
        mz = ref["mz"][i]
        lL = mz - mz * mzTol / 1e6
        uL = mz + mz * mzTol / 1e6
        rt = ref["RT"][i]
        intensity = ref["intensity"][i]
        if z == 0:
            # For a reference feature with undetermined charge, consider all possible charges in compared feature
            rowInd = np.where((comp["mz"] >= lL) & (comp["mz"] < uL))[0]
        else:
            rowInd = np.where((comp["mz"] >= lL) & (comp["mz"] < uL) & (comp['z'] == z))[0]

        if len(rowInd) > 0:
            rowInd = rowInd[0]
            rtShifts.append(comp["RT"][rowInd] - rt)
            mzShifts.append((comp["mz"][rowInd] - mz) / mz * 1e6)
            comp = np.delete(comp, rowInd, 0)
            j += 1

        i += 1

    # For more robust calculation, top and bottom 10% values are trimmed
    rtShifts = np.array(rtShifts)
    rtShifts = rtShifts[(rtShifts >= np.percentile(rtShifts, 10)) &
                     (rtShifts <= np.percentile(rtShifts, 90))]
    mzShifts = np.array(mzShifts)
    mzShifts = mzShifts[(mzShifts >= np.percentile(mzShifts, 10)) &
                     (mzShifts <= np.percentile(mzShifts, 90))]
    return rtShifts, mzShifts


def localCalibration(ref, comp, rtSd, mzSd, params, type):
    '''
    n = ref.shape[0]    # Number of rows in "ref"
    if not isinstance(rtSd, (list, np.ndarray)):
        rtSd = np.repeat(rtSd, n)
    if not isinstance(mzSd, (list, np.ndarray)):
        mzSd = np.repeat(mzSd, n)
    sdWidth = float(params['sd_width'])
    rtTol = rtSd * sdWidth
    mzTol = mzSd * sdWidth

    # Look for comparable features between "ref" and "comp"
    # and prepare LOESS modeling
    # refRt, refMz, compRt, compMz = [], [], [], []
    refInd, compInd = [], []
    for i in range(n):   # For each feature in "ref", look for a matching one in "comp"
        z = ref["z"][i]
        mz = ref["mz"][i]
        rt = ref["RT"][i]
        intensity = ref["intensity"][i]
        rtDev = comp["RT"] - rt
        mzDev = (comp["mz"] - mz) / mz * 1e6    # Unit of PPM
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
            if rowInd in compInd:
                continue
            else:
                refInd.append(i)
                compInd.append(rowInd)
    '''
    refInd, compInd = matchFeatures(ref, comp, rtSd, mzSd, params)

    # Perform LOESS regression to calibrated RT and m/z
    ro.r.source("./R/loess.R")
    rLoess = ro.globalenv["loess.as"]   # Note that control = loess.control(surface = "direct") should be set for prediction
    rPredict = ro.r("predict")
    refRt = ref["RT"][refInd]
    compRt = comp["RT"][compInd]
    refMz = ref["mz"][refInd]
    compMz = comp["mz"][compInd]

    if type == "RT":    # RT-calibration
        if (refRt == compRt).all():
            rtShifts = 1e-6 * np.random.normal(len(refRt))
        else:
            rtShifts = compRt - refRt
        mod = rLoess(FloatVector(compRt), FloatVector(rtShifts))
        compRt = compRt - np.array(mod.rx2("fitted"))   # Calibrated RT based on the model

        # Calculate a new (dynamic) RT-tolerance
        if (refRt == compRt).all():
            rtShifts = 1e-6 * np.random.normal(len(refRt))
        else:
            rtShifts = compRt - refRt   # Calibrated RT-shifts
        ind = np.where((rtShifts > np.percentile(rtShifts, 10)) &
                       (rtShifts < np.percentile(rtShifts, 90)))[0]
        modRtSd = rLoess(FloatVector(compRt[ind]), FloatVector(rtShifts[ind] ** 2))
        rtSd = np.sqrt(np.maximum(0, rPredict(modRtSd, FloatVector(ref["RT"]))))

        # Calculate a new (dynamic) m/z-tolerance
        if (refMz == compMz).all():
            mzShifts = 1e-6 * np.random.normal(len(refMz))
        else:
            mzShifts = (compMz - refMz) / refMz * 1e6
        # Sometimes, the variation of mzShifts cannot be captured when trimming is applied
        # so, the trimming is not used for mzShifts
        modMzSd = rLoess(FloatVector(compMz), FloatVector(mzShifts ** 2), 1, "aicc", "gaussian")
        mzSd = np.sqrt(np.maximum(0, rPredict(modMzSd, FloatVector(ref["mz"]))))

        # Calibration of the entire comp["RT"]
        comp["RT"] = comp["RT"] - rPredict(mod, FloatVector(comp["RT"]))
    elif type == "mz":  # m/z-calibration
        if (refMz == compMz).all():
            mzShifts = 1e-6 * np.random.normal(len(refMz))
        else:
            mzShifts = (compMz - refMz) / refMz * 1e6
        mod = rLoess(FloatVector(compMz), FloatVector(mzShifts), 1, "aicc", "gaussian")
        compMz = compMz * (1 + np.array(mod.rx2("fitted")) / 1e6)   # Calibrated m/z basd on the model

        # Calculate a new (dynamic) m/z-tolerance
        if (refMz == compMz).all():
            mzShifts = 1e-6 * np.random.normal(len(refMz))
        else:
            mzShifts = (compMz - refMz) / refMz * 1e6   # Calibrated m/z-shifts
        modMzSd = rLoess(FloatVector(compMz), FloatVector(mzShifts ** 2), 1, "aicc", "gaussian")
        mzSd = np.sqrt(np.maximum(0, rPredict(modMzSd, FloatVector(ref["mz"]))))

        # Calculate a new (dynamic) RT-tolerance
        if (refRt == compRt).all():
            rtShifts = 1e-6 * np.random.normal(len(refRt))
        else:
            rtShifts = compRt - refRt
        ind = np.where((rtShifts > np.percentile(rtShifts, 10)) &
                       (rtShifts < np.percentile(rtShifts, 90)))[0]
        modRtSd = rLoess(FloatVector(compRt[ind]), FloatVector(rtShifts[ind] ** 2), 1, "aicc", "gaussian")
        rtSd = np.sqrt(np.maximum(0, rPredict(modRtSd, FloatVector(ref["RT"]))))

        # Calibration of the entire comp["mz"]
        comp["mz"] = comp["mz"] / (1 + np.array(rPredict(mod, FloatVector(comp["mz"]))) / 1e6)

    return ref, comp, rtSd, mzSd


def matchFeatures(ref, comp, rtSd, mzSd, params):
    ref = np.sort(ref, order = "intensity")[::-1]   # Sort features in descending order of intensity
    comp = np.sort(comp, order = "intensity")[::-1]

    n = ref.shape[0]
    if not isinstance(rtSd, (list, np.ndarray)):
        rtSd = np.repeat(rtSd, n)
    if not isinstance(mzSd, (list, np.ndarray)):
        mzSd = np.repeat(mzSd, n)
    rtSd[rtSd == 0] = min(rtSd[rtSd > 0])
    mzSd[mzSd == 0] = min(mzSd[mzSd > 0])
    sdWidth = float(params['sd_width'])
    rtTol = rtSd * sdWidth
    mzTol = mzSd * sdWidth

    # Look for matching features between "ref" and "comp"
    refInd, compInd = [], []
    for i in range(n):   # For each feature in "ref", look for a matching one in "comp"
        z = ref["z"][i]
        mz = ref["mz"][i]
        rt = ref["RT"][i]
        intensity = ref["intensity"][i]
        rtDev = comp["RT"] - rt
        mzDev = (comp["mz"] - mz) / mz * 1e6    # Unit of PPM
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
            if rowInd in compInd:
                continue
            else:
                refInd.append(i)
                compInd.append(rowInd)

    return refInd, compInd


###########################################
################ Main part ################
###########################################

paramFile = r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\jumpm_negative.params"
featureFiles = [r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.1.dev.feature",
                r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.1.dev.feature",
                r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.1.dev.feature"]

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
for i in range(nFiles):
    data = np.genfromtxt(featureFiles[i], delimiter = "\t", dtype = None, names = True)
    fArray.append(data)

if nFiles > 1: # Multiple feature files -> alignment is required
    print ("  Feature calibration")
    print ("  ===================")

    ###################################
    # Selection of a reference sample #
    ###################################
    if params["reference_feature"] == "0":
        # A sample with the largest median of top 100 intensities is set to a reference run
        refNo = 0
        refIntensity = 0
        for i in range(nFiles):
            tmpIntensity = np.median(sorted(fArray[i]["intensity"], reverse = True)[0: 100])
            if tmpIntensity >= refIntensity:
                refNo = i
                refIntensity = tmpIntensity
    else:
        try:
            refNo = featureFiles.index(params["reference_feature"])
        except:
            sys.exit("  'reference_feature' parameter should be correctly specified")
    print ("  %s is chosen as the reference run" % os.path.basename(featureFiles[refNo]))

    ############################################################
    # Calibration of features against those in a reference run #
    ############################################################
    rtSdArray, mzSdArray = [], []
    for i in range(nFiles):
        if i != refNo:
            print ("  " + os.path.basename(featureFiles[i]) + " is being aligned against the reference run (it may take a while)")
            fArray[i], rtSd, mzSd = calibrateFeatures(fArray[refNo], fArray[i], params)
            rtSdArray.append(rtSd)
            mzSdArray.append(mzSd)
        else:
            rtSdArray.append("NA")
            mzSdArray.append("NA")

    print ("  Calibration summary")
    print ("  After calibration, RT- and m/z-shifts of each run (against the reference run) are centered to zero")
    print ("  Variations (i.e. standard deviation) of RT- and m/z-shifts are as follows,")
    print ("  Filename\t\t\t#features\tSD of RT-shifts [second]\tSD of m/z-shifts [ppm]")
    for i in range(nFiles):
        featureName = os.path.basename(featureFiles[i])
        nFeatures = str(fArray[i].shape[0])
        if i != refNo:
            meanRtSd = "%.6f" % np.mean(rtSdArray[i])
            meanMzSd = "%.6f" % np.mean(mzSdArray[i])
        else:
            meanRtSd = "NA"
            meanMzSd = "NA"
        print ("  " + featureName + "\t\t\t" + nFeatures + "\t" + meanRtSd + "\t" + meanMzSd)
    print ()

    #################################################################
    # Identification of fully-aligned features for further analysis #
    #################################################################
    print ("  Feature alignment")
    print ("  =================")
    for i in range(nFiles):
        nFeatures = str(fArray[i].shape[0])
        if i != refNo:
            refName = os.path.basename(featureFiles[refNo])
            compName = os.path.basename(featureFiles[i])
            print ("  " + refName + ": %d features (reference run)" % fArray[refNo].shape[0])
            print ("  " + compName + ": %d features (compared run)" % fArray[i].shape[0])
            # matchFeatures(fArray[refNo], fArray[i], rtSdArray[i], mzSdArray[i], params)







