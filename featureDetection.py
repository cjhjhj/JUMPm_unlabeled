#!/usr/bin/python

import sys, re, numpy as np
from typing import List, Any

from pyteomics import mzxml


def getParams(paramFile):
    parameters = dict()
    with open(paramFile, 'r') as file:
        for line in file:
            if re.search(r'^#', line) or re.search(r'^\s', line):
                continue
            line = re.sub(r'#.*', '', line)  # Remove comments (start from '#')
            line = re.sub(r'\s*', '', line)  # Remove all whitespaces
            key = line.split('=')[0]
            val = line.split('=')[1]
            parameters[key] = val
    return parameters


def detectPeaks(spec, params):
    # Parameters
    if params['data_acquisition_mode'] == "1":
        isCentroided = 1
    elif params['data_acquisition_mode'] == "2":
        isCentroided = 0
    else:
        print("Please set the proper 'data_acquisition_mode' parameter")
        sys.exit("")
    intensityThreshold = 0  # May come from a parameter file

    # m/z and intensity arrays from a spectrum object
    mzArray = spec['m/z array']
    intensityArray = spec['intensity array']
    nPeaks = len(mzArray)
    newMzArray = []
    newIntensityArray = []

    # Detect peaks (i.e. centroidization of MS1 spectrum)
    if isCentroided == 0:  # i.e. Profile mode MS1
        for i in range(2, nPeaks - 1):
            if intensityArray[i] > 0:
                # Consider 2 points before and after the point of interest x, i.e. 5 point window
                b2, b1, x, a1, a2 = intensityArray[(i - 2):(i + 3)]
                if x >= intensityThreshold:
                    if isMax(b2, b1, x, a1, a2):
                        # If x is the local maximum in a 5-point window, lower and upper bounds for a peak will be explored
                        # Refer Figure 1a and b in the paper, Cox and Mann, Nature Biotech. 2008; 26: 1367-22
                        minInd = findMinPeakIndex(i, intensityArray)
                        maxInd = findMaxPeakIndex(i, intensityArray)
                        if (maxInd - minInd) > 2:
                            newMz, newIntensity = findPeakCenter(minInd, i, maxInd, mzArray, intensityArray)
                            newMzArray.append(newMz)
                            newIntensityArray.append(newIntensity)

        # Update 'spec' object
        spec['m/z array'] = newMzArray
        spec['intensity array'] = newIntensityArray

    # Do nothing for centroid mode MS1
    return spec


def isMax(b2, b1, x, a1, a2):
    if x > b1 and x > a1:
        return True
    if x > b2 and x == b1 and x > a1:
        return True
    if x > b1 and x == a1 and x > a2:
        return True
    return False


def findMinPeakIndex(ind, array):
    while ind > 0 and array[ind] != 0 and array[ind - 1] <= array[ind]:
        ind -= 1
    return ind + 1


def findMaxPeakIndex(ind, array):
    count = len(array)
    while ind < count and array[ind] != 0 and array[ind + 1] <= array[ind]:
        ind += 1
    return ind - 1


def findPeakCenter(minInd, centerInd, maxInd, mz, intensity):
    # Find the center of a peak composed of five data points
    centerMz = 0
    centerIntensity = 0
    for i in range(minInd, maxInd + 1):
        if intensity[i] >= centerIntensity:
            centerIntensity = intensity[i]  # Take the maximum intensity as a peak intensity

    # There's a plateau, bu others are zeros
    if minInd == maxInd:
        centerMz = mz[maxInd]
        return centerMz, centerIntensity

    # Right-angled triangle-shaped peak
    if minInd == centerInd:
        centerMz = estimate2(mz[centerInd], mz[centerInd + 1], intensity[centerInd], intensity[centerInd + 1])
        return centerMz, centerIntensity

    # Left-angled triangle-shaped peak
    if maxInd == centerInd:
        centerMz = estimate2(mz[centerInd - 1], mz[centerInd], intensity[centerInd - 1], intensity[centerInd])
        return centerMz, centerIntensity

    # Typical bell(triangle)-shaped peak
    centerMz = estimate3(mz[centerInd - 1], mz[centerInd], mz[centerInd + 1], intensity[centerInd - 1], intensity[centerInd], intensity[centerInd + 1])
    return centerMz, centerIntensity


def estimate2(m1, m2, i1, i2):
    centerVal = (m1 * i1 + m2 * i2) / (i1 + i2) # Intensity-weighted average of m/z
    return centerVal


def estimate3(m1, m2, m3, i1, i2, i3):
    l1 = np.log(i1)
    l2 = np.log(i2)
    l3 = np.log(i3)
    centerVal = 0.5 * ((l1 - l2) * (m3 ** 2 - m1 ** 2) - (l1 - l3) * (m2 ** 2 - m1 ** 2)) / ((l1 - l2) * (m3 - m1) - (l1 - l3) * (m2 - m1));
    return centerVal


def findPeakMatch(array, value, tolerance):
    mzArray = np.asarray(array['m/z array'])
    ind = np.abs(mzArray - value).argmin()    # Index of the closest element to "value" in the array
    difference = np.abs(mzArray[ind] - value) / value * 1e6
    if difference <= tolerance:
        return True, ind
    else:
        return False, 0


def reduceMS1(spec, array):
    # Input
    # spec: spectrum object read by pyteomics
    # array: index of 'm/z array' (and 'intensity array') to be retained
    redMzArray = [spec['m/z array'][i] for i in array]
    redIntensityArray = [spec['intensity array'][i] for i in array]

    # Replace m/z array and intensity array of spec with reduced ones
    spec['m/z array'] = redMzArray
    spec['intensity array'] = redIntensityArray

    # Extract intensity information of "noisy peaks"
    noiseIndex = np.setdiff1d(range(0, len(spec['intensity array'])), array)
    noiseIntensityArray = [spec['intensity array'][i] for i in noiseIndex]
    noiseIntensityLevel = np.percentile(noiseIntensityArray, 75)    # Third quartile = median of top 50%

    # Return "reduced" spectrum and noise intensity level of the input spectrum (scan)
    return spec, noiseIntensityLevel

# Input: mzXML file

# To-do: expand to mzML

inputFile = "/Research/Projects/7Metabolomics/JUMPm/IROAsamples/IROA_IS_NEG_1.mzXML"
reader = mzxml.read(inputFile)

##############
# Parameters #
##############
paramFile = "jumpm_negative_desktop.params"
params = getParams(paramFile)
# firstScan = int(params['first_scan_extraction'])
# lastScan = int(params['last_scan_extraction'])
firstScan = 1000
lastScan = 1200
gap = int(params['skipping_scans'])
scanWindow = gap + 1
matchTolerance = float(params['mass_tolerance_peak_matching'])

noiseInfo = {}  # Empty dictionary for noise level information
features = []
nFeatures = -1

with reader:
    ############################
    # Get MS1 scan information #
    ############################
    nMS1 = 0
    # ms1ScanNumArray: List[Any] = []
    ms1ScanNumArray = []
    for spec in reader:
        msLevel = int(spec['msLevel'])
        scanNum = int(spec['num'])
        if msLevel == 1 and firstScan <= scanNum <= lastScan:
            nMS1 += 1
            ms1ScanNumArray.append(spec['num'])
        elif scanNum > lastScan:
            break

    ##########################
    # 3D-peak identification #
    ##########################
    cacheArray = []
    oldMinInd = -1
    oldMaxInd = -1
    for i in range(0, nMS1):
        ############################################################
        # Put former and latter scans to a cache (temporary) array #
        # Current scan index = i                                   #
        # Former and latter scan indexes are from the variables of #
        # minInd, maxInd, oldMinInd and oldMaxInd                  #
        ############################################################
        minInd = max(0, i - gap - 1)
        maxInd = min(nMS1 - 1, i + gap + 1)
        if i == 0:
            for j in range(0, maxInd + 1):
                spec = reader[ms1ScanNumArray[j]]
                spec = detectPeaks(spec, params)
                cacheArray.append(spec)
        else:
            for j in range(oldMinInd, minInd):
                cacheArray.pop(0)   # Equivalent to "shift" in Perl, remove the first element from an array
            for j in range(oldMinInd + 1, maxInd + 1):
                spec = reader[ms1ScanNumArray[j]]
                spec = detectPeaks(spec, params)
                cacheArray.append(spec)

        ##################################################################################################
        # Reduction step                                                                                 #
        # For each scan, retain MS1 peaks which can form a 3D-peak with other peaks in neighboring scans #
        # Reduction of MS1 peaks is performed in forward/backward direction                              #
        ##################################################################################################
        p = cacheArray[i - minInd]  # MS1 spectrum under consideration
        retainIndexArray = []
        for j in range(0, len(p['m/z array'])):
            # For each j-th peak in "p" (spec object), look for matched peaks in "q" (spec object of neighboring scans)
            cm = p['m/z array'][j]  # m/z of j-th peak
            match = 0
            nTries = 0
            # Backward search
            for k in range(i - 1, minInd + 1, -1):  # Decreasing index, k
                q = cacheArray[k - minInd]
                match, ind = findPeakMatch(q, cm, matchTolerance)
                if match:
                    break
                nTries += 1
                if nTries > scanWindow:
                    break

            # Forward search
            if not match:
                for k in range(i + 1, maxInd + 1):
                    q = cacheArray[k - minInd]
                    match, ind = findPeakMatch(q, cm, matchTolerance)
                    if match:
                        break
                    nTries += 1
                    if nTries > scanWindow:
                        break

            if match:
                retainIndexArray.append(j)

        # Reduction of peaks in "p" and extract "noisy spectrum/scan" information
        p, noiseLevel = reduceMS1(p, retainIndexArray)
        noiseInfo[p['num']] = noiseLevel

        #####################################################################################################
        # Generation of 3D-peaks by combining peaks within the specified mass tolerance in the former scans #
        # 3D-peak generation is performed in backward direction                                             #
        #####################################################################################################
        cacheArray[i - minInd] = p  # Replace cacheArray with the reduced spectrum
        for j in range(0, len(p['m/z array'])):
            cm = p['m/z array'][j]
            match = 0
            nTries = 0
            qIndexArray = []
            for k in range(i - 1, minInd + 1, -1):  # Backward search to merge peaks
                q = cacheArray[k - minInd]
                match, ind = findPeakMatch(q, cm, matchTolerance)
                if match:
                    qIndexArray.append(ind)

            if match:   # i.e. If there's at least one peak that can be merged with p['m/z array'][j] in q
                qIndexArray = np.unique(qIndexArray)
            else:
                # There's no matchable peak (with the current peak, p['m/z array'][j]) in q
                # In this case, the current peak is added to a new feature
                if i < nMS1:
                    nFeatures += 1
                    cacheArray[i - minInd]['featureIndex'].append(nFeatures)
                    features[nFeatures]['mz'].append(p['m/z array'][j])
                    features[nFeatures]['intensity'].append(p['intensity array'][j])
                    features[nFeatures]['scanNum'].append(p['num'])
                    features[nFeatures]['rt'].append(p['retentionTime'])

        oldMinInd = minInd
        oldMaxInd = maxInd

    # Re-numbering of 3D-peaks (some 3D-peaks are "singleton")


















