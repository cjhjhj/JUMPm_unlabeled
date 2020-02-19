#!/usr/bin/python

import numpy as np
import re
import sys

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
    centerMz = estimate3(mz[centerInd - 1], mz[centerInd], mz[centerInd + 1], intensity[centerInd - 1],
                         intensity[centerInd], intensity[centerInd + 1])
    return centerMz, centerIntensity


def estimate2(m1, m2, i1, i2):
    centerVal = (m1 * i1 + m2 * i2) / (i1 + i2)  # Intensity-weighted average of m/z
    return centerVal


def estimate3(m1, m2, m3, i1, i2, i3):
    l1 = np.log(i1)
    l2 = np.log(i2)
    l3 = np.log(i3)
    centerVal = 0.5 * ((l1 - l2) * (m3 ** 2 - m1 ** 2) - (l1 - l3) * (m2 ** 2 - m1 ** 2)) / (
                (l1 - l2) * (m3 - m1) - (l1 - l3) * (m2 - m1))
    return centerVal


def findPeakMatch(array, value, tolerance):
    mzArray = np.asarray(array['m/z array'])
    ind = np.abs(mzArray - value).argmin()  # Index of the closest element to "value" in the array
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
    noiseIntensityLevel = np.percentile(noiseIntensityArray, 75)  # Third quartile = median of top 50%

    # Return "reduced" spectrum and noise intensity level of the input spectrum (scan)
    return spec, noiseIntensityLevel


def getScans(reader, array, scanNums, currIndex, lastIndex, windowSize, params):
    # When currIndex = 0 (i.e. first scan)
    # array = [currIndex-th scan,
    #          (currIndex + 1)-th scan,
    #          ...
    #          (currIndex + windowSize)-th scan

    # When 0 < currIndex < lastIndex (i.e. between first and last scan)
    # 1. Remove the first element from array
    # 2. Centroid (detectPeak) (currIndex + windowSize)-th scan
    # 3. Append the "centroided" (currIndex + windowSize)-th scan to the array

    # Example
    # currIndex = 0, windowSize = 3 (i.e. gap = 2, allow "skip" up to 2 scans), then
    # array = [centroided 0-th scan,
    #          centroided 1-th scan,
    #          centroided 2-th scan,
    #          centroided 3-th scan]
    # currIndex = 1, then,
    # 1. Remove the first element from array
    #    array = [centroided 1-th scan,
    #             centroided 2-th scan,
    #             centroided 3-th scan]
    # 2. Centroid (currIndex + windowSize)-th scan, i.e. 4-th scan
    # 3. Append the above centroided 4-th scan
    #    array = [centroided 1-th scan,
    #             centroided 2-th scan,
    #             centroided 3-th scan,
    #             centroided 4-th scan]

    if len(array) == 0:  # when currentIndex = 0
        array = []
        for i in range(currIndex, currIndex + windowSize + 1):
            spec = reader[scanNums[i]]
            spec = detectPeaks(spec, params)
            array.append(spec)
    elif (currIndex + windowSize) < lastIndex:
        array.pop(0)
        spec = reader[scanNums[currIndex + windowSize]]
        spec = detectPeaks(spec, params)
        array.append(spec)
    else:
        array.pop(0)

    return array


def insertPeak(mz, intensity, f, fMz, specs, params):
    # Input arguments
    # mz = m/z of the query peak
    # intensity = intensity of the query peak
    # specs = array of spectra, 0-th element is the query spectrum
    # f = array of features
    # fMz = array of features' m/z values (center m/z)
    # params = dictionary of parameters
    tol = float(params['mass_tolerance_peak_matching'])
    scanNumber = int(specs[0]['num'])  # MS1 scan number of the query spectrum
    rt = float(specs[0]['retentionTime'])  # RT of the query spectrum

    if len(f) > 0:  # If any feature exists
        # Check whether an existing feature can be extended
        doInsert = False
        fInd = np.where((abs(fMz - mz) / mz * 1e6) <= tol)[0]

        if len(fInd) == 0:  # A new feature will be created for the query peak
            # Check next scan(s) whether the query peak can form a feature
            for i in range(1, len(specs)):
                mzArray = np.array(specs[i]['m/z array'])
                mzInd = np.argmin(abs(mzArray - mz))
                diff = abs(mz - mzArray[mzInd]) / mz * 1e6
                if diff < tol:
                    doInsert = True
                    break
        else:  # An existing feature will grow
            if len(fInd) > 1:
                # When the query peak is matched to multiple features,
                # the feature having the lowest index (representative feature) will grow
                # Other matched features will be merged to the representative feature
                for i in range(1, len(fInd)):
                    f[fInd[0]]['mz'].extend(f[fInd[i]]['mz'])
                    f[fInd[0]]['intensity'].extend(f[fInd[i]]['intensity'])
                    f[fInd[0]]['scanNumber'].extend(f[fInd[i]]['scanNumber'])
                    f[fInd[0]]['rt'].extend(f[fInd[i]]['rt'])

                # Then, remove "merged" features (list comprehension looks slow)
                for i in sorted(fInd[1:], reverse = True):
                    del f[i]
                    del fMz[i]

                # f['centerMz'] = [x for i, x in enumerate(f['centerMz']) if i not in fInd[1:]]
                # f['mz'] = [x for i, x in enumerate(f['mz']) if i not in fInd[1:]]
                # f['intensity'] = [x for i, x in enumerate(f['intensity']) if i not in fInd[1:]]
                # f['scanNumber'] = [x for i, x in enumerate(f['scanNumber']) if i not in fInd[1:]]
                # f['rt'] = [x for i, x in enumerate(f['rt']) if i not in fInd[1:]]

            # Update the feature
            fInd = fInd[0]
            f[fInd]['mz'].append(mz)
            f[fInd]['intensity'].append(intensity)
            f[fInd]['scanNumber'].append(scanNumber)
            f[fInd]['rt'].append(rt)
            centerMz = np.dot(np.array(f[fInd]['mz']), np.array(f[fInd]['intensity'])) / np.sum(
                np.array(f[fInd]['intensity']))
            fMz[fInd] = centerMz
    else:
        # When there's no feature yet, a new feature will be created with the following information
        doInsert = True

    if doInsert is True:
        if len(f) > 0:
            # New feature is appended to the feature array
            ind = np.argmin(abs(fMz - mz))
            if mz > fMz[ind]:
                ind += 1
            f.insert(ind, {'mz': [mz], 'intensity': [intensity], 'scanNumber': [scanNumber], 'rt': [rt]})
            fMz.insert(ind, mz)
        else:
            # Initial feature array is created (it should be executed only once)
            f.append({'mz': [mz], 'intensity': [intensity], 'scanNumber': [scanNumber], 'rt': [rt]})
            fMz.append(mz)

    return f, fMz, specs

'''
class progressBar():
    def __init__(self, nTot):
        # Initialization
        self.nTot = nTot
        self.minorTicks = 5
        self.majorTicks = 25
        self.cycleInterval = 50
        self.cycleCharArray = ["-", "/", '\\', '|']
        self.cycleI = 0
        self.curPercent = 0
        self.counter = 0
        self.backupCursor = 0

    def increment(self):
        self.counter += 1
        if self.counter % self.cycleInterval == 0:
            if self.backupCursor == 1:
                print("\b", end='')
            print(self.cycleCharArray[self.cycleI], end='')
            self.cycleI = (self.cycleI + 1) % len(self.cycleCharArray)
            self.backupCursor = 1
            sys.stdout.flush()
        if (self.counter / self.nTot * 100) >= self.curPercent:
            sys.stdout.flush()
            if self.curPercent % self.majorTicks == 0:
                if self.backupCursor == 1:
                    print("\b", end='')
                print("%d%%" % self.curPercent, end='')
                self.backupCursor = 0
                if self.curPercent == 100:
                    print()
            else:
                if self.backupCursor == 1:
                    print("\b", end='')
                print(".", end='')
            self.curPercent += self.minorTicks
'''

class progressBar:
    def __init__(self, total):
        self.total = total
        self.barLength = 20
        self.count = 0
        self.progress = 0
        self.block = 0
        self.status = ""

    def increment(self):
        self.count += 1
        self.progress = self.count / self.total
        self.block = int(round(self.barLength * self.progress))
        if self.progress == 1:
            self.status = "Done...\r\n"
        else:
            self.status = ""
        #         self.status = str(self.count) + "/" + str(self.total)
        text = "\rProgress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block),
                                                   int(self.progress * 100), self.status)
        sys.stdout.write(text)
        sys.stdout.flush()


######################
##### Main part ######
######################
# Input: mzXML file

# To-do: expand to mzML

inputFile = r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.mzXML"
reader = mzxml.read(inputFile)
# reader = mzxml.read("U:/Research/Projects/7Metabolomics/JUMPm/IROAsamples/IROA_IS_NEG_1.mzXML")

##############
# Parameters #
##############
paramFile = "jumpm_negative_desktop.params"
params = getParams(paramFile)
firstScan = int(params['first_scan_extraction'])
lastScan = int(params['last_scan_extraction'])
gap = int(params['skipping_scans'])
# firstScan = 1000
# lastScan = 1200
# gap = 1
scanWindow = gap + 1
matchTolerance = float(params['mass_tolerance_peak_matching'])

##################
# Initialization #
##################
noiseInfo = {}  # Empty dictionary for noise level information
features = []
ms1ToFeatures = {}  # It is necessary for the labeled approach
featuresCenterMzArray = []
nFeatures = -1

with reader:
    ############################
    # Get MS1 scan information #
    ############################
    nScans = 0
    # ms1ScanNumArray: List[Any] = []
    scanNumArray = []
    print ("Extraction of MS1 spectra from input file(s)")
    for spec in reader:
        msLevel = int(spec['msLevel'])
        scanNum = int(spec['num'])
        if msLevel == 1 and firstScan <= scanNum <= lastScan:
            nScans += 1
            scanNumArray.append(spec['num'])  # Scan numbers of MS1 scans
        elif scanNum > lastScan:
            break

    print ("Done\n")

    ################################
    # Feature (3D-peak) generation #
    ################################
    print ("Feature detection ")
    progress = progressBar(nScans)
    specArray = []
    for i in range(0, nScans):
        progress.increment()
        # Make an array containing spectra; specArray
        # specArray[0] = current MS1 (i.e i-th MS1 scan)
        # specArray[1] = (i + 1)-th MS1 scan
        # ...
        # specArray[scanWindow] = (i + scanWindow)-th MS1 scan
        specArray = getScans(reader, specArray, scanNumArray, i, nScans, scanWindow, params)

        # Insert peaks
        # 3D-peak generation step by inserting peaks in neighboring scan(s)
        for j in range(0, len(specArray[0]['m/z array'])):
            mz_j = specArray[0]['m/z array'][j]
            intensity_j = specArray[0]['intensity array'][j]
            features, featuresCenterMzArray, specArray = insertPeak(mz_j, intensity_j, features, featuresCenterMzArray,
                                                                    specArray, params)

    for i in range(0, len(features)):
        features[i]['centerMz'] = featuresCenterMzArray[i]
        for j in range(0, len(features[i]['scanNumber'])):
            scanNum = features[i]['scanNumber'][j]
            if scanNum not in ms1ToFeatures:
                ms1ToFeatures[scanNum] = {'mz': [features[i]['mz'][j]],
                                          'intensity':[features[i]['intensity'][j]]}
            else:
                ms1ToFeatures[scanNum]['mz'].append(features[i]['mz'][j])
                ms1ToFeatures[scanNum]['intensity'].append(features[i]['intensity'][j])




    print ("Finished the feature detection")