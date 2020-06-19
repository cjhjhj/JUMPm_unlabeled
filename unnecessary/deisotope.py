#!/usr/bin/python

import sys, os, re, numpy as np
from pyteomics import mzxml

def findCharge(spec, tol, queryMz, queryIntensity, maxCharge = 6):
    delC = 1.00335 # Difference of mass between C12 and C13
    lL = queryMz - delC - queryMz * tol / 1e6  # Lower limit
    uL = queryMz + delC + queryMz * tol / 1e6  # Upper limit
    subSpec = spec[(spec[:, 0] >= lL) & (spec[:, 0] <= uL)]
    if subSpec.shape[0] > 0:
        strongestMz = queryMz
        strongestIntensity = queryIntensity
        for mz, intensity in subSpec:
            if intensity >= strongestIntensity:
                continue
            charge = np.around(1 / abs(mz - strongestMz))
            if (charge > maxCharge) or (charge == 0):
                continue
            err = abs(abs(mz - strongestMz) - (delC / charge))
            if err > tol * strongestMz / 1e6:
                continue





'''
def findCharge(spec, tol, queryMz, queryIntensity, maxCharge = 6):
    # Parameter/constants
    delC = 1.003355 # Mass difference between C12 and C13

    # Output variables
    charge = 0
    intRatio = 0
    isotopeMz = 0
    peakC = ""

    lL = queryMz - delC - queryMz * tol / 1e6 # Lower limit
    uL = queryMz + delC + queryMz * tol / 1e6 # Upper limit
    subSpec = spec[(spec[:, 0] >= lL) & (spec[:, 0] <= uL)]
    if subSpec.shape[0] > 0:
        strongestMz = queryMz
        strongestIntensity = queryIntensity
        for mz, intensity in subSpec:
            if intensity >= strongestIntensity:

                # To-do
                # As of now, nothing would happen if a potential isotopic peak (M+1) is higher than the current peak (M)
                # However, it may happen (?) -> Check DOPEGAL feature detection part

                continue
            val = np.around(1 / abs(mz - strongestMz))
            if (val > maxCharge) or (val == 0):
                continue
            err = abs(abs(mz - strongestMz) - (delC / val))
            if err > tol * strongestMz / 1e6:
                continue
            intRatio = intensity / strongestIntensity * 100
            charge = val
            isotopeMz = mz
            if mz < strongestMz:
                peakC = "C13"
            else:
                peakC = "C12"
            break
    return (charge, intRatio, isotopeMz, peakC)               


def deisotopeMS1(inputFile):
    spec = np.loadtxt(inputFile, skiprows = 1)
    spec = spec[(-spec[:, 1]).argsort()]    # Sort the array in descending order of intensity

    # Parameters and/or constants
    tol = 10 # Decharge ppm
    maxCharge = 6 # Maximum charge state allowed

    # Output variables
    specDict = {}
    intRatioDict = {}
    cPeakDict = {}
    chargeDict = {}

    for mz, intensity in spec:
        specDict[mz] = intensity
        charge, intRatio, isotopeMz, cPeak = findCharge(spec, tol, mz, intensity, maxCharge)
        intRatioDict[mz] = intRatio
        if charge == 0:
            continue
        else:
            chargeDict[mz] = charge
            chargeDict[isotopeMz] = charge
            cPeakDict[mz] = cPeak

    # The following step is a decharging of peaks ("changeMH" function in Perl script)
    # It is used to decharge multiply-charged peaks to singly-charged peaks
    # The information of singly-charged peaks is stored in /peaks directory, and
    # old peaks (i.e. mulitply-charged) are removed from "specDict"
    for mz, _ in spec:
        if mz in chargeDict and chargeDict[mz] > 1:
            newMz = (mz - H) * chargeDict[mz] + H # Singly-charged m/z = neutral mass + H
            if newMz > 1500:
                del specDict[mz], intRatioDict[mz], cPeakDict[mz], chargeDict[mz]
            specDict[newMz] = specDict[mz]
            print ("%f deleted" % mz)
            del specDict[mz]

    return (specDict, intRatioDict, cPeakDict, chargeDict)
'''

# Read a dta-format file containing MS1 spectra (*.MS1 file)
inputFile = "/Research/Projects/7Metabolomics/JUMPm_v1.8.0_XWang/Yeast_4Plex_RP_Pos_1/Yeast_4Plex_RP_Pos_1.1/13002.MS1"

# 1. Deisotope of MS1 peaks
# If there's no feature with charge state > 1, "specDict" is the same as inputFile
specDict, intRatioDict, cPeakDict, chargeDict = deisotopeMS1(inputFile)
print("%d features are detected" % len(specDict))

# 2. Clustering of isotopically labeled peak pair candidates
spec = np.loadtxt(inputFile, skiprows=1)
spec = spec[(-spec[:, 1]).argsort()]  # Sort the array in descending order of intensity

# Clustering
clusterDict = {}
for queryMz, queryIntensity in spec:
    # Initialization of "clusterDict" at "mz"
    clusterDict[queryMz] = {}
    clusterDict[queryMz][0] = {}
    clusterDict[queryMz][0]['mz'] = queryMz
    clusterDict[queryMz][0]['int'] = queryIntensity

    # Constants
    tolPpm = 10  # Clustering ppm
    delC = 1.003355

    # Variables
    previousIntensity = queryIntensity
    selectedIntensity = queryIntensity

    # Backward search
    searchLoop = 1
    flag = 1
    while (searchLoop and flag):
        lL = queryMz - delC * searchLoop - tolPpm * queryMz / 1e6
        uL = queryMz - delC * searchLoop + tolPpm * queryMz / 1e6
        flag = 0
        subSpec = spec[(spec[:, 0] >= lL) & (spec[:, 0] <= uL)]
        if np.size(subSpec, 0) > 0:
            flag = 1
            for mz, intensity in subSpec:
                if intensity > previousIntensity:  # C12
                    clusterDict[queryMz][-searchLoop] = {}
                    clusterDict[queryMz][-searchLoop]['mz'] = queryMz
                    if clusterDict[queryMz][-searchLoop].get('int') is not None:
                        clusterDict[queryMz][-searchLoop]['int'] += intensity
                    else:
                        clusterDict[queryMz][-searchLoop]['int'] = intensity
                    searchLoop += 1
                    ind = np.where(spec[:, 0] == queryMz)
                    spec = np.delete(spec, ind, 0)
                elif intensity < previousIntensity:  # C13
                    clusterDict[queryMz][-searchLoop] = {}
                    clusterDict[queryMz][-searchLoop]['mz'] = queryMz
                    if clusterDict[queryMz][-searchLoop].get('int') is not None:
                        clusterDict[queryMz][-searchLoop]['int'] += intensity
                    else:
                        clusterDict[queryMz][-searchLoop]['int'] = intensity
                    searchLoop += 1
                    ind = np.where(spec[:, 0] == mz)
                    spec = np.delete(spec, ind, 0)
                else:
                    searchLoop = 0

    # Forward search
    searchLoop = 1
    flag = 1
    while (searchLoop and flag):
        lL = queryMz + delC * searchLoop - tolPpm * queryMz / 1e6
        uL = queryMz + delC * searchLoop + tolPpm * queryMz / 1e6
        flag = 0
        subSpec = spec[(spec[:, 0] >= lL) & (spec[:, 0] <= uL)]
        if np.size(subSpec, 0) > 0:
            flag = 1
            for mz, intensity in subSpec:
                if intensity > previousIntensity:  # C12
                    clusterDict[queryMz][searchLoop] = {}
                    clusterDict[queryMz][searchLoop]['mz'] = queryMz
                    if clusterDict[queryMz][searchLoop].get('int') is not None:
                        clusterDict[queryMz][searchLoop]['int'] += intensity
                    else:
                        clusterDict[queryMz][searchLoop]['int'] = intensity
                    searchLoop += 1
                    ind = np.where(spec[:, 0] == mz)
                    spec = np.delete(spec, ind, 0)
                elif intensity < previousIntensity:  # C13
                    clusterDict[queryMz][searchLoop] = {}
                    clusterDict[queryMz][searchLoop]['mz'] = queryMz
                    if clusterDict[queryMz][searchLoop].get('int') is not None:
                        clusterDict[queryMz][searchLoop]['int'] += intensity
                    else:
                        clusterDict[queryMz][searchLoop]['int'] = intensity
                    searchLoop += 1
                    ind = np.where(spec[:, 0] == queryMz)
                    spec = np.delete(spec, ind, 0)
                else:
                    searchLoop = 0

print(clusterDict)
