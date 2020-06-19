#!/usr/bin/python

import os, sys, re, numpy as np, pandas as pd, pickle
from featureDetection import detectFeatures
from featureAlignment import alignFeatures
from featureToMS2 import ms2ForFeatures
from librarySearch import searchLibrary
from pyteomics import mzxml

inputFiles = [r"/Research/Projects/7Metabolomics/htan_IROA/2020/hilic/IROA_neg_target/IROA_neg_target.mzXML"]
paramFile = r"/Research/Projects/7Metabolomics/htan_IROA/2020/hilic/align_test/jumpm_negative.params"
nFiles = len(inputFiles)

#####################
# Feature detection #
#####################
featureArray = []
for i in range(nFiles):
    f = detectFeatures(inputFiles[i], paramFile)
    featureArray.append(f)

#####################
# Feature alignment #
#####################
featureFiles = [os.path.basename(i) for i in inputFiles]
fullFeatures, partialFeatures, unalignedFeatures = alignFeatures(featureArray, featureFiles, paramFile)

############################################
# MS2 spectrum generation for each feature #
############################################
fullFeatures, featureToScan = ms2ForFeatures(fullFeatures, inputFiles, paramFile)

##################
# Library search #
##################
libFile = "library.db"
searchLibrary(fullFeatures, libFile, paramFile)