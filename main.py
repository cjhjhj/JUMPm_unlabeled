#!/usr/bin/python

import os, re, utils, numpy as np, pandas as pd
from featureDetection import detectFeatures
from featureAlignment import alignFeatures
from featureToMS2 import ms2ForFeatures
from librarySearch import searchLibrary

print("  Jump -m started")
print()

# paramFile = r"/Research/Projects/7Metabolomics/htan_IROA/2020/hilic/align_test/jumpm_negative.params"
paramFile = r"/Research/Projects/7Metabolomics/Dev/JUMPm_unlabel_python/jumpm_positive.params"
params = utils.getParams(paramFile)

#####################
# Feature detection #
#####################
print("  #####################")
print("  # Feature detection #")
print("  #####################")
featureArray = []
if params["skip_feature_detection"] == "0":
    # inputFiles = [r"/Research/Projects/7Metabolomics/htan_IROA/2020/hilic/IROA_neg_target/IROA_neg_target.mzXML"]
    inputFiles = [r"/Research/Projects/7Metabolomics/htan_IROA/2019/c18/test/IROA_c18_target1/IROA_c18_target1.mzXML",
                  r"/Research/Projects/7Metabolomics/htan_IROA/2019/c18/test/IROA_c18_target2/IROA_c18_target2.mzXML"]
    nFiles = len(inputFiles)
    for i in range(nFiles):
        f = detectFeatures(inputFiles[i], paramFile)
        featureArray.append(f)
elif params["skip_feature_detection"] == "1":
    # Note that .feature file(s) and .mzXML file(s) should be in the same directory
    inputFiles = []
    for file in params["feature_files"]:
        # Read .feature file and append to featureArray
        f = pd.read_csv(file, sep = "\t")
        featureArray.append(f)

        # Define inputFiles array containing mzXML file(s)
        inputFiles.append(re.sub(".feature$", ".mzXML", file))

print()

#####################
# Feature alignment #
#####################
# When there are multiple runs (samples), features from them will be aligned against a reference run
# (the reference run is selected by (1) number of features, (2) intensity level of most strongest features,
#  or (3) user-specified run)
# Output variables
# 1. fullFeatures: pandas DataFrame of fully-aligned features (with run-specific information)
# 2. partialFeatures: pandas DataFrame of partially-aligned features (with run-specific information)
# 3. unalignedFeatures: array of pandas DataFrame of unaligned features
# In addition, those features are written to files
print("  #####################")
print("  # Feature alignment #")
print("  #####################")
fileNames = [os.path.basename(i) for i in inputFiles]
fullFeatures, partialFeatures, unalignedFeatures = alignFeatures(featureArray, fileNames, paramFile)
print()

############################################
# MS2 spectrum generation for each feature #
############################################
print("  #######################################")
print("  # Processing MS2 spectra for features #")
print("  #######################################")
fullFeatures, featureToScan = ms2ForFeatures(fullFeatures, inputFiles, paramFile)
print()

##################
# Library search #
##################
print("  ##################")
print("  # Library search #")
print("  ##################")
libFile = r"/Research/Projects/7Metabolomics/library/StJude/stjude_library_c18p.db"
res = searchLibrary(fullFeatures, libFile, paramFile)
print()