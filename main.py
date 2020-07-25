#!/usr/bin/python

import os, re, utils, numpy as np, pandas as pd
from featureDetection import detectFeatures
from featureAlignment import alignFeatures
from featureToMS2 import ms2ForFeatures
from librarySearch import searchLibrary


# paramFile = r"/Research/Projects/7Metabolomics/htan_IROA/2020/hilic/align_test/jumpm_negative.params"
paramFile = r"/Research/Projects/7Metabolomics/Dev/JUMPm_unlabel_python/jumpm_positive.params"
params = utils.getParams(paramFile)

#####################
# Feature detection #
#####################
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
fileNames = [os.path.basename(i) for i in inputFiles]
fullFeatures, partialFeatures, unalignedFeatures = alignFeatures(featureArray, fileNames, paramFile)
print()
# Print out features (fully-, partially and un-aligned features)

############################################
# MS2 spectrum generation for each feature #
############################################
fullFeatures, featureToScan = ms2ForFeatures(fullFeatures, inputFiles, paramFile)
print()

##################
# Library search #
##################
# libFile = "library.db"
# searchLibrary(fullFeatures, libFile, paramFile)
# print()