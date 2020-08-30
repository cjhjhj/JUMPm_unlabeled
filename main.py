#!/usr/bin/python

import sys, os, re, utils, numpy as np, pandas as pd
from featureDetection import detectFeatures
from featureAlignment import alignFeatures
from featureToMS2 import ms2ForFeatures
from librarySearch import searchLibrary
from datetime import datetime

##################
# Initialization #
##################
# args = sys.argv
# del args[0]
# paramFile = args[0]
# inputFiles = args[1:]

# For desktop debugging,
# paramFile = r"/Research/Projects/7Metabolomics/Dev/JUMPm_unlabel_python/comparison_test/python/jumpm_positive.params"
# inputFiles = [r"/Research/Projects/7Metabolomics/Dev/JUMPm_unlabel_python/comparison_test/python/IROA_c18_target1.mzXML",
#               r"/Research/Projects/7Metabolomics/Dev/JUMPm_unlabel_python/comparison_test/python/IROA_c18_target2.mzXML"]

paramFile = r"/Research/Projects/7Metabolomics/Dev/JUMPm_unlabel_python/koa_wt/library_jumpm_negative.params"
inputFiles = [r"/Research/Projects/7Metabolomics/Dev/JUMPm_unlabel_python/koa_wt/neg_ko_a1/neg_ko_a1.feature",
              r"/Research/Projects/7Metabolomics/Dev/JUMPm_unlabel_python/koa_wt/neg_ko_a2/neg_ko_a2.feature",
              r"/Research/Projects/7Metabolomics/Dev/JUMPm_unlabel_python/koa_wt/neg_ko_a3/neg_ko_a3.feature"]


print()
print("  Jump -m started")
now = datetime.now()
nowString = now.strftime("%Y/%m/%d %H:%M:%S")
print("  " + nowString)
params = utils.getParams(paramFile)

try:
    #####################
    # Feature detection #
    #####################
    print()
    print("  #####################")
    print("  # Feature detection #")
    print("  #####################")
    featureArray = []
    if params["skip_feature_detection"] == "0":
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
            mzxmlFile = re.sub(".feature$", ".mzXML", file)
            if os.path.isfile(mzxmlFile):
                inputFiles.append(mzxmlFile)
            else:
                sys.exit("  {} should be in the directory where {} is located".format(os.path.basename(mzxmlFile), os.path.basename(file)))
        print("  According to the parameter setting, the feature detection is skipped")
        print("  Feature(s) is/are obtained from .feature file(s)")

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
    res = searchLibrary(fullFeatures, paramFile)
    print()

    print("  Jump -m finished")
    now = datetime.now()
    nowString = now.strftime("%Y/%m/%d %H:%M:%S")
    print("  " + nowString)
except KeyboardInterrupt:
    sys.exit()