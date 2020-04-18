#!/usr/bin/python

import os, sys, re, numpy as np, pandas as pd, pickle
from featureDetection import detectFeatures
from featureAlignment import alignFeatures
from featureToMS2 import ms2ForFeatures
from pyteomics import mzxml


inputFiles = [r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.mzXML",
              r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.mzXML",
              r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.mzXML"]
nFiles = len(inputFiles)
paramFile = r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\jumpm_negative_desktop.params"

#####################
# Feature detection #
#####################
featureArray = []
for i in range(nFiles):
    f = detectFeatures(inputFiles[i], paramFile)


    # Optional part (to export features to files) ###########
    df = pd.DataFrame.from_records(f)
    outputFile = os.path.basename(inputFiles[i])
    outputFile = os.path.splitext(outputFile)[0] + ".feature"
    df.to_csv(outputFile, sep = "\t", index = False)
    #########################################################

    featureArray.append(f)


# # Optional part (to load features from files) ##############################
# for i in range(nFiles):
#     file = os.path.splitext(os.path.basename(inputFiles[i]))[0] + ".feature"
#     f = pd.read_csv(file, sep = "\t")
#     featureArray.append(f)
# #############################################################################

#####################
# Feature alignment #
#####################
featureFiles = [os.path.basename(i) for i in inputFiles]
full, partial, unaligned = alignFeatures(featureArray, featureFiles, paramFile)

# Optional part (to save fully-aliged features to a file) ###################
file = "IROA_IS_NEG_Fully_Aligned.feature"
full.to_csv(file, sep = "\t", index = False)
#############################################################################


############################################
# MS2 spectrum generation for each feature #
############################################

# # Optional part (to load features from files) ###############################
# full = pd.read_csv("IROA_IS_NEG_Fully_Aligned.feature", sep = "\t")
# #############################################################################

ms2Array = ms2ForFeatures(full, inputFiles, paramFile)

# Optional part (to save fully-aligned features and their MS2 spectra) ########
pickle.dump([full, ms2Array], open("featuresAndMS2_full_ms2Array.pickle", "wb"))

# # Optional part (to load fully-aligned features and their MS2 spectra) ########
# data = pickle.load(open("featuresAndMS2_full_ms2Array.pickle", "rb"))
# full = data[0]
# ms2Array = data[1]
# ###############################################################################


# Library search
