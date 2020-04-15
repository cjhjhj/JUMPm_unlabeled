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
# for i in range(nFiles):
#     f = detectFeatures(inputFiles[i], paramFile)
#
#
#     # # Optional part (to export features to files) ###########
#     # df = pd.DataFrame.from_records(f)
#     # outputFile = os.path.basename(inputFiles[i])
#     # outputFile = os.path.splitext(outputFile)[0] + ".feature"
#     # df.to_csv(outputFile, sep = "\t", index = False)
#     # #########################################################
#
#     featureArray.append(f)


# Optional part (to load features from files) ##############################
for i in range(nFiles):
    file = os.path.splitext(os.path.basename(inputFiles[i]))[0] + ".feature"
    f = pd.read_csv(file, sep = "\t")
    f = f.to_records(index = False)
    f.dtype.names = ("mz", "intensity", "z", "RT", "minRT", "maxRT", "MS1", "minMS1", "maxMS1", "SNratio", "PercentageTF", "isotope")
    featureArray.append(f)
#############################################################################

#####################
# Feature alignment #
#####################
featureFiles = [os.path.basename(i) for i in inputFiles]
full, partial, unaligned = alignFeatures(featureArray, featureFiles, paramFile)

############################################
# MS2 spectrum generation for each feature #
############################################
ms2Array = ms2ForFeatures(full, inputFiles, paramFile)
print()
