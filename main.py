#!/usr/bin/python

import os, sys, re, subprocess, numpy as np, pandas as pd, pickle, utils
from pathlib import Path

from featureDetection import detectFeatures
from featureAlignment import alignFeatures
from featureToMS2 import ms2ForFeatures
from pyteomics import mzxml

# # Input arguments
# n = len(sys.argv)
# paramFile = sys.argv[1]
# inputFiles = sys.argv[2:n]
# nFiles = len(inputFiles)

# Test in desktop
paramFile = r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\jumpm_negative_desktop.params"
inputFiles = [r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.mzXML",
              r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.mzXML",
              r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.mzXML"]
nFiles = len(inputFiles)

###############################
# File and directory handling #
###############################
params = utils.getParams(paramFile)
# if params["skip_feature_detection"] == "0":
#     # Perform feature detection using .mzXML file(s)
#     # Prepare directories according to the names of input .mzXML file(s)
#     print("  Using the following files")
#     for file in inputFiles:
#         print("  %s" % os.path.basename(file))
#         # Each mzXML file and its feature file will be located at "dataDir"
#         dataDir = os.path.splitext(file)[0]
#         Path(dataDir).mkdir(exist_ok=True)
# elif params["skip_feature_detection"] == "1":
#     # Skip feature detection, but directly use .feature files
#     # Check the existence of .feature file(s) and corresponding mzXML file(s)


featureFiles = [r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.1.feature",
                r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.1.feature",
                r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.1.feature"]
mzxmlFiles = [r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.mzXML",
              r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.mzXML",
              r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.mzXML"]

######################################
# Feature alignment (using R script) #
######################################
binDir = os.path.abspath(os.path.dirname(__file__))
command = [r"C:\Program Files\R\R-3.6.2\bin\Rscript", binDir + r"\R\alignment.R",
           paramFile, ",".join(featureFiles), "./test/jumpm.log", "./test/"]
subprocess.call(command)



#####################
# Feature detection #
#####################
# featureArray = []
# for i in range(nFiles):
#     f = detectFeatures(inputFiles[i], paramFile)
#
#
#     # Optional part (to export features to files) ###########
#     df = pd.DataFrame.from_records(f)
#     outputFile = os.path.basename(inputFiles[i])
#     outputFile = os.path.splitext(outputFile)[0] + ".feature"
#     df.to_csv(outputFile, sep = "\t", index = False)
#     #########################################################
#
#     # featureArray.append(f)


# # Optional part (to load features from files) ##############################
# for i in range(nFiles):
#     file = os.path.splitext(os.path.basename(inputFiles[i]))[0] + ".feature"
#     f = pd.read_csv(file, sep = "\t")
#     featureArray.append(f)
# #############################################################################


#######################################
# Feature alignment (using R scripts) #
#######################################
# featureFiles = [r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.1.feature",
#               r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.1.feature",
#               r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.1.feature"]
#
# # cmd = r"C:\Program Files\R\R-3.6.2\bin\Rscript --vanilla ./R/test.R 10 'abc'"
# # cmd = r"C:\Program Files\R\R-3.6.2\bin\R --vanilla test.R 10 'abc'"
# # out = subprocess.Popen(cmd, shell = True)
# subprocess.call([r'C:\Program Files\R\R-3.6.2\bin\Rscript', "test.R", "10", "abc"])
# print()





'''
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

'''