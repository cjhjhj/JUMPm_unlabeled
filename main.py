#!/usr/bin/python

import os, sys, re, numpy as np
from featureDetection2 import detectFeatures
# from featureAlignment import alignFeatures
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
    featureArray.append(f)


print ()
