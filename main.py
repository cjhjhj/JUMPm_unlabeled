#!/usr/bin/python

import os, sys, re, numpy as np
from featureDetection import detectFeatures
from featureAlignment import alignFeatures
from pyteomics import mzxml

# inputFile = r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.mzXML"
# reader = mzxml.read(inputFile)

inputFiles = [r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.mzXML",
              r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.mzXML",
              r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.mzXML"]
nFiles = len(inputFiles)
paramFile = r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\jumpm_negative_desktop.params"

#####################
# Feature detection #
#####################
featureArray = []
for i in range(nFiles):
    f = featureDetection.featureDetection(inputFiles[i], paramFile)
    featureArray.append(f)


print ()
