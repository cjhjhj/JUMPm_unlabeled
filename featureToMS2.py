#!/usr/bin/python

import sys, re, numpy as np
import statsmodels as stat
import utils
from pyteomics import mzxml

paramFile = r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\jumpm_negative.params"
featureFiles = [r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.1.feature",
                r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.1.feature",
                r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.1.feature"]
mzXMLFiles = [r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.mzXML",
              r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.mzXML",
              r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.mzXML"]

# Check the length of featureFiles and mzXMLFiles
if len(featureFiles) != len(mzXMLFiles):
    sys.exit("Input feature files and mzXML files do not match")
nFiles = len(featureFiles)

################################
# Load parameters and features #
################################
params = utils.readParams(paramFile)
ppiThreshold = "max"    # Hard-coded
pctTfThreshold = 50 # Hard-coded
tolIsolation = float(params['isolation_window'])
tolPrecursor = float(params['tol_precursor'])
tolIntraMS2Consolidation = float(params['tol_intra_ms2_consolidation'])
tolInterMS2Consolidation = float(params['tol_inter_ms2_consolidation'])

# Features from .feature files are stored in fArray. For example,
# featureFiles = [file1, file2, file3]
# fArray[0] = features from file1 (which has column names like 'index', 'mz', etc.)
# fArray[1] = features from file2
# ...
# The array of m/z values from the first feature file can be accessed by fArray[0]['mz']
# The header of .feature file is used as column names of the array
# Note that "/" (slash) is ignored when the file is loaded through "genfromtxt"
fArray = []
for i in range(0, len(featureFiles)):
    data = np.genfromtxt(featureFiles[i], delimiter = "\t", dtype = None, names = True)
    fArray.append(data)

########################################################
# Find MS2 scans responsible for features within a run #
########################################################
featureToScan = np.empty(shape = [, nFiles])
for i in range(0, nFiles):
    print ("Reading %s" % mzXMLFiles[i])
    reader = mzxml.read(mzXMLFiles[i])
    with reader:
        for spec in reader:
            if int(spec['num']) > int(params['last_scan_extraction']):
                break

            # Find MS2 scans which satisfy
            # 1. Their precursor m/z values are within feature's m/z +/- isolation_window
            # 2. Their presumable survey scans are between feature's min. and max. MS1 scan numbers
            # 3. Feature width should be less than a threshold
            msLevel = int(spec['msLevel'])
            if msLevel == 1:
                surveyScanNumber = int(spec['num'])
            elif msLevel == 2:
                precMz = float(spec['precursorMz'][0]['precursorMz'])
                fInd = np.where((fArray[i]['minMS1ScanNumber'] < surveyScanNumber) &
                               (fArray[i]['maxMS1ScanNumber'] > surveyScanNumber) &
                               (fArray[i]['mz'] >= (precMz - tolIsolation)) &
                               (fArray[i]['mz'] <= (precMz + tolIsolation)) &
                               (fArray[i]['PercentageTF'] < pctTfThreshold))[0]
                if len(fInd) == 0:
                    continue

                # Check the intensities of candidate features at the very preceding MS1 scan
                # For example, let's assume that candidate features are as follows for a MS2 scan #141
                # index  mz        z  MS1  minMS1 maxMS1 Intensity
                # 3      218.1498  0  136  1      951    37544
                # 20     220.0705  0  126  1      1446   91709
                # 40     218.8597  6  18   1      764    91745
                # 65     220.1052  0  1    1      1248   355843
                # Also, suppose that the very preceding MS1 scan (i.e. survey scan) is scan#140
                # Then, we need to check the intensities of those candidate features at scan#140
                # example) For the 1st feature whose representative m/z = 218.1498,
                #          1. Get MS1 spectrum of scan#140
                #          2. Open a m/z window with a tolerance; [218.1498 - 10ppm, 218.1498 + 10ppm]
                #          3. Look for the MS1 peak with the highest intensity within the window, and record the intensity
                #          4. If the intensity is the highest among candidate features, choose it for the MS2 scan
                #          5. Otherwise, check the next candidate feature
                #          6. Instead of choosing one feature with the highest intensity,
                #             PPI can be used and multiple features may be chosen for each MS2
                surveyScan = reader[str(surveyScanNumber)]
                ppi = []
                for j in range(0, len(fInd)):
                    mz = fArray[i]['mz'][fInd[j]]
                    lL = mz - mz * tolPrecursor * 1e6
                    uL = mz + mz * tolPrecursor * 1e6
                    jInd = np.where((surveyScan['m/z array'] >= lL) & (surveyScan['m/z array'] <= uL))[0]
                    if len(jInd) > 0:
                        ppi.append(np.max(surveyScan['intensity array'][jInd]))
                    else:
                        ppi.append(0)

                if sum(ppi) == 0:
                    continue
                ppi = ppi / np.sum(ppi) * 100  # Convert intensities to percentage values
                if ppiThreshold == "max":
                    fInd = fInd[np.argmax(ppi)]
                else:
                    # ppiThreshold should be a numeric value
                    fInd = fInd[np.where(ppi > ppiThreshold)]
                if len(fInd) == 0:  # Last check of candidate feature indexes
                    continue

                for j in range(0, len(fInd)):





print ()
