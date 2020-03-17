#!/usr/bin/python

import os, sys, re, pickle, utils, numpy as np
from pyteomics import mzxml




# For development purpose #################################################################
with open("featureToMS2.pickle", "rb") as f:
    vars = pickle.load(f)

featureFiles = vars[0]
full = vars[1]
params = vars[2]

mzXMLFiles = [r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.mzXML",
              r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.mzXML",
              r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.mzXML"]

# Check the length of featureFiles and mzXMLFiles
if len(featureFiles) != len(mzXMLFiles):
    sys.exit("Input feature files and mzXML files do not match")
nFiles = len(featureFiles)


# readers = []
# for i in range(nFiles):
#     reader = mzxml.read(mzXMLFiles[i])
#     readers.append(reader)

############################################################################################





######################################
# Load parameters and initialization #
######################################
ppiThreshold = "max"    # Hard-coded
pctTfThreshold = 50 # Hard-coded
tolIsolation = float(params['isolation_window'])
tolPrecursor = float(params['tol_precursor'])
tolIntraMS2Consolidation = float(params['tol_intra_ms2_consolidation'])
tolInterMS2Consolidation = float(params['tol_inter_ms2_consolidation'])
featureToScan = np.empty((full.shape[0], nFiles))
featureToScan[:] = np.nan

########################################################
# Find MS2 scans responsible for features within a run #
########################################################
m = -1  # Index for input (mzXML) files
for file in mzXMLFiles:


    # File handling may need to be revised according to the format of wrapper ######

    reader = mzxml.read(file)
    m += 1

    ################################################################################

    # Column name should be well organized and managed #############################

    # Handling column names
    fileBasename, _ = os.path.splitext(os.path.basename(file))
    r = re.compile(fileBasename + ".*minMS1.*")
    colMinMS1 = list(filter(r.match, full.dtype.names))[0]
    r = re.compile(fileBasename + ".*maxMS1.*")
    colMaxMS1 = list(filter(r.match, full.dtype.names))[0]
    r = re.compile(fileBasename + ".*mz.*")
    colMz = list(filter(r.match, full.dtype.names))[0]
    r = re.compile(fileBasename + ".*TF.*")
    colTF = list(filter(r.match, full.dtype.names))[0]

    #################################################################################

    with reader:
        for spec in reader:
            msLevel = int(spec['msLevel'])
            if msLevel == 1:
                surveyScanNumber = int(spec['num'])
            elif msLevel == 2:
                # Find MS2 scans which satisfy
                # 1. Their precursor m/z values are within feature's m/z +/- isolation_window
                # 2. Their presumable survey scans are between feature's min. and max. MS1 scan numbers
                # 3. Feature width should be less than a threshold
                precMz = float(spec['precursorMz'][0]['precursorMz'])
                fInd = np.where((full[colMinMS1] < surveyScanNumber) &
                                (full[colMaxMS1] > surveyScanNumber) &
                                (full[colMz] >= (precMz - tolIsolation)) &
                                (full[colMz] <= (precMz + tolIsolation)) &
                                (full[colTF] < pctTfThreshold))[0]
                if len(fInd) == 0:
                    continue
                else:
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
                        mz = full[colMz][fInd[j]]
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
                        fInd = np.array([fInd[np.argmax(ppi)]])
                    else:
                        # ppiThreshold should be a numeric value
                        fInd = fInd[np.where(ppi > ppiThreshold)]
                    if len(fInd) == 0:  # Last check of candidate feature indexes
                        continue
                    else:
                        # Mapping between features and MS2 scan numbers
                        for j in range(len(fInd)):
                            if np.isnan(featureToScan[fInd[j], m]):
                                featureToScan[fInd[j], m] = int(spec["num"])
                            else:
                                featureToScan[fInd[j], m].append(int(spec["num"]))

        print()