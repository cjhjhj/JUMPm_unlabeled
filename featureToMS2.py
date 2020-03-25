#!/usr/bin/python

import os, sys, re, pickle, utils, numpy as np, pandas as pd
from numpy.lib.recfunctions import append_fields
from pyteomics import mzxml


def intraConsolidation(ms2, scans, tol):
    # Sort MS2 spectra according to their total ion current (descending order)
    scans = scans.split(";")
    tic = [sum(ms2[key]["intensity"]) for key in scans]
    ind = np.argmax(tic)
    spec = ms2[scans[ind]]  # MS2 spectrum with the highest total ion current
    del scans[ind]
    if len(scans) > 1:
        for i in range(len(scans)):
            p = ms2[scans[i]]
            for j in range(len(spec["mz"])):
                if len(p["mz"]) == 0:
                    break
                mz = spec["mz"][j]
                lL = mz - mz * tol / 1e6
                uL = mz + mz * tol / 1e6
                ind = np.where((p["mz"] >= lL) & (p["mz"] <= uL))[0]
                if len(ind) > 0:
                    ind = ind[np.argmax(p["intensity"][ind])]
                    spec["intensity"][j] += p["intensity"][ind]
                    p["mz"] = np.delete(p["mz"], ind)
                    p["intensity"] = np.delete(p["intensity"], ind)
            spec["mz"] = np.append(spec["mz"], p["mz"])
            spec["intensity"] = np.append(spec["intensity"], p["intensity"])
    # Sort the spectrum in ascending order of m/z
    ind = np.argsort(spec["mz"])
    spec["mz"] = spec["mz"][ind]
    spec["intensity"] = spec["intensity"][ind]
    return spec


def interConsolidation(specs, tol):
    specs = [i for i in specs if i is not None] # Skip "None"
    tic = [sum(s["intensity"]) for s in specs]
    if len(tic) > 1:
        ind = np.argmax(tic)
        spec = specs[ind]  # Reference MS2 spectrum for merging others for a feature
        del specs[ind]
    else:
        spec = specs[0]
        specs = []

    if len(specs) > 0:
        for i in range(len(specs)):
            p = specs[i]
            for j in range(len(spec["mz"])):
                if len(p["mz"]) == 0:
                    break
                mz = spec["mz"][j]
                lL = mz - mz * tol / 1e6
                uL = mz + mz * tol / 1e6
                ind = np.where((p["mz"] >= lL) & (p["mz"] <= uL))[0]
                if len(ind) > 0:
                    ind = ind[np.argmax(p["intensity"][ind])]
                    spec["intensity"][j] += p["intensity"][ind]
                    p["mz"] = np.delete(p["mz"], ind)
                    p["intensity"] = np.delete(p["intensity"], ind)
            spec["mz"] = np.append(spec["mz"], p["mz"])
            spec["intensity"] = np.append(spec["intensity"], p["intensity"])
    # Sort the spectrum in ascending order of m/z
    ind = np.argsort(spec["mz"])
    spec["mz"] = spec["mz"][ind]
    spec["intensity"] = spec["intensity"][ind]

    # Simplification of the merged spectrum
    # 1. Limit the number of peaks in each .dta file to 100
    # 2. Divide m/z-range into 10 bins (e.g. 0~100, 100~200, etc.) and retain the 10 largest peaks in each bin
    if len(spec["mz"]) > 100:
        nBins = 10
        bins = np.linspace(min(spec["mz"]), max(spec["mz"]), (nBins + 1))
        filteredSpec = {"mz": [], "intensity": []}
        for i in range(nBins):
            ind = np.where((spec["mz"] >= bins[i]) & (spec["mz"] < bins[i + 1]))[0]
            if len(ind) > 10:
                # Select 10 highest intensity peaks
                ind10 = sorted(range(len(spec["intensity"][ind])), key = lambda j: spec["intensity"][ind][j], reverse = True)[:10]
                ind = ind[ind10]
            filteredSpec["mz"] = np.append(filteredSpec["mz"], spec["mz"][ind])
            filteredSpec["intensity"] = np.append(filteredSpec["intensity"], spec["intensity"][ind])
        spec = filteredSpec
    # Sort the spectrum in ascending order of m/z
    ind = np.argsort(spec["mz"])
    spec["mz"] = spec["mz"][ind]
    spec["intensity"] = spec["intensity"][ind]
    return spec


# For development purpose #################################################################
with open("featureToMS2_3.pickle", "rb") as f:
    vars = pickle.load(f)

featureFiles = vars[0]
full = vars[1]
full.dtype.names = [name.replace("ScanNumber", "") if name.endswith("ScanNumber") else name for name in full.dtype.names]
params = vars[2]
featureToScan = vars[3]
featureToSpec = vars[4]

mzXMLFiles = [r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.mzXML",
              r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.mzXML",
              r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.mzXML"]

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
tolIsolation = float(params["isolation_window"])
tolPrecursor = float(params["tol_precursor"])
tolIntraMS2Consolidation = float(params["tol_intra_ms2_consolidation"])
tolInterMS2Consolidation = float(params["tol_inter_ms2_consolidation"])

df = pd.DataFrame(data = full)
nFeatures = df.shape[0]
# featureToScan = np.empty((nFeatures, nFiles), dtype = object)
# featureToSpec = np.empty((nFeatures, nFiles), dtype = object)
#
# ########################################################
# # Find MS2 scans responsible for features within a run #
# ########################################################
# m = -1  # Index for input (mzXML) files
# for file in mzXMLFiles:
#
#
#     # File handling may need to be revised according to the format of wrapper ######
#     m += 1
#
#
#     print ("  Reading %s" % os.path.basename(file))
#     reader = mzxml.MzXML(file)
#     nScans = len(reader)
#
#     fileBasename, _ = os.path.splitext(os.path.basename(file))
#     colInd = [i for i, item in enumerate(df.columns) if re.search(fileBasename, item)]
#     subDf = df.iloc[:, colInd]
#     subDf.columns = [s.split("_")[-1] for s in subDf.columns]
#     ms2Dict = {}
#     progress = utils.progressBar(nScans)
#     print("  Finding MS2 spectra of %s responsible for features" % os.path.basename(file))
#     for i in range(nScans):
#         progress.increment()
#         if i < min(subDf["minMS1"]) or i > max(subDf["maxMS1"]):
#             continue
#
#
#             # In fact, when i > max(subDf["maxMS1"]), we can use "break",i.e. no need to iterate
#             # However, progressBar may stop
#
#
#         else:
#             spec = reader[str(i)]
#             msLevel = int(spec["msLevel"])
#             if msLevel == 1:
#                survey = spec
#             elif msLevel == 2:
#                 # Find MS2 scans which satisfy
#                 # 1. Their precursor m/z values are within feature's m/z +/- isolation_window
#                 # 2. Their presumable survey scans are between feature's min. and max. MS1 scan numbers
#                 # 3. Feature width should be less than a threshold
#                 precMz = float(spec["precursorMz"][0]["precursorMz"])
#                 fInd = np.where((subDf["minMS1"] < int(survey["num"])) &
#                                 (subDf["maxMS1"] > int(survey["num"])) &
#                                 (subDf["mz"] >= (precMz - tolIsolation)) &
#                                 (subDf["mz"] <= (precMz + tolIsolation)) &
#                                 (subDf["PercentageTF"] < pctTfThreshold))[0]
#                 if len(fInd) == 0:
#                     continue
#                 else:
#                     # Check the intensities of candidate features at the very preceding MS1 scan
#                     # For example, let's assume that candidate features are as follows for a MS2 scan #141
#                     # index  mz        z  MS1  minMS1 maxMS1 Intensity
#                     # 3      218.1498  0  136  1      951    37544
#                     # 20     220.0705  0  126  1      1446   91709
#                     # 40     218.8597  6  18   1      764    91745
#                     # 65     220.1052  0  1    1      1248   355843
#                     # Also, suppose that the very preceding MS1 scan (i.e. survey scan) is scan#140
#                     # Then, we need to check the intensities of those candidate features at scan#140
#                     # example) For the 1st feature whose representative m/z = 218.1498,
#                     #          1. Get MS1 spectrum of scan#140
#                     #          2. Open a m/z window with a tolerance; [218.1498 - 10ppm, 218.1498 + 10ppm]
#                     #          3. Look for the MS1 peak with the highest intensity within the window, and record the intensity
#                     #          4. If the intensity is the highest among candidate features, choose it for the MS2 scan
#                     #          5. Otherwise, check the next candidate feature
#                     #          6. Instead of choosing one feature with the highest intensity,
#                     #             PPI can be used and multiple features may be chosen for each MS2
#                     ppi = []
#                     for i in range(len(fInd)):
#                         mz = subDf["mz"][fInd[i]]
#                         lL = mz - mz * tolPrecursor / 1e6
#                         uL = mz + mz * tolPrecursor / 1e6
#                         ind = np.where((survey["m/z array"] >= lL) & (survey["m/z array"] <= uL))[0]
#                         if len(ind) > 0:
#                             ppi.append(np.max(survey["intensity array"][ind]))
#                         else:
#                             ppi.append(0)
#
#                     if sum(ppi) == 0:
#                         continue
#                     ppi = ppi / np.sum(ppi) * 100  # Convert intensities to percentage values
#                     if ppiThreshold == "max":
#                         fInd = np.array([fInd[np.argmax(ppi)]])
#                     else:
#                         # ppiThreshold should be a numeric value
#                         fInd = fInd[np.where(ppi > ppiThreshold)]
#                     if len(fInd) == 0:  # Last check of candidate feature indexes
#                         continue
#                     else:
#                         # Add this MS2 scan information to ms2Dict
#                         ms2Dict[spec["num"]] = {}
#                         ms2Dict[spec["num"]]["mz"] = spec["m/z array"]
#                         ms2Dict[spec["num"]]["intensity"] = spec["intensity array"]
#                         # Mapping between features and MS2 scan numbers
#                         for i in range(len(fInd)):
#                             if featureToScan[fInd[i], m] is None:
#                                 featureToScan[fInd[i], m] = spec["num"]
#                             else:
#                                 featureToScan[fInd[i], m] += ";" + spec["num"]
#
#     ##############################################################
#     # Consolidation of MS2 spectra for each feature within a run #
#     ##############################################################
#     print("  Merging MS2 spectra within each feature")
#     progress = utils.progressBar(nFeatures)
#     for i in range(nFeatures):
#         progress.increment()
#         if featureToScan[i, m] is None:
#             continue
#         else:
#             spec = intraConsolidation(ms2Dict, featureToScan[i, m], tolIntraMS2Consolidation)
#             featureToSpec[i, m] = spec

#############################################
# Consolidation of MS2 spectra between runs #
#############################################
print ("  Merging MS2 spectra between runs for each feature")
progress = utils.progressBar(nFeatures)
for i in range(nFeatures):
    progress.increment()
    if sum(featureToSpec[i, :] == None) == nFiles:
        continue
    spec = interConsolidation(featureToSpec[i, :], tolInterMS2Consolidation)
print ()





