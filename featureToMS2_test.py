#!/usr/bin/python

import os, sys, re, pickle, utils, numpy as np, pandas as pd
from pyteomics import mzxml
from featureToMS2 import interConsolidation, intraConsolidation

# f = open('fully_aligned_features.pickle', 'rb')
# full = pickle.load(f)
# f.close()
# mzxmlFiles = [r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.mzXML",
#               r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.mzXML",
#               r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.mzXML"]
# paramFile = r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\jumpm_negative_desktop.params"


full = pd.read_csv("./IROA_IS_NEG/.IROA_IS_NEG_fully_aligned.feature", sep="\t", index_col=False)
for col in full.columns:
    if col.endswith("minMS1ScanNumber"):
        full.rename({col: re.sub("minMS1ScanNumber", "minMS1", col)}, axis=1, inplace=True)
    elif col.endswith("maxMS1ScanNumber"):
        full.rename({col: re.sub("maxMS1ScanNumber", "maxMS1", col)}, axis=1, inplace=True)
    elif col.endswith("Intensity"):
        full.rename({col: re.sub("Intensity", "intensity", col)}, axis=1, inplace=True)
    elif col.endswith("PercentageofTF"):
        full.rename({col: re.sub("PercentageofTF", "PercentageTF", col)}, axis=1, inplace=True)

mzxmlFiles = [r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.mzXML",
              r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_2.mzXML",
              r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_3.mzXML"]
paramFile = r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\jumpm_negative_desktop.params"
full = full.to_records(index=False)  # Change pd.dataframe to np.recarray for internal computation




######################################
# Load parameters and initialization #
######################################
ppiThreshold = "max"  # Hard-coded
pctTfThreshold = 100  # Hard-coded
params = utils.getParams(paramFile)
tolIsolation = float(params["isolation_window"])
tolPrecursor = float(params["tol_precursor"])
tolIntraMS2Consolidation = float(params["tol_intra_ms2_consolidation"])
tolInterMS2Consolidation = float(params["tol_inter_ms2_consolidation"])
nFeatures = len(full)
nFiles = len(mzxmlFiles)
featureToScan = np.empty((nFeatures, nFiles), dtype=object)
featureToSpec = np.empty((nFeatures, nFiles), dtype=object)

m = -1  # Index for input files
for file in mzxmlFiles:
    m += 1
    reader = mzxml.MzXML(file)
    fileBasename, _ = os.path.splitext(os.path.basename(file))
    colNames = [item for item in full.dtype.names if re.search(fileBasename, item)]
    subset = full[colNames]
    subset.dtype.names = [s.split("_")[-1] for s in subset.dtype.names]
    ms2Dict = {}
    minScan, maxScan = int(min(subset["minMS1"])), int(max(subset["maxMS1"]))
    progress = utils.progressBar(maxScan - minScan + 1)
    print("  %s is being processed" % os.path.basename(file))
    for i in range(minScan, maxScan + 1):
        progress.increment()
        spec = reader[str(i)]
        msLevel = spec["msLevel"]
        if msLevel == 1:
            surveyNum = i
        elif msLevel == 2:
            # Find MS2 scans which satisfy the following conditions
            precMz = spec["precursorMz"][0]["precursorMz"]
            survey = reader[str(surveyNum)]
            fInd = np.where((surveyNum >= subset["minMS1"]) &
                            (surveyNum <= subset["maxMS1"]) &
                            (subset["mz"] >= (precMz - tolIsolation)) &
                            (subset["mz"] <= (precMz + tolIsolation)) &
                            (subset["PercentageTF"] <= pctTfThreshold))[0]
            if len(fInd) > 0:
                ppi = []
                for i in range(len(fInd)):
                    mz = subset["mz"][fInd[i]]
                    lL = mz - mz * tolPrecursor / 1e6
                    uL = mz + mz * tolPrecursor / 1e6
                    ind = np.where((survey["m/z array"] >= lL) & (survey["m/z array"] <= uL))[0]
                    if len(ind) > 0:
                        ppi.append(np.max(survey["intensity array"][ind]))
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
                    # Add this MS2 scan information to ms2Dict
                    ms2Dict[spec["num"]] = {}
                    ms2Dict[spec["num"]]["mz"] = spec["m/z array"]
                    ms2Dict[spec["num"]]["intensity"] = spec["intensity array"]
                    # Mapping between features and MS2 scan numbers
                    for i in range(len(fInd)):
                        if featureToScan[fInd[i], m] is None:
                            featureToScan[fInd[i], m] = spec["num"]
                        else:
                            featureToScan[fInd[i], m] += ";" + spec["num"]

    print("  Merging MS2 spectra within a run for each feature")
    progress = utils.progressBar(nFeatures)
    for i in range(nFeatures):
        progress.increment()
        if featureToScan[i, m] is not None:
            spec = intraConsolidation(ms2Dict, featureToScan[i, m], tolIntraMS2Consolidation)
            featureToSpec[i, m] = spec
    print()

print("  Merging MS2 spectra between runs for each feature")
specArray = np.array([])
progress = utils.progressBar(nFeatures)
for i in range(nFeatures):
    progress.increment()
    if np.sum(featureToSpec[i] == None) == nFiles:
        specArray = np.append(specArray, None)
    else:
        # Heuristic charge determination across runs
        # 1. Charge state other than 0
        # 2. More frequent charge state
        # 3. If the same frequency, choose the lower one
        colNames = [col for col in full.dtype.names if col.endswith("_z")]
        charges = [c for c in full[colNames][i] if c > 0]
        if len(charges) == 0:
            charge = 1
        else:
            charge = max(set(charges), key=charges.count)

        spec = interConsolidation(featureToSpec[i, :], tolInterMS2Consolidation)
        spec["charge"] = charge
        specArray = np.append(specArray, spec)

# "specArray" is the list of (consolidated) MS2 spectra
# specArray[i] is the MS2 spectrum corresponding to the i-th feature
# If there's no MS2 spectrum, then specArray[i] is None
print()