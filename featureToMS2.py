#!/usr/bin/python

import os, sys, re, pickle, time, utils, numpy as np, pandas as pd
from numpy.lib.recfunctions import append_fields
from pyteomics import mzxml


def intraConsolidation(ms2, scans, tol):
    # Sort MS2 spectra according to their total ion current (descending order)
    scans = scans.split(";")
    tic = np.array([sum(ms2[key]["intensity"]) for key in scans])
    ind = (-tic).argsort()  # Sort "tic" in descending order
    scans = [scans[i] for i in ind]
    spec = ms2[scans[0]]  # MS2 spectrum with the highest total ion current
    del scans[0]
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
    specs = [i for i in specs if i is not None]  # Skip "None"
    tic = [sum(s["intensity"]) for s in specs]
    if len(tic) > 1:
        tic = np.array(tic)
        ind = (-tic).argsort()  # Sort "tic" in descending order
        specs = [specs[i] for i in ind]
        spec = specs[0]
        del specs[0]
    else:
        spec = specs[0]
        specs = []

    # tic = [sum(s["intensity"]) for s in specs]
    # if len(tic) > 1:
    #     ind = np.argmax(tic)
    #     spec = specs[ind]  # Reference MS2 spectrum for merging others for a feature
    #     del specs[ind]
    # else:
    #     spec = specs[0]
    #     specs = []

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
                ind10 = sorted(range(len(spec["intensity"][ind])), key=lambda j: spec["intensity"][ind][j],
                               reverse=True)[:10]
                ind = ind[ind10]
            filteredSpec["mz"] = np.append(filteredSpec["mz"], spec["mz"][ind])
            filteredSpec["intensity"] = np.append(filteredSpec["intensity"], spec["intensity"][ind])
        spec = filteredSpec
    # Sort the spectrum in ascending order of m/z
    ind = np.argsort(spec["mz"])
    spec["mz"] = spec["mz"][ind]
    spec["intensity"] = spec["intensity"][ind]
    return spec


def ms2ForFeatures(full, mzxmlFiles, paramFile):
    print("  Identification of MS2 spectra for the features")
    full = full.to_records(index=False)  # Change pd.dataframe to np.recarray for internal computation

    ######################################
    # Load parameters and initialization #
    ######################################
    ppiThreshold = "max"  # Hard-coded
    params = utils.getParams(paramFile)
    pctTfThreshold = float(params["max_percentage_RT_range"])
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
    return specArray
