#!/usr/bin/python

import os, sys, re, time, shutil, utils, numpy as np, pandas as pd
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
    if len(scans) > 0:
        for i in range(len(scans)):
            p = ms2[scans[i]]
            for j in range(len(spec["mz"])):
                if len(p["mz"]) == 0:
                    break
                mz = spec["mz"][j]
                intensity = spec["intensity"][j]
                lL = mz - mz * tol / 1e6
                uL = mz + mz * tol / 1e6
                ind = np.where((p["mz"] >= lL) & (p["mz"] <= uL))[0]
                if len(ind) > 0:
                    ind = ind[np.argmax(p["intensity"][ind])]
                    spec["mz"][j] = (mz * intensity + p["mz"][ind] * p["intensity"][ind]) \
                                    / (intensity + p["intensity"][ind])  # New m/z = weighted average
                    spec["intensity"][j] += p["intensity"][ind] # New intensity = sum of intensities
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
                intensity = spec["intensity"][j]
                lL = mz - mz * tol / 1e6
                uL = mz + mz * tol / 1e6
                ind = np.where((p["mz"] >= lL) & (p["mz"] <= uL))[0]
                if len(ind) > 0:
                    ind = ind[np.argmax(p["intensity"][ind])]
                    spec["mz"][j] = (mz * intensity + p["mz"][ind] * p["intensity"][ind]) \
                                    / (intensity + p["intensity"][ind])  # New m/z = weighted average
                    spec["intensity"][j] += p["intensity"][ind]  # New intensity = sum of intensities
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
    print("  ==============================================")
    full = full.to_records(index = False)  # Change pd.DataFrame to np.RecArray for internal computation (speed issue)

    ######################################
    # Load parameters and initialization #
    ######################################
    params = utils.getParams(paramFile)
    # ppiThreshold = "max"  # Hard-coded
    ppiThreshold = params["ppi_threshold_of_features"]
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
        print("  Looking for MS2 scan(s) responsible for each feature")
        for i in range(minScan, maxScan + 1):
            progress.increment()
            spec = reader[str(i)]
            msLevel = spec["msLevel"]
            if msLevel == 1:
                surveyNum = i
            elif msLevel == 2:
                # Find MS2 scans which satisfy the following conditions

                # From the discussion around June 2020,
                # 1. In ReAdW-derived mzXML files, precursor m/z values are in two tags: "precursorMz" and "filterLine"
                # 2. Through Haiyan's manual inspection, the real precursor m/z value is closer to one in "filterLine" tag
                # 3. So, in this script, precursor m/z of MS2 scan is obtained from "filterLine" tag
                # 4. Note that it may be specific to ReAdW-derived mzXML files since MSConvert-derived mzXML files do not have "filterLine" tag
                # 4.1. In this case, maybe the use of mzML (instead of mzXML) would be a solution (to-do later)

                # precMz = spec["precursorMz"][0]["precursorMz"]  # Precursor m/z from "precursorMz" tag
                p = re.search("([0-9.]+)\\@", spec["filterLine"])
                precMz = float(p.group(1))
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
                        ppiThreshold = float(ppiThreshold)
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

        print("  Merging MS2 spectra for each feature within a run (it may take a while)")
        progress = utils.progressBar(nFeatures)
        for i in range(nFeatures):
            progress.increment()
            if featureToScan[i, m] is not None:
                spec = intraConsolidation(ms2Dict, featureToScan[i, m], tolIntraMS2Consolidation)
                featureToSpec[i, m] = spec

    print("  Merging MS2 spectra for each feature between runs when there are multiple runs")
    print("  Simplification of MS2 spectrum for each feature by retaining the most strongest 100 peaks")
    specArray = np.array([])
    progress = utils.progressBar(nFeatures)
    for i in range(nFeatures):
        progress.increment()
        if np.sum(featureToSpec[i] == None) == nFiles:
            specArray = np.append(specArray, None)
        else:
            spec = interConsolidation(featureToSpec[i, :], tolInterMS2Consolidation)
            specArray = np.append(specArray, spec)

    # "specArray" is the list of (consolidated) MS2 spectra
    # specArray[i] is the MS2 spectrum corresponding to the i-th feature
    # If there's no MS2 spectrum, then specArray[i] is None
    df = utils.generateSummarizedFeatureFile(nFeatures, full, specArray, params)

    # Move mzXML files to the directory(ies) where individual .feature files are located
    for file in mzxmlFiles:
        baseFilename = os.path.basename(file)
        featureDirectory = os.path.join(os.getcwd(), os.path.splitext(baseFilename)[0])
        os.rename(file, os.path.join(featureDirectory, baseFilename))

    return df, featureToScan
