#!/usr/bin/python

import os, re, pickle, logging, utils, numpy as np, pandas as pd
from pyteomics import mzxml


def groupMzValues(mz, ppm):
    mzdiff = np.diff(mz)
    if ppm > 0:
        res = mzdiff >= (mz[:-1] * ppm / 1e6)
    res = np.insert(res, 0, 0)
    res = np.cumsum(res)
    return res


# def estimateMzScattering(x):
#     # Kernel density estimation of "mzdiff" using Silverman's rule-of-thumb
#     sigma = np.std(x)
#     q75, q25 = np.percentile(x, [75 ,25])
#     iqr = q75 - q25
#     n = len(x)
#     bw = 0.9 * min(sigma, iqr / 1.34) * (n ** (-1 / 5)) # Bandwidth for KDE using rule-of-thumb
#     kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(x.reshape(-1, 1))
#     xScores = np.linspace(min(x) - 3 * bw, max(x) + 3 * bw, 512)
#     scores = np.exp(kde.score_samples(xScores.reshape(-1, 1)))
#     idx = np.where(np.diff(np.sign(np.diff(scores))) == 2)[0] + 1
#     if len(idx) > 1:
#         idx = idx[0]
#     return xScores[idx].item()


def mergeMs2(mzs, ints, mzGroups):
    mzArray, intArray = [], []
    for i in np.unique(mzGroups):
        idx = np.where(mzGroups == i)[0]
        mz = np.average(mzs[idx], weights=ints[idx])
        intensity = sum(ints[idx])
        mzArray.append(mz)
        intArray.append(intensity)
    spec = {"mz": np.array(mzArray), "intensity": np.array(intArray)}
    return spec


def intraConsolidation(ms2, scans, tol):
    # input arguments
    # ms2: dictionary of MS2 spectra; key = scan number, val = {"mz": np.array, "intensity": np.array}
    # scans: array of "feature-to-scanNumber" [# features x # runs]
    #        scans[i, j] = list of MS2 scan(number)s of j-th run corresponding to i-th feature
    # tol: m/z tolerance for merging MS2 peaks

    scans = scans.split(";")
    if len(scans) > 1:
        mzs = np.array([])
        ints = np.array([])
        # Extract m/z (and intensity) values from individual MS2 spectrum and merge them to "mzs" (and "ints")
        for s in scans:
            mzs = np.append(mzs, ms2[s]["mz"])
            ints = np.append(ints, ms2[s]["intensity"])
        # Sort m/z values (and intensities accordingly)
        idx = np.argsort(mzs)
        mzs = mzs[idx]
        ints = ints[idx]
        mzGroups = groupMzValues(mzs, tol)
        spec = mergeMs2(mzs, ints, mzGroups)
    else:
        spec = ms2[scans[0]]
    return spec



def interConsolidation(specs, tol):
    # input arguments
    # specs: array of "feature-to-spectrum" [# features x # runs]
    #        specs[i, j] = (intra-consolidated) MS2 spectrum of j-th run corresponding to i-th feature
    # tol: m/z tolerance for merging MS2 peaks
    specs = [i for i in specs if i is not None]  # Skip "None"
    if len(specs) > 1:
        mzs = np.array([])
        ints = np.array([])
        # Extract m/z (and intensity) values from individual MS2 spectrum and merge them to "mzs" (and "ints")
        for s in specs:
            mzs = np.append(mzs, s["mz"])
            ints = np.append(ints, s["intensity"])
        # Sort m/z values (and intensities accordingly)
        idx = np.argsort(mzs)
        mzs = mzs[idx]
        ints = ints[idx]
        mzGroups = groupMzValues(mzs, tol)
        spec = mergeMs2(mzs, ints, mzGroups)
    else:
        spec = specs[0]
    spec = simplifyMs2(spec)
    return spec


def simplifyMs2(spec):
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
    logging.info("  Identification of MS2 spectra for the features")
    logging.info("  ==============================================")
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

    #################################################
    # Assignment of MS2 spectra to features         #
    # Consolidation of MS2 spectra for each feature #
    #################################################
    m = -1  # Index for input files
    for file in mzxmlFiles:
        m += 1
        reader = mzxml.MzXML(file)
        fileBasename, _ = os.path.splitext(os.path.basename(file))
        colNames = [item for item in full.dtype.names if re.search(fileBasename + "_", item)]
        subset = full[colNames]
        subset.dtype.names = [s.split("_")[-1] for s in subset.dtype.names]
        ms2Dict = {}
        minScan, maxScan = int(np.nanmin(subset["minMS1"])), int(np.nanmax(subset["maxMS1"]))
        progress = utils.progressBar(maxScan - minScan + 1)
        print("  %s is being processed" % os.path.basename(file))
        print("  Looking for MS2 scan(s) responsible for each feature")
        logging.info("  %s is being processed" % os.path.basename(file))
        logging.info("  Looking for MS2 scan(s) responsible for each feature")
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
        logging.info("  Merging MS2 spectra for each feature within a run (it may take a while)")
        progress = utils.progressBar(nFeatures)
        for i in range(nFeatures):
            progress.increment()
            if featureToScan[i, m] is not None:
                spec = intraConsolidation(ms2Dict, featureToScan[i, m], tolIntraMS2Consolidation)
                featureToSpec[i, m] = spec

    print("  Merging MS2 spectra for each feature between runs when there are multiple runs")
    print("  Simplification of MS2 spectrum for each feature by retaining the most strongest 100 peaks")
    logging.info("  Merging MS2 spectra for each feature between runs when there are multiple runs")
    logging.info("  Simplification of MS2 spectrum for each feature by retaining the most strongest 100 peaks")
    specArray = np.array([])
    progress = utils.progressBar(nFeatures)
    for i in range(nFeatures):
        progress.increment()
        if np.sum(featureToSpec[i] == None) == nFiles:
            specArray = np.append(specArray, None)
        else:
            spec = interConsolidation(featureToSpec[i, :], tolInterMS2Consolidation)
            specArray = np.append(specArray, spec)

    ###############################
    # MS2 processing for features #
    ###############################
    # "specArray" is the list of (consolidated) MS2 spectra
    # specArray[i] is the MS2 spectrum corresponding to the i-th feature
    # If there's no MS2 spectrum, then specArray[i] is None
    df = utils.summarizeFeatures(full, params)
    # Add the mean m/z of feature and its charge state to the beginning of MS2 spectrum (similar to .dta file)
    for i in range(nFeatures):
        if specArray[i] is not None:
            specArray[i]["mz"] = np.insert(specArray[i]["mz"], 0, df["feature_m/z"].iloc[i])
            specArray[i]["intensity"] = np.insert(specArray[i]["intensity"], 0, df["feature_z"].iloc[i])
    df["MS2"] = specArray
    df = df.sort_values(by="feature_m/z", ignore_index=True)  # Features are sorted by "feature_m/z"
    df.insert(loc=0, column="feature_num", value=df.index + 1)
    # df["feature_num"] = df.index + 1  # Update "feature_num" according to the ascending order of "feature_m/z" (as sorted)

    # Write MS2 spectra to files
    filePath = os.path.join(os.getcwd(), "align_" + params["output_name"])
    ms2Path = os.path.join(filePath, "MS2")
    if not os.path.exists(ms2Path):
        os.mkdir(ms2Path)
    for i in range(df.shape[0]):
        if df["MS2"].iloc[i] is not None:
            fileName = os.path.join(ms2Path, "f" + str(i + 1) + ".MS2")
            dfMS2 = pd.DataFrame.from_dict(df["MS2"].iloc[i])
            dfMS2.to_csv(fileName, index=False, header=False, sep="\t")

    # Save fully-aligned features with their MS2 spectra (i.e. res) for debugging purpose
    # When the pipeline gets mature, this part needs to be removed
    pickle.dump(df, open(os.path.join(filePath, ".fully_aligned_feature.pickle"), "wb"))    # Make the file be hidden

    ##########################
    # Handling mzXML file(s) #
    ##########################
    # Move mzXML files to the directory(ies) where individual .feature files are located
    if params["skip_feature_detection"] == "0":
        for file in mzxmlFiles:
            baseFilename = os.path.basename(file)
            featureDirectory = os.path.join(os.getcwd(), os.path.splitext(baseFilename)[0])
            os.rename(file, os.path.join(featureDirectory, baseFilename))

    return df, featureToScan
