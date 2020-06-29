#!/usr/bin/python

import sys, os, sqlite3, utils, numpy as np, pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.vectors import FloatVector
from featureAlignment import loess
from statsmodels.distributions.empirical_distribution import ECDF


def calcMS2Similarity(featSpec, libSpec):
    # Calculation of MS2 similarity between a feature and a library compound
    # Reference: Clustering millions of tandem mass spectra, J Proteome Res. 2008; 7: 113-22

    # Input arguments
    # featSpec (dictionary): MS2 spectrum of a feature (key = "mz", "intensity")
    # libSpec (dictionary): MS2 spectrum of a library compound (key = "mz", "intensity", "index" (ignorable))

    k = min(30, min(len(featSpec["mz"]), len(libSpec["mz"])))

    # Keep $k strongest peaks in both spectra
    # featDict[mz] = intensity
    # libDict[mz] = intensity
    featDict, libDict = {}, {}
    ind = np.argsort([-i for i in featSpec["intensity"]])
    for i in ind[0:k]:
        featDict[featSpec["mz"][i]] = featSpec["intensity"][i]
    ind = np.argsort([-i for i in libSpec["intensity"]])
    for i in ind[0:k]:
        libDict[libSpec["mz"][i]] = libSpec["intensity"][i]

    # Join two sets of m/z values and make a new set of unique m/z values
    # Duplicate masses are removed as follows
    # - We consider two peaks to have a similar mass if they are within 0.5 Da from each other)
    # - For those two peaks having similar mass, the lower one will be the unique one
    #   (e.g. One peak with m/z = 100 and the other peak with m/z = 100.4 -> they will be merged to m/z = 100)
    mzArray = list(featDict.keys()) + list(libDict.keys())
    mzArray = sorted(mzArray)
    mzDict = {}
    val = 0
    for mz in mzArray:
        if abs(mz - val) <= 0.5:
            mzDict[mz] = val
        else:
            mzDict[mz] = mz
            val = mz

    # Reduction of spectrum to a vector by assigning to each intensity to the unique m/z bins
    # And then, calculate the similarity; normalized dot-product
    s = {}
    for key, val in mzDict.items():
        s[val] = {}
        s[val]["feat"] = 0
        s[val]["lib"] = 0

    for key, val in mzDict.items():
        if key in featDict:
            s[val]["feat"] += np.sqrt(featDict[key])
        if key in libDict:
            s[val]["lib"] += np.sqrt(libDict[key])

    den, num1, num2 = 0, 0, 0
    for mz in s.keys():
        den += s[mz]["feat"] * s[mz]["lib"]
        num1 += s[mz]["feat"] ** 2
        num2 += s[mz]["lib"] ** 2

    if num1 * num2 == 0:
        normDotProduct = 0
    else:
        normDotProduct = den / np.sqrt(num1 * num2)
    return normDotProduct


def searchLibrary(full, libFile, paramFile):
    ###################
    # Input arguments #
    ###################
    # full = fully-aligned features and their MS2 (consolidated) spectra (Pandas DataFrame)
    # libFile = SQLite-version of a library file (.db file)
    #           It should include entities and MS2 spectra of library compounds
    #           - Main table: entities
    #           - Sub tables: individual MS2 spectra
    # paramFile = parameter file (a text file)

    ##################################
    # Load parameters and initialize #
    ##################################
    try:
        params = utils.getParams(paramFile)
    except:
        sys.exit("Parameter file cannot be found or cannot be loaded")

    condition = params["LC_column"].lower()
    if params["mode"] == "1":
        condition = condition + "p"
    elif params["mode"] == "-1":
        condition = condition + "n"
    else:
        sys.exit("'mode' parameter should be either 1 or -1")
    proton = 1.007276466812
    matchMzTol = 10  # Unit of ppm

    #############################
    # Open sqlite-based library #
    #############################
    try:
        conn = sqlite3.connect(libFile)
    except:
        sys.exit("Library file cannot be found or cannot be loaded")

    command = "select * from library"   # Load a master table containing entities of library compounds
    library = pd.read_sql_query(command, conn)
    nCompounds = library.shape[0]
    library = library.replace("na", np.nan)  # If there's "na", replace it with np.nan
    library.columns = library.columns.str.replace(" ", "")
    libRt = library[condition + "_rt"].to_numpy(dtype="float")
    libRt = [rt * 60 for rt in libRt] # Note that library RT is using "minute"
    libM = library["monoisotopic_mass"].to_numpy(dtype="float")  # This "monoisotopic mass" is a neutral mass
    libZ = library[condition + "_charge"].to_numpy(dtype="float")
    if params["mode"] == "1":
        libMz = (libM + libZ * proton) / libZ
    elif params["mode"] == "-1":
        libMz = (libM - libZ * proton) / libZ

    ############################
    # Load feature information #
    ############################
    # Organize feature information
    nFeatures = full.shape[0]
    fMS2 = full["MS2"].tolist()
    mzCol = [col for col in full.columns if col.lower().endswith("_mz")]
    fMeanMz = (full[mzCol].mean(axis=1)).to_numpy()
    rtCol = [col for col in full.columns if col.lower().endswith("_rt")]
    fMeanRt = (full[rtCol].mean(axis=1)).to_numpy()
    intCol = [col for col in full.columns if col.lower().endswith("_intensity")]
    fMeanIntensity = (full[intCol].mean(axis=1)).to_numpy()
    fZ = np.zeros(len(fMS2))
    fZ[:] = np.nan
    for i in range(len(fMS2)):
        if fMS2[i] is not None:
            fZ[i] = fMS2[i]["charge"]

    #####################################################
    # RT-alignment between features and library entries #
    #####################################################
    # Preparation of LOESS-based RT alignment
    x, y = np.array([]), np.array([])  # To be used for RT-alignment
    for i in range(len(libMz)):
        if np.isnan(libRt[i]):
            continue
        refMz = libMz[i]
        refRt = libRt[i]
        refZ = libZ[i]
        compInd = np.nan

        # Previously, when multiple features correspond to a library compound, one with the largest intensity was chosen
        # In the revision, one with the closest m/z to the library compound is chosen
        minRtDiff = 600  # Arbitrary threshold of minimum RT-difference between a feature and library compounds
        for j in range(nFeatures):
            compMz = fMeanMz[j]
            compRt = fMeanRt[j]
            compZ = fZ[j]
            if np.isnan(compZ):
                continue
            mzDiff = abs(refMz - compMz) / refMz * 1e6
            rtDiff = abs(refRt - compRt)
            if mzDiff <= matchMzTol and rtDiff <= minRtDiff:
                if refZ == 0:
                    compInd = j
                    minRtDiff = rtDiff
                else:
                    if compZ == 0 or compZ == refZ:
                        compInd = j
                        minRtDiff = rtDiff
        if not np.isnan(compInd):
            x = np.append(x, fMeanRt[compInd])  # RT of features
            y = np.append(y, (fMeanRt[compInd] - refRt))  # RTshift

    # LOESS modeling
    rLoess = loess()
    rPredict = ro.r("predict")

    # Truncation of y (RT-shift)
    truncatedMean = np.mean(y[(y >= np.quantile(y, 0.1)) & (y <= np.quantile(y, 0.9))])
    truncatedSd = np.std(y[(y >= np.quantile(y, 0.1)) & (y <= np.quantile(y, 0.9))])
    lL = truncatedMean - 3 * truncatedSd
    uL = truncatedMean + 3 * truncatedSd
    ind = np.where((y >= lL) & (y <= uL))[0]

    # LOESS fitting
    mod = rLoess(FloatVector(x[ind]), FloatVector(y[ind]))  # LOESS between featureRT vs. RTshift
    fMeanRt = fMeanRt - rPredict(mod, FloatVector(fMeanRt))

    # Empirical CDF of alignment (absolute) residuals (will be used to calculate RT shift-based scores)
    ecdfRT = ECDF(abs(np.array(mod.rx2("residuals"))))


    ########################################
    # Match features and library compounds #
    ########################################
    # For library MS2 spectra, the path itself cannot be used as a table name in sqlite because of slashes
    # So, in library building, the table name is set to 'path' (bracketed by single quotations)
    # When retrieving, the command should be as follows,
    # SELECT * FROM '''/research/.../sjm00001p.MS2'''

    # Match features and library compounds
    res = np.zeros((nFeatures, nCompounds))
    res[:] = np.nan
    for i in range(nFeatures):
        if np.isnan(fZ[i]):  # When MS2 spectrum of the feature is not defined, skip it
            continue
        for j in range(nCompounds):
            if np.isnan(libRt[j]):  # When there's no experimental information of a compound, skip it
                continue
            if 0 < fZ[i] != libZ[j] > 0:
                continue
            mzDiff = abs(fMeanMz[i] - libMz[j]) / libMz[j] * 1e6
            if mzDiff < matchMzTol:
                # Calculate the similarity score between feature and library MS2 spectra
                featSpec = fMS2[i]
                sqlQuery = r"select * from " + "'''" + library[condition + "_linkms2"][j] + "'''"
                libSpec = pd.read_sql_query(sqlQuery, conn).to_dict(orient="list")
                pMS2 = 1 - calcMS2Similarity(featSpec, libSpec)   # p-value-like score
                pMS2 = max(np.finfo(float).eps, pMS2)

                # Calculate the (similarity?) score based on RT-difference
                pRT = ecdfRT(abs(fMeanRt[i] - libRt[j]))  # Also, p-value-like score
                pRT = max(np.finfo(float).eps, pRT)

                # Harmonic mean of two p(like)-value
                p = 1 / (0.5 / pMS2 + 0.5 / pRT)    # Harmonic mean with equal weights
                print(i, j, pMS2, pRT, p)
    conn.close()
    print()
