#!/usr/bin/python

import sys, os, sqlite3, utils, numpy as np, pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.vectors import FloatVector
from featureAlignment import loess


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
    # libFile = SQLite-version of library file (.db file)
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
    matchMzTol = 20  # Unit of ppm
    matchRtTol = 10  # Unit of second

    #############################
    # Open sqlite-based library #
    #############################
    try:
        conn = sqlite3.connect(libFile)
    except:
        sys.exit("Library file cannot be found or cannot be loaded")

    command = "select * from library"
    library = pd.read_sql_query(command, conn)
    nCompounds = library.shape[0]
    library = library.replace("na", np.nan)  # If there's "na", replace it with np.nan
    libRt = library[condition + "_rt"].to_numpy(dtype="float")
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
    ms2Array = full["MS2"].tolist()
    mzCol = [col for col in full.columns if col.endswith("_mz")]
    fMeanMz = (full[mzCol].mean(axis=1)).to_numpy()
    rtCol = [col for col in full.columns if col.endswith("_RT")]
    fMeanRt = (full[rtCol].mean(axis=1)).to_numpy()
    intCol = [col for col in full.columns if col.endswith("_Intensity")]
    fMeanIntensity = (full[intCol].mean(axis=1)).to_numpy()
    fZ = np.zeros(len(ms2Array))
    fZ[:] = np.nan
    for i in range(len(ms2Array)):
        if ms2Array[i] is not None:
            fZ[i] = ms2Array[i]["charge"]

    #####################################################
    # RT-alignment between features and library entries #
    #####################################################
    # Preparation of LOESS-based RT alignment
    x, y = np.array([]), np.array([])  # To be used for RT-alignment
    for i in range(len(libMz)):
        if np.isnan(libRt[i]):
            continue
        refMz = libMz[i]
        refRt = libRt[i] * 60  # Convert minute to second
        refZ = libZ[i]
        mzDiff = abs(fMeanMz - refMz) / refMz * 1e6
        if refZ == 0:
            ind = np.where((mzDiff <= matchMzTol) & ~(np.isnan(fZ)))[0]
        else:
            ind = np.where((mzDiff <= matchMzTol) & (fZ == refZ) & ~(np.isnan(fZ)))[0]

        if len(ind) == 0:
            continue
        else:
            if len(ind) > 1:
                maxInd = np.argmax(fMeanIntensity[ind])
                ind = ind[maxInd]
            else:
                ind = ind[0]
            x = np.append(x, fMeanRt[ind])  # RT of features
            y = np.append(y, (fMeanRt[ind] - refRt))  # RTshift

    # LOESS modeling
    rLoess = loess()
    rPredict = ro.r("predict")

    # Truncation of y (RT-shift)
    truncatedMean = np.mean(y[(y >= np.quantile(y, 0.1)) & (y <= np.quantile(y, 0.9))])
    truncatedSd = np.std(y[(y >= np.quantile(y, 0.1)) & (y <= np.quantile(y, 0.9))])
    lL = truncatedMean - 3 * truncatedSd
    uL = truncatedMean + 3 * truncatedSd
    ind = np.where((y >= lL) & (y <= uL))[0]
    mod = rLoess(FloatVector(x[ind]), FloatVector(y[ind]))  # LOESS between featureRT vs. RTshift
    fNewRt = fMeanRt - rPredict(mod, FloatVector(fMeanRt))

    ########################################
    # Match features and library compounds #
    ########################################
    # For library MS2 spectra, the path itself cannot be used as a table name in sqlite because of slashes
    # So, in library building, the table name is set to 'path' (bracketed by single quotations)
    # Example)
    # SELECT * FROM /research/.../sjm00001p.MS2 -> Not okay
    # SELECT * FROM '/research/.../sjm00001p.MS2' -> Okay

    # Match features and library compounds
    res = np.zeros((nFeatures, nCompounds))
    res[:] = np.nan
    rtDiffArray = []
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
                # Calculate the similarity between feature and library MS2 spectra
                featSpec = ms2Array[i]
                sqlQuery = r"select * from " + "'" + library[condition + "_linkms2"][j] + "'"
                libSpec = pd.read_sql_query(sqlQuery, conn).to_dict(orient="list")
                sim = calcMS2Similarity(featSpec, libSpec)
                print(i, j, sim)

                # # If the feature has already "well" matched to a library compound, then skip
                # if ~np.isnan(res[i, j]) and res[i, j] >= sim:
                #     continue
    conn.close()
    print()


'''
##################################
# Load parameters and initialize #
##################################
paramFile = r"C:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\jumpm_negative_desktop.params"
# paramFile = r"C:\Research\Projects\7Metabolomics\Dev\JUMPm_unlabel_python\librarySearch\jumpm_positive.params"
params = utils.getParams(paramFile)
condition = params["LC_column"].lower()
if params["mode"] == "1":
    condition = condition + "p"
elif params["mode"] == "-1":
    condition = condition + "n"
else:
    sys.exit("'mode' parameter should be either 1 or -1")
proton = 1.007276466812
H = 1.0078250321
matchMzTol = 10  # Unit of ppm
matchRtTol = 10  # Unit of second

#############################
# Open sqlite-based library #
#############################
conn = sqlite3.connect("library.db")
command = "select * from library"
library = pd.read_sql_query(command, conn)
nCompounds = library.shape[0]
library = library.replace("na", np.nan)  # If there's "na", replace it with np.nan
libRt = library[condition + "_rt"].to_numpy(dtype = "float")
libM = library["monoisotopic_mass"].to_numpy(dtype = "float")  # This "monoisotopic mass" is a neutral mass
libZ = library[condition + "_charge"].to_numpy(dtype = "float")
if params["mode"] == "1":
    libMz = (libM + libZ * proton) / libZ
elif params["mode"] == "-1":
    libMz = (libM - libZ * proton) / libZ

############################
# Load feature information #
############################
# Features and their MS2 spectra need to be loaded
data = pickle.load(open("featuresAndMS2_full_ms2Array.pickle", "rb"))
full = data[0]
ms2Array = data[1]

# full = pd.read_csv("./librarySearch/test_fully_aligned.feature", sep = "\t")

# Organize feature information
nFeatures = full.shape[0]
mzCol = [col for col in full.columns if col.endswith("_mz")]
fMeanMz = (full[mzCol].mean(axis = 1)).to_numpy()
rtCol = [col for col in full.columns if col.endswith("_RT")]
fMeanRt = (full[rtCol].mean(axis = 1)).to_numpy()
intCol = [col for col in full.columns if col.endswith("_Intensity")]
fMeanIntensity = (full[intCol].mean(axis = 1)).to_numpy()
fZ = np.ones(full.shape[0])
fZ = np.zeros(len(ms2Array))
fZ[:] = np.nan
for i in range(len(ms2Array)):
    if ms2Array[i] is not None:
        fZ[i] = ms2Array[i]["charge"]

#####################################################
# RT-alignment between features and library entries #
#####################################################
# Preparation of LOESS-based RT alignment
x, y = np.array([]), np.array([])  # To be used for RT-alignment
for i in range(len(libMz)):
    if np.isnan(libRt[i]):
        continue
    refMz = libMz[i]
    refRt = libRt[i]* 60    # Convert minute to second
    refZ = libZ[i]
    mzDiff = abs(fMeanMz - refMz) / refMz * 1e6
    if refZ == 0:
        ind = np.where((mzDiff <= matchMzTol) & ~(np.isnan(fZ)))[0]
    else:
        ind = np.where((mzDiff <= matchMzTol) & (fZ == refZ) & ~(np.isnan(fZ)))[0]

    if len(ind) == 0:
        continue
    else:
        if len(ind) > 1:
            maxInd = np.argmax(fMeanIntensity[ind])
            ind = ind[maxInd]
        else:
            ind = ind[0]
        x = np.append(x, fMeanRt[ind])  # RT of features
        y = np.append(y, (fMeanRt[ind] - refRt))  # RTshift

# LOESS modeling
rLoess = loess()
rPredict = ro.r("predict")

# Truncation of y (RT-shift)
truncatedMean = np.mean(y[(y >= np.quantile(y, 0.1)) & (y <= np.quantile(y, 0.9))])
truncatedSd = np.std(y[(y >= np.quantile(y, 0.1)) & (y <= np.quantile(y, 0.9))])
lL = truncatedMean - 3 * truncatedSd
uL = truncatedMean + 3 * truncatedSd
ind = np.where((y >= lL) & (y <= uL))[0]
mod = rLoess(FloatVector(x[ind]), FloatVector(y[ind]))  # LOESS between featureRT vs. RTshift
fNewRt = fMeanRt - rPredict(mod, FloatVector(fMeanRt))

########################################
# Match features and library compounds #
########################################
# For library MS2 spectra, the path itself cannot be used as a table name in sqlite because of slashes
# So, in library building, the table name is set to 'path' (bracketed by single quotations)
# Example)
# SELECT * FROM /research/.../sjm00001p.MS2 -> Not okay
# SELECT * FROM '/research/.../sjm00001p.MS2' -> Okay

# Match features and library compounds
res = np.zeros((nFeatures, nCompounds))
res[:] = np.nan
rtDiffArray = []
for i in range(nFeatures):
    if np.isnan(fZ[i]): # When MS2 spectrum of the feature is not defined, skip it
        continue
    for j in range(nCompounds):
        if np.isnan(libRt[j]):  # When there's no experimental information of a compound, skip it
            continue
        if 0 < fZ[i] != libZ[j] > 0:
            continue
        mzDiff = abs(fMeanMz[i] - libMz[j]) / libMz[j] * 1e6
        if mzDiff < matchMzTol:
            # Calculate the similarity between feature and library MS2 spectra
            featSpec = ms2Array[i]
            sqlQuery = r"select * from " + "'" + library[condition + "_linkms2"][j] + "'"
            libSpec = pd.read_sql_query(sqlQuery, conn).to_dict(orient = "list")
            sim = calcMS2Similarity(featSpec, libSpec)
            print(i, j, sim)

            # # If the feature has already "well" matched to a library compound, then skip
            # if ~np.isnan(res[i, j]) and res[i, j] >= sim:
            #     continue
conn.close()
print()
'''
