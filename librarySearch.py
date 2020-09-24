#!/usr/bin/python

import sys, os, re, sqlite3, utils, numpy as np, pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.vectors import FloatVector
from featureAlignment import loess
from statsmodels.distributions.empirical_distribution import ECDF
from scipy import stats


def calcMS2Similarity(featSpec, libSpec, params):
    # Calculation of MS2 similarity between a feature and a library compound
    # Reference: Clustering millions of tandem mass spectra, J Proteome Res. 2008; 7: 113-22

    # Input arguments
    # featSpec (dictionary): MS2 spectrum of a feature (key = "mz", "intensity")
    # libSpec (dictionary): MS2 spectrum of a library compound (key = "mz", "intensity", "index" (ignorable))
    nPeaks = int(params["num_peaks_ms2_similarity"])    # Default = 30 according to the above reference
    k = min(nPeaks, min(len(featSpec["mz"]), len(libSpec["mz"])))

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


def searchLibrary(full, paramFile):
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
    nFeatures = full.shape[0]

    #############################
    # Open sqlite-based library #
    #############################
    print("  Library is being loaded")
    libFile = params["library"]
    try:
        conn = sqlite3.connect(libFile)
    except:
        sys.exit("Library file cannot be found or cannot be loaded.")

    #####################################################
    # RT-alignment between features and library entries #
    #####################################################
    # Check whether 'rt' column of the library is numeric value or not
    if params["library_rt_alignment"] == "1":
        preQuery = r"SELECT rt FROM library ORDER BY ROWID ASC LIMIT 1"
        preDf = pd.read_sql_query(preQuery, conn)
        if "float" not in str(preDf["rt"].dtype):
            params["library_rt_alignment"] = "2"

    if params["library_rt_alignment"] == "0":
        print("  According to the parameter, RT-alignment is not performed between features and library compounds")
    elif params["library_rt_alignment"] == "2":
        print("  Although the parameter is set to perform RT-alignment against the library, there is/are non-numeric value(s) in the library")
        print("  Therefore, RT-alignment is not performed")
        params["library_rt_alignment"] = "0"
    elif params["library_rt_alignment"] == "1":
        print("  RT-alignment is being performed between features and library compounds")
        # Preparation of LOESS-based RT alignment
        x, y = np.array([]), np.array([])  # To be used for RT-alignment
        for i in range(nFeatures):
            compMz = full["feature_m/z"].iloc[i]
            compRt = full["feature_RT"].iloc[i]
            compZ = full["feature_z"].iloc[i]
            if params["mode"] == "1":  # Positive mode
                compMass = compZ * (compMz - proton)
            elif params["mode"] == "-1":  # Negative mode
                compMass = compZ * (compMz + proton)
            if np.isnan(compZ):
                continue
            if compZ == 0:
                sqlQuery = r"SELECT rt FROM library WHERE abs((?) - mass) / mass * 1e6 < (?) AND abs((?) - rt) < 600"
                df = pd.read_sql_query(sqlQuery, conn, params=(compMass, matchMzTol, compRt))
            else:
                sqlQuery = r"SELECT rt FROM library WHERE abs((?) - mass) / mass * 1e6 < (?) AND abs((?) - rt) < 600 AND charge = (?)"
                # type(compZ) = numpy int64
                # For some reasons, "numpy int64" cannot be directly used to SQL statement
                # (I don't understand why. It seems that "numpy float64" works in SQL statement)
                # To avoid this problem, in the following SQL query, compZ is converted into "int" type again
                df = pd.read_sql_query(sqlQuery, conn, params=(compMass, matchMzTol, compRt, int(compZ)))
            if not df.empty:
                refRt = min(df["rt"])
                x = np.append(x, compRt)  # RT of features
                y = np.append(y, compRt - refRt)  # RT-shift

        # LOESS modeling
        rLoess = loess()
        rPredict = ro.r("predict")

        # Truncation of y (RT-shift)
        truncatedMean = np.mean(y[(y >= np.quantile(y, 0.1)) & (y <= np.quantile(y, 0.9))])
        truncatedSd = np.std(y[(y >= np.quantile(y, 0.1)) & (y <= np.quantile(y, 0.9))])
        lL = truncatedMean - 3 * truncatedSd
        uL = truncatedMean + 3 * truncatedSd
        ind = np.where((y >= lL) & (y <= uL))[0]

        # LOESS fitting and calibrate feature RT
        mod = rLoess(FloatVector(x[ind]), FloatVector(y[ind]))  # LOESS between featureRT vs. RTshift
        full["feature_RT"] = full["feature_RT"] - rPredict(mod, FloatVector(full["feature_RT"]))

        # Empirical CDF of alignment (absolute) residuals (will be used to calculate RT shift-based scores)
        ecdfRt = ECDF(abs(np.array(mod.rx2("residuals"))))

    ########################################
    # Match features and library compounds #
    ########################################
    # Match features and library compounds
    print("  Features are being compared with library compounds")
    res = {"no": [], "feature_index": [], "feature_m/z": [], "feature_RT": [], "feature_intensity": [],
           "id": [], "formula": [], "name": [], "SMILES": [], "InchiKey": [], "collision_energy": [],
           "RT_shift": [], "RT_score": [], "MS2_score": [], "combined_score": []}
    n = 0
    progress = utils.progressBar(nFeatures)
    for i in range(nFeatures):
        progress.increment()
        # Feature information
        fZ = full["feature_z"].iloc[i]
        fSpec = full["MS2"].iloc[i]
        if np.isnan(fZ) or fSpec is None:  # When MS2 spectrum of the feature is not defined, skip it
            continue
        fMz = full["feature_m/z"].iloc[i]
        fRt = full["feature_RT"].iloc[i]
        fIntensity = full["feature_intensity"].iloc[i]
        if params["mode"] == "1":  # Positive mode
            fMass = fZ * (fMz - proton)
        elif params["mode"] == "-1":  # Negative mode
            fMass = fZ * (fMz + proton)

        # Retrieve library compounds satisfying conditions and calculate MS2-based and RT-based similarity (if exist)
        sqlQuery = r"SELECT * FROM library WHERE abs(((?) - mass) / mass * 1e6) < (?)"
        df = pd.read_sql_query(sqlQuery, conn, params=(fMass, matchMzTol))
        if not df.empty:
            for j in range(df.shape[0]):
                uid = df["id"].iloc[j]
                uid = uid.replace("##Decoy_", "")
                sqlQuery = r"SELECT * FROM {}".format(uid)
                try:
                    libSpec = pd.read_sql_query(sqlQuery, conn)
                except:
                    continue
                if not libSpec.empty:
                    n += 1
                    # Calculate the score based on MS2 spectrum
                    libSpec = libSpec.to_dict(orient="list")
                    simMs2 = calcMS2Similarity(fSpec, libSpec, params)
                    pMs2 = 1 - simMs2  # p-value-like score (the smaller, the better)
                    pMs2 = max(np.finfo(float).eps, pMs2)   # Prevent the underflow caused by 0

                    # Calculate the (similarity?) score based on RT-shift
                    if params["library_rt_alignment"] == "1":
                        rtShift = fRt - df["rt"].iloc[j]
                        pRt = ecdfRt(abs(rtShift))  # Also, p-value-like score (the smaller, the better)
                        pRt = max(np.finfo(float).eps, pRt)
                        simRt = 1 - pRt
                        # p = 1 / (0.5 / pMS2 + 0.5 / pRt)  # Combined p-value using harmonic mean with equal weights
                        p = 1 - stats.chi2.cdf(-2 * (np.log(pMs2) + np.log(pRt)), 4)    # Fisher's method
                        # p = -2 * (np.log(pMs2) + np.log(pRt))   # Fisher's method used in Perl pipeline (the smaller, the better)
                    else:
                        rtShift = None
                        pRt = 1
                        simRt = 1 - pRt
                        p = pMs2

                    # Output
                    libId = df["id"].iloc[j]
                    libFormula = df["formula"].iloc[j]
                    libName = df["name"].iloc[j]
                    libSmiles = df["smiles"].iloc[j]
                    libInchiKey = df["inchikey"].iloc[j]
                    libEnergy = df["collision_energy"].iloc[j]

                    res["no"].append(n)
                    res["feature_index"].append(i + 1)
                    res["feature_m/z"].append(fMz)
                    res["feature_RT"].append(fRt / 60)
                    res["feature_intensity"].append(fIntensity)
                    res["id"].append(libId)
                    res["formula"].append(libFormula)
                    res["name"].append(libName)
                    res["SMILES"].append(libSmiles)
                    res["InchiKey"].append(libInchiKey)
                    res["collision_energy"].append(libEnergy)
                    res["RT_shift"].append(rtShift)
                    # res["RT_score"].append(abs(-np.log10(pRt)))  # Scores are transformed by -log10
                    # res["MS2_score"].append(abs(-np.log10(pMS2)))

                    # Haiyan's preference
                    # RT_score and MS2_score: 0 ~ 1 (bad to good)
                    res["RT_score"].append(simRt)
                    res["MS2_score"].append(simMs2)
                    res["combined_score"].append(abs(-np.log10(p)))

    conn.close()
    res = pd.DataFrame.from_dict(res)
    filePath = os.path.join(os.getcwd(), "align_" + params["output_name"])
    outputFile = os.path.join(filePath, "align_" + params["output_name"] + ".library_matches")
    res.to_csv(outputFile, sep = "\t", index = False)

    return res
