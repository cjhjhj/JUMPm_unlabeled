import re, sys, os, numpy as np, pandas as pd


def getParams(paramFile):
    parameters = dict()
    with open(paramFile, 'r') as file:
        for line in file:
            if re.search(r'^#', line) or re.search(r'^\s', line):
                continue
            line = re.sub(r'#.*', '', line)  # Remove comments (start from '#')
            line = re.sub(r'\s*', '', line)  # Remove all whitespaces

            # Exception for "feature_files" parameter
            if "feature_files" in parameters and line.endswith("feature"):
                parameters["feature_files"].append(line)
            else:
                key = line.split('=')[0]
                val = line.split('=')[1]
                if key == "feature_files":
                    parameters[key] = [val]
                else:
                    parameters[key] = val
    return parameters


def readFeatures(featureFile):
    with open(featureFile, 'r') as file:
        features = file.readlines()
    return features


class progressBar:
    def __init__(self, total):
        self.total = total
        self.barLength = 20
        self.count = 0
        self.progress = 0
        self.block = 0
        self.status = ""

    def increment(self):
        self.count += 1
        self.progress = self.count / self.total
        self.block = int(round(self.barLength * self.progress))
        if self.progress == 1:
            self.status = "Done...\r\n"
        else:
            self.status = ""
        #         self.status = str(self.count) + "/" + str(self.total)
        text = "\r  Progress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block),
                                                     int(self.progress * 100), self.status)
        sys.stdout.write(text)
        sys.stdout.flush()


def generateSummarizedFeatureFile(nFeatures, full, ms2, params):
    filePath = os.path.join(os.getcwd(), "align_" + params["output_name"])

    #############################################################
    # Summarize fully-aligned features and write them to a file #
    #############################################################
    # This file contains "summarized" information of fully-aligned features
    # e.g. mean m/z, mean intensity, mean RT of fully-aligned features and so on
    #      width and SNratio are from the reference run
    mzCol = [col for col in full.dtype.names if col.lower().endswith("_mz")]
    rtCol = [col for col in full.dtype.names if col.lower().endswith("_rt")]
    intensityCol = [col for col in full.dtype.names if col.lower().endswith("_intensity")]
    chargeCol = [col for col in full.dtype.names if col.lower().endswith("_z")]
    minRtCol = [col for col in full.dtype.names if col.lower().endswith("_minrt")]
    maxRtCol = [col for col in full.dtype.names if col.lower().endswith("_maxrt")]
    snCol = [col for col in full.dtype.names if col.lower().endswith("_snratio")]

    df = pd.DataFrame.from_records(full)
    res = pd.DataFrame()
    res["feature_m/z"] = df[mzCol].mean(axis=1)
    res["feature_RT"] = df[rtCol].mean(axis=1)
    res["feature_intensity"] = df[intensityCol].mean(axis=1)
    res["feature_SNratio"] = df[snCol].mean(axis=1)
    res["feature_width"] = pd.DataFrame((df[maxRtCol].values - df[minRtCol].values) / 60).mean(axis=1)
    # Handling charges
    res["feature_z"] = df[chargeCol].mode(axis=1)[0]
    res["feature_z"][res["feature_z"] == 0] = 1
    res["feature_z"] = res["feature_z"].astype(int)
    if params["mode"] == "-1":
        res["feature_ion"] = "[M-" + res["feature_z"].astype(str) + "H]" + res["feature_z"].astype(str) + "-"
        res["feature_ion"] = res["feature_ion"].replace("[M-1H]1-", "[M-H]-")
    elif params["mode"] == "1":
        res["feature_ion"] = "[M+" + res["feature_z"].astype(str) + "H]" + res["feature_z"].astype(str) + "+"
        res["feature_ion"] = res["feature_ion"].replace("[M+1H]1+", "[M+H]+")

    # Add the mean m/z of feature and its charge state to the beginning of MS2 spectrum (similar to .dta file)
    for i in range(nFeatures):
        if ms2[i] is not None:
            ms2[i]["mz"] = np.insert(ms2[i]["mz"], 0, res["feature_m/z"].iloc[i])
            ms2[i]["intensity"] = np.insert(ms2[i]["intensity"], 0, res["feature_z"].iloc[i])
    res["MS2"] = ms2

    # Write the summarized fully-aligned features to a file
    fullName = os.path.join(filePath, params["output_name"] + "_summarized_fully_aligned.feature")
    res = res.sort_values(by="feature_m/z", ignore_index=True)
    res["feature_num"] = res.index + 1  # Update "feature_num" according to the ascending order of "feature_m/z" (as sorted)
    resColumns = ["feature_num", "feature_ion", "feature_z", "feature_m/z", "feature_RT",
                 "feature_width", "feature_SNratio", "feature_intensity"]
    res.to_csv(fullName, columns=resColumns, index=False, sep="\t")

    # Write MS2 spectra to files
    filePath = os.path.join(filePath, "MS2")
    if not os.path.exists(filePath):
        os.mkdir(filePath)
    for i in range(res.shape[0]):
        if res["MS2"].loc[i] is not None:
            fileName = os.path.join(filePath, "f" + str(i + 1) + ".MS2")
            dfMS2 = pd.DataFrame.from_dict(res["MS2"].loc[i])
            dfMS2.to_csv(fileName, index=False, header=False, sep="\t")

    '''
    featureDict = {"feature_num": [], "feature_ion": [], "feature_z":[], "feature_m/z": [], "feature_RT": [],
                   "feature_width": [], "feature_SNratio": [], "feature_intensity": []}
    mzCol = [col for col in full.dtype.names if col.lower().endswith("_mz")]
    rtCol = [col for col in full.dtype.names if col.lower().endswith("_rt")]
    intensityCol = [col for col in full.dtype.names if col.lower().endswith("_intensity")]
    chargeCol = [col for col in full.dtype.names if col.lower().endswith("_z")]
    minRtCol = [col for col in full.dtype.names if col.lower().endswith("_minrt")]
    maxRtCol = [col for col in full.dtype.names if col.lower().endswith("_maxrt")]
    snCol = [col for col in full.dtype.names if col.lower().endswith("_snratio")]
    nRuns = len(mzCol)
    for i in range(nFeatures):
        # m/z, RT, intensity, width and SNratio of each feature are obtained by averaging information across runs
        nRuns =
        fMz = sum(full[mzCol][i]) / nRuns
        fRt = sum(full[rtCol][i]) / nRuns
        fIntensity = sum(full[intensityCol][i]) / nRuns
        fWidth = 0
        for j in range(nRuns):
            fWidth += full[maxRtCol[j]][i] - full[minRtCol[j]][i]
        fWidth = (fWidth / 60) / nRuns  # Unit of minute
        fSNratio = sum(full[snCol][i]) / nRuns

        # Heuristic charge determination across runs
        # 1. Charge state other than 0
        # 2. More frequent charge state
        # 3. If the same frequency, choose the lower one
        charges = [int(c) for c in full[chargeCol][i] if c > 0]
        if len(charges) == 0:
            charge = 1
        else:
            charge = max(set(charges), key=charges.count)
        fZ = charge
        if params["mode"] == "-1":  # Negative mode
            if charge > 1:
                fIon = "[M-" + str(charge) + "H]" + str(charge) + "-"
            else:
                fIon = "[M-H]-"
        if params["mode"] == "1":  # Positive mode
            if charge > 1:
                fIon = "[M+" + str(charge) + "H]" + str(charge) + "+"
            else:
                fIon = "[M+H]+"

        # Add the mean m/z of feature and its charge state to the beginning of MS2 spectrum (similar to .dta file)
        if ms2[i] is not None:
            ms2[i]["mz"] = np.insert(ms2[i]["mz"], 0, fMz)
            ms2[i]["intensity"] = np.insert(ms2[i]["intensity"], 0, charge)

        # Summarize featureDict (it is going to be used to make a pandas DataFrame)
        featureDict["feature_num"].append(i + 1)
        featureDict["feature_ion"].append(fIon)
        featureDict["feature_z"].append(charge)
        featureDict["feature_m/z"].append(fMz)
        featureDict["feature_RT"].append(fRt)
        featureDict["feature_width"].append(fWidth)
        featureDict["feature_SNratio"].append(fSNratio)
        featureDict["feature_intensity"].append(fIntensity)

    # Write the summarized fully-aligned features to a file
    fullName = os.path.join(filePath, params["output_name"] + "_summarized_fully_aligned.feature")
    df = pd.DataFrame.from_dict(featureDict)
    df["MS2"] = ms2
    df = df.sort_values(by = "feature_m/z", ignore_index = True)
    df["feature_num"] = df.index + 1    # Update "feature_num" according to the ascending order of "feature_m/z" (as sorted)
    dfColumns = ["feature_num", "feature_ion", "feature_m/z", "feature_RT",
                 "feature_width", "feature_SNratio", "feature_intensity"]
    df.to_csv(fullName, columns = dfColumns, index = False, sep = "\t")
    
    # Write MS2 spectra to files
    filePath = os.path.join(filePath, "MS2")
    if not os.path.exists(filePath):
        os.mkdir(filePath)

    for i in range(df.shape[0]):
        if df["MS2"].loc[i] is not None:
            fileName = os.path.join(filePath, "f" + str(i + 1) + ".MS2")
            dfMS2 = pd.DataFrame.from_dict(df["MS2"].loc[i])
            dfMS2.to_csv(fileName, index = False, header = False, sep = "\t")
    '''
    return df


def generateFeatureFile(full, partial, unaligned, params):
    # Input arguments
    # full: fully-aligned features (numpy recarray)
    # partial: partially-aligned features (numpy recarray, or None)
    # unaligned: un-aligned features (list of numpy recarray, or None)
    # params: dictionary of parameters
    filePath = os.path.join(os.getcwd(), "align_" + params["output_name"])
    if not os.path.exists(filePath):
        os.mkdir(filePath)

    ################################
    # Write fully-aligned features #
    ################################
    # This file contains fully-aligned features with run-specific information
    # Since the run-specific information is required for MS2 processing, it should be kept here
    fullName = os.path.join(filePath, params["output_name"] + "_fully_aligned.feature")
    dfFull = pd.DataFrame(full)
    dfFull["meanMz"] = dfFull.filter(regex=(".*mz$")).mean(axis=1)
    dfFull = dfFull.sort_values(by = "meanMz")  # To make it consistent with the "summarized feature" file, features are sorted by mean m/z
    colNames = dfFull.columns.tolist()
    colNames = colNames[-1:] + colNames[:-1]
    dfFull = dfFull[colNames]
    dfFull.to_csv(fullName, index = False, sep = "\t")

    ##############################################
    # Write partially-aligned features, if exist #
    ##############################################
    if partial is not None:
        dfPartial = pd.DataFrame(partial)
        if dfPartial.shape[0] > 0:
            partialName = os.path.join(filePath, params["output_name"] + "_partially_aligned.feature")
            dfPartial.to_csv(partialName, index = False, sep = "\t")
    else:
        dfPartial = None

    #######################################
    # Write un-aligned features, if exist #
    #######################################
    dfArrayUnaligned = None
    if unaligned is not None:
        dfArrayUnaligned = []
        for un in unaligned:
            unName = [col for col in un.dtype.names if col.endswith('_mz')][0]
            unName = os.path.join(filePath, re.sub("_mz", "", unName) + "_unaligned.feature")
            dfUn = pd.DataFrame(un)
            dfUn.to_csv(unName, index = False, sep = "\t")
            dfArrayUnaligned.append(dfUn)
    else:
        dfArrayUnaligned = None

    return dfFull, dfPartial, dfArrayUnaligned
