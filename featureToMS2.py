#!/usr/bin/python

import os, sys, re, pickle, utils, numpy as np
from pyteomics import mzxml

def ms2Consolidation(ms2, scans, tol, type):
    # inputReader = pyteomics reader object for a mzXML file
    # scans = MS2 scan number(s) for the feature
    # tol = tolerance for merging MS2 spectra
    # type = either "intra" (within a run) or "inter" (between runs)
    if type == "intra":
        # Sort MS2 spectra according to their total ion current (descending order)
        scans = scans.split(";")
        totIonCurrent = [ms2[key]["totIonCurrent"] for key in scans]
        ind = np.argmax(totIonCurrent)
        spec = ms2[scans[ind]] # MS2 spectrum with the highest total ion current
        del scans[ind]
    elif type == "inter":
        print ()
    else:
        sys.exit ("  MS2 consolidation type should be either 'intra' or 'inter'")

    if len(scans) > 1:
        for i in range(len(scans)):
            if type == "intra":
                p = ms2Dict[scans[i]]
            elif type == "inter":
                print ()
            else:
                sys.exit("  MS2 consolidation type should be either 'intra' or 'inter'")

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

    return spec


# For development purpose #################################################################
with open("featureToMS2_2.pickle", "rb") as f:
    vars = pickle.load(f)

featureFiles = vars[0]
full = vars[1]
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

nFeatures = full.shape[0]


'''
featureToScan = np.empty((nFeatures, nFiles), dtype = object)
featureToSpec = np.empty((nFeatures, nFiles), dtype = object)

########################################################
# Find MS2 scans responsible for features within a run #
########################################################
m = -1  # Index for input (mzXML) files
for file in mzXMLFiles:


    # File handling may need to be revised according to the format of wrapper ######

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

    # print ("  Reading %s" % os.path.basename(file))
    # reader = mzxml.read(file)
    # nScans = len(list(reader))  # Traverse the iterator (i.e. reader) to count the number of total spectra

    nScans = 45000

    ms2Dict = {}
    with mzxml.read(file) as reader:
        print ("  Finding MS2 spectra of %s responsible for features" % os.path.basename(file))
        progress = utils.progressBar(nScans)
        for spec in reader:
            progress.increment()
            msLevel = int(spec["msLevel"])
            if msLevel == 1:
                # surveyScanNumber = spec['num']
                survey = spec
            elif msLevel == 2:
                # Find MS2 scans which satisfy
                # 1. Their precursor m/z values are within feature's m/z +/- isolation_window
                # 2. Their presumable survey scans are between feature's min. and max. MS1 scan numbers
                # 3. Feature width should be less than a threshold
                precMz = float(spec["precursorMz"][0]["precursorMz"])
                fInd = np.where((full[colMinMS1] < int(survey["num"])) &
                                (full[colMaxMS1] > int(survey["num"])) &
                                (full[colMz] >= (precMz - tolIsolation)) &
                                (full[colMz] <= (precMz + tolIsolation)) &
                                (full[colTF] < pctTfThreshold))[0]
                # fInd = np.where((full[colMinMS1] < int(surveyScanNumber)) &
                #                 (full[colMaxMS1] > int(surveyScanNumber)) &
                #                 (full[colMz] >= (precMz - tolIsolation)) &
                #                 (full[colMz] <= (precMz + tolIsolation)) &
                #                 (full[colTF] < pctTfThreshold))[0]
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

                    # surveyMzArray = survey["m/z array"]
                    # surveyIntArray = survey["intensity array"]
                    # surveyScan = reader[surveyScanNumber]
                    ppi = []
                    for i in range(len(fInd)):
                        mz = full[colMz][fInd[i]]
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
                        ms2Dict[spec["num"]]["totIonCurrent"] = spec["totIonCurrent"]
                        # Mapping between features and MS2 scan numbers
                        for i in range(len(fInd)):
                            if featureToScan[fInd[i], m] is None:
                                featureToScan[fInd[i], m] = spec["num"]
                            else:
                                featureToScan[fInd[i], m] += ";" + spec["num"]

        ##############################################################
        # Consolidation of MS2 spectra for each feature within a run #
        ##############################################################
        print  ("  Merging MS2 spectra within each feature")
        progress = utils.progressBar(nFeatures)
        for i in range(nFeatures):
            progress.increment()
            if featureToScan[i, m] is None:
                continue
            else:
                spec = ms2Consolidation(ms2Dict, featureToScan[i, m], tolIntraMS2Consolidation, "intra")
                featureToSpec[i, m] = spec

'''
#############################################
# Consolidation of MS2 spectra between runs #
#############################################
print ("  Merging MS2 spectra between runs for each feature")

print ()



# # Development purpose ##################################
# with open("featureToMS2_2.pickle", "wb") as f:
#     pickle.dump([featureFiles, full, params, featureToScan], f)

