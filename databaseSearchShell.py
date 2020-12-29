import sys, os, re, pickle, subprocess, utils, pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
ro.r['options'](warn=-1)


def adductDictionary(params):
    adduct = {}
    for key, val in params.items():
        if key.startswith("adduct"):
            key = re.sub(r'adduct_', '', key)
            if key == "NH3":
                key = "NH4" # In MetFrag, NH3 is used as NH4
            adduct[key] = float(val)

    # In MetFrag, the mass-shift is integerized
    if params["mode"] == "1":
        for key, val in adduct.items():
            adduct[key] = int(val + 0.5) + 1
        adduct["H"] = 1
    elif params["mode"] == "-1":
        for key, val in adduct.items():
            adduct[key] = int(val + 0.5) - 1
        adduct["-H"] = -1

    return adduct


def generateFiles(feature, params, precursorIonMode):
    num = feature["feature_num"]
    paramFile = "metfrag_params_" + str(num) + ".txt"
    ms2File = "metfrag_data_" + str(num) + ".txt"
    outputName = "metfrag_result_" + str(num)
    outputFile = "metfrag_result_" + str(num) + ".csv"
    proton = 1.007276466812
    if params["mode"] == "1":
        mass = feature["feature_z"] * (feature["feature_m/z"] - proton) # Neutral (monoisotopic) mass
    elif params["mode"] == "-1":
        mass = feature["feature_z"] * (feature["feature_m/z"] + proton) # Neutral (monoisotopic) mass

    # Parameter file for MetFrag
    f = open(paramFile, "w")
    f.write("PeakListPath = {}\n".format(ms2File))

    # Even though this code is designed to use a database including PubChem, PubChem is not available in practice due
    # to the following reasons. 1) Many queries to PubChem through MetFrag are not stable -> problem of MetFrag (
    # according to authors) 2) An alternative of using a local CSV file is possible, but it takes too long because of
    # the size of file (> 30GB)
    if params["database"].lower() == "pubchem":
        f.write("MetFragDatabaseType = PubChem\n")
    else:
        if os.path.isfile(params["database"]):
            f.write("MetFragDatabaseType = LocalCSV\n")
            f.write("LocalDatabasePath = {}\n".format(params["database"]))
        else:
            sys.exit("Please check the path of a database file (.csv)")
    f.write("DatabaseSearchRelativeMassDeviation = {}\n".format(params["mass_tolerance_formula_search"]))
    f.write("FragmentPeakMatchRelativeMassDeviation = {}\n".format(params["mass_tolerance_ms2_peaks"]))
    f.write("NeutralPrecursorMass = {}\n".format(mass))
    f.write("PrecursorIonMode = {}\n".format(precursorIonMode))   # It may contain adduct information. Refer https://ipb-halle.github.io/MetFrag/projects/metfragcl/
    if params["mode"] == "1":
        f.write("IsPositiveIonMode = True\n")
    elif params["mode"] == "-1":
        f.write("IsPositiveIonMode = False\n")
    f.write("MetFragScoreTypes = FragmenterScore\n")
    f.write("MetFragScoreWeights = 1.0\n")
    f.write("MetFragCandidateWriter = CSV\n")
    f.write("SampleName = {}\n".format(outputName))
    f.write("ResultsPath = .\n")
    f.write("MaximumTreeDepth = 2\n")
    f.write("MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter\n")
    f.write("MetFragPostProcessingCandidateFilter = InChIKeyFilter\n")
    f.close()

    # MS2 data file for MetFrag
    ms2Dict = feature["MS2"]
    df = pd.DataFrame.from_dict(ms2Dict, orient = "columns")
    df = df.drop([0])
    df.to_csv(ms2File, sep = "\t", index = False, header = False)

    return paramFile, ms2File, outputFile


def lipidFrag():
    rstring = """
    rLipidFrag = function(inputFile, mode) {
        # generate models
        cwd = getwd()  
        setwd("/hpcf/authorized_apps/proteomics_apps/jumpm/python/lipidfrag/R")
        source("lipidfrag_train.r")
        source("lipidfrag_predict.r")

        # initialise the prediction models
        if (mode == "1") {
            models <- generate.model("/hpcf/authorized_apps/proteomics_apps/jumpm/python/lipidfrag/data/model_scores_pos.txt")
        } else {
            models <- generate.model("/hpcf/authorized_apps/proteomics_apps/jumpm/python/lipidfrag/data/model_scores_neg.txt")
        }

        # predict the lipid (main-/sub-)class for a given MetFrag annotation
        predicted = predict.lipidmaps.class(inputFile, models, sep=",")
        setwd(cwd)    
        return(predicted)
    }
    """
    return ro.r(rstring)


def runMetFrag(feature, params):
    if feature["MS2"] is not None:
        dfAll = pd.DataFrame()
        adducts = adductDictionary(params)
        for k, v in adducts.items():
            paramFile, ms2File, outputFile = generateFiles(feature, params, v)

            # MetFrag should be installed first and its path should be put to the following command
            cmd = "java -jar " + params["metfrag"] + " " + paramFile + "> /dev/null 2>&1" # "> /dev/null 2>&1" is for linux only
            subprocess.call(cmd, shell=True)    # Subprocess is recommended instead of os.system
            df = pd.read_csv(outputFile, keep_default_na=False)
            if not df.empty:
                # Run LipidFrag, if necessary
                if "lipidfrag" in params and params["lipidfrag"] == "1":
                    rLF = lipidFrag()
                    pred = rLF(os.path.abspath(outputFile), 1)
                    with localconverter(ro.default_converter + pandas2ri.converter):  # Conversio from rpy2 object to pandas dataframe
                        pred = ro.conversion.rpy2py(pred)
                    if not pred.empty:
                        df = pd.merge(df, pred[["Identifier", "LipidMapsClass"]], how="left", on=["Identifier", "Identifier"])

                # Organize the output dataframe
                if params["database"].lower() == "pubchem":
                    df = df.rename(columns = {"IUPACName": "CompoundName"})
                df["feature_index"] = feature["feature_num"]
                df["feature_m/z"] = feature["feature_m/z"]
                df["feature_RT"] = feature["feature_RT"]
                if k == "-H":
                    df["Ion"] = "[M" + str(k) + "]"
                else:
                    df["Ion"] = "[M+" + str(k) + "]"
                if params["mode"] == "1":
                    df["Ion"] = df["Ion"] + "+"
                elif params["mode"] == "-1":
                    df["Ion"] = df["Ion"] + "-"
                intensityCols = [col for col in feature.keys() if col.lower().endswith("_intensity")]
                for c in intensityCols:
                    df[c] = feature[c]

                # Formatting of output dataframe
                if params["lipidfrag"] == "1":
                    columns = ["feature_index", "feature_m/z", "feature_RT"] + intensityCols + \
                              ["Identifier", "OtherIDs(PubChem;ChEBI;KEGG;HMDB;SwissLipid;LipidBank;PlantFA)", "MolecularFormula",
                               "CompoundName", "SystematicName", "Synonyms", "Abbreviation", "Category", "MainClass", "SubClass",
                               "Ion", "SMILES", "InChIKey", "FragmenterScore", "LipidMapsClass"]
                else:
                    columns = ["feature_index", "feature_m/z", "feature_RT"] + intensityCols + \
                              ["Identifier", "MolecularFormula", "CompoundName", "Ion", "SMILES", "InChIKey", "FragmenterScore"]
                df = df[columns]
                dfAll = dfAll.append(df, ignore_index=True)
            os.remove(paramFile)
            os.remove(ms2File)
            os.remove(outputFile)
        return dfAll
    else:
        return None


########
# Main #
########
featureFile = sys.argv[1]
paramFile = sys.argv[2]
try:
    params = utils.getParams(paramFile)
    features = pickle.load(open(featureFile, "rb"))
except:
    sys.exit("Parameter file cannot be found or cannot be loaded")
res = pd.DataFrame()
for idx, row in features.iterrows():
    df = runMetFrag(row, params)
    res = res.append(df, ignore_index=True)
outputFile = os.path.splitext(featureFile)[0] + ".csv"
res.to_csv(outputFile, sep="\t", index=False, na_rep="NA")

