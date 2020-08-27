import os, sys, time, pandas as pd, numpy as np


def generateFiles(feature):
    num = feature["feature_num"]
    paramFile = "metfrag_params_" + str(num) + ".txt"
    ms2File = "metfrag_data_" + str(num) + ".txt"
    outputName = "metfrag_result_" + str(num)
    outputFile = "metfrag_result_" + str(num) + ".csv"


    proton = 1.007276466812
    mass = feature["feature_z"] * (feature["feature_m/z"] - proton) # Neutral (monoisotopic) mass

    # Parameter file for MetFrag
    f = open(paramFile, "w")
    f.write("PeakListPath = {}\n".format(ms2File))
    f.write("MetFragDatabaseType = LocalCSV\n")
    f.write("LocalDatabasePath = /Research/Projects/7Metabolomics/Database/HMDB/hmdb_metabolites.csv\n")
    f.write("NeutralPrecursorMass = {}\n".format(mass))
    f.write("DatabaseSearchRelativeMassDeviation = 10\n")
    f.write("FragmentPeakMatchRelativeMassDeviation = 5\n")
    f.write("PrecursorIonMode = 1\n")
    f.write("IsPositiveIonMode = True\n")
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
    df = pd.DataFrame.from_dict(ms2Dict, orient="columns")
    df = df.drop([0])
    df.to_csv(ms2File, sep="\t", index=False, header=False)

    return paramFile, ms2File, outputFile


def generateParams(feature):
    num = feature["feature_num"]
    ms2File = "metfrag_data_" + num + ".txt"
    f = open("metfrag_params.txt", "w")
    f.write("PeakListPath = metfrag_data.txt\n")
    f.write("LocalDatabasePath = /Research/Projects/7Metabolomics/Database/HMDB/hmdb_metabolites.csv\n")
    f.write("NeutralPrecursorMass = {}\n".format(mass))
    f.write("DatabaseSearchRelativeMassDeviation = 10\n")
    f.write("FragmentPeakMatchRelativeMassDeviation = 5\n")
    f.write("PrecursorIonMode = 1\n")
    f.write("IsPositiveIonMode = True\n")
    f.write("MetFragScoreTypes = FragmenterScore\n")
    f.write("MetFragScoreWeights = 1.0\n")
    f.write("MetFragCandidateWriter = CSV\n")
    f.write("SampleName = metfrag_result\n")
    f.write("ResultsPath = .\n")
    f.write("MaximumTreeDepth = 2\n")
    f.write("MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter\n")
    f.write("MetFragPostProcessingCandidateFilter = InChIKeyFilter\n")
    f.close()


def generateMS2(feature):
    ms2Dict = feature["MS2"]
    num = feature["feature_num"]
    ms2File = "metfrag_data_" + num + ".txt"
    df = pd.DataFrame.from_dict(ms2Dict, orient="columns")
    df = df.drop([0])
    df.to_csv(ms2File, sep="\t", index=False, header=False)
    return ms2File


def run(feature):
    if feature["MS2"] is not None:
        # # Initialization for MetFrag
        # csvFile = r"/Research/Projects/7Metabolomics/Database/HMDB/hmdb_metabolites.csv"
        # matchMzTol = 10
        # proton = 1.007276466812
        #
        # # Extract feature information and generate files for MetFrag
        # fZ = feature["feature_z"]
        # fMz = feature["feature_m/z"]
        # fMass = fZ * (fMz - proton)  # Neutral (monoisotopic) mass
        # ms2File = generateMS2(feature)
        # paramFile = generateParams(fMass)
        #
        # # Execute MetFrag through system command
        # cmd = "java -jar ../../MetFrag/MetFrag2.4.5-CL.jar metfrag_params.txt"
        # os.system(cmd)
        # df = pd.read_csv("metfrag_result.csv")

        paramFile, ms2File, outputFile = generateFiles(feature)
        cmd = "java -jar ../../MetFrag/MetFrag2.4.5-CL.jar " + paramFile
        os.system(cmd)
        time.sleep(0.1)
        df = pd.read_csv(outputFile)
        os.remove(paramFile)
        os.remove(ms2File)
        os.remove(outputFile)

        return df
