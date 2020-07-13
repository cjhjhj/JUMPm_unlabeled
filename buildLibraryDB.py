#!/usr/bin/python

import glob, os, sqlite3, utils, numpy as np, pandas as pd

# Path information
libPath = r"/Research/Projects/7Metabolomics/library/ours/"

# Read a text file containing the information of metabolomes
txtFile = libPath + r"Metabolome_library_v3.1.1.txt"
df = pd.read_csv(txtFile, sep = "\t", engine = "python")
df["type"] = pd.Series(["target"] * df.shape[0])

# Addition of decoys (by adding 3 * proton to the neutral mass)
proton = 1.007276466812
dfDecoy = df.copy()
for i in range(dfDecoy.shape[0]):
    dfDecoy.loc[i, "idstjude"] = dfDecoy.loc[i, "idstjude"] + "_decoy"
    dfDecoy.loc[i, "monoisotopic_mass"] += 3 * proton # This way prevents 'SettingwithCopyWarning'
    dfDecoy.loc[i, "type"] = "decoy"

# Merge target and decoy data frames into one
df = df.append(dfDecoy, ignore_index = True)

# Open a sqlite database and create tables
conn = sqlite3.connect("library.db")

# Master table for library compounds
df.to_sql("library", conn, if_exists = "replace")

# # Tables for column conditions
# colInd = np.where(df.columns == "InChIKey")[0][0]
# entryDf = df.iloc[:, 0:(colInd + 1)]
# conds = [i for i in df.columns if i.endswith("linkms2")]
# conds = list(set([i.split("_")[0] for i in conds]))   # Unique column conditions
# for cond in conds:
#     subDf = df.filter(regex = cond)
#     subDf = pd.concat([entryDf, subDf], axis = 1)
#     subDf.to_sql(cond, conn, if_exists = "replace")

# Individual MS2 spectrum
pathArray = df["c18p_linkms2"]
filePath = libPath + r"/c18p/*.MS2"
files = glob.glob(filePath)
progress = utils.progressBar(len(files))
for i in range(len(files)):
    progress.increment()
    file = files[i]
    dfMs2 = pd.read_csv(file, sep = "\t", engine = "python")
    dfMs2.columns = ["mz", "intensity"]
    tableName = pathArray.values[pathArray.str.contains(os.path.basename(file), regex = False)][0]
    tableName = str("'") + tableName + str("'")
    dfMs2.to_sql(tableName, conn, if_exists = "replace")

# Individual MS2 spectrum
pathArray = df["hilicn_linkms2"]
filePath = libPath + r"/hilicn/*.MS2"
files = glob.glob(filePath)
progress = utils.progressBar(len(files))
for i in range(len(files)):
    progress.increment()
    file = files[i]
    dfMs2 = pd.read_csv(file, sep = "\t", engine = "python")
    dfMs2.columns = ["mz", "intensity"]
    tableName = pathArray.values[pathArray.str.contains(os.path.basename(file), regex = False)][0]
    tableName = str("'") + tableName + str("'")
    dfMs2.to_sql(tableName, conn, if_exists = "replace")

# Close the sqlite database
conn.close()
