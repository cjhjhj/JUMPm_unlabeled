#!/usr/bin/python

import glob, os, sqlite3, utils, numpy as np, pandas as pd

# Read a text file containing the information of metabolomes
txtFile = r"/Research/Projects/7Metabolomics/Library/Metabolome_library_v3.1.1.txt"
df = pd.read_csv(txtFile, sep = "\t", engine = "python")

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
files = glob.glob(r"/Research/Projects/7Metabolomics/Library/c18p/*.MS2")
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
files = glob.glob(r"/Research/Projects/7Metabolomics/Library/hilicn/*.MS2")
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
