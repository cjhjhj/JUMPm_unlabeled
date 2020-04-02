#!/usr/bin/python

import glob, os, sqlite3, numpy as np, pandas as pd

# Read a text file containing the information of metabolomes
txtFile = r"C:\Research\Projects\7Metabolomics\Library\Metabolome_library_v2.1.2.txt"
df = pd.read_csv(txtFile, sep = "\t", engine = "python")

# Open a sqlite database and create tables
conn = sqlite3.connect("library.db")

# Library entry table
colInd = np.where(df.columns == "InChIKey")[0][0]
entryDf = df.iloc[:, 0:(colInd + 1)]
# subDf.to_sql("libEntry", conn, if_exists = "replace")

# Tables for column conditions
conds = [i for i in df.columns if i.endswith("linkms2")]
conds = list(set([i.split("_")[0] for i in conds]))   # Unique column conditions
for cond in conds:
    subDf = df.filter(regex = cond)
    subDf = pd.concat([entryDf, subDf], axis = 1)
    subDf.to_sql(cond, conn, if_exists = "replace")

# Individual MS2 spectrum
files = glob.glob(r"C:\Research\Projects\7Metabolomics\Library\c18p\*.MS2")
for file in files:
    df = pd.read_csv(file, sep = "\t", engine = "python")
    df.columns = ["mz", "intensity"]
    df.to_sql(os.path.basename(file), conn, if_exists = "replace")

# Individual MS2 spectrum
files = glob.glob(r"C:\Research\Projects\7Metabolomics\Library\hilicn\*.MS2")
for file in files:
    df = pd.read_csv(file, sep = "\t", engine = "python")
    df.columns = ["mz", "intensity"]
    df.to_sql(os.path.basename(file), conn, if_exists = "replace")
# Close the sqlite database
conn.close()
