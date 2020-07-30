#!/usr/bin/python

import os, sqlite3, pandas as pd

# Custom script to generate our own metabolite library
# Template is a text file containing compound metadata
# MS2 spectrum of each metabolite (when available) is stored in a separate text file (.MS2)

# Initialization (path of library template file and experimental condition)
templateFile = r"/Research/Projects/7Metabolomics/library/StJude/Metabolome_library_v3.1.1.txt"
condition = "c18p"    # Column name and ion mode, e.g. hilicn = HILIC column with negative ion mode

# Read a text file containing the information of metabolomes
df = pd.read_csv(templateFile, sep = "\t", engine = "python")
df["type"] = pd.Series(["target"] * df.shape[0])

# Addition of decoys (by adding 3 * proton to the neutral mass)
proton = 1.007276466812
dfDecoy = df.copy()
for i in range(dfDecoy.shape[0]):
    dfDecoy.loc[i, "idstjude"] = "##Decoy_" + dfDecoy.loc[i, "idstjude"]
    dfDecoy.loc[i, "monoisotopic_mass"] += 3 * proton # This way prevents 'SettingwithCopyWarning'
    dfDecoy.loc[i, "type"] = "decoy"

# Merge target and decoy DataFrames into one
df = df.append(dfDecoy, ignore_index = True)

# Library of the input condition
compoundCols = ['idstjude', 'idkegg', 'idhmdb', 'PC_CID', 'PC_SID', 'CHEBI', 'idmetlin',
                 'name', 'synonym', 'formula', 'monoisotopic_mass', 'CAS', 'SMILES', 'InChIKey']
spectrumCols = [col for col in df.columns if col.lower().startswith(condition)]
dfLib = pd.concat([df[compoundCols], df[spectrumCols]], axis = 1)    # Similar to cbind in R

# Save the library DataFrame to a file
dbName = "stjude_library_" + condition + ".db"
dbName = os.path.join(os.path.dirname(templateFile), dbName)
conn = sqlite3.connect(dbName)
dfLib.to_sql("compoundTable", conn, if_exists = "replace")

# Processing of MS2 spectrum
ms2Array = list()
for i in range(dfLib.shape[0]):
    colName = condition + "_linkms2"
    ms2Path = df[colName].iloc[i]
    if ms2Path != "na":
        # Open .MS2 file, write peaks to a DataFrame and save the DataFrame to the database
        # Although the path of .MS2 file contains a letter representing an ion mode (either "p" or "n"),
        # "tableName" does not have the letter to make it consistent with "idstjude"
        tableName = tableName = os.path.splitext(os.path.basename(ms2Path))[0][:-1]
        dfMs2 = pd.read_csv(ms2Path, sep = "\t", engine = "python")
        dfMs2.columns = ["mz", "intensity"]
        dfMs2.to_sql(tableName, conn, if_exists = "replace")

conn.close()
