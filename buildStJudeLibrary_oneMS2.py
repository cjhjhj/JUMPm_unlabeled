#!/usr/bin/python

import os, sqlite3, pandas as pd

# Custom script to generate our own metabolite library
# Template is a text file containing compound metadata
# MS2 spectrum of each metabolite (when available) is stored in a separate text file (.MS2)

################################
# Handling of library template #
################################
# Initialization (path of library template file and experimental condition)
templateFile = r"/Research/Projects/7Metabolomics/library/StJude/Metabolome_library_v3.1.1.txt"
condition = "c18p"    # Column name and ion mode, e.g. hilicn = HILIC column with negative ion mode
dbName = "stjude_library_" + condition + ".db"
dbName = os.path.join(os.path.dirname(templateFile), dbName)
conn = sqlite3.connect(dbName)

# Read a text file containing the information of metabolomes
df = pd.read_csv(templateFile, sep = "\t", engine = "python")

##################################
# Preparation of a library table #
##################################
# Columns will be as follows
# Compound metadata
# 1. id: unique key of a compound/metabolite
# 2. other_ids: IDs of the compound/metabolite from other public databases
# 3. name: official name of the compound/metabolite
# 4. synonym: synonym of the compound/metabolite
# 5. formula: chemical formula
# 6. mass: monoisotopic neutral mass of the compound/metabolite
# 7. smiles: SMILES of the compound/metabolite (if any)
# 8. inchikey: InChIKey of the compound/metabolite (if any)
# Spectrum metadata
# 9. rt: RT of the compound/metabolite (if any)
# 10. charge: precursor charge state of the compound/metabolite (if any)
# 11. energy: collision energy of MS2

# Compound metadata
colNameOtherIds = "other_ids(KEGG;HMDB;PubChem_CID;PubChem_SID;ChEBI;METLIN;CAS)"    # Handling "other_ids"
otherIds = df['idkegg'].map(str) + ';' + df['idhmdb'].map(str) + ';' + df['PC_CID'].map(str) + ";" + \
df["PC_SID"].map(str) + ";" + df["CHEBI"].map(str) + ";" + df["idmetlin"].map(str) + df["CAS"].map(str)
dfLib = df[["idstjude", "name", "synonym", "formula", "monoisotopic_mass", "SMILES", "InChIKey"]]
dfLib = dfLib.rename(columns = {"idstjude": "id", "monoisotopic_mass": "mass"})
dfLib.columns = dfLib.columns.str.lower()    # Column names are all lowercases

# Spectral metadata table
colNames = [col for col in df.columns if col.lower().startswith(condition)]
dfLib["rt"] = df[condition + "_rt"]
dfLib["rt"] = [float(val) * 60 if val != "na" else None for val in dfLib["rt"]] # Convert to "second" unit
dfLib["charge"] = df[condition + "_charge"]

##############################
# Processing of MS2 spectrum #
##############################
dfAllMs2 = pd.DataFrame()
ms2Array = list()
for i in range(dfLib.shape[0]):
    colName = condition + "_linkms2"
    ms2Path = df[colName].iloc[i]
    if ms2Path != "na":
        # Open .MS2 file, write peaks to a DataFrame and save the DataFrame to the database
        # Although the path of .MS2 file contains a letter representing an ion mode (either "p" or "n"),
        # "uid" does not have the letter to make it consistent with "idstjude"
        uid = os.path.splitext(os.path.basename(ms2Path))[0][:-1]
        dfMs2 = pd.read_csv(ms2Path, sep = "\t", engine = "python")    # Header (precursor m/z and charge) is ignored
        dfMs2.columns = ["mz", "intensity"]
        dfMs2["id"] = dfLib["id"].iloc[i]
        dfAllMs2 = dfAllMs2.append(dfMs2, ignore_index = True)
dfAllMs2.to_sql("ms2", conn, if_exists = "replace", index = False)   # Table name is "ms2"

#################################################################
# Addition of decoys (by adding 3 * proton to the neutral mass) #
#################################################################
proton = 1.007276466812
dfDecoy = dfLib.copy()
for i in range(dfDecoy.shape[0]):
    dfDecoy.loc[i, "id"] = "##Decoy_" + dfDecoy.loc[i, "id"]
    dfDecoy.loc[i, "mass"] += 3 * proton # This way prevents 'SettingwithCopyWarning'

# Merge target and decoy DataFrames into one
dfLib = dfLib.append(dfDecoy, ignore_index = True)
dfLib.to_sql("library", conn, if_exists = "replace", index = False)    # Table name is "library"

conn.close()

'''
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
'''