#!/usr/bin/python

import sys, os, sqlite3, numpy as np, pandas as pd

# Custom script to generate our own lipid library
# Template is a text file containing compound metadata
# MS2 spectrum of each compound (when available) is stored in a separate text file (.MS2)

# Usage
# python buildStJudeLipidLibrary.py <template file (full path)> <column and mode information (e.g. hilicn, c18p, etc.)>

################################
# Handling of library template #
################################
# Initialization (path of library template file and experimental condition)
# templateFile = r"/Users/jcho/OneDrive - St. Jude Children's Research Hospital/UDrive/Research/Projects/7Metabolomics/Library/Lipid/lipid_library_v0.03_desktop.txt"
# condition = "c18n"    # Column name and ion mode, e.g. hilicn = HILIC column with negative ion mode
templateFile = sys.argv[1]
condition = sys.argv[2]
dbName = "stjude_library_lipid_" + condition + ".db"
dbName = os.path.join(os.path.dirname(templateFile), dbName)
conn = sqlite3.connect(dbName)
proton = 1.007276466812

# Read a text file containing the information of metabolomes
# df = pd.read_csv(templateFile, sep = "\t", engine = "python")
df = pd.read_csv(templateFile, sep="\t", encoding="latin1")
ind = df[condition + "_linkms2"] != "na"
df = df[ind].reset_index(drop = True)

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
# 12. precursor_mz: precursor m/z value
# 13. ion_type: ion type of precursor (e.g. [M+H]+, [M+NH4]+, etc.)

# Compound metadata
# Note that the lipid table does not have some identifiers such as InChIKey
colNameOtherIds = "other_ids(KEGG;HMDB;PubChem_CID;CAS)"    # Handling "other_ids"
otherIds = df['idkegg'].map(str) + ';' + df['idhmdb'].map(str) + ';' + df['PC_CID'].map(str) + ";" + df["CAS"].map(str)
dfLib = df[["idstjude", "name", "synonym", "formula", "monoisotopic_mass", "SMILES"]]
dfLib = dfLib.rename(columns = {"idstjude": "id", "monoisotopic_mass": "mass"})
dfLib.columns = dfLib.columns.str.lower()    # Column names are all lowercases
dfLib[colNameOtherIds] = otherIds

# Spectral metadata table
dfLib["collision_energy"] = df[condition + "_ms2setting"]
dfLib["rt"] = df[condition + "_rt"]
dfLib["rt"] = [float(val) * 60 if val != "na" else None for val in dfLib["rt"]] # Convert to "second" unit
dfLib["charge"] = pd.to_numeric(df[condition + "_charge"])
dfLib["precursor_mz"] = np.nan
sign = ""
if condition[-1] == "p":
    sign = "+"
elif condition[-1] == "n":
    sign = "-"
dfLib["ion_type"] = df[condition + "_adduct"]
for i in range(dfLib.shape[0]):
    if dfLib["ion_type"].loc[i] == "na":
        if dfLib["charge"].loc[i] == 1:
            dfLib["ion_type"].loc[i] = "[M" + sign + "H]" + sign
        else:
            dfLib["ion_type"].loc[i] = "[M" + sign + str(int(dfLib["charge"].loc[i])) + "H]" + sign
dfLib = dfLib.assign(inchikey="NA")

##############################
# Processing of MS2 spectrum #
##############################
for i in range(dfLib.shape[0]):
    colName = condition + "_linkms2"
    ms2Path = df[colName].iloc[i]
    if ms2Path != "na":
        # Open .MS2 file, write peaks to a DataFrame and save the DataFrame to the database
        # Although the path of .MS2 file contains a letter representing an ion mode (either "p" or "n"),
        # "uid" does not have the letter to make it consistent with "idstjude"
        uid = os.path.splitext(os.path.basename(ms2Path))[0][:-1]
        dfMs2 = pd.read_csv(ms2Path, sep="\s+|\t", engine="python")
        dfLib.iloc[i, dfLib.columns.get_loc("precursor_mz")] = float(dfMs2.columns[0])  # To avoid pandas warning
        dfMs2.columns = ["mz", "intensity"]
        dfMs2.to_sql(uid, conn, if_exists = "replace", index = False)   # Table name is the same as compound id (e.g. sjm00001)

#################################################################
# Addition of decoys (by adding 3 * proton to the neutral mass) #
#################################################################
dfDecoy = dfLib.copy()
for i in range(dfDecoy.shape[0]):
    dfDecoy.loc[i, "id"] = "##Decoy_" + dfDecoy.loc[i, "id"]
    dfDecoy.loc[i, "mass"] += 3 * proton # This way prevents 'SettingwithCopyWarning'
if condition[-1] == "p":
    dfDecoy["precursor_mz"] += 3 * proton / dfDecoy["charge"]
elif condition[-1] == "n":
    dfDecoy["precursor_mz"] -= 3 * proton / dfDecoy["charge"]

# Merge target and decoy DataFrames into one
dfLib = dfLib.append(dfDecoy, ignore_index = True)
dfLib.to_sql("library", conn, if_exists = "replace", index = False)    # Table name is "library"

conn.close()
