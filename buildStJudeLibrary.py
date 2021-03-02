#!/usr/bin/python

import sys, os, sqlite3, pandas as pd, numpy as np

# Custom script to generate our own metabolite library
# Template is a text file containing compound metadata
# MS2 spectrum of each metabolite (when available) is stored in a separate text file (.MS2)

# Usage
# python buildStJudeLibrary.py <template file (full path)> <column and mode information (e.g. hilicn, c18p, etc.)>

################################
# Handling of library template #
################################
# Initialization (path of library template file and experimental condition)
# templateFile = r"/Research/Projects/7Metabolomics/Library/StJude/Metabolome_library_v4.1.1.txt"
# condition = "c18p"    # Column name and ion mode, e.g. hilicn = HILIC column with negative ion mode
templateFile = sys.argv[1]
condition = sys.argv[2]
dbName = "stjude_library_" + condition + ".db"
dbName = os.path.join(os.path.dirname(templateFile), dbName)
proton = 1.007276466812
conn = sqlite3.connect(dbName)

# Read a text file containing the information of metabolomes
df = pd.read_csv(templateFile, sep="\t", engine="python")
ind = df[condition + "_linkms2"] != "na"
df = df[ind].reset_index(drop=True)

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
colNameOtherIds = "other_ids(KEGG;HMDB;PubChem_CID;PubChem_SID;ChEBI;METLIN;CAS)"    # Handling "other_ids"
otherIds = df['idkegg'].map(str) + ';' + df['idhmdb'].map(str) + ';' + df['PC_CID'].map(str) + ";" + \
df["PC_SID"].map(str) + ";" + df["CHEBI"].map(str) + ";" + df["idmetlin"].map(str) + df["CAS"].map(str)
dfLib = df[["idstjude", "name", "synonym", "formula", "monoisotopic_mass", "SMILES", "InChIKey"]]
dfLib = dfLib.rename(columns={"idstjude": "id", "monoisotopic_mass": "mass"})
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
        dfMs2 = pd.read_csv(ms2Path, sep="\t")
        dfLib["precursor_mz"].iloc[i] = float(dfMs2.columns[0])
        dfMs2.columns = ["mz", "intensity"]
        dfMs2.to_sql(uid, conn, if_exists="replace", index=False)   # Table name is the same as compound id (e.g. sjm00001)

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
dfLib = dfLib.append(dfDecoy, ignore_index=True)
dfLib.to_sql("library", conn, if_exists="replace", index=False)    # Table name is "library"

conn.close()
