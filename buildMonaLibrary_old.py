#!/usr/bin/python

import sys, os, sqlite3, pyodbc, re, numpy as np, pandas as pd
from sqlalchemy import create_engine, event

sdfFile = r"/Research/Projects/7Metabolomics/library/MoNA/MoNA-export-LC-MS-MS_Positive_Mode.sdf"

# Open a SQLite database and write the library DataFrame to the database
dbName = os.path.splitext(sdfFile)[0] + ".db"
conn = sqlite3.connect(dbName)
# engine = create_engine('sqlite:///' + dbName, echo=False)


# @event.listens_for(engine, "before_cursor_execute")
# def receive_before_cursor_execute(conn, cursor, statement, params, context, executemany):
#     if executemany:
#         cursor.fast_executemany = True


########################
# Initialize variables #
########################
n = 0
proton = 1.007276466812
flagSynonym, flagComment, flagMS2 = 0, 0, 0
uid, otherIds, name, synonym, formula, mass, energy, inchikey, smiles, rt, charge = \
    "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"
smile, kegg, hmdb, pcid, psid, chebi, chemspider, cas = \
    "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"

# Dictionaries for library compounds and MS2 spectra
# They will be changed to pandas DataFrames and saved to a SQLite file
dictLib = {"id": [], "other_ids": [], "name": [], "synonym": [], "formula": [], "mass": [],
           "energy": [], "smiles": [], "inchikey": [], "rt": [], "charge": []}
dictAllMs2 = {"id": [], "mz": [], "intensity": []}

####################
# Parse a SDF file #
####################
with open(sdfFile, encoding="utf-8") as f:
    for line in f:  # Read one line at a time to save memory usage
        line = line.strip()
        if line.startswith(">"):
            flagSynonym, flagComment, flagMS2 = 0, 0, 0
            if line.endswith("<ID>"):
                uid = f.readline().strip()
            elif line.endswith("<NAME>"):
                name = f.readline().strip()
            elif line.endswith("<SYNONYMS>"):
                flagSynonym = 1
                synonym = f.readline().strip()
            elif line.endswith("<FORMULA>"):
                formula = f.readline().strip()
            elif line.endswith("<EXACT MASS>"):
                mass = float(f.readline().strip())
            elif line.endswith("<INCHIKEY>"):
                inchikey = f.readline().strip()
            elif line.endswith("<COMMENT>"):
                # SMILES and Other IDs (KEGG, HMDB, PubChem_CID, PubChem_SID, ChEBI, ChemSpider, CAS) need to be parsed
                flagComment = 1
            elif line.endswith("<PRECURSOR TYPE>"):
                ionType = f.readline().strip()
                # Charge states of compounds were manually inspected by comparing <PRECURSOR TYPE> field and
                # adduct information (https://fiehnlab.ucdavis.edu/images/files/software/ESI-MS-adducts-2020.xls)
                # For the positive mode file
                if ionType == "[M+2H]+" or ionType == "[M+H+K]+" or ionType == "[M+2H]++" or ionType == "[M]++" or ionType == "[M+H+Na]+" or ionType == "[M+H]2+":
                    charge = 2
                # For the negative mode file
                elif ionType == "[M-2H]-" or ionType == "[M-2H]--":
                    charge = 2
                else:
                    charge = 1
            elif line.endswith("<COLLISION ENERGY>"):
                energy = f.readline().strip()
            elif line.endswith("<MASS SPECTRAL PEAKS>"):
                flagMS2 = 1
                dictMs2 = {"id": [], "mz": [], "intensity": []}
        elif line.endswith("$$$$"):
            # Output formatting: "$$$$" is a separator of each compound
            dictLib["id"].append(uid)
            otherIds = kegg + ";" + hmdb + ";" + pcid + ";" + psid + ";" + chebi + ";" + chemspider + ";" + cas
            dictLib["other_ids"].append(otherIds)
            dictLib["name"].append(name)
            dictLib["synonym"].append(synonym)
            dictLib["formula"].append(formula)
            dictLib["mass"].append(mass)
            dictLib["energy"].append(energy)
            dictLib["inchikey"].append(inchikey)
            dictLib["smiles"].append(smiles)
            dictLib["rt"].append(rt)
            # Some entries do not have <PRECURSOR TYPE> field information. In this cae, charge is set to 1
            if charge is "NA":
                charge = 1
            dictLib["charge"].append(charge)

            # If MS2 spectrum exists (i.e. dictMs2 is not None), it is processed (intensity normalization) and merged to dictAllMs2
            try:
                dictMs2
            except NameError:
                dictMs2 = None
            if dictMs2 is not None:
                # Normalize intensities of MS2 spectrum (base peak = 100)
                dictMs2["intensity"] = [val / max(dictMs2["intensity"]) * 100 for val in dictMs2["intensity"]]
                dfMs2 = pd.DataFrame.from_dict(dictMs2, orient="columns")
                dfDecoyMs2 = dfMs2.copy()
                for i in range(dfDecoyMs2.shape[0]):
                    dfDecoyMs2.loc[i, "id"] = "##Decoy_" + dfDecoyMs2.loc[i, "id"]
                dfMs2 = dfMs2.append(dfDecoyMs2, ignore_index=True)
                if n == 0:
                    dfMs2.to_sql("ms2", con=conn, index=False, if_exists="replace")
                else:
                    dfMs2.to_sql("ms2", con=conn, index=False, if_exists="append")

            # # Generate a decoy compound corresponding to the current compound
            # uid = "##Decoy_" + uid  # Change unique id (with the prefix of ##Decoy)
            # dictLib["id"].append(uid)
            # dictLib["other_ids"].append(otherIds)
            # dictLib["name"].append(name)
            # dictLib["synonym"].append(synonym)
            # dictLib["formula"].append(formula)
            # dictLib["mass"].append(mass + 3 * proton)    # Change mass by adding 3 protons
            # dictLib["energy"].append(energy)
            # dictLib["inchikey"].append(inchikey)
            # dictLib["smiles"].append(smiles)
            # dictLib["rt"].append(rt)
            # dictLib["charge"].append(charge)
            # if dictMs2 is not None:
            #     dictMs2["id"] = ["##Decoy_" + val for val in dictMs2["id"]]
            #     for k, v in dictAllMs2.items():
            #         v.extend(dictMs2[k])

            # Re-initialize variables (for next compound)
            uid, otherIds, name, synonym, formula, mass, energy, inchikey, smiles, rt, charge = \
                "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"
            smile, kegg, hmdb, pcid, psid, chebi, chemspider, cas = \
                "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"
            flagSynonym, flagComment, flagMS2 = 0, 0, 0

            n += 1
            if n % 1000 == 0:
                text = "\r  {0} entries are parsed and written to a database".format(n)
                sys.stdout.write(text)
                sys.stdout.flush()
        else:
            if flagSynonym == 1 and line != "":
                synonym = synonym + ";" + f.readline().strip()
            elif flagComment == 1 and line != "":
                if line.startswith("SMILES"):
                    smiles = line.replace("SMILES=", "")
                elif line.startswith("retention time"):
                    rt = line.replace("retention time=", "")
                    if rt.endswith("min") or rt.endswith("minutes") or rt.endswith("minute") or rt.endswith("m"):
                        try:
                            rt = float(re.search("[0-9.]+",
                                                 rt).group()) * 60  # Convert to numeric value and the unit of second
                        except:
                            rt = "NA"
                    elif rt.endswith("sec") or rt.endswith("seconds") or rt.endswith("second") or rt.endswith("s"):
                        try:
                            rt = float(re.search("[0-9.]+", rt).group())
                        except:
                            rt = "NA"
                else:
                    key = line.split("=")[0]
                    if key == "kegg":
                        kegg = line.replace(key + "=", "")
                    elif key == "hmdb":
                        hmdb = line.replace(key + "=", "")
                    elif key == "pubchem cid":
                        pcid = line.replace(key + "=", "")
                    elif key == "pubchem sid":
                        psid = line.replace(key + "=", "")
                    elif key == "chebi":
                        chebi = line.replace(key + "=", "")
                    elif key == "chemspider":
                        chemspider = line.replace(key + "=", "")
                    elif key == "cas":
                        cas = line.replace(key + "=", "")
            elif flagMS2 == 1 and line != "":
                mz, intensity = line.split(" ")
                dictMs2["id"].append(uid)
                dictMs2["mz"].append(float(mz))
                dictMs2["intensity"].append(float(intensity))

###########################################################################
# Organize the DataFrame of library table (with the generation of decoys) #
###########################################################################
dfLib = pd.DataFrame.from_dict(dictLib, orient = "columns")
dfLib = dfLib.rename(columns = {"other_ids": "other_ids(KEGG;HMDB;PubChem_CID;PubChem_SID;ChEBI;ChemSpider;CAS)",
                                "energy": "collision_energy"})
dfDecoy = dfLib.copy()
proton = 1.007276466812
for i in range(dfDecoy.shape[0]):
    dfDecoy.loc[i, "id"] = "##Decoy_" + dfDecoy.loc[i, "id"]
    dfDecoy.loc[i, "mass"] += 3 * proton # This way prevents 'SettingwithCopyWarning'
dfLib = dfLib.append(dfDecoy, ignore_index = True)
dfLib.to_sql("library", conn, if_exists = "replace")    # Table name is "library"

# ############################################################
# # Organize the DataFrame of MS2 spectra (including decoys) #
# ############################################################
# dfMs2 = pd.DataFrame.from_dict(dictMs2, orient = "columns")
# dfDecoyMs2 = dfMs2.copy()
# for i in range(dfDecoyMs2.shape[0]):
#     dfDecoyMs2.loc[i, "id"] = "##Decoy_" + dfDecoyMs2.loc[i, "id"]
# dfMs2 = dfMs2.append(dfDecoyMs2, ignore_index = True)
# dfMs2.to_sql("ms2", conn, if_exists = "replace")

conn.close()
'''
import sys, os, sqlite3, re, numpy as np, pandas as pd

sdfFile = r"/Research/Projects/7Metabolomics/library/MoNA/MoNA-export-LC-MS-MS_Positive_Mode.sdf"


########################
# Initialize variables #
########################
n = 0
proton = 1.007276466812
flagSynonym, flagComment, flagMS2 = 0, 0, 0
uid, otherIds, name, synonym, formula, mass, energy, inchikey, smiles, rt, charge = \
    "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"
smile, kegg, hmdb, pcid, psid, chebi, chemspider, cas = \
    "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"

# Dictionaries for library compounds and MS2 spectra
# They will be changed to pandas DataFrames and saved to a SQLite file
dictLib = {"id": [], "other_ids": [], "name": [], "synonym": [], "formula": [], "mass": [],
           "energy": [], "smiles": [], "inchikey": [], "rt": [], "charge": []}
dictAllMs2 = {"id": [], "mz":[], "intensity": []}

# Open a SQLite database and write the library DataFrame to the database
dbName = os.path.splitext(sdfFile)[0] + ".db"
conn = sqlite3.connect(dbName)

####################
# Parse a SDF file #
####################
with open(sdfFile, encoding = "utf-8") as f:
    for line in f:  # Read one line at a time to save memory usage
        line = line.strip()
        if line.startswith(">"):
            flagSynonym, flagComment, flagMS2 = 0, 0, 0
            if line.endswith("<ID>"):
                uid = f.readline().strip()
            elif line.endswith("<NAME>"):
                name = f.readline().strip()
            elif line.endswith("<SYNONYMS>"):
                flagSynonym = 1
                synonym = f.readline().strip()
            elif line.endswith("<FORMULA>"):
                formula = f.readline().strip()
            elif line.endswith("<EXACT MASS>"):
                mass = float(f.readline().strip())
            elif line.endswith("<INCHIKEY>"):
                inchikey = f.readline().strip()
            elif line.endswith("<COMMENT>"):
                # SMILES and Other IDs (KEGG, HMDB, PubChem_CID, PubChem_SID, ChEBI, ChemSpider, CAS) need to be parsed
                flagComment = 1
            elif line.endswith("<PRECURSOR TYPE>"):
                ionType = f.readline().strip()
                # Charge states of compounds were manually inspected by comparing <PRECURSOR TYPE> field and
                # adduct information (https://fiehnlab.ucdavis.edu/images/files/software/ESI-MS-adducts-2020.xls)
                # For the positive mode file
                if ionType == "[M+2H]+" or ionType == "[M+H+K]+" or ionType == "[M+2H]++" or ionType == "[M]++" or ionType == "[M+H+Na]+" or ionType == "[M+H]2+":
                    charge = 2
                # For the negative mode file
                elif ionType == "[M-2H]-" or ionType == "[M-2H]--":
                    charge = 2
                else:
                    charge = 1
            elif line.endswith("<COLLISION ENERGY>"):
                energy = f.readline().strip()
            elif line.endswith("<MASS SPECTRAL PEAKS>"):
                flagMS2 = 1
                dictMs2 = {"id": [], "mz": [], "intensity": []}
        elif line.endswith("$$$$"):
            # Output formatting: "$$$$" is a separator of each compound
            dictLib["id"].append(uid)
            otherIds = kegg + ";" + hmdb + ";" + pcid + ";" + psid + ";" + chebi + ";" + chemspider + ";" + cas
            dictLib["other_ids"].append(otherIds)
            dictLib["name"].append(name)
            dictLib["synonym"].append(synonym)
            dictLib["formula"].append(formula)
            dictLib["mass"].append(mass)
            dictLib["energy"].append(energy)
            dictLib["inchikey"].append(inchikey)
            dictLib["smiles"].append(smiles)
            dictLib["rt"].append(rt)
            # Some entries do not have <PRECURSOR TYPE> field information. In this cae, charge is set to 1
            if charge is "NA":
                charge = 1
            dictLib["charge"].append(charge)

            # If MS2 spectrum exists (i.e. dictMs2 is not None), it is processed (intensity normalization) and merged to dictAllMs2
            try:
                dictMs2
            except NameError:
                dictMs2 = None
            if dictMs2 is not None:
                # Normalize intensities of MS2 spectrum (base peak = 100)
                dictMs2["intensity"] = [val / max(dictMs2["intensity"]) * 100 for val in dictMs2["intensity"]]
                for k, v in dictAllMs2.items():
                    v.extend(dictMs2[k])

            # # Generate a decoy compound corresponding to the current compound
            # uid = "##Decoy_" + uid  # Change unique id (with the prefix of ##Decoy)
            # dictLib["id"].append(uid)
            # dictLib["other_ids"].append(otherIds)
            # dictLib["name"].append(name)
            # dictLib["synonym"].append(synonym)
            # dictLib["formula"].append(formula)
            # dictLib["mass"].append(mass + 3 * proton)    # Change mass by adding 3 protons
            # dictLib["energy"].append(energy)
            # dictLib["inchikey"].append(inchikey)
            # dictLib["smiles"].append(smiles)
            # dictLib["rt"].append(rt)
            # dictLib["charge"].append(charge)
            # if dictMs2 is not None:
            #     dictMs2["id"] = ["##Decoy_" + val for val in dictMs2["id"]]
            #     for k, v in dictAllMs2.items():
            #         v.extend(dictMs2[k])

            # Re-initialize variables (for next compound)
            uid, otherIds, name, synonym, formula, mass, energy, inchikey, smiles, rt, charge = \
                "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"
            smile, kegg, hmdb, pcid, psid, chebi, chemspider, cas = \
                "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"
            flagSynonym, flagComment, flagMS2 = 0, 0, 0

            n += 1
            if n % 1000 == 0:
                text = "\r  {0} entries are parsed and written to a database".format(n)
                sys.stdout.write(text)
                sys.stdout.flush()
        else:
            if flagSynonym == 1 and line != "":
                synonym = synonym + ";" + f.readline().strip()
            elif flagComment == 1 and line != "":
                if line.startswith("SMILES"):
                    smiles = line.replace("SMILES=", "")
                elif line.startswith("retention time"):
                    rt = line.replace("retention time=", "")
                    if rt.endswith("min") or rt.endswith("minutes") or rt.endswith("minute") or rt.endswith("m"):
                        try:
                            rt = float(re.search("[0-9.]+", rt).group()) * 60   # Convert to numeric value and the unit of second
                        except:
                            rt = "NA"
                    elif rt.endswith("sec") or rt.endswith("seconds") or rt.endswith("second") or rt.endswith("s"):
                        try:
                            rt = float(re.search("[0-9.]+", rt).group())
                        except:
                            rt = "NA"
                else:
                    key = line.split("=")[0]
                    if key == "kegg":
                        kegg = line.replace(key + "=", "")
                    elif key == "hmdb":
                        hmdb = line.replace(key + "=", "")
                    elif key == "pubchem cid":
                        pcid = line.replace(key + "=", "")
                    elif key == "pubchem sid":
                        psid = line.replace(key + "=", "")
                    elif key == "chebi":
                        chebi = line.replace(key + "=", "")
                    elif key == "chemspider":
                        chemspider = line.replace(key + "=", "")
                    elif key == "cas":
                        cas = line.replace(key + "=", "")
            elif flagMS2 == 1 and line != "":
                mz, intensity = line.split(" ")
                dictMs2["id"].append(uid)
                dictMs2["mz"].append(float(mz))
                dictMs2["intensity"].append(float(intensity))

###########################################################################
# Organize the DataFrame of library table (with the generation of decoys) #
###########################################################################
dfLib = pd.DataFrame.from_dict(dictLib, orient = "columns")
dfLib = dfLib.rename(columns = {"other_ids": "other_ids(KEGG;HMDB;PubChem_CID;PubChem_SID;ChEBI;ChemSpider;CAS)",
                                "energy": "collision_energy"})
dfDecoy = dfLib.copy()
proton = 1.007276466812
for i in range(dfDecoy.shape[0]):
    dfDecoy.loc[i, "id"] = "##Decoy_" + dfDecoy.loc[i, "id"]
    dfDecoy.loc[i, "mass"] += 3 * proton # This way prevents 'SettingwithCopyWarning'
dfLib = dfLib.append(dfDecoy, ignore_index = True)
dfLib.to_sql("library", conn, if_exists = "replace")    # Table name is "library"

############################################################
# Organize the DataFrame of MS2 spectra (including decoys) #
############################################################
dfMs2 = pd.DataFrame.from_dict(dictMs2, orient = "columns")
dfDecoyMs2 = dfMs2.copy()
for i in range(dfDecoyMs2.shape[0]):
    dfDecoyMs2.loc[i, "id"] = "##Decoy_" + dfDecoyMs2.loc[i, "id"]
dfMs2 = dfMs2.append(dfDecoyMs2, ignore_index = True)
dfMs2.to_sql("ms2", conn, if_exists = "replace")

conn.close()
'''
