#!/usr/bin/python

import sys, os, sqlite3, re, numpy as np, pandas as pd

# Custom script to generate MoNA metabolite library
# Template is a SDF file downloaded from MoNA website
# MS2 spectrum of each metabolite (when available) is also in the SDF file (needs to be parsed)

# Usage
# python buildMonaLibrary.py <template file (full path)> <column and mode information (e.g. hilicn, c18p, etc.)>

# Open a SQLite database and write the library DataFrame to the database
# sdfFile = r"/Research/Projects/7Metabolomics/library/MoNA/MoNA-export-LC-MS-MS_Positive_Mode.sdf"
# sdfFile = r"/Research/Projects/7Metabolomics/library/MoNA/MoNA-export-LipidBlast.sdf"
sdfFile = sys.argv[1]
condition = sys.argv[2]
dbName = os.path.splitext(sdfFile)[0] + ".db"
conn = sqlite3.connect(dbName)

########################
# Initialize variables #
########################
n = 0
proton = 1.007276466812
flagSynonym, flagComment, flagMS2, nPeaks = 0, 0, 0, 0
uid, otherIds, name, synonym, formula, energy, inchikey, smiles, iontype, rt, mass, precmz, charge = \
    "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", None, None, None, None
smile, kegg, hmdb, pcid, psid, chebi, chemspider, cas = \
    "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"

# Dictionaries for library compounds and MS2 spectra
# They will be changed to pandas DataFrames and saved to a SQLite file
dictLib = {"id": [], "other_ids": [], "name": [], "synonym": [], "formula": [], "precursor_mz": [], "mass": [],
           "ion_type": [], "energy": [], "smiles": [], "inchikey": [], "rt": [], "charge": []}

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
                uid = uid.replace("-", "_") # Convert hyphen(s) (i.e. "-") to underbar(s) (i.e. "_")
            elif line.endswith("<NAME>"):
                name = f.readline().strip()
            elif line.endswith("<SYNONYMS>"):
                flagSynonym = 1
                synonym = f.readline().strip()
            elif line.endswith("<FORMULA>"):
                formula = f.readline().strip()
            elif line.endswith("<PRECURSOR M/Z>"):
                precmz = float(f.readline().strip())
            elif line.endswith("<EXACT MASS>"):
                mass = float(f.readline().strip())
            elif line.endswith("<INCHIKEY>"):
                inchikey = f.readline().strip()
            elif line.endswith("<COMMENT>"):
                # SMILES and Other IDs (KEGG, HMDB, PubChem_CID, PubChem_SID, ChEBI, ChemSpider, CAS) need to be parsed
                flagComment = 1
            elif line.endswith("<PRECURSOR TYPE>"):
                iontype = f.readline().strip()
                # Charge states of compounds were manually inspected by comparing <PRECURSOR TYPE> field and
                # adduct information (https://fiehnlab.ucdavis.edu/images/files/software/ESI-MS-adducts-2020.xls)
                # For the positive mode file
                if iontype == "[M+2H]+" or iontype == "[M+H+K]+" or iontype == "[M+2H]++" or iontype == "[M]++" or iontype == "[M+H+Na]+" or iontype == "[M+H]2+":
                    charge = 2
                # For the negative mode file
                elif iontype == "[M-2H]-" or iontype == "[M-2H]--":
                    charge = 2
                else:
                    charge = 1
            elif line.endswith("<COLLISION ENERGY>"):
                energy = f.readline().strip()
            elif line.endswith("<NUM PEAKS>"):
                nPeaks = int(f.readline().strip())
            elif line.endswith("<MASS SPECTRAL PEAKS>") and 0 < nPeaks < 1000:
                flagMS2 = 1
                dictMs2 = {"mz": [], "intensity": []}
        elif line.endswith("$$$$"):
            # Output formatting: "$$$$" is a separator of each compound
            if mass != "NA":
                dictLib["id"].append(uid)
                otherIds = kegg + ";" + hmdb + ";" + pcid + ";" + psid + ";" + chebi + ";" + chemspider + ";" + cas
                dictLib["other_ids"].append(otherIds)
                dictLib["name"].append(name)
                dictLib["synonym"].append(synonym)
                dictLib["formula"].append(formula)
                dictLib["precursor_mz"].append(precmz)
                dictLib["mass"].append(mass)
                dictLib["ion_type"].append(iontype)
                dictLib["energy"].append(energy)
                dictLib["inchikey"].append(inchikey)
                dictLib["smiles"].append(smiles)
                dictLib["rt"].append(rt)
                # Some entries do not have <PRECURSOR TYPE> field information. In this cae, charge is set to 1
                if charge is None:
                    charge = 1
                dictLib["charge"].append(charge)

                # If MS2 spectrum exists (i.e. dictMs2 is not None), it is processed (intensity normalization) and written to the database
                try:
                    dictMs2
                except NameError:
                    dictMs2 = None
                if dictMs2 is not None:
                    # Normalize intensities of MS2 spectrum (base peak = 100)
                    dictMs2["intensity"] = [val / max(dictMs2["intensity"]) * 100 for val in dictMs2["intensity"]]
                    listMs2 = []
                    for i in range(len(dictMs2["mz"])):
                        listMs2.append({"mz": str(dictMs2["mz"][i]), "intensity": str(dictMs2["intensity"][i])})
                    createQuery = "CREATE TABLE " + uid + " (mz REAL, intensity REAL)"
                    conn.execute(createQuery)
                    insertQuery = "INSERT INTO " + uid + " (mz, intensity) VALUES (:mz, :intensity)"
                    conn.executemany(insertQuery, listMs2)
                    # conn.executemany('''INSERT INTO ms2 (mz, intensity) VALUES (:mz, :intensity);''', listMs2)
                    conn.commit()

            # Re-initialize variables (for next compound)
            uid, otherIds, name, synonym, formula, energy, inchikey, smiles, iontype, rt, mass, precmz, charge = \
                "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", None, None, None, None
            smile, kegg, hmdb, pcid, psid, chebi, chemspider, cas = \
                "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"
            flagSynonym, flagComment, flagMS2, nPeaks = 0, 0, 0, 0
            dictMs2 = None

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
                            rt = float(re.search("[0-9.]+", rt).group()) * 60  # Convert to numeric value and the unit of second
                        except (NameError, AttributeError):
                            rt = None
                    elif rt.endswith("sec") or rt.endswith("seconds") or rt.endswith("second") or rt.endswith("s"):
                        try:
                            rt = float(re.search("[0-9.]+", rt).group())
                        except (NameError, AttributeError):
                            rt = None
                elif line.startswith("exact mass") and mass == "NA":
                    mass = float(line.replace("exact mass=", ""))
                elif line.startswith("ion type") and iontype == "NA":
                    iontype = line.replace("ion type=", "")
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
                dictMs2["mz"].append(float(mz))
                dictMs2["intensity"].append(float(intensity))

print("\n  Finished inserting MS2 spectra to the database")

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
if condition[-1] == "p":
    dfDecoy["precursor_mz"] = (dfDecoy["mass"] + dfDecoy["charge"] * proton) / dfDecoy["charge"]
elif condition[-1] == "n":
    dfDecoy["precursor_mz"] = (dfDecoy["mass"] - dfDecoy["charge"] * proton) / dfDecoy["charge"]
dfLib = dfLib.append(dfDecoy, ignore_index = True)
dfLib.to_sql("library", conn, if_exists = "replace")    # Table name is "library"
conn.close()
print("  Finished inserting a library table to the database")
print()
