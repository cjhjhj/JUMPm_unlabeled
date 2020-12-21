#!/usr/bin/python

# Generate a CSV file for LipidMAPS database
# This file is for using MetFragCL. Even though it is said that MetFragCL can accept SDF file format, it does not work well
# Currently, only CSV file format works well
import sys, os, pandas as pd


sdfFile = sys.argv[1]
# sdfFile = "/Research/Projects/7Metabolomics/Database/LipidMaps/structures.sdf"

########################
# Initialize variables #
########################
uid, otherIds, name, sysname, synonyms, abbreviation, formula, smiles, inchi, inchikey, category, mainclass, subclass, mass = \
    "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", None
pubchem, chebi, kegg, hmdb, swisslipid, lipidbank, plantfa = "NA", "NA", "NA", "NA", "NA", "NA", "NA"
lmapsDict = {"Identifier": [], "CompoundName": [], "MolecularFormula": [], "MonoisotopicMass": [],
             "SystematicName": [], "Synonyms": [], "Abbreviation": [], "Category": [], "MainClass": [], "SubClass": [],
             "SMILES": [], "InChI": [], "InChIKey": [], "OtherIDs(PubChem;ChEBI;KEGG;HMDB;SwissLipid;LipidBank;PlantFA)": []}

####################
# Parse a SDF file #
####################
n = 0
with open(sdfFile, encoding="utf-8") as f:
    for line in f:  # Read one line at a time to save memory usage
        line = line.strip()
        if line.startswith(">"):
            flagSynonym, flagComment, flagMS2 = 0, 0, 0
            if line.endswith("<LM_ID>"):
                uid = f.readline().strip()
                uid = uid.replace("-", "_") # Convert hyphen(s) (i.e. "-") to underbar(s) (i.e. "_")
            elif line.endswith("<NAME>"):
                name = f.readline().strip()
            elif line.endswith("<SYSTEMATIC_NAME>"):
                sysname = f.readline().strip()
            elif line.endswith("<ABBREVIATION>"):
                abbreviation = f.readline().strip()
            elif line.endswith("<SYNONYMS>"):
                synonyms = f.readline().strip()
            elif line.endswith("<CATEGORY>"):
                category = f.readline().strip()
            elif line.endswith("<MAIN_CLASS>"):
                mainclass = f.readline().strip()
            elif line.endswith("<SUB_CLASS>"):
                subclass = f.readline().strip()
            elif line.endswith("<FORMULA>"):
                formula = f.readline().strip()
            elif line.endswith("<EXACT_MASS>"):
                mass = float(f.readline().strip())
            elif line.endswith("<INCHI>"):
                inchi = f.readline().strip()
            elif line.endswith("<INCHI_KEY>"):
                inchikey = f.readline().strip()
            elif line.endswith("<SMILES>"):
                smiles = f.readline().strip()
            elif line.endswith("<PUBCHEM_ID>"):
                pubchem = f.readline().strip()
            elif line.endswith("<CHEBI_ID>"):
                chebi = f.readline().strip()
            elif line.endswith("<KEGG_ID>"):
                kegg = f.readline().strip()
            elif line.endswith("<HMDB_ID>"):
                hmdb = f.readline().strip()
            elif line.endswith("<SWISSLIPIDS_ID>"):
                swisslipid = f.readline().strip()
            elif line.endswith("<LIPIDBANK_ID>"):
                lipidbank = f.readline().strip()
            elif line.endswith("<PLANTFA_ID>"):
                plantfa = f.readline().strip()
        elif line.endswith(r"$$$$"):
            # Output formatting: "$$$$" is a separator of each compound
            lmapsDict["Identifier"].append(uid)
            lmapsDict["CompoundName"].append(name)
            lmapsDict["MolecularFormula"].append(formula)
            lmapsDict["MonoisotopicMass"].append(mass)
            lmapsDict["SystematicName"].append(sysname)
            lmapsDict["Synonyms"].append(synonyms)
            lmapsDict["Abbreviation"].append(abbreviation)
            lmapsDict["Category"].append(category)
            lmapsDict["MainClass"].append(mainclass)
            lmapsDict["SubClass"].append(subclass)
            lmapsDict["SMILES"].append(smiles)
            lmapsDict["InChI"].append(inchi)
            lmapsDict["InChIKey"].append(inchikey)
            otherIds = pubchem + ";" + chebi + ";" + kegg + ";" + hmdb + ";" + swisslipid + ";" + lipidbank + ";" + plantfa
            lmapsDict["OtherIDs(PubChem;ChEBI;KEGG;HMDB;SwissLipid;LipidBank;PlantFA)"].append(otherIds)

            # Re-initialize variables (for next compound)
            uid, otherIds, name, sysname, synonyms, abbreviation, formula, smiles, inchi, inchikey, category, mainclass, subclass, mass = \
                "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", None
            pubchem, chebi, kegg, hmdb, swisslipid, lipidbank, plantfa = "NA", "NA", "NA", "NA", "NA", "NA", "NA"

            n += 1
            if n % 1000 == 0:
                text = "\r  {0} entries are parsed and written to a database".format(n)
                sys.stdout.write(text)
                sys.stdout.flush()

# Create a DataFrame for database entries and their decoys
df = pd.DataFrame.from_dict(lmapsDict, orient="columns")
df = df[~df["MonoisotopicMass"].isna()]
proton = 1.007276466812
dfDecoy = df.copy()
dfDecoy["Identifier"] = "##Decoy_" + dfDecoy["Identifier"]
dfDecoy["MonoisotopicMass"] = pd.to_numeric(dfDecoy["MonoisotopicMass"]) + 3 * proton
df = df.append(dfDecoy, ignore_index=True)
csvFile = os.path.splitext(sdfFile)[0] + ".csv"
df.to_csv(csvFile, index=False)
print("  Done")
