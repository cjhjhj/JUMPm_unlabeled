#!/usr/bin/python

import sys, os, sqlite3, pyodbc, re, numpy as np, pandas as pd
from sqlalchemy import create_engine, event

sdfFile = r"/Research/Projects/7Metabolomics/library/MoNA/MoNA-export-LipidBlast.sdf"

########################
# Initialize variables #
########################
n = 0
proton = 1.007276466812
flagSynonym, flagComment, flagMS2, nPeaks = 0, 0, 0, 0
uid, otherIds, name, synonym, formula, energy, inchikey, smiles, rt, mass, charge = \
    "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", None, None
smile, kegg, hmdb, pcid, psid, chebi, chemspider, cas = \
    "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"

listNumPeaks = []
dictCharges = {}
nTotalCompounds, nValidCompounds, nMS2Compounds = 0, 0, 0

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
            elif line.endswith("<EXACT MASS>"):
                mass = float(f.readline().strip())
            elif line.endswith("<INCHIKEY>"):
                inchikey = f.readline().strip()
            elif line.endswith("<COMMENT>"):
                # SMILES and Other IDs (KEGG, HMDB, PubChem_CID, PubChem_SID, ChEBI, ChemSpider, CAS) need to be parsed
                flagComment = 1
            elif line.endswith("<PRECURSOR TYPE>"):
                ionType = f.readline().strip()
                # Investigate unique precursor types
                if ionType in dictCharges:
                    dictCharges[ionType] += 1
                else:
                    dictCharges[ionType] = 1
            elif line.endswith("<COLLISION ENERGY>"):
                energy = f.readline().strip()
            elif line.endswith("<NUM PEAKS>"):
                nPeaks = int(f.readline().strip())
                listNumPeaks.append(nPeaks)
            elif line.endswith("<MASS SPECTRAL PEAKS>"):
                flagMS2 = 1
                dictMs2 = 1
        elif line.endswith("$$$$"):
            nTotalCompounds += 1
            # Output formatting: "$$$$" is a separator of each compound
            if mass != "NA":
                nValidCompounds += 1
                # If MS2 spectrum exists (i.e. dictMs2 is not None), it is processed (intensity normalization) and written to the database
                try:
                    dictMs2
                except NameError:
                    dictMs2 = None
                if dictMs2 is not None:
                    nMS2Compounds += 1

            # Re-initialize variables (for next compound)
            uid, otherIds, name, synonym, formula, energy, inchikey, smiles, rt, mass, charge = \
                "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", None, None
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
                if line.startswith("exact mass") and mass == "NA":
                    mass = float(line.replace("exact mass=", ""))
            elif flagMS2 == 1 and line != "":
                dictMs2 = 1

print(nTotalCompounds, nValidCompounds, nMS2Compounds)
print(dictCharges)
# listNumPeaks