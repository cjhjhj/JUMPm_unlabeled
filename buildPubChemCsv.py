#!/usr/bin/python

# Build a "LocalCSV" for PubChem database
import sys, os, gzip, shutil, pandas as pd, numpy as np
from ftplib import FTP

# Initialization
H = 1.0078250321    # Mass of hydrogen (global variable)
outputFile = "pubchem.csv"

# Connect to FTP server #
ftp = FTP("ftp.ncbi.nlm.nih.gov")
ftp.login()
ftp.cwd('pubchem/Compound/CURRENT-Full/SDF/')
files = ftp.nlst()
ftp.quit()
files = [file for file in files if file.endswith("gz")]

for file in files:
    # Downloading .gz file from NCBI FTP server
    filename = os.path.basename(file)
    with open(filename, "wb") as f:
        print("  Downloading {} ...".format(filename))
        ftp = FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login()
        ftp.cwd('pubchem/Compound/CURRENT-Full/SDF/')
        ftp.retrbinary("RETR %s" % file, f.write)
        ftp.quit()

    # Decompress the downloaded file
    print("  Decompressing {} ...".format(filename))
    gzFile = filename # e.g. Compound_xxx_yyy.xml.gz
    sdfFile = gzFile.rsplit(".", 1)[0]  # Compound_xxx_yyy.xml
    with gzip.open(gzFile, "r") as fIn, open(sdfFile, "wb") as fOut:
        shutil.copyfileobj(fIn, fOut)

    # Read a SDF file
    pcDict = {"Identifier": [], "CompoundName": [], "MolecularFormula": [], "MonoisotopicMass": [],
              "SMILES": [], "InChI": [], "InChIKey": []}
    uid, name, inchi, inchikey, formula, smiles, mass = "NA", "NA", "NA", "NA", "NA", "NA", None
    n = 0
    with open(sdfFile, encoding="utf-8") as f:
        for line in f:  # not using readlines(), as this consumes the memory
            line = line.strip()
            if line.endswith("<PUBCHEM_COMPOUND_CID>"):
                uid = f.readline().strip()
            elif line.endswith("<PUBCHEM_IUPAC_NAME>"):
                name = f.readline().strip()
            elif line.endswith("<PUBCHEM_IUPAC_INCHI>"):
                inchi = f.readline().strip()
            elif line.endswith("<PUBCHEM_IUPAC_INCHIKEY>"):
                inchikey = f.readline().strip()
            elif line.endswith("<PUBCHEM_MOLECULAR_FORMULA>"):
                formula = f.readline().strip()
            elif line.endswith("<PUBCHEM_OPENEYE_CAN_SMILES>"):
                smiles = f.readline().strip()
            elif line.endswith("<PUBCHEM_MONOISOTOPIC_WEIGHT>"):
                mass = f.readline().strip()
            elif line.endswith("$$$$"):
                pcDict["Identifier"].append(uid)
                pcDict["CompoundName"].append(name)
                pcDict["InChI"].append(inchi)
                pcDict["InChIKey"].append(inchikey)
                pcDict["MolecularFormula"].append(formula)
                pcDict["SMILES"].append(smiles)
                pcDict["MonoisotopicMass"].append(mass)
                uid, name, inchi, inchikey, formula, smiles, mass = "NA", "NA", "NA", "NA", "NA", "NA", None
                n += 1
                if n % 1000 == 0:
                    text = "\r  {0} entries are parsed and written to a database".format(n)
                    sys.stdout.write(text)
                    sys.stdout.flush()
    f.close()
    print("\n")

    # Write the parsed information to a CSV file
    df = pd.DataFrame.from_dict(pcDict, orient = "columns")
    df = df[~df["MonoisotopicMass"].isna()]
    df.to_csv(outputFile, mode = "a", index = False, header = not os.path.exists(outputFile))

    os.remove(gzFile)
    os.remove(sdfFile)
