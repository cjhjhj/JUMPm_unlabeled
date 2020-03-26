#!/usr/bin/python

# Build PubChem database using SQLite
import sys, os, sqlite3, gzip, shutil
from lxml import etree
from ftplib import FTP

##################
# Create a table #
##################
conn = sqlite3.connect('pubchem.db')   # SQLite database for PubChem
conn.execute('''CREATE TABLE PUBCHEM
                (ID INTEGER PRIMARY KEY NOT NULL,
                PUBCHEMID VARCHAR(255) NOT NULL,
                FORMULA VARCHAR(255) NOT NULL,
                SMILES VARCHAR(255),
                INCHIKEY VARCHAR(255),
                NAME VARCHAR(255),
                MONOMASS REAL NOT NULL,
                TYPE VARCHAR(255));''')

#########################
# Connect to FTP server #
#########################
ftp = FTP("ftp.ncbi.nlm.nih.gov")
ftp.login()
ftp.cwd('pubchem/Compound/CURRENT-Full/SDF/')
files = ftp.nlst()
ftp.quit()
files = [file for file in files if file.endswith("gz")]

H = 1.0078250321    # Mass of hydrogen (global variable)
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

    # Read a SDF file and write to sqlite database
    n = 0
    with open(sdfFile) as f:
        for line in f:  # not using readlines(), as this consumes the memory
            line = line.strip()
            if line.endswith("<PUBCHEM_COMPOUND_CID>"):
                pubchemId = f.readline().strip()
            elif line.endswith("<PUBCHEM_IUPAC_NAME>"):
                name = f.readline().strip()
            elif line.endswith("<PUBCHEM_IUPAC_INCHIKEY>"):
                inchiKey = f.readline().strip()
            elif line.endswith("<PUBCHEM_MOLECULAR_FORMULA>"):
                formula = f.readline().strip()
            elif line.endswith("<PUBCHEM_OPENEYE_CAN_SMILES>"):
                smiles = f.readline().strip()
            elif line.endswith("<PUBCHEM_MONOISOTOPIC_WEIGHT>"):
                monoMass = f.readline().strip()
            elif line.endswith("$$$$"):
                decoyMass = str(float(monoMass) + H)
                conn.execute(
                    'INSERT INTO PUBCHEM (PUBCHEMID, FORMULA, SMILES, INCHIKEY, NAME, MONOMASS, TYPE) VALUES (?, ?, ?, ?, ?, ?, ?)',
                    (pubchemId, formula, smiles, inchiKey, name, monoMass, "target"))
                conn.execute(
                    'INSERT INTO PUBCHEM (PUBCHEMID, FORMULA, SMILES, INCHIKEY, NAME, MONOMASS, TYPE) VALUES (?, ?, ?, ?, ?, ?, ?)',
                    (pubchemId, formula, smiles, inchiKey, name, decoyMass, "decoy"))
                n += 1
                if n % 1000 == 0:
                    text = "\r  {0} entries are parsed and written to a database".format(n)
                    sys.stdout.write(text)
                    sys.stdout.flush()
    conn.commit()
    f.close()
    print("\n")
    os.remove(gzFile)
    os.remove(sdfFile)

conn.close()
