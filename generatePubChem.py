#!/usr/bin/python

# Build HMDB database using SQLite
import sys, os, sqlite3, gzip, shutil
from lxml import etree
from ftplib import FTP

# ftp = FTP("ftp.ncbi.nlm.nih.gov")
# ftp.login()
# ftp.cwd('pubchem/Compound/CURRENT-Full/XML/')
# files = ftp.nlst()
# # files = [file for file in files if file.endswith("gz")]
#
#
# files = [file for file in files if "Compound_050000001_050500000" in file]
#
# for file in files:
#     # Downloading .gz file from NCBI FTP server
#     filename = os.path.basename(file)
#     with open(filename, "wb") as f:
#         print("Downloading {} ...".format(filename))
#         ftp.retrbinary("RETR %s" % file, f.write)
#
#         print("Decompressing {} ...".format(filename))
#
#     # Decompress the downloaded file
#     infile = filename # e.g. Compound_xxx_yyy.xml.gz
#     outfile = infile.rsplit(".", 1)[0]  # Compound_xxx_yyy.xml
#     with gzip.open(infile, "r") as fIn, open(outfile, "wb") as fOut:
#         shutil.copyfileobj(fIn, fOut)




##################
# Create a table #
##################
# conn = sqlite3.connect('pubchem.db')   # SQLite database for HMDB
# conn.execute('''CREATE TABLE HMDB
#                 (ID INTEGER PRIMARY KEY NOT NULL,
#                 HMDBID VARCHAR(255) NOT NULL,
#                 FORMULA VARCHAR(255) NOT NULL,
#                 SMILES VARCHAR(255),
#                 INCHIKEY VARCHAR(255),
#                 NAME VARCHAR(255),
#                 MONOMASS REAL NOT NULL,
#                 TYPE VARCHAR(255));''')

################################################
# Read a XML file and write to a database file #
################################################
inputFile = "Compound_050000001_050500000.xml"
n = 0
H = 1.0078250321    # Mass of hydrogen
try:
    # Reading a XML file
    for event, element in etree.iterparse(inputFile, tag = "{*}PC-Compound"):
        for child in element:
            tag = etree.QName(child).localname  # Remove namespace from the tag
            val = child.text
            if tag == "accession" and val is not None:
                hmdbId = val
                n += 1
                if n % 1000 == 0:
                    text = "\r  {0} entries are parsed and written to a database".format(n)
                    sys.stdout.write(text)
                    sys.stdout.flush()
            elif tag == "chemical_formula" and val is not None:
                formula = val
            elif tag == "smiles" and val is not None:
                smiles = val
            elif tag == "inchi" and val is not None:
                inchi = val
            elif tag == "inchikey" and val is not None:
                inchiKey = val
            elif tag == "monisotopic_molecular_weight" and val is not None:
                monoMass = val
                decoyMass = str(float(monoMass) + H)
            elif tag == "name" and val is not None:
                name = val
        # conn.execute('INSERT INTO HMDB (HMDBID, FORMULA, SMILES, INCHIKEY, NAME, MONOMASS, TYPE) VALUES (?, ?, ?, ?, ?, ?, ?)',
        #              (hmdbId, formula, smiles, inchiKey, name, monoMass, "target"))
        # conn.execute('INSERT INTO HMDB (HMDBID, FORMULA, SMILES, INCHIKEY, NAME, MONOMASS, TYPE) VALUES (?, ?, ?, ?, ?, ?, ?)',
        #              (hmdbId, formula, smiles, inchiKey, name, decoyMass, "decoy"))
        element.clear()
    # conn.commit()
    # conn.close()
except ImportError:
    print ("Error in import lxml package")
