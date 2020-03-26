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
ftp.cwd('pubchem/Compound/CURRENT-Full/XML/')
files = ftp.nlst()
ftp.quit()
files = [file for file in files if file.endswith("gz")]

H = 1.0078250321    # Mass of hydrogen (global variable)


# For test
files = ["Compound_050000001_050500000.xml.gz",
         "Compound_099000001_099500000.xml.gz",
         "Compound_048500001_049000000.xml.gz"]


for file in files:
    # Downloading .gz file from NCBI FTP server
    filename = os.path.basename(file)
    with open(filename, "wb") as f:
        print("  Downloading {} ...".format(filename))
        ftp = FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login()
        ftp.cwd('pubchem/Compound/CURRENT-Full/XML/')
        ftp.retrbinary("RETR %s" % file, f.write)
        ftp.quit()

    # Decompress the downloaded file
    print("  Decompressing {} ...".format(filename))
    gzFile = filename # e.g. Compound_xxx_yyy.xml.gz
    xmlFile = gzFile.rsplit(".", 1)[0]  # Compound_xxx_yyy.xml
    with gzip.open(gzFile, "r") as fIn, open(xmlFile, "wb") as fOut:
        shutil.copyfileobj(fIn, fOut)

    # Read a XML file and write to sqlite database
    n = 0
    flagId, flagFormula, flagSmiles, flagCanonical = 0, 0, 0, 0
    flagInchikey, flagMonomass, flagName, flagPreferred = 0, 0, 0, 0
    try:
        context = etree.iterparse(xmlFile)
        for action, elem in context:
            elem.tag = etree.QName(elem).localname  # Remove namespace from the tag
            tag = elem.tag
            val = elem.text
            if tag == "PC-CompoundType_id_cid" and val is not None:
                if n > 0:
                    decoyMass = str(float(monoMass) + H)
                    conn.execute(
                        'INSERT INTO PUBCHEM (PUBCHEMID, FORMULA, SMILES, INCHIKEY, NAME, MONOMASS, TYPE) VALUES (?, ?, ?, ?, ?, ?, ?)',
                        (pubchemId, formula, smiles, inchiKey, name, monoMass, "target"))
                    conn.execute(
                        'INSERT INTO PUBCHEM (PUBCHEMID, FORMULA, SMILES, INCHIKEY, NAME, MONOMASS, TYPE) VALUES (?, ?, ?, ?, ?, ?, ?)',
                        (pubchemId, formula, smiles, inchiKey, name, decoyMass, "decoy"))
                    pubchemId = val
                else:
                    pubchemId = val
                n += 1
                if n % 1000 == 0:
                    text = "\r  {0} entries are parsed and written to a database".format(n)
                    sys.stdout.write(text)
                    sys.stdout.flush()
            if val == "Molecular Formula":
                flagFormula = 1
            if flagFormula == 1 and tag == "PC-InfoData_value_sval" and val is not None:
                formula = val
                flagFormula = 0
            if val == "SMILES":
                flagSmiles = 1
            if val == "Canonical":  # Right after "SMILES"
                flagCanonical = 1
            if flagSmiles == 1 and flagCanonical == 1 and tag == "PC-InfoData_value_sval" and val is not None:
                smiles = val
                flagSmiles, flagCanonical = 0, 0
            if val == "InChIKey":
                flagInchikey = 1
            if flagInchikey == 1 and tag == "PC-InfoData_value_sval" and val is not None:
                inchiKey = val
                flagInchikey = 0
            if val == "MonoIsotopic":
                flagMonomass = 1
            if flagMonomass == 1 and tag == "PC-InfoData_value_fval" and val is not None:
                monoMass = val
                flagMonomass = 0
            if val == "IUPAC Name":
                flagName = 1
            if val == "Preferred":  # Right after "IUPAC Name". PubChem website uses "preferred" name
                flagPreferred = 1
            if flagName == 1 and flagPreferred == 1 and tag == "PC-InfoData_value_sval" and val is not None:
                name = val
                flagName, flagPreferred = 0, 0
            elem.clear()
        conn.commit()
        print("\n")
        os.remove(gzFile)
        os.remove(xmlFile)
    except ImportError:
        print("Error in import lxml package")

conn.close()

