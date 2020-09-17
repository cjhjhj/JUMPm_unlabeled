#!/usr/bin/python

# Build PubChem database using SQLite
import sys, os, sqlite3, gzip, shutil
from ftplib import FTP

##################
# Create a table #
##################
conn = sqlite3.connect('pubchem.db')   # SQLite database for PubChem
createQuery = "CREATE TABLE pubchem (id TEXT, other_ids TEXT, name TEXT, synonym TEXT," \
              "formula TEXT, mass REAL, smiles TEXT, inchikey TEXT)"
conn.execute(createQuery)

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
    recordList = []
    n = 0
    with open(sdfFile, encoding="utf-8") as f:
        for line in f:  # not using readlines(), as this consumes the memory
            line = line.strip()
            if line.endswith("<PUBCHEM_COMPOUND_CID>"):
                uid = f.readline().strip()
            elif line.endswith("<PUBCHEM_IUPAC_NAME>"):
                name = f.readline().strip()
            elif line.endswith("<PUBCHEM_IUPAC_SYSTEMATIC_NAME>"):
                synonym = f.readline().strip()
            elif line.endswith("<PUBCHEM_IUPAC_INCHIKEY>"):
                inchikey = f.readline().strip()
            elif line.endswith("<PUBCHEM_MOLECULAR_FORMULA>"):
                formula = f.readline().strip()
            elif line.endswith("<PUBCHEM_OPENEYE_CAN_SMILES>"):
                smiles = f.readline().strip()
            elif line.endswith("<PUBCHEM_MONOISOTOPIC_WEIGHT>"):
                mass = f.readline().strip()
            elif line.endswith("$$$$"):
                recordList.append((uid, "NA", name, synonym, formula, mass, smiles, inchikey))
                n += 1
                if n % 1000 == 0:
                    text = "\r  {0} entries are parsed and written to a database".format(n)
                    sys.stdout.write(text)
                    sys.stdout.flush()
    insertQuery = "INSERT INTO pubchem (id, other_ids, name, synonym, formula, mass, smiles, inchikey) VALUES " \
                  "(?, ?, ?, ?, ?, ?, ?, ?)"
    conn.executemany(insertQuery, recordList)
    conn.commit()
    f.close()
    print("\n")
    os.remove(gzFile)
    os.remove(sdfFile)

conn.close()
