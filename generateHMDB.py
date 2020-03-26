#!/usr/bin/python

# Build HMDB database using SQLite
import sys, os, sqlite3
from lxml import etree

##################
# Create a table #
##################
conn = sqlite3.connect('hmdb.db')   # SQLite database for HMDB
conn.execute('''CREATE TABLE HMDB
                (ID INTEGER PRIMARY KEY NOT NULL,
                HMDBID VARCHAR(255) NOT NULL,
                FORMULA VARCHAR(255) NOT NULL,
                SMILES VARCHAR(255),
                INCHIKEY VARCHAR(255),
                NAME VARCHAR(255),
                MONOMASS REAL NOT NULL,
                TYPE VARCHAR(255));''')

################################################
# Read a XML file and write to a database file #
################################################
inputFile = "hmdb_metabolites.xml"
n = 0
H = 1.0078250321    # Mass of hydrogen
try:
    # Reading a XML file
    for event, element in etree.iterparse(inputFile, tag = "{*}metabolite"):
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
        conn.execute('INSERT INTO HMDB (HMDBID, FORMULA, SMILES, INCHIKEY, NAME, MONOMASS, TYPE) VALUES (?, ?, ?, ?, ?, ?, ?)',
                     (hmdbId, formula, smiles, inchiKey, name, monoMass, "target"))
        conn.execute('INSERT INTO HMDB (HMDBID, FORMULA, SMILES, INCHIKEY, NAME, MONOMASS, TYPE) VALUES (?, ?, ?, ?, ?, ?, ?)',
                     (hmdbId, formula, smiles, inchiKey, name, decoyMass, "decoy"))
        element.clear()
    conn.commit()
    conn.close()
except ImportError:
    print ("Error in import lxml package")
