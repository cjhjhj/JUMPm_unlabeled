#!/usr/bin/python

# Build HMDB database using SQLite
import sqlite3
conn = sqlite3.connect('hmdb.db')   # SQLite database for HMDB

##################
# Create a table #
##################
conn.execute('''CREATE TABLE HMDB
                (ID INTEGER PRIMARY KEY NOT NULL,
                HMDBID VARCHAR(255) NOT NULL,
                FORMULA VARCHAR(255) NOT NULL,
                SMILES VARCHAR(255),
                INCHIKEY VARCHAR(255),
                NAME VARCHAR(255),
                MONOMASS REAL NOT NULL);''')

################################################
# Read a XML file and write to a database file #
################################################
inputFile = "hmdb_metabolites.xml"
n = 0
try:
    # Reading a XML  file
    from lxml import etree
    for event, element in etree.iterparse(inputFile, tag = "{*}metabolite"):
        for child in element:
            tag = etree.QName(child).localname  # Remove namespace from the tag
            val = child.text
            if tag == "accession" and val is not None:
                hmdbId = val
                n += 1
                if n % 1000 == 0:
                    print ("%d entries are parsed and written to a database" % n)
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
            elif tag == "name" and val is not None:
                name = val
        conn.execute('INSERT INTO HMDB (HMDBID, FORMULA, SMILES, INCHIKEY, NAME, MONOMASS) VALUES (?, ?, ?, ?, ?, ?)',
                     (hmdbId, formula, smiles, inchiKey, name, monoMass))
        element.clear()
    conn.commit()
    conn.close()
except ImportError:
    print ("Error in import lxml package")
