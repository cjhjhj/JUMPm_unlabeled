#!/usr/bin/python

# Build HMDB database using SQLite
import sys, os, sqlite3

##################
# Create a table #
##################
# xmlFile = sys.argv[1]
# XML is from /https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip (All metabolites)
xmlFile = r"/Research/Projects/7Metabolomics/Database/HMDB/hmdb_metabolites.xml"
dbName = os.path.splitext(xmlFile)[0] + ".db"
conn = sqlite3.connect(dbName)
createQuery = "CREATE TABLE hmdb (id TEXT, other_ids TEXT, name TEXT, synonym TEXT," \
              "formula TEXT, mass REAL, smiles TEXT, inchikey TEXT)"
conn.execute(createQuery)

################################################
# Read a XML file and write to a database file #
################################################
recordList = []
n = 0
try:
    # Reading a XML  file
    from lxml import etree
    for event, element in etree.iterparse(xmlFile, tag = "{*}metabolite"):
        for child in element:
            tag = etree.QName(child).localname  # Remove namespace from the tag
            val = child.text
            if tag == "accession" and val is not None:
                uid = val   # HMDB ID (unique identifier)
                n += 1
                if n % 1000 == 0:
                    text = "\r  {0} entries are parsed and written to a database".format(n)
                    sys.stdout.write(text)
                    sys.stdout.flush()
            elif tag == "synonyms":
                synonym = [c.text for c in child]
                synonym = ";".join(synonym)
            elif tag == "chemical_formula":
                if val is None:
                    formula = "NA"
                else:
                    formula = val
            elif tag == "smiles":
                if val is None:
                    smiles = "NA"
                else:
                    smiles = val
            elif tag == "inchikey":
                if val is None:
                    inchikey = "NA"
                else:
                    inchikey = val
            elif tag == "monisotopic_molecular_weight":
                mass = val  # Monoisotopic (neutral) mass
            elif tag == "name":
                if val is None:
                    name = "NA"
                else:
                    name = val
            elif tag == "kegg_id":
                if val is None:
                    kegg = "NA"
                else:
                    kegg = val
            elif tag == "chemspider_id":
                if val is None:
                    chemspider= "NA"
                else:
                    chemspider = val
            elif tag == "metlin_id":
                if val is None:
                    metlin = "NA"
                else:
                    metlin = val
            elif tag == "pubchem_compound_id":
                if val is None:
                    pcid = "NA"
                else:
                    pcid = val
        otherId = ";".join([kegg, chemspider, metlin, pcid])
        recordList.append((uid, otherId, name, synonym, formula, mass, smiles, inchikey))
        element.clear()

    insertQuery = "INSERT INTO hmdb (id, other_ids, name, synonym, formula, mass, smiles, inchikey) VALUES " \
                  "(?, ?, ?, ?, ?, ?, ?, ?)"
    conn.executemany(insertQuery, recordList)
    conn.commit()
    conn.close()
except ImportError:
    print ("Error in import lxml package")
