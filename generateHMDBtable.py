#!/usr/bin/python

# Build HMDB database using SQLite
import sys, os
from lxml import etree

################################################
# Read a XML file and write to a database file #
################################################
inputFile = "hmdb_small.xml"
n = 0
outputFile = "hmdb.txt"
with open(outputFile, "w") as f:
    f.write("HMDBID\tName\tFormula\tMonoisotopicMass\tSMILES\tInChIKey\n")
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
                elif tag == "name" and val is not None:
                    name = val
                elif tag == "chemical_formula" and val is not None:
                    formula = val
                elif tag == "monisotopic_molecular_weight" and val is not None:
                    monoMass = val
                elif tag == "smiles" and val is not None:
                    smiles = val
                elif tag == "inchi" and val is not None:
                    inchi = val
                elif tag == "inchikey" and val is not None:
                    inchiKey = val
            line = u'\t'.join((hmdbId, name, formula, monoMass, smiles, inchiKey)).encode('utf-8')
            f.write(line)
            # f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (hmdbId, name, formula, monoMass, smiles, inchiKey))
            element.clear()
    except ImportError:
        print ("Error in import lxml package")
