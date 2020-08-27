#!/usr/bin/python

# Generate a CSV file for HMDB database
# This file is for using MetFragCL. Even though it is said that MetFragCL can accept SDF file format, it does not work well
# Currently, only CSV file format works well
import sys, os, pandas as pd

##################
# Create a table #
##################
# xmlFile = sys.argv[1]
# XML is from /https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip (All metabolites)
xmlFile = r"/Research/Projects/7Metabolomics/Database/HMDB/hmdb_metabolites.xml"

###############################################################################
# Read a XML file and create a pandas DataFrame (to be written to a CSV file) #
###############################################################################
hmdbDict = {"Identifier": [], "CompoundName": [], "MolecularFormula": [], "MonoisotopicMass": [],
            "IUPAC_Name": [], "Traditional_IUPAC_Name": [], "CAS_RN": [], "SMILES": [],
            "InChI": [], "InChIKey": []}
n = 0
try:
    # Reading a XML  file
    from lxml import etree
    for event, element in etree.iterparse(xmlFile, tag = "{*}metabolite"):
        for child in element:
            tag = etree.QName(child).localname  # Remove namespace from the tag
            val = child.text
            if tag == "accession":
                hmdbDict["Identifier"].append(val)
                n += 1
                if n % 1000 == 0:
                    text = "\r  {0} entries are parsed and written to a database".format(n)
                    sys.stdout.write(text)
                    sys.stdout.flush()
            elif tag == "name":
                hmdbDict["CompoundName"].append(val)
            elif tag == "chemical_formula":
                hmdbDict["MolecularFormula"].append(val)
            elif tag == "monisotopic_molecular_weight":
                hmdbDict["MonoisotopicMass"].append(val)
            elif tag == "iupac_name":
                hmdbDict["IUPAC_Name"].append(val)
            elif tag == "traditional_iupac":
                hmdbDict["Traditional_IUPAC_Name"].append(val)
            elif tag == "cas_registry_number":
                hmdbDict["CAS_RN"].append(val)
            elif tag == "smiles":
                hmdbDict["SMILES"].append(val)
            elif tag == "inchi":
                hmdbDict["InChI"].append(val)
            elif tag == "inchikey":
                hmdbDict["InChIKey"].append(val)
        element.clear()
except ImportError:
    print("Error in import lxml package")

df = pd.DataFrame.from_dict(hmdbDict, orient="columns")
df = df[~df["MonoisotopicMass"].isna()]
csvFile = os.path.splitext(xmlFile)[0] + ".csv"
df.to_csv(csvFile, index = False)