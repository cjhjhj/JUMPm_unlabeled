#!/usr/bin/python

# Initiation of hdf5 format database
# import h5py
# db = h5py.File("PubChem.h5", "w")
# dt = h5py.special_dtype(vlen = str)

inputFile = "Compound_050000001_050500000.xml"
n = 0
try:
    from lxml import etree
    for _, element in etree.iterparse(inputFile, tag = "{*}PC-Compound"):
        children = list(element);
        print (children);
        
#         for child in element:
#             tag = etree.QName(child).localname  # Remove namespace from the tag
#             val = child.text
#             print (tag)
#             if tag == "PC-CompoundType_id_cid" and val is not None:
#                 pubChemId = val
#                 n += 1
#                 if n % 1000 == 0:
#                     print ("%d entries are parsed" % n)
#             elif tag == "chemical_formula" and val is not None:
#                 formula = val
#             elif tag == "smiles" and val is not None:
#                 smiles = val
#             elif tag == "inchi" and val is not None:
#                 inchi = val
#             elif tag == "inchikey" and val is not None:
#                 inchiKey = val
#             elif tag == "monisotopic_molecular_weight" and val is not None:
#                 monoMass = val
#             elif tag == "name" and val is not None:
#                 name = val
        
#         # Write to a hdf5 file
#         intMonoMass = str(int(float(monoMass)))
#         if intMonoMass not in db:
#             # Create a new group and add a metabolite entry
#             db.create_group(intMonoMass)
#             db[intMonoMass].create_dataset('hmdbId', (1, ), dtype = dt, maxshape = (None, ))
#             db[intMonoMass]['hmdbId'][0] = hmdbId
#             db[intMonoMass].create_dataset('formula', (1, ), dtype = dt, maxshape = (None, ))
#             db[intMonoMass]['formula'][0] = formula
#             db[intMonoMass].create_dataset('smiles', (1, ), dtype = dt, maxshape = (None, ))
#             db[intMonoMass]['smiles'][0] = smiles
#             db[intMonoMass].create_dataset('inchi', (1, ), dtype = dt, maxshape = (None, ))
#             db[intMonoMass]['inchi'][0] = inchi
#             db[intMonoMass].create_dataset('inchiKey', (1, ), dtype = dt, maxshape = (None, ))
#             db[intMonoMass]['inchiKey'][0] = inchiKey
#             db[intMonoMass].create_dataset('monoMass', (1, ), dtype = "float64", maxshape = (None, ))
#             db[intMonoMass]['monoMass'][0] = float(monoMass)
#             db[intMonoMass].create_dataset('name', (1, ), dtype = dt, maxshape = (None, ))
#             db[intMonoMass]['name'][0] = name
#         else:
#             # Add entry to existing datasets
#             m = len(db[intMonoMass]['hmdbId'])
#             db[intMonoMass]['hmdbId'].resize(db[intMonoMass]['hmdbId'].shape[0] + 1, axis = 0)
#             db[intMonoMass]['hmdbId'][m] = hmdbId
#             db[intMonoMass]['formula'].resize(db[intMonoMass]['formula'].shape[0] + 1, axis = 0)
#             db[intMonoMass]['formula'][m] = formula
#             db[intMonoMass]['smiles'].resize(db[intMonoMass]['smiles'].shape[0] + 1, axis = 0)
#             db[intMonoMass]['smiles'][m] = smiles
#             db[intMonoMass]['inchi'].resize(db[intMonoMass]['inchi'].shape[0] + 1, axis = 0)
#             db[intMonoMass]['inchi'][m] = inchi
#             db[intMonoMass]['inchiKey'].resize(db[intMonoMass]['inchiKey'].shape[0] + 1, axis = 0)
#             db[intMonoMass]['inchiKey'][m] = inchiKey
#             db[intMonoMass]['monoMass'].resize(db[intMonoMass]['monoMass'].shape[0] + 1, axis = 0)
#             db[intMonoMass]['monoMass'][m] = float(monoMass)
#             db[intMonoMass]['name'].resize(db[intMonoMass]['name'].shape[0] + 1, axis = 0)
#             db[intMonoMass]['name'][m] = name
        element.clear()
except ImportError:
    print ("Error in import lxml package")
    
