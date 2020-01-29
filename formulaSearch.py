#!/usr/bin/python

import h5py

# Read the database file
db = h5py.File('HMDB.h5', 'r')

# Query conditions
queryMass = 110.392 # Should be a neutral mass
queryMassTol = 10 # Unit of ppm
lowerLimit = queryMass - queryMass * queryMassTol / 1e6
upperLimit = queryMass + queryMass * queryMassTol / 1e6
lowerIndex = int(lowerLimit)
upperIndex = int(upperLimit)
