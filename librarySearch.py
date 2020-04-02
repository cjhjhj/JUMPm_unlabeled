#!/usr/bin/python

import sys, os, sqlite3, numpy as np, pandas as pd
import utils

##################################
# Load parameters and initialize #
##################################
paramFile = r"C:\Research\Projects\7Metabolomics\IROA_HILIC_NEG_Target\jumpm_negative.params"
params = utils.getParams(paramFile)
condition = params["LC_column"].lower()
if params["mode"] == "1":
    condition = condition + "p"
elif params["mode"] == "-1":
    condition = condition + "n"
else:
    sys.exit("'mode' parameter should be either 1 or -1")
proton = 1.007276466812
H = 1.0078250321
matchMzTol = 10 # Unit of ppm
matchRtTol = 10 # Unit of second

#############################
# Open sqlite-based library #
#############################
# conn = sqlite3.connect("flights.db")

############################
# Load feature information #
############################


print ()