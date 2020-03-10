#!/usr/bin/python

import subprocess, numpy as np
import pandas as pd

# rPath = "C:\\Program Files\\R\\R-3.6.2\\bin\\Rscript.exe"
# script = "test.R"
# cmd = [rPath, script, 'a', 'b', 'C', 'd']
# output = subprocess.call(cmd, shell=False)
# print ("")

# import rpy2.robjects as ro
# from rpy2.robjects.vectors import IntVector, FloatVector
# import numpy as np
#
# # ro.r('x = c()')
# # ro.r('x[1] = 22')
# # ro.r('x[2] = 44')
# # print (ro.r('x'))
#
# # r = ro.r
# # ro.globalenv['args'] = ["abc"]
# # r.source("test.R")
#
# # ro.globalenv["args"] = ['abc', 'cde']
# # ro.r.source('test.R')
#
# ro.r.source("./R/featureAlignment.R")
# rLowess = ro.globalenv['loess.as']
# x = np.random.random(100)
# xnew = np.random.random(50)
# y = np.sin(x)
# # x = x.tolist()
# # y = y.tolist()
# rPredict = ro.r('predict')
# model = rLowess(FloatVector(x), FloatVector(y))
# yhat = np.array(model.rx2("fitted"))
# yhat_new = np.array(rPredict(model, FloatVector(xnew)))
# print (yhat_new)


from numpy.lib.recfunctions import append_fields

file = r"U:\Research\Projects\7Metabolomics\JUMPm\IROAsamples\IROA_IS_NEG_1.1.dev.feature"
# data = pd.read_csv(file, sep = "\t")
# data = np.loadtxt(file, delimiter = "\t", skiprows = 1)
data = np.genfromtxt(file, delimiter = "\t", dtype = None, names = True)

ind = np.zeros((data.shape[0], ))
ind[1] = 1
ind[3] = 1

data = append_fields(data, 'Aligned', data = ind)

print()