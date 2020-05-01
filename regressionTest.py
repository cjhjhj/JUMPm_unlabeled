#!/usr/bin/python

import os, sys, numpy as np, pandas as pd
from scipy.interpolate import interp1d
import statsmodels.api as sm


df = pd.read_csv("rt.txt", sep = "\t", header = None)
ind = np.where(~np.isnan(df[1]))
x = df.iloc[ind][0].to_numpy()
y = df.iloc[ind][1].to_numpy()
xnew = df[0].to_numpy()

# Kernel regression
model1 = sm.nonparametric.KernelReg(y, x, "c", bw='aic')
res1 = model1.fit()

# # LOWESS regression
# model2 = sm.nonparametric.lowess(y, x)
# res2 = model2.fit()
print()