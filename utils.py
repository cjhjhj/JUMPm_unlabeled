import re, sys, os, numpy as np


def readParams(paramFile):
    parameters = dict()
    with open(paramFile, 'r') as file:
        for line in file:
            if re.search(r'^#', line) or re.search(r'^\s', line):
                continue
            line = re.sub(r'#.*', '', line)  # Remove comments (start from '#')
            line = re.sub(r'\s*', '', line)  # Remove all whitespaces
            key = line.split('=')[0]
            val = line.split('=')[1]
            parameters[key] = val
    return parameters

def readFeatures(featureFile):
    with open(featureFile, 'r') as file:
        features = file.readlines()
    return features



# def writeFeatures(featureFile):
#     return