import re, sys, os, numpy as np


def getParams(paramFile):
    parameters = dict()
    with open(paramFile, 'r') as file:
        for line in file:
            if re.search(r'^#', line) or re.search(r'^\s', line):
                continue
            line = re.sub(r'#.*', '', line)  # Remove comments (start from '#')
            line = re.sub(r'\s*', '', line)  # Remove all whitespaces

            # Exception for "feature_files" parameter
            if "feature_files" in parameters and line.endswith("feature"):
                parameters["feature_files"].append(line)
            else:
                key = line.split('=')[0]
                val = line.split('=')[1]
                if key == "feature_files":
                    parameters[key] = [val]
                else:
                    parameters[key] = val
    return parameters


def readFeatures(featureFile):
    with open(featureFile, 'r') as file:
        features = file.readlines()
    return features


class progressBar:
    def __init__(self, total):
        self.total = total
        self.barLength = 20
        self.count = 0
        self.progress = 0
        self.block = 0
        self.status = ""

    def increment(self):
        self.count += 1
        self.progress = self.count / self.total
        self.block = int(round(self.barLength * self.progress))
        if self.progress == 1:
            self.status = "Done...\r\n"
        else:
            self.status = ""
        #         self.status = str(self.count) + "/" + str(self.total)
        text = "\r  Progress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block),
                                                   int(self.progress * 100), self.status)
        sys.stdout.write(text)
        sys.stdout.flush()

