import FWCore.ParameterSet.Config as cms
import os

def parseJobSectionOption():
    usage = """Needed environment variables:
    DATASET  : Name of dataset
    SECTION  : Section number
    MAXFILES : Number of files per section
"""

    if 'DATASET' not in os.environ or \
       'SECTION' not in os.environ or \
       'MAXFILES' not in os.environ:
        print usage

    dataset  = os.environ["DATASET"]
    section  = int(os.environ["SECTION"])
    maxFiles = int(os.environ["MAXFILES"])

    return dataset, section, maxFiles

def loadDataset(sample):
    lines = open("%s/src/TopAnalysis/TTbarDilepton/data/dataset-%s.txt" % (os.environ["CMSSW_BASE"], sample)).readlines()
    files = []
    for line in lines:
        line = line.strip()
        if len(line) == 0 or line[0] == '#': continue
        if '.root' not in line: continue
        files.append(line)

    return files

def calculateRange(files, section, maxFiles):
    begin = section*maxFiles
    end = min(len(files), (section+1)*maxFiles)

    if begin > end:
        print "Warning in loadDataset_cff/calculateRange: range exceeds number of flies in the dataset"
        begin = end

    return begin, end

def isRealData(dataset):
    if 'Run20' in dataset:
        return True
    else:
        return False
