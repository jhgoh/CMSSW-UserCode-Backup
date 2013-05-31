#!/usr/bin/env python

from TopAnalysis.TTbarDilepton.dataset_cff import *
import sys, os
from datetime import datetime

doSubmit = True
cmgVersion = "V5_13_0"
#cmgVersion = "V5_12_0_44X"
maxFiles = 50

datasets = [
    "DoubleElectron-Run2012A", "DoubleElectron-Run2012B", "DoubleElectron-Run2012C", "DoubleElectron-Run2012D",
    "DoubleMu-Run2012A", "DoubleMu-Run2012B", "DoubleMu-Run2012C", "DoubleMu-Run2012D",
    "MuEG-Run2012A", "MuEG-Run2012B", "MuEG-Run2012C", "MuEG-Run2012D",

    "QCDMu-Summer12", 
    "QCDEM20To30-Summer12", "QCDEM30To80-Summer12", "QCDEM80To170-Summer12", 

    "SingleTop-Summer12", 
    "WJets-Summer12", 
    "ZJets-Summer12", "ZJets10To50-Summer12", 
    "WW-Summer12", "WZ-Summer12", "ZZ-Summer12", 

    "TTJets-Summer12", 
    "TTJetsLL-Summer12", 

    "TTJetsM161-Summer12", "TTJetsM163-Summer12", "TTJetsM166-Summer12", "TTJetsM169-Summer12", 
    "TTJetsM175-Summer12", "TTJetsM178-Summer12", "TTJetsM181-Summer12", "TTJetsM184-Summer12", 

    "TTJetsMatchingdown-Summer12", "TTJetsMatchingup-Summer12", 
    "TTJetsScaledown-Summer12", "TTJetsScaleup-Summer12", 

    "TTJetsMCNLO-Summer12", "TTJetsPowheg-Summer12", 

    "TTH-Summer12", 
] 

runScript = open("run.sh", "w")
print>>runScript, """#!/bin/bash

cd %s
eval `scram runtime -sh`

export CMGVERSION=%s

export DATASET=$1
export SECTION=$2
export MAXFILES=$3

cmsRun ntuple_cfg.py

""" % (os.getcwd(), cmgVersion)
runScript = None

os.system("chmod +x run.sh")

if not os.path.isdir("log"): os.mkdir("log")
if not os.path.isdir("ntuple"): os.mkdir("ntuple")
if not os.path.isdir("ntuple/unmerged"): os.mkdir("ntuple/unmerged")

submitLog = open("log/submit.log", "w")
print>>submitLog, "Submitting jobs", datetime.now()
for dataset in datasets:
    files = loadDataset(cmgVersion, dataset)
    nFiles = len(files)
    nJobs = nFiles/maxFiles
    if nJobs*maxFiles < nFiles: nJobs += 1

    print "Submitting", dataset
    print "    nFiles   =", nFiles
    print "    nJobs    =", nJobs
    print "    maxFiles =", maxFiles

    print>>submitLog, "Submitting", dataset
    print>>submitLog, "    nFiles   =", nFiles
    print>>submitLog, "    nJobs    =", nJobs
    print>>submitLog, "    maxFiles =", maxFiles

    for section in range(nJobs):
        cmd = "bsub -q 8nh -oo log/%s_%s.log run.sh %s %s %s" % (dataset, section, dataset, section, maxFiles)
        print cmd
        print>>submitLog, cmd
        if doSubmit: os.system(cmd)
print>>submitLog, "Submit done", datetime.now()
submitLog = None
