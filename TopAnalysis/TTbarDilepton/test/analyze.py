#!/usr/bin/env python

import sys, os
from multiprocessing import Process

sys.argv.append('-b')
from ROOT import *

sys.path.append("../python")
gROOT.ProcessLine(".x rootlogon.C")
gSystem.CompileMacro("../src/TTbarDileptonNtupleAnalyzer.cc")

#from sampleInfo import *
srcDir = "ntuple/merged"
samples = {}

samples["mm"] = [
    "Run2012A-DoubleMu", "Run2012B-DoubleMu", "Run2012C-DoubleMu", 
    "Summer12-ZJets", "Summer12-ZJets10To50", "Summer12-WJets",
    "Summer12-TTJets",
    "Summer12-SingleTop",
    "Summer12-WW", "Summer12-WZ", "Summer12-ZZ",
]

samples["ee"] = [
    "Run2012A-DoubleElectron", "Run2012B-DoubleElectron", "Run2012C-DoubleElectron", 
    "Summer12-ZJets", "Summer12-ZJets10To50", "Summer12-WJets",
    "Summer12-TTJets",
    "Summer12-SingleTop",
    "Summer12-WW", "Summer12-WZ", "Summer12-ZZ",
]

samples["me"] = [
    "Run2012A-MuEG", "Run2012B-MuEG", "Run2012C-MuEG", 
    "Summer12-ZJets", "Summer12-ZJets10To50", "Summer12-WJets",
    "Summer12-TTJets",
    "Summer12-SingleTop",
    "Summer12-WW", "Summer12-WZ", "Summer12-ZZ",
]

def runModule(srcDir, sample, mode, verbose):
    module = TTbarDileptonNtupleAnalyzer("%s/%s.root" % (srcDir, sample), mode, "hist/%s_%s.root" % (sample, mode))
    module.analyze(verbose)
    module.endJob(verbose)

procs = []
for mode in ["mm", "ee", "me"]:
    for sample in samples[mode]:
#        runModule(srcDir, sample, mode, 2)
        proc = Process(target=runModule, args=(srcDir, sample, mode, 1))
        procs.append(proc)

        proc.start()

for proc in procs:
    proc.join()

sys.path.append(".")
from draw import *

