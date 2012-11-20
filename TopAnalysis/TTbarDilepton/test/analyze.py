#!/usr/bin/env python

from ROOT import *
import sys, os
from multiprocessing import Process

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

def runModule(module):
    module.analyze(1)
    module.endJob(1)

procs = []
for mode in ["mm", "ee", "me"]:
    for sample in samples[mode]:
        module = TTbarDileptonNtupleAnalyzer("%s/%s.root" % (srcDir, sample), mode, "hist/%s_%s.root" % (sample, mode))

        runModule(module)
#        proc = Process(target=runModule, args=(module,))
#        procs.append(proc)

#        proc.start()

#for proc in procs:
#    proc.join()

