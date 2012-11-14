#!/usr/bin/env python

from ROOT import *
import sys, os
from multiprocessing import Process

sys.path.append("../python")
gROOT.ProcessLine(".x rootlogon.C")
gSystem.CompileMacro("../src/TTbarDileptonNtupleAnalyzer.cc")

#from sampleInfo import *
srcDir = "/users/jhgoh/data/CMS/ntuple/Generic/20121112_1"
samples = [
    "Run2012A-DoubleElectron", "Run2012A-DoubleMu", "Run2012A-MuEG",
    "Run2012B-DoubleElectron", "Run2012B-DoubleMu", "Run2012B-MuEG",
    "Run2012C-DoubleElectron", "Run2012C-DoubleMu", "Run2012C-MuEG",
    "Summer12-ZJets",
    "Summer12-ZJets10To50",
    "Summer12-WJets",
    "Summer12-TTJets",
    "Summer12-SingleTop",
    "Summer12-WW",
    "Summer12-WZ",
    "Summer12-ZZ",
]

modules = []
for sample in samples:
    module = TTbarDileptonNtupleAnalyzer("%s/ntuple_%s.root" % (srcDir, sample), "hist/hist_%s.root" % sample)
    modules.append(module)

def runModule(module):
    module.analyze(1)
    module.endJob(1)

procs = []
for module in modules:
    proc = Process(target=runModule, args=(module,))
    proc.start()

    procs.append(proc)

#while procs:
#    procs.pop().join()

#for module in modules:
#    module.analyze(1)
#    module.endJob(1)
