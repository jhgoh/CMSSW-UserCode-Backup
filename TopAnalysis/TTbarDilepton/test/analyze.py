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
    "ZJets", "ZJets10To50", "WJets",
    "TTJets",
    "SingleTop",
    "WW", "WZ", "ZZ",
]

samples["ee"] = [
    "Run2012A-DoubleElectron", "Run2012B-DoubleElectron", "Run2012C-DoubleElectron", 
    "ZJets", "ZJets10To50", "WJets",
    "TTJets",
    "SingleTop",
    "WW", "WZ", "ZZ",
]

samples["me"] = [
    "Run2012A-MuEG", "Run2012B-MuEG", "Run2012C-MuEG", 
    "ZJets", "ZJets10To50", "WJets",
    "TTJets",
    "SingleTop",
    "WW", "WZ", "ZZ",
]

from TopAnalysis.TTbarDilepton.PlotTool import *

def runModule(srcDir, sample, mode, verbose):
    #module = TTbarDileptonNtupleAnalyzer("%s/%s.root" % (srcDir, sample), mode, "hist/%s_%s.root" % (sample, mode))
    #module.analyze(verbose)
    #module.endJob(verbose)

    ana = NtupleAnalyzerLite("Run2012", "Summer12", mode, sample, "ntuple/merged", "hEventCounter", "event")
    if mode == "ee":
        ana.addCategory("Step 0", "@electrons.size() >= 2 && electrons_Q[0]+electrons_Q[1] == 0")
        ana.addCategory("Step 1", "@electrons.size() >= 2 && electrons_Q[0]+electrons_Q[1] == 0 && @jets.size() >= 2")
        #ana.addHistogram("(electrons[0]+electrons[1]).mass()", "hLL_m", "Dilepton mass;Dilepton mass (GeV/c^{2});Events", 100, 0, 500)
        ana.addHistogram("electrons[0].pt()", "hElectron1_pt", "Electron1 p_{T};Transverse momentum p_{T} (GeV/c);Events", 100, 0, 500)
    elif mode == "mm":
        ana.addCategory("Step 0", "@muons.size() >= 2 && muons_Q[0]+muons_Q[1] == 0")
        ana.addCategory("Step 1", "@muons.size() >= 2 && muons_Q[0]+muons_Q[1] == 0 && @jets.size() >= 2")
        ana.addHistogram("muons[0].pt()", "hMuon1_pt", "Muon1 p_{T};Transverse momentum p_{T} (GeV/c);Events", 100, 0, 500)
    elif mode == "me":
        ana.addCategory("Step 0", "@muons.size() >= 1 && @electrons.size() >= 1 && electrons_Q[0]+muons_Q[1] == 0")
        ana.addCategory("Step 1", "@muons.size() >= 1 && @electrons.size() >= 1 && electrons_Q[0]+muons_Q[1] == 0 && @jets.size() >= 2")
        ana.addHistogram("muons[0].pt()", "hMuon1_pt", "Muon1 p_{T};Transverse momentum p_{T} (GeV/c);Events", 100, 0, 500)
        ana.addHistogram("electrons[0].pt()", "hElectron1_pt", "Electron1 p_{T};Transverse momentum p_{T} (GeV/c);Events", 100, 0, 500)
        #ana.addHistogram("(electrons[0]+muons[0]).mass()", "hLL_m", "Dilepton mass;Dilepton mass (GeV/c^{2});Events", 100, 0, 500)
        #ana.addHistogram("(muons[0]+muons[1]).mass()", "hLL_m", "Dilepton mass;Dilepton mass (GeV/c^{2});Events", 100, 0, 500)
    ana.addHistogram("@jets.size()", "hNJets", "Jet multiplicity;Jet multiplicity;Events", 10, 0, 10)
    ana.scanCutFlow()

procs = []
for mode in ["mm", "ee", "me"]:
    for sample in samples[mode]:
#        runModule(srcDir, sample, mode, 2)
        proc = Process(target=runModule, args=(srcDir, sample, mode, 1))
        procs.append(proc)

        proc.start()

for proc in procs:
    proc.join()

#sys.path.append(".")
#from draw import *

