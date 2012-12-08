#!/usr/bin/env python

import sys, os
from multiprocessing import Process

sys.argv.append('-b')
from ROOT import *

sys.path.append("../python")
gROOT.ProcessLine(".x rootlogon.C")
#gSystem.CompileMacro("../src/TTbarDileptonNtupleAnalyzer.cc")

#from sampleInfo import *
srcDir = "ntuple/merged"
samples = {}

samples["mm"] = [
    "Run2012A-DoubleMu", "Run2012B-DoubleMu", "Run2012C-DoubleMu", 
    "ZJets", "ZJets10To50", "WJets",
    "TTJets",
    "SingleToptW", "SingleTopBartW",
    "WW", "WZ", "ZZ",
]

samples["ee"] = [
    "Run2012A-DoubleElectron", "Run2012B-DoubleElectron", "Run2012C-DoubleElectron", 
    "ZJets", "ZJets10To50", "WJets",
    "TTJets",
    "SingleToptW", "SingleTopBartW",
    "WW", "WZ", "ZZ",
]

samples["me"] = [
    "Run2012A-MuEG", "Run2012B-MuEG", "Run2012C-MuEG", 
    "ZJets", "ZJets10To50", "WJets",
    "TTJets",
    "SingleToptW", "SingleTopBartW",
    "WW", "WZ", "ZZ",
]

from TopAnalysis.TTbarDilepton.PlotTool import *

def runModule(srcDir, sample, mode, verbose):
    if verbose > 0:
        print "Starting to analyze %s-%s" % (sample, mode)
    #module = TTbarDileptonNtupleAnalyzer("%s/%s.root" % (srcDir, sample), mode, "hist/%s_%s.root" % (sample, mode))
    #module.analyze(verbose)
    #module.endJob(verbose)

    vallotModeMap = {"ee":"ElEl", "mm":"MuMu", "me":"MuEl"}
    ana = NtupleAnalyzerLite("Run2012", "Summer12", mode, sample, "vallot", "%s/EventSummary" % vallotModeMap[mode], "%s/tree" % vallotModeMap[mode])

    cuts = {
        "ee":["ZMass > 12 && PairSign < 0", "isIso > 0", "abs(ZMass-91.2) > 15", "nJet30 >= 2", "MET > 30", "nbjets30_CSVM > 0"],
        "mm":["ZMass > 12 && PairSign < 0", "isIso > 0", "abs(ZMass-91.2) > 15", "nJet30 >= 2", "MET > 30", "nbjets30_CSVM > 0"],
        "me":["ZMass > 12 && PairSign < 0", "isIso > 0", "abs(ZMass-91.2) > -1", "nJet30 >= 2", "MET > 20", "nbjets30_CSVM > 0"],
    }                                                                                      
    cutsQCD = {                                                                            
        "ee":["ZMass > 12 && PairSign > 0", "isIso < 0", "abs(ZMass-91.2) > 15", "nJet30 >= 2", "MET > 30", "nbjets30_CSVM > 0"],
        "mm":["ZMass > 12 && PairSign > 0", "isIso < 0", "abs(ZMass-91.2) > 15", "nJet30 >= 2", "MET > 30", "nbjets30_CSVM > 0"],
        "me":["ZMass > 12 && PairSign > 0", "isIso < 0", "abs(ZMass-91.2) > -1", "nJet30 >= 2", "MET > 20", "nbjets30_CSVM > 0"],
    }                                                                                      
    cutsDY = {                                                                             
        "ee":["ZMass > 12 && PairSign < 0", "isIso > 0", "abs(ZMass-91.2) < 15", "nJet30 >= 2", "MET > 30", "nbjets30_CSVM > 0"],
        "mm":["ZMass > 12 && PairSign < 0", "isIso > 0", "abs(ZMass-91.2) < 15", "nJet30 >= 2", "MET > 30", "nbjets30_CSVM > 0"],
        "me":["ZMass > 12 && PairSign < 0", "isIso > 0", "abs(ZMass-91.2) < 15", "nJet30 >= 2", "MET > 20", "nbjets30_CSVM > 0"],
    }
    for i in range(len(cuts[mode])):
        ana.addCategory("Step %d" % i, " && ".join(cuts[mode][:i+1]))
        ana.addCategory("QCD Step %d" % i, " && ".join(cutsQCD[mode][:i+1]))
        ana.addCategory("DY Step %d" % i, " && ".join(cutsDY[mode][:i+1])) 

    ana.addHistogram("ZMass", "hLL_m", "Dilepton mass;Dilepton mass (GeV/c^{2});Events", 100, 0, 500)
    ana.addHistogram("kinttbarM", "hKinsTT_m", "Kinematic t#bar{t} mass;t#bar{t} mass (GeV/c^{2});Events", 100, 0, 2000)
    ana.addHistogram("ttbar.M_", "hVsumTT_m", "Vector sum t#bar{t} mass;t#bar{t} mass (GeV/c^{2});Events", 100, 0, 2000)
    ana.addHistogram("nJet30", "hNJets", "Jet multiplicity;Jet multiplicity;Events", 10, 0, 10)
    ana.addHistogram("MET", "hMET_pt", "Missing transverse momentum;Missing transverse momentum (GeV/c);Events", 100, 0, 500)
    ana.addHistogram("nvertex", "hNVertex", "Number of primary vertices;Number of vertices;Events", 30, 0, 30)
    ana.addHistogram("nbjets30_CSVM", "hNBJet", "Number of b tag;Number of b tagged jets;Events", 10, 0, 10)

    ana.scanCutFlow()

    if verbose > 0:
        print "Finished analysis of %s-%s" % (sample, mode)

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

