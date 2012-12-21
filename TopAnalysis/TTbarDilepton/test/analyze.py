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
    "Summer12-ZJets", "Summer12-ZJets10To50", "Summer12-WJets",
    "Summer12-TTJets", "Summer12-TTbarTuneZ2",
    "Summer12-SingleToptW", "Summer12-SingleTopBartW",
    "Summer12-WW", "Summer12-WZ", "Summer12-ZZ",
]

samples["ee"] = [
    "Run2012A-DoubleElectron", "Run2012B-DoubleElectron", "Run2012C-DoubleElectron",
    "Summer12-ZJets", "Summer12-ZJets10To50", "Summer12-WJets",
    "Summer12-TTJets", "Summer12-TTbarTuneZ2",
    "Summer12-SingleToptW", "Summer12-SingleTopBartW",
    "Summer12-WW", "Summer12-WZ", "Summer12-ZZ",
]

samples["me"] = [
    "Run2012A-MuEG", "Run2012B-MuEG", "Run2012C-MuEG",
    "Summer12-ZJets", "Summer12-ZJets10To50", "Summer12-WJets",
    "Summer12-TTJets", "Summer12-TTbarTuneZ2",
    "Summer12-SingleToptW", "Summer12-SingleTopBartW",
    "Summer12-WW", "Summer12-WZ", "Summer12-ZZ",
]

from TopAnalysis.TTbarDilepton.PlotTool import *

def runModule(srcDir, sample, mode, verbose):
    module = TTbarDileptonNtupleAnalyzer("%s/%s.root" % (srcDir, sample), mode, "hist/%s_%s.root" % (sample, mode))
    module.analyze(verbose)
    module.endJob(verbose)

def runVallot(sample, mode, verbose):
    if verbose > 0: print "Starting to analyze %s-%s" % (sample, mode)

    vallotModeMap = {"ee":"ElEl", "mm":"MuMu", "me":"MuEl"}
    inputTreePath = "vallot/%s.root:%s/tree" % (sample, vallotModeMap[mode])
    outputFilePath = "hist/%s-%s.root" % (sample, mode)
    hNEventPath = "%s/EventSummary" % vallotModeMap[mode]
    if 'Run' in sample:
        ana = NtupleAnalyzerLite(inputTreePath, outputFilePath, hNEventPath, 1)
    else:
        precut = ""
        if sample == "Summer12-TTbarTuneZ2": precut = "ttbarGen.diLeptonic_ == 0"
        ana = NtupleAnalyzerLite(inputTreePath, outputFilePath, hNEventPath, hNEventPath+":1", weightVar="puweight", precut=precut)

    cuts = {
        "ee":["ZMass > 12 && PairSign < 0", "lep1.relIso03() < 0.15 && lep2.relIso03() < 0.15", "abs(ZMass-91.2) > 15", "MET > 30", "nJet30>=2", "nbjets30_CSVM >= 1"],
        "mm":["ZMass > 12 && PairSign < 0", "lep1.relIso03() < 0.15 && lep2.relIso03() < 0.15", "abs(ZMass-91.2) > 15", "MET > 30", "nJet30>=2", "nbjets30_CSVM >= 1"],
        "me":["ZMass > 12 && PairSign < 0", "lep1.relIso03() < 0.15 && lep2.relIso03() < 0.15", "abs(ZMass-91.2) > -1", "MET >  0", "nJet30>=2", "nbjets30_CSVM >= 1"],
    }
    cutsQCD = {
        "ee":["ZMass > 12 && PairSign > 0", "lep1.relIso03() > 0.15 && lep2.relIso03() > 0.15", "abs(ZMass-91.2) > 15", "MET > 30", "nbjets30_CSVM > 0"],
        "mm":["ZMass > 12 && PairSign > 0", "lep1.relIso03() > 0.15 && lep2.relIso03() > 0.15", "abs(ZMass-91.2) > 15", "MET > 30", "nbjets30_CSVM > 0"],
        "me":["ZMass > 12 && PairSign > 0", "lep1.relIso03() > 0.15 && lep2.relIso03() > 0.15", "abs(ZMass-91.2) > -1", "MET > 0", "nbjets30_CSVM > 0"],
    }
    cutsDY = {
        "ee":["ZMass > 12 && PairSign < 0", "lep1.relIso03() < 0.15 && lep2.relIso03() < 0.15", "abs(ZMass-91.2) < 15", "MET > 30", "nbjets30_CSVM > 0"],
        "mm":["ZMass > 12 && PairSign < 0", "lep1.relIso03() < 0.15 && lep2.relIso03() < 0.15", "abs(ZMass-91.2) < 15", "MET > 30", "nbjets30_CSVM > 0"],
        "me":["ZMass > 12 && PairSign < 0", "lep1.relIso03() < 0.15 && lep2.relIso03() < 0.15", "abs(ZMass-91.2) < 15", "MET > 0", "nbjets30_CSVM > 0"],
    }

    hists = [
      ["hNVertex", "hLL_m", "hNJets", "hMET_pt", "hNBJet", "hGenTT_m", "hVsumTT_m", "hKinsTT_m",],
      ["hNVertex", "hLL_m", "hNJets", "hMET_pt", "hNBJet", "hGenTT_m", "hVsumTT_m", "hKinsTT_m",],
      ["hNVertex", "hLL_m", "hNJets", "hMET_pt", "hNBJet", "hGenTT_m", "hVsumTT_m", "hKinsTT_m",],
      ["hNVertex", "hLL_m", "hNJets", "hMET_pt", "hNBJet", "hGenTT_m", "hVsumTT_m", "hKinsTT_m",],
      ["hNVertex", "hLL_m", "hNJets", "hMET_pt", "hNBJet", "hGenTT_m", "hVsumTT_m", "hKinsTT_m", "hGenTT_m__KinsTT_m", "hGenTT_m__VsumTT_m",],
      ["hNVertex", "hLL_m", "hNJets", "hMET_pt", "hNBJet", "hGenTT_m", "hVsumTT_m", "hKinsTT_m", "hGenTT_m__KinsTT_m", "hGenTT_m__VsumTT_m"],
      ["hNVertex", "hLL_m", "hNJets", "hMET_pt", "hNBJet", "hGenTT_m", "hVsumTT_m", "hKinsTT_m", "hGenTT_m__KinsTT_m", "hGenTT_m__VsumTT_m"],
    ]

    ana.addH1("ZMass", "hLL_m", "Dilepton mass;Dilepton mass (GeV/c^{2});Events", 100, 0, 500)
    ana.addH1("kinttbarM", "hKinsTT_m", "Kinematic t#bar{t} mass;t#bar{t} mass (GeV/c^{2});Events", 400, 0, 2000)
    ana.addH1("ttbar.M_", "hVsumTT_m", "Vector sum t#bar{t} mass;t#bar{t} mass (GeV/c^{2});Events", 400, 0, 2000)
    ana.addH1("nJet30", "hNJets", "Jet multiplicity;Jet multiplicity;Events", 10, 0, 10)
    ana.addH1("MET", "hMET_pt", "Missing transverse momentum;Missing transverse momentum (GeV/c);Events", 100, 0, 500)
    ana.addH1("nvertex", "hNVertex", "Number of primary vertices;Number of vertices;Events", 30, 0, 30)
    ana.addH1("nbjets30_CSVM", "hNBJet", "Number of b tag;Number of b tagged jets;Events", 10, 0, 10)
    ana.addH1("ttbarGen.mass_", "hGenTT_m", "Generator level t#bar{t} mass;t#bar{t} mass (GeV/c^{2});Events", 400, 0, 2000)
    ana.addH2("ttbarGen.mass_:kinttbarM", "hGenTT_m__KinsTT_m", "Response matrix;Generator level t#bar{t} mass;t#bar{t} mass (GeV/c^{2});Kinematic solution t#bar{t} mass", 400, 0, 2000, 400, 0, 2000)
    ana.addH2("ttbarGen.mass_:ttbar.M_", "hGenTT_m__VsumTT_m", "Response matrix;Generator level t#bar{t} mass;t#bar{t} mass (GeV/c^{2});Vector sum t#bar{t} mass", 400, 0, 2000, 400, 0, 2000)

    for i in range(len(cuts[mode])):
        ana.addCategory("Step %d" % i, cuts[mode][i], hists[i])
        #ana.addCategory("QCD %d" % i, " && ".join(cutsQCD[mode][:i+1]))
        #ana.addCategory("DY %d" % i, " && ".join(cutsDY[mode][:i+1]))
    #ana.addCategory("Step 7", "", hists[0], ["renewCut"])

    ana.scanCutFlow()

    if verbose > 0:
        print "Finished analysis of %s-%s" % (sample, mode)

procs = []
for mode in ["mm", "ee", "me"]:
    for sample in samples[mode]:
#        proc = Process(target=runModule, args=(srcDir, sample, mode, 1))
        if 'Run' in sample:
            proc = Process(target=runVallot, args=(sample, mode, 1))
        else:
            proc = Process(target=runVallot, args=(sample, mode, 1))
        procs.append(proc)

        proc.start()

for proc in procs:
    proc.join()

#sys.path.append(".")
#from draw import *

