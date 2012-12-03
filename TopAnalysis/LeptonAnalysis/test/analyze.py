#!/usr/bin/env python

from ROOT import *
from TopAnalysis.TTbarDilepton.PlotTool import *

gROOT.ProcessLine(".x rootlogon.C")

gROOT.cd()

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

for dsName in samples["mm"]:
    anaMM = NtupleAnalyzerLite("Run2012", "Summer12", "mm", dsName, "ntuple", "hEvents", "tree")
    anaMM.addCategory("m_OS", "q == 0 && isoDbeta1 < 0.15 && met > 30 && m > 12 && abs(m-91.2) > 15")
    anaMM.addCategory("m_SS", "q != 0 && isoDbeta1 > 0.15 && met < 30 && m > 12 && abs(m-91.2) > 15")
    anaMM.addHistogram("isoDbeta1", "hMuon1_relIso", "Muon isolation;Relative isolation;Entries", 100, 0, 2)
    anaMM.addHistogram("isoDbeta2", "hMuon2_relIso", "Muon isolation;Relative isolation;Entries", 100, 0, 2)
    anaMM.scanCutFlow()

for dsName in samples["ee"]:
    anaEE = NtupleAnalyzerLite("Run2012", "Summer12", "ee", dsName, "ntuple", "hEvents", "tree")
    anaEE.addCategory("e_OS", "q == 0 && isoRho1 < 0.15 && id1_mvaTrigV0 > 0.5 && met > 30 && m > 12 && abs(m-91.2) > 15")
    anaEE.addCategory("e_SS", "q != 0 && isoRho1 > 0.15 && id1_mvaTrigV0 > 0.5 && met < 30 && m > 12 && abs(m-91.2) > 15")
    anaEE.addHistogram("id1_mvaTrigV0", "hElectron1_mva", "Electron ID;MVA discriminator;Entries", 100, -1, 1)
    anaEE.addHistogram("id2_mvaTrigV0", "hElectron2_mva", "Electron ID;MVA discriminator;Entries", 100, -1, 1)
    anaEE.addHistogram("isoRho1", "hElectron1_relIso", "Electron isolation;Relative isolation;Entries", 100, 0, 2)
    anaEE.addHistogram("isoRho2", "hElectron2_relIso", "Electron isolation;Relative isolation;Entries", 100, 0, 2)

    anaEE.scanCutFlow()

for dsName in samples["me"]:
    anaME = NtupleAnalyzerLite("Run2012", "Summer12", "me", dsName, "ntuple", "hEvents", "tree")
    anaME.addCategory("e_OS", "q == 0 && isoDbeta1 < 0.15 && met > 20 && m > 12")
    anaME.addCategory("e_SS", "q != 0 && isoDbeta1 > 0.15 && met < 30 && m > 12")
    anaME.addCategory("m_OS", "q == 0 && isoRho2 < 0.15 && id2_mvaTrigV0 > 0.5 && met > 20")
    anaME.addCategory("m_SS", "q != 0 && isoRho2 > 0.15 && id2_mvaTrigV0 > 0.5 && met < 30")
    anaME.addHistogram("id2_mvaTrigV0", "hElectron1_mva", "Electron ID;MVA discriminator;Entries", 100, -1, 1)
    anaME.addHistogram("isoRho2", "hElectron1_relIso", "Electron isolation;Relative isolation;Entries", 100, 0, 2)
    anaME.addHistogram("isoDbeta1", "hMuon1_relIso", "Muon isolation;Relative isolation;Entries", 100, 0, 2)
    anaME.scanCutFlow()

