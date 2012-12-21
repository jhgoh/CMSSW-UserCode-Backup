import FWCore.ParameterSet.Config as cms
import sys, os

sample = os.environ["SAMPLE"]
section = int(os.environ["SECTION"])
nSection = int(os.environ["NSECTION"])

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("KoPFA.TopAnalyzer.topAnalysis_cff")
process.load("KoPFA.CommonTools.eventFilter_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
for line in open("/afs/cern.ch/user/j/jhgoh/public/sources/CMG/V5_10_0/%s.txt" % sample).readlines():
    line = line.strip()
    if 'root' not in line: continue
    if '#' == line[0]: continue
    process.source.fileNames.append(line)

if 'Run2012' in sample:
    isMC = False  
    from CMGTools.Common.Tools.applyJSON_cff import *
    json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
    applyJSON(process, json)
#process.load("KoPFA.CommonTools.Sources.CMG.V5_10_0.Run2012.cmgTuple_%s_cff" % sample)
else:
    isMC = True
    #process.load("KoPFA.CommonTools.Sources.CMG.V5_10_0.Summer12.cmgTuple_%s_cff" % sample)
    #process.source.fileNames = process.source.fileNames[:4000]
#process.maxEvents.input = 1000

nFile = len(process.source.fileNames)
begin = section*int(nFile/nSection)
end = min(nFile, (section+1)*int(nFile/nSection))
process.source.fileNames = process.source.fileNames[begin:end]

#from Configuration.AlCa.autoCond import autoCond
if isMC:
    process.GlobalTag.globaltag = "START53_V7A::All"
else:
	process.GlobalTag.globaltag = "GR_R_53_V14::All"
	process.source.lumisToProcess = cms.untracked.vstring()
	from CMGTools.Common.Tools.applyJSON_cff import applyJSON
	#json = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-201678_8TeV_PromptReco_Collisions12_JSON.txt"
	json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-203002_8TeV_PromptReco_Collisions12_JSON_v2.txt'
	applyJSON(process, json)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple/ntuple_%s_%04d.root" % (sample, section)),
)

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltHighLevel.throw = False
process.hltEE = process.hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",))
process.hltMM = process.hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*",))
process.hltME = process.hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",))

from KoPFA.CommonTools.PileUpWeight_cff import *
process.PUweight.PileUpRD   = PileUpRD2012
process.PUweight.PileUpRDup = PileUpRD2012UP
process.PUweight.PileUpRDdn = PileUpRD2012DN
process.PUweight.PileUpMC   = Summer12

process.load("TopAnalysis.LeptonAnalysis.patRelIsoLepton_cfi")
process.patElectronsWithRelIso.coneSize = cms.double(0.3)
process.patElectronsWithRelIso.cut = "isPF && passConversionVeto && gsfTrack.trackerExpectedHitsInner.numberOfHits <= 0 && abs(eta) < 2.5 && pt > 5"
process.patMuonsWithRelIso.coneSize = cms.double(0.3)
process.patMuonsWithRelIso.cut = "isPFMuon && (isGlobalMuon || isTrackerMuon) && abs(eta) < 2.5 && pt > 5"

process.ee = cms.EDAnalyzer("DoubleElectronAnalyzer",
    lepton1 = cms.InputTag("patElectronsWithRelIso"),
    lepton2 = cms.InputTag("patElectronsWithRelIso"),
    weight = cms.InputTag("PUweight", "weight"),
    met = cms.InputTag("cmgPFMET"),
    eventCounter = cms.InputTag("prePathCounter"),
    genPariclesLabel = cms.InputTag("genParticlesPruned"),
    idNames1 = cms.vstring(
        "mvaNonTrigV0", "mvaTrigV0", 
        "eidLoose",  "eidMedium", "eidSuperTight", "eidTight", "eidVeryLoose",
        "simpleEleId60relIso", "simpleEleId70relIso", "simpleEleId80relIso", 
        "simpleEleId85relIso", "simpleEleId90relIso", "simpleEleId95relIso",
    ),
    idNames2 = cms.vstring(
        "mvaNonTrigV0", "mvaTrigV0", 
        "eidLoose",  "eidMedium", "eidSuperTight", "eidTight", "eidVeryLoose",
        "simpleEleId60relIso", "simpleEleId70relIso", "simpleEleId80relIso", 
        "simpleEleId85relIso", "simpleEleId90relIso", "simpleEleId95relIso",
    ),
)

process.mm = cms.EDAnalyzer("DoubleMuonAnalyzer",
    lepton1 = cms.InputTag("patMuonsWithRelIso"),
    lepton2 = cms.InputTag("patMuonsWithRelIso"),
    weight = cms.InputTag("PUweight", "weight"),
    met = cms.InputTag("cmgPFMET"),
    eventCounter = cms.InputTag("prePathCounter"),
    genPariclesLabel = cms.InputTag("genParticlesPruned"),
    idNames1 = cms.vstring(),
    idNames2 = cms.vstring(),
)

process.me = cms.EDAnalyzer("MuEGAnalyzer",
    lepton1 = cms.InputTag("patMuonsWithRelIso"),
    lepton2 = cms.InputTag("patElectronsWithRelIso"),
    weight = cms.InputTag("PUweight", "weight"),
    met = cms.InputTag("cmgPFMET"),
    eventCounter = cms.InputTag("prePathCounter"),
    genPariclesLabel = cms.InputTag("genParticlesPruned"),
    idNames1 = cms.vstring(),
    idNames2 = cms.vstring(
        "mvaNonTrigV0", "mvaTrigV0",
        "eidLoose",  "eidMedium", "eidSuperTight", "eidTight", "eidVeryLoose",
        "simpleEleId60relIso", "simpleEleId70relIso", "simpleEleId80relIso",
        "simpleEleId85relIso", "simpleEleId90relIso", "simpleEleId95relIso",
    ),
)

if isMC:
    process.pEE = cms.Path(
        process.hltEE +
        process.PUweight+
        process.patElectronsWithRelIso +
        process.ee
    )
    process.pMM = cms.Path(
        process.hltMM +
        process.PUweight+
        process.patMuonsWithRelIso +
        process.mm
    )
    process.pME = cms.Path(
        process.hltME +
        process.PUweight+
        process.patElectronsWithRelIso + process.patMuonsWithRelIso +
        process.me
    )
else:
    process.pEE = cms.Path(
        process.hltEE +
        process.patElectronsWithRelIso +
        process.ee
    )
    process.pMM = cms.Path(
        process.hltMM +
        process.patMuonsWithRelIso +
        process.mm
    )
    process.pME = cms.Path(
        process.hltME +
        process.patElectronsWithRelIso + process.patMuonsWithRelIso +
        process.me
    )

#print process.source.fileNames[0]
#print process.source.fileNames[-1]
