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
if 'Run2012' in sample:
    isMC = False  
    process.load("KoPFA.CommonTools.Sources.CMG.V5_10_0.Run2012.cmgTuple_%s_cff" % sample)
else:
    isMC = True
    process.load("KoPFA.CommonTools.Sources.CMG.V5_10_0.Summer12.cmgTuple_%s_cff" % sample)
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
    json = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-201678_8TeV_PromptReco_Collisions12_JSON.txt"
    applyJSON(process, json)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple/ntuple_%s_%04d.root" % (sample, section)),
)

from KoPFA.CommonTools.PileUpWeight_cff import *
process.PUweight.PileUpRD   = PileUpRD2012
process.PUweight.PileUpRDup = PileUpRD2012UP
process.PUweight.PileUpRDdn = PileUpRD2012DN
process.PUweight.PileUpMC   = Summer12

process.load("TopAnalysis.LeptonAnalysis.patRelIsoLepton_cfi")
process.patElectronsWithRelIso04 = process.patElectronsWithRelIso.clone(coneSize = cms.double(0.4))
process.patMuonsWithRelIso04 = process.patMuonsWithRelIso.clone(coneSize = cms.double(0.4))

process.genElectron = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticlesPruned"),
    cut = cms.string("abs(pdgId) == 11 && status == 3"),
)

process.genMuon = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticlesPruned"),
    cut = cms.string("abs(pdgId) == 13 && status == 3"),
)

process.genElectronFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("genElectron"),
    minNumber = cms.uint32(2),
)

process.genMuonFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("genMuon"),
    minNumber = cms.uint32(2),
)

process.genMuEGMuonFilter = process.genMuonFilter.clone(minNumber = cms.uint32(1))
process.genMuEGElectronFilter = process.genElectronFilter.clone(minNumber = cms.uint32(1))
process.genMuEGFilter = cms.Sequence(process.genMuEGMuonFilter+process.genMuEGElectronFilter)

process.ee = cms.EDAnalyzer("DoubleElectronAnalyzer",
    lepton1 = cms.InputTag("patElectronsWithRelIso"),
    lepton2 = cms.InputTag("patElectronsWithRelIso"),
    weight = cms.InputTag("PUweight", "weight"),
    met = cms.InputTag("patMETs"),
    eventCounter = cms.InputTag("prePathCounter"),
    isoDR = cms.double(0.3),
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
    met = cms.InputTag("patMETs"),
    eventCounter = cms.InputTag("prePathCounter"),
    isoDR = cms.double(0.3),
    idNames1 = cms.vstring(),
    idNames2 = cms.vstring(),
)

process.em = cms.EDAnalyzer("MuEGAnalyzer",
    lepton1 = cms.InputTag("patMuonsWithRelIso"),
    lepton2 = cms.InputTag("patElectronsWithRelIso"),
    weight = cms.InputTag("PUweight", "weight"),
    met = cms.InputTag("patMETs"),
    eventCounter = cms.InputTag("prePathCounter"),
    isoDR = cms.double(0.3),
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
        process.genElectron* process.genElectronFilter*
        process.PUweight+
        process.ee
    )
    process.eeOthers = process.ee.clone()
    process.pEEOthers = cms.Path(
        process.genElectron* ~process.genElectronFilter*
        process.PUweight+
        process.eeOthers
    )
    process.pMM = cms.Path(
        process.genMuon* process.genMuonFilter*
        process.PUweight+
        process.mm
    )
    process.mmOthers = process.mm.clone()
    process.pMMOthers = cms.Path(
        process.genMuon* ~process.genMuonFilter*
        process.PUweight+
        process.mmOthers
    )
    process.pEM = cms.Path(
        process.genMuon* process.genMuEGMuonFilter*
        process.genElectron* process.genMuEGElectronFilter*
        process.PUweight+
        process.em
    )
    process.emOthers = process.em.clone()
    process.pEMOthers = cms.Path(
        process.genMuon* ~process.genMuEGMuonFilter*
        process.genElectron* ~process.genMuEGElectronFilter*
        process.PUweight+
        process.emOthers
    )
else:
    process.pEE = cms.Path(
        process.ee
    )
    process.pMM = cms.Path(
        process.mm
    )
    process.pEM = cms.Path(
        process.em
    )

#print process.source.fileNames[0]
#print process.source.fileNames[-1]
