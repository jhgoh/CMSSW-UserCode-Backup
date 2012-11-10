import FWCore.ParameterSet.Config as cms
import sys, os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
#process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.GlobalTag.globaltag = "START52_V12::All"

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
from TopAnalysis.TTbarDilepton.dataset_cff import *
dataset, section, nFiles = parseJobSectionOption()
files = loadDataset(dataset)
begin, end = calculateRange(files, section, nFiles)
process.source.fileNames = files[begin:end]

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple/unmerged/ntuple_%s_%03d.root" % (dataset, section)),
)

process.genParticleCount = cms.EDFilter("GenParticleCountFilter",
    src = cms.InputTag("genParticlesPruned"),
    cut = cms.string("status == 3"),
    ids = cms.vint32(13, -13, 11, -11),
    minNumber = cms.uint32(1), 
    maxNumber = cms.uint32(1),
)

process.genParticleTauVeto =cms.EDFilter("GenParticleCountFilter",
    src = cms.InputTag("genParticlesPruned"),
    cut = cms.string("status == 3"),
    ids = cms.vint32(15, -15),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(0),
)

process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
    src = cms.InputTag("genParticlesPruned"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False)
)

process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25),
)

process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter", 
    src = cms.InputTag('offlinePrimaryVertices'),
    filterParams =  cms.PSet(
        minNdof = cms.double(4.),
        maxZ    = cms.double(24.), 
        maxRho  = cms.double(2.)
    ),
    filter = cms.bool(True),
)

process.load("Configuration.StandardSequences.Generator_cff")
process.genParticlesForJetsNoNu.src = "genParticlesPruned"

process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
from TopAnalysis.TTbarDilepton.trigger_cff import *

process.load("TopAnalysis.GeneratorTools.genJetAssociation_cff")
process.load("TopAnalysis.GeneratorTools.lumiWeight_cff")

process.event = cms.EDAnalyzer("EventTupleProducer",
    eventCounters = cms.vstring(
        'prePathCounter',
        'postPathCounter',
    ),

    doMCMatch = cms.bool(not isRealData(dataset)),
    gen = cms.InputTag("genParticlesPruned"),
    genJetToPartonsMap = cms.InputTag("genJetToPartons"),
    recoToGenJetMap = cms.InputTag("recoToGenJet"),

    weight = cms.string("lumiWeight"),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    met = cms.InputTag("cmgPFMET"),

    electron = cms.PSet(
        src = cms.InputTag("cmgElectronSel"),
        cut = cms.string(
            "pt>20 && abs(eta) < 2.5"
            " && sourcePtr.get.gsfTrack.isNonnull && sourcePtr.get.gsfTrack.trackerExpectedHitsInner.numberOfLostHits<2"
            " && sourcePtr.get.gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1 "
            #" && !(1.4442 < abs(sourcePtr.get.superCluster.eta) && abs(sourcePtr.get.superCluster.eta) < 1.5660)"
            " && relIso(0.5, 0, 0.3) < 0.15 && sourcePtr.get.dB < 0.04"
            ' && passConversionVeto && sourcePtr.get.electronID("mvaTrigV0") >= 0'
        ),
    ),
    muon = cms.PSet(
        src = cms.InputTag("cmgMuonSel"),
        cut = cms.string(
            "abs(eta) < 2.4 && pt > 20" 
            #" && sourcePtr.get.dB < 0.2"
            " && sourcePtr.get.isPFMuon && (sourcePtr.get.isGlobalMuon || sourcePtr.get.isTrackerMuon)"
            #" && sourcePtr.get.normChi2 < 10"
            #" && sourcePtr.get.track.hitPattern.trackerLayersWithMeasurement > 5"
            #" && sourcePtr.get.globalTrack.hitPattern.numberOfValidMuonHits > 0"
            #" && sourcePtr.get.innerTrack.hitPattern.numberOfValidPixelHits > 0"
            #" && sourcePtr.get.numberOfMatchedStations > 1"
            " && relIso(0.5, 0, 0.3) < 0.2"
        ),
    ),
    jet = cms.PSet(
        src = cms.InputTag("cmgPFJetSelCHS"),
        cut = cms.string(
            " abs(eta) < 2.5 && pt > 30 && nConstituents > 1"
            " && component(5).fraction < 0.99 && component(4).fraction < 0.99"
            " && (abs(eta) >= 2.4 || component(2).fraction < 0.99 )"
            " && (abs(eta) >= 2.4 || component(1).fraction > 0 )"
            " && (abs(eta) >= 2.4 || component(1).number > 0 ) "
        ),
        bTagType = cms.string("combinedSecondaryVertexBJetTags"),
    ),
)

if isRealData(dataset):
    process.p = cms.Path(
        process.goodOfflinePrimaryVertices
    #  + process.hltHighLevel
      * process.event
    )
else:
    process.p = cms.Path(
        process.goodOfflinePrimaryVertices
    #  + process.genParticleCount + process.genParticleTauVeto
    #  + process.hltHighLevel
      + process.genParticlesForJetsNoNu * process.ak5GenJetsNoNu
      + process.recoToGenJet + process.genJetToPartons
      + process.lumiWeight
      * process.event
    #   process.printDecay
    )

