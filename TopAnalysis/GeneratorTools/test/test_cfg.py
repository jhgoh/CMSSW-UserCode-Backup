import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames.append("/store/cmst3/user/cmgtools/CMG/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_10_0/patTuple_0.root")

process.recoToGenJet = cms.EDProducer("RecoToGenJetAssociator",
    genJets = cms.InputTag("ak5GenJetsNoNu"),
    recoJets = cms.InputTag("cmgPFJetSelCHS"),
    cuts = cms.PSet(
        maxDR = cms.double(0.5),
    ),
)

process.genJetToPartons = cms.EDProducer("GenJetPartonAssociator",
    genJets = cms.InputTag("ak5GenJetsNoNu"),
#    genParticles = cms.InputTag("genParticlesPruned"),
    genParticles = cms.InputTag("genParticles"),
    cuts = cms.PSet(
        minNConstituent = cms.int32(1),
        minFracConstituent = cms.double(0.5),
    ),
)

process.out = cms.OutputModule("PoolOutputModule", 
    fileName = cms.untracked.string("out.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_*_*_Ana",
    )
)

process.outPath = cms.EndPath(process.out)

process.p = cms.Path(
    process.recoToGenJet
  + process.genJetToPartons
)

