import FWCore.ParameterSet.Config as cms

globalTag = 'IDEAL_V12'

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = cms.string('%s::All' % globalTag)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# PAT sequences
process.load("PhysicsTools/PatAlgos/patSequences_cff")
process.load("PhysicsTools/PatAlgos/patEventContent_cff")
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

# Muon selection
process.goodPatMuons = cms.EDProducer("PATMuonSelector",
    src = cms.InputTag("cleanLayer1Muons"),
    cut = cms.string("pt > 10 & abs(eta) < 2.5 & isGlobalMuon()")
)

process.goodPatElectrons = cms.EDProducer("PATElectronSelector",
    src = cms.InputTag("cleanLayer1Electrons"),
    cut = cms.string("pt > 10 & abs(eta) < 2.5")
)

process.goodPatMuonCountFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodPatMuons"),
    minNumber = cms.uint32(2)
)

process.goodPatElectronCountFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodPatElectrons"),
    minNumber = cms.uint32(2)
)

# Candidate production
process.negDeltaToEMu = cms.EDProducer("DileptonProducer",
    lepton1 = cms.PSet(
        src = cms.InputTag("goodPatElectrons"),
        type = cms.string("electron"),
        charge = cms.int32(-1)
    ),
    lepton2 = cms.PSet(
        src = cms.InputTag("goodPatMuons"),
        type = cms.string("muon"),
        charge = cms.int32(-1)
    ),
    chargeConjugation = cms.bool(False)
)

process.posDeltaToEMu = cms.EDProducer("DileptonProducer",
    lepton1 = cms.PSet(
        src = cms.InputTag("goodPatElectrons"),
        type = cms.string("electron"),
        charge = cms.int32(1)
    ),
    lepton2 = cms.PSet(
        src = cms.InputTag("goodPatMuons"),
        type = cms.string("muon"),
        charge = cms.int32(1)
    ),
    chargeConjugation = cms.bool(False)
)

process.posDeltaCountFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("posDeltaToEMu"),
    minNumber = cms.uint32(1)
)

process.negDeltaCountFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("negDeltaToEMu"),
    minNumber = cms.uint32(1)
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('HppEMu.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *', 'keep *_*_*_Ana')
)
process.outpath = cms.EndPath(process.out)

# Module sequences and Paths
process.muonFilterSeq = cms.Sequence(process.goodPatMuons*process.goodPatMuonCountFilter)
process.electronFilterSeq = cms.Sequence(process.goodPatElectrons*process.goodPatElectronCountFilter)
process.candCombineSeq = cms.Sequence(process.posDeltaToEMu+process.negDeltaToEMu)
process.candFilterSeq = cms.Sequence(process.posDeltaCountFilter + process.negDeltaCountFilter)
#process.deltaAnalysisSeq = cms.Sequence(process.fourMuonAnalyzer)

process.p = cms.Path(process.muonFilterSeq*process.electronFilterSeq*process.candCombineSeq*process.candFilterSeq)
