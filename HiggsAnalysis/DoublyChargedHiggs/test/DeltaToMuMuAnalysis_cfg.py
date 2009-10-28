import FWCore.ParameterSet.Config as cms

# Set variables from the os environment
globalTag = 'MC_31X_V8'

# Load Standard CMSSW process initial configurations
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

# Set datafiles
process.source.fileNames.append('file:PAT.root')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# PAT sequences
process.load("PhysicsTools/PatAlgos/patSequences_cff")
process.load("PhysicsTools/PatAlgos/patEventContent_cff")
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

# Muon selection
process.goodPatMuons = cms.EDProducer("PATMuonSelector",
    src = cms.InputTag("cleanLayer1Muons"),
    cut = cms.string('pt > 10 & abs(eta) < 2.5 & isGlobalMuon() & charge < 0')
)

# Candidate production
process.deltaToMuMu = cms.EDProducer("DileptonProducer",
    lepton1 = cms.PSet(
        src = cms.InputTag("goodPatMuons"),
        type = cms.string("muon"),
        charge = cms.int32(1)
    ),
    lepton2 = cms.PSet(
        src = cms.InputTag("goodPatMuons"),
        type = cms.string("muon"),
        charge = cms.int32(1)
    ),
    chargeConjugation = cms.bool(True)
)

# File output
#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('hHppMuMu.root')
#)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('HppMuMu.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *', 'keep *_*_*_Ana')
)
process.outpath = cms.EndPath(process.out)

# Module sequences and Paths
process.deltaCombineSeq = cms.Sequence(process.goodPatMuons * process.deltaToMuMu)
#process.deltaAnalysisSeq = cms.Sequence(process.fourMuonAnalyzer)

process.p = cms.Path(process.deltaCombineSeq)
#process.p = cms.Path(process.deltaCombineSeq*
#                     process.deltaAnalysisSeq)
