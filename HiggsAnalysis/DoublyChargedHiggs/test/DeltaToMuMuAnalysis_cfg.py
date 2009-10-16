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

# Candidate production
process.posDeltaToMuMu = cms.EDProducer("DileptonProducer",
    lepton1 = cms.PSet(
        src = cms.InputTag("cleanLayer1Muons"),
        charge = cms.int32(1),
        type = cms.string("muon")
    ),
    lepton2 = cms.PSet(
        src = cms.InputTag("cleanLayer1Muons"),
        charge = cms.int32(1),
        type = cms.string("muon")
    ),
    chargeConjugation = cms.bool(False)
)

process.negDeltaToMuMu = cms.EDProducer("DileptonProducer",
    lepton1 = cms.PSet(
        src = cms.InputTag("cleanLayer1Muons"),
        charge = cms.int32(1),
        type = cms.string("muon")
    ),
    lepton2 = cms.PSet(
        src = cms.InputTag("cleanLayer1Muons"),
        charge = cms.int32(1),
        type = cms.string("muon")
    ),
    chargeConjugation = cms.bool(False)
)

# User analysis block

# File output
#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('hHppMuMu.root')
#)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('HppMuMu.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output
    outputCommands = cms.untracked.vstring('drop *', 'keep *_*_*_Ana') # you need a '*' to unpack the list of commands 'patEventContent'
)
process.outpath = cms.EndPath(process.out)

# Module sequences and Paths
process.deltaCombineSeq = cms.Sequence(process.posDeltaToMuMu+process.negDeltaToMuMu)
#process.deltaAnalysisSeq = cms.Sequence(process.fourMuonAnalyzer)

process.p = cms.Path(process.deltaCombineSeq)
#process.p = cms.Path(process.deltaCombineSeq*
#                     process.deltaAnalysisSeq)
