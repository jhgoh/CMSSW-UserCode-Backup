import FWCore.ParameterSet.Config as cms

# Set variables from the os environment
globalTag = 'IDEAL_V12'

# Load Standard CMSSW process initial configurations
process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = cms.string('%s::All' % globalTag)

# Set datafiles
#process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
process.source = cms.Source("PoolSource")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# User analysis block
## Analyzer to have a overall look 
process.hToEMuPromptMu = cms.EDAnalyzer("HiggsToEMuAnalyzer",
    muon1 = cms.InputTag("cleanLayer1Muons"),
    e1 = cms.InputTag("cleanLayer1Electrons"),
    muonSelection = cms.string("GlobalMuonPromptTight")
)

process.hToEMuAllMu = cms.EDAnalyzer("HiggsToEMuAnalyzer",
    muon1 = cms.InputTag("cleanLayer1Muons"),
    e1 = cms.InputTag("cleanLayer1Electrons"),
    muonSelection = cms.string("AllGlobalMuons")
)

# File output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hist_HiggsToEMu.root')
)

# Module sequences and Paths
process.emuAnaSeq = cms.Sequence(process.hToEMuPromptMu+process.hToEMuAllMu)

process.p = cms.Path(process.emuAnaSeq)

