import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.GlobalTag.globaltag = cms.string('GR10_P_V7::All')

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(),
    lumisToProcess = cms.untracked.VLuminosityBlockRange()
)

# Set input files and LumiSections
#sys.path.append("samples")
#from muonJetMETSkim_cff import *
#process.source.fileNames = fileNames
#process.source.lumisToProcess = lumisToProcess

process.load("HLTrigger.TPGAnalysis.muonAnalysis_cff")
process.load("HLTrigger.TPGAnalysis.dqmFourVector_cff")

process.p = cms.Path(
    process.muonCommonFilters * 
    process.muonHLTAnalyzer + 
    process.hltResults * process.hltFourVectorClient
)

process.out = cms.EndPath(process.dqmSaver)

