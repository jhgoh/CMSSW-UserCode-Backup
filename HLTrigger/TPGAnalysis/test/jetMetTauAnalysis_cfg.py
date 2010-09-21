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
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(),
    lumisToProcess = cms.untracked.VLuminosityBlockRange()
)

# Set input files and LumiSections
#import sys, os
#sys.path.append("samples")
#from TPGSkim_goldenSample_139407_cff import *
#process.source.fileNames = fileNames
#process.source.lumisToProcess = lumisToProcess
#process.GlobalTag.globaltag = globalTag

# Analyzer modules
process.load("HLTrigger.TPGAnalysis.jetMetTauAnalysis_cff")
process.load("HLTrigger.TPGAnalysis.dqmFourVector_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("jetMetTauAnalysis.root")
)

process.p = cms.Path(
    #process.jetMetTauCommonFilters *
    process.jetMetTauMinBiasCommonFilters *
    process.jetHLTAnalyzer #+
    #process.hltResults * process.hltFourVectorClient
)

#process.outPath = cms.EndPath(process.dqmSaver)
