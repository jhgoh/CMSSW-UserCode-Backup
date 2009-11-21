import FWCore.ParameterSet.Config as cms

# Set variables from the os environment
globalTag = 'GR09_P_V4'

# Load Standard CMSSW process initial configurations
process = cms.Process('Ana')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = cms.string('%s::All' % globalTag)

# Read input sample list
import sys
sys.path.append('test')
sys.path.append('.')

from source_beamSplash_976_993_cfg import *
#from source_beamSplash_120026_cfg import *
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.MessageLogger.cerr.FwkReport.reportEvery = 3

# File output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hist_beamSplash_120026.root')
)

# Modules
process.muonTimingAnalyzer = cms.EDAnalyzer("MuonTimingAnalyzer",
    digisLabel = cms.InputTag('muonRPCDigis'),
    recHitsLabel = cms.InputTag('rpcRecHits')
)

process.cosmicAnalysisSeq = cms.Sequence(process.muonTimingAnalyzer)

process.p = cms.Path(process.cosmicAnalysisSeq)
