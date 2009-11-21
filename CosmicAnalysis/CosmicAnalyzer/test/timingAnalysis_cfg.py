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

#from source_beamSplash_121943_cfg import *
from source_beamSplash_121964_cfg import *
#from source_beamSplash_121993_cfg import *
#from source_beamSplash_976_993_cfg import *
#from source_beamSplash_120026_cfg import *
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.MessageLogger.cerr.FwkReport.reportEvery = 100

# File output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hist_beamSplash_121964.root')
)

# Analysis modules
process.muonTimingAnalyzer = cms.EDAnalyzer("MuonTimingAnalyzer",
    digiLabel = cms.InputTag('muonRPCDigis'),
    recHitLabel = cms.InputTag('rpcRecHits')
)
process.cosmicAnalysisSeq = cms.Sequence(process.muonTimingAnalyzer)

# Unpack RPC Digis
import EventFilter.RPCRawToDigi.rpcUnpacker_cfi
process.muonRPCDigis = EventFilter.RPCRawToDigi.rpcUnpacker_cfi.rpcunpacker.clone()
process.rpcUnpackerSeq = cms.Sequence(process.muonRPCDigis)

# Define run path
process.p = cms.Path(process.rpcUnpackerSeq*process.cosmicAnalysisSeq)
