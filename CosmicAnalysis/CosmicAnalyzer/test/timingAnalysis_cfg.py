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

import os
if 'RUNNUMBER' not in os.environ:
  runNumber = "121943"
else :
  runNumber = os.environ['RUNNUMBER']

if runNumber == "120020" :
  from source_beamSplash_120020_cfg import *
elif runNumber == "120026" :
  from source_beamSplash_120026_cfg import *
elif runNumber == "120042" :
  from source_beamSplash_120042_cfg import *
elif runNumber == "121943" :
  from source_beamSplash_121943_cfg import *
elif runNumber == "121964" :
  from source_beamSplash_121964_cfg import *
elif runNumber == "121993" :
  from source_beamSplash_121993_cfg import *
elif runNumber == "errStream" :
  from source_beamSplash_errStream_cfg import *
else :
  print "Invalid run number input"
  sys.exit()

process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.MessageLogger.cerr.FwkReport.reportEvery = 100

# File output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hist_beamSplash_%s.root' % runNumber)
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
