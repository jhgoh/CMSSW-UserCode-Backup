import FWCore.ParameterSet.Config as cms

# Set variables from the os environment
globalTag = 'CRAFT09_R_V4'

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

from source_cfg import *
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

process.MessageLogger.cerr.FwkReport.reportEvery = 100

# CondDB
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.CondDBCommon.connect = 'sqlite_file:/afs/cern.ch/user/d/dpagano/public/dati.db'
#process.CondDBCommon.DBParameter.authenticationPath='./'

process.rpcMon = cms.ESSource("PoolDBESSource",
  process.CondDBCommon,
  timetype = cms.string('timestamp'),
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('RPCObImonRcd'),
      tag = cms.string('Imon_v3')
    ),
    cms.PSet(
      record = cms.string('RPCObVmonRcd'),
      tag = cms.string('Vmon_v3')
    ),
    cms.PSet(
      record = cms.string('RPCObTempRcd'),
      tag = cms.string('Temp_v3')
    ),
    cms.PSet(
      record = cms.string('RPCObPVSSmapRcd'),
      tag = cms.string('PVSS_v3')
    )
  )
)

# File output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hist.root')
)

# Modules
import time
process.muonRPCAnalyzer = cms.EDAnalyzer("MuonRPCAnalyzer",
    digiLabel = cms.InputTag('muonRPCDigis'),
    histoDimensions = cms.PSet(
        minUTime = cms.untracked.uint32(4294915000),
        maxUTime = cms.untracked.uint32(4294915000+3600*4),
#        minUTime = cms.untracked.uint32( int(time.mktime( (2008, 11,  9, 0, 0, 0, 0, 0, 0) )) ),
#        maxUTime = cms.untracked.uint32( int(time.mktime( (2008, 11, 10, 0, 0, 0, 0, 0, 0) )) ),
        dTime = cms.untracked.uint32(30*60)
    )
)

# Sequences and Paths

process.cosmicAnalysisSeq = cms.Sequence(process.muonRPCAnalyzer)

process.p = cms.Path(process.cosmicAnalysisSeq)
