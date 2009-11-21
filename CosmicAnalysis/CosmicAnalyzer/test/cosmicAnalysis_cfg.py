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
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.MessageLogger.cerr.FwkReport.reportEvery = 100

# CondDB
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.CondDBCommon.connect = 'frontier://cmsfrontier.cern.ch:8000/FrontierProd/CMS_COND_31X_RPC'
"""
cmscond_list_iov -c frontier://cmsfrontier.cern.ch:8000/FrontierProd/CMS_COND_31X_RPC
"""
#process.CondDBCommon.connect = 'sqlite_file:/afs/cern.ch/user/d/dpagano/public/dati.db'
#process.CondDBCommon.DBParameter.authenticationPath='./'

process.rpcMon = cms.ESSource("PoolDBESSource",
  process.CondDBCommon,
  timetype = cms.string('timestamp'),
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('RPCObImonRcd'),
      tag = cms.string('Imon_STD')
    ),
    cms.PSet(
      record = cms.string('RPCObVmonRcd'),
      tag = cms.string('Vmon_STD')
    ),
    cms.PSet(
      record = cms.string('RPCObTempRcd'),
      tag = cms.string('Temp_STD')
    ),
    cms.PSet(
      record = cms.string('RPCObPVSSmapRcd'),
      tag = cms.string('RPCPVSSmap_STD')
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
        minUTime = cms.untracked.uint64(long(time.mktime( (2009, 11,  3, 0, 0, 0, 0, 0, 0) ))),
        maxUTime = cms.untracked.uint64(long(time.mktime( (2009, 11, 17, 0, 0, 0, 0, 0, 0) ))),
        dTime = cms.untracked.uint64(4*3600)
    )
)

process.cosmicAnalysisSeq = cms.Sequence(process.muonRPCAnalyzer)

process.p = cms.Path(process.cosmicAnalysisSeq)
