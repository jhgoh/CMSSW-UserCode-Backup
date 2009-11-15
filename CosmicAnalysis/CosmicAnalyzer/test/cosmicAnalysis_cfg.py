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
process.muonRPCAnalyzer = cms.EDAnalyzer("MuonRPCAnalyzer",
    digiLabel = cms.InputTag('muonRPCDigis'),
    histoDimensions = cms.PSet(
        minDate = cms.untracked.uint32(91108),
        maxDate = cms.untracked.uint32(101108),
        minTime = cms.untracked.uint32(0),
        maxTime = cms.untracked.uint32(235959),
        dDate = cms.untracked.uint32(0),
        dTime = cms.untracked.uint32(40000)
    )
)

# Sequences and Paths

process.cosmicAnalysisSeq = cms.Sequence(process.muonRPCAnalyzer)

process.p = cms.Path(process.cosmicAnalysisSeq)
