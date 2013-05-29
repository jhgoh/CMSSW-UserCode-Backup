import FWCore.ParameterSet.Config as cms
import sys, os

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
    fileNames = cms.untracked.vstring('file:/home/jhgoh/data/TPGAnalysis/6E9ECE95-3992-DF11-813C-001D09F27067.root',),
    lumisToProcess = cms.untracked.VLuminosityBlockRange()
)

# Set input files and LumiSections
sys.path.append("samples")
from muonJetMETTau_ReReco_v2_cff import *
#from muMonitor_cff import *
process.source.fileNames = fileNames
process.source.lumisToProcess = lumisToProcess
#jobSection = int(os.environ["JOBSECTION"])
#process.source.lumisToProcess = selectLumi(jobSection, lumisToProcess)
process.GlobalTag.globaltag = globalTag

# Analyzer modules
process.load("HLTrigger.TPGAnalysis.jetMetTauAnalysis_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("result.root")
)

process.p = cms.Path(
    #process.jetMetTauCommonFilters *
    process.jetMetTauMinBiasCommonFilters *
    process.jetHLTAnalyzer
)

