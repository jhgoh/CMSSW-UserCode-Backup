import FWCore.ParameterSet.Config as cms
import sys
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.TrackingTools.MuonTrackLoader_cff import *

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

# configure HLT
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('NOT (36 OR 37 OR 38 OR 39)')

process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

process.oneGoodVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 20 && position.Rho <= 2"),
    filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

process.L1HLT =cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    HLTPaths = cms.vstring('HLT_MinBiasBSC'),           # provide list of HLT paths (or patterns) you want
    eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
    andOr = cms.bool(True),  # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(False)    # throw exception on unknown path names
)

process.muonCounterFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('muons'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(99999)
)

process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagator_cfi")
process.load("HLTrigger.TPGAnalysis.muonHLTAnalyzer_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("muonHLTAnalysis.root")
)

process.p = cms.Path(
    process.hltLevel1GTSeed * 
    process.noscraping * process.oneGoodVertexFilter * process.muonCounterFilter + 
    process.muonHLTAnalyzer
)
