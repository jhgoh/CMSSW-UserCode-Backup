import FWCore.ParameterSet.Config as cms

# configure HLT
from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff import *
from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff import *
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import *

hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
#hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('NOT (36 OR 37 OR 38 OR 39)')

muonOrthogonalTriggers = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
    HLTPaths = cms.vstring("HLT_Jet30U", "HLT_Jet50U", "HLT_Jet70U"),
    eventSetupPathsKey = cms.string(''),
    andOr = cms.bool(True), # Set true to activate with "OR"
    throw = cms.bool(False)
)

noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(15),
    maxd0 = cms.double(2)
)

muonCounterFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('muons'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(99999)
)

muonMinBiasCommonFilters = cms.Sequence(
    hltLevel1GTSeed * noscraping * 
    primaryVertexFilter * muonCounterFilter
)

muonCommonFilters = cms.Sequence(
    hltLevel1GTSeed * muonOrthogonalTriggers * noscraping * 
    primaryVertexFilter * muonCounterFilter
)

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.TrackingTools.MuonTrackLoader_cff import *

from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagator_cfi import *

#TFileService = cms.Service("TFileService",
#    fileName = cms.string("TF.root")
#)

from HLTrigger.TPGAnalysis.muonHLTAnalyzer_cfi import *
