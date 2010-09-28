import FWCore.ParameterSet.Config as cms

# configure HLT
from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff import *
from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff import *
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import *

hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
#hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('NOT (36 OR 37 OR 38 OR 39)')

jetMetTauOrthogonalTriggers = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
    HLTPaths = cms.vstring("HLT_Mu3", "HLT_Mu5", "HLT_Mu7", "HLT_Mu9"),
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

jetCounterFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("ak5CaloJets"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(99999)
)

jetMetTauCommonFilters = cms.Sequence(
    hltLevel1GTSeed * jetMetTauOrthogonalTriggers * noscraping *
    primaryVertexFilter * jetCounterFilter
)

jetMetTauMinBiasCommonFilters = cms.Sequence(
    hltLevel1GTSeed * noscraping *
    primaryVertexFilter * jetCounterFilter
)

from HLTrigger.TPGAnalysis.jetHLTAnalyzer_cfi import *

