import FWCore.ParameterSet.Config as cms

noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

dhTriggerFilter = cms.EDFilter("HLTHighLevel",
    #TriggerResultsTag = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = cms.vstring('HLT_Mu9', 'HLT_DoubleMu3', 'HLT_Mu15', 'HLT_Ele15_SW_EleId_L1R'),           # provide list of HLT paths (or patterns) you want
    eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
    andOr = cms.bool(True),  # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(False)    # throw exception on unknown path names
)

# Vertex selection
oneGoodVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 20 && position.Rho <= 2"),
    filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

eventFilters = cms.Sequence(
    noscraping * dhTriggerFilter * oneGoodVertexFilter
)

