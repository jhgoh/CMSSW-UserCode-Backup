import FWCore.ParameterSet.Config as cms

# List of 6 separate analysis
dhGenEventFilterToEEEE = cms.EDFilter("DHGenEventFilter",
    genLabel = cms.InputTag("genParticles"),
    decay1 = cms.string("ee"),
    decay2 = cms.string("ee")
)

dhGenEventFilterToMMMM = cms.EDFilter("DHGenEventFilter",
    genLabel = cms.InputTag("genParticles"),
    decay1 = cms.string("mm"),
    decay2 = cms.string("mm")
)

dhGenEventFilterToEMEM = cms.EDFilter("DHGenEventFilter",
    genLabel = cms.InputTag("genParticles"),
    decay1 = cms.string("em"),
    decay2 = cms.string("em")
)

dhGenEventFilterToEEMM = cms.EDFilter("DHGenEventFilter",
    genLabel = cms.InputTag("genParticles"),
    decay1 = cms.string("ee"),
    decay2 = cms.string("mm")
)

dhGenEventFilterToEEEM = cms.EDFilter("DHGenEventFilter",
    genLabel = cms.InputTag("genParticles"),
    decay1 = cms.string("ee"),
    decay2 = cms.string("em")
)

dhGenEventFilterToEMMM = cms.EDFilter("DHGenEventFilter",
    genLabel = cms.InputTag("genParticles"),
    decay1 = cms.string("em"),
    decay2 = cms.string("mm")
)

