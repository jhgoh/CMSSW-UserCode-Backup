import FWCore.ParameterSet.Config as cms

jetHLTAnalyzer = cms.EDAnalyzer("JetHLTAnalyzer",
    interestedFilterName = cms.string("hlt1jet50U"),
    recoJet = cms.InputTag("ak5CaloJets"),
    cut = cms.PSet(
        l1MinEt = cms.double(30),
        workingPointEt = cms.double(50),
        maxL1DeltaR = cms.double(0.4),
        maxHLTDeltaR = cms.double(0.3)
    ),
    JetIDParams = cms.PSet(
        useRecHits      = cms.bool(True),
        hbheRecHitsColl = cms.InputTag("hbhereco"),
        hoRecHitsColl   = cms.InputTag("horeco"),
        hfRecHitsColl   = cms.InputTag("hfreco"),
        ebRecHitsColl   = cms.InputTag("ecalRecHit", "EcalRecHitsEB"),
        eeRecHitsColl   = cms.InputTag("ecalRecHit", "EcalRecHitsEE")
    )
)
