import FWCore.ParameterSet.Config as cms

jetHLTAnalyzer = cms.EDAnalyzer("JetHLTAnalyzer",
    interestedFilterName = cms.string("hlt1jet50U"),
    recoJet = cms.InputTag("ak5CaloJets"),
    cut = cms.PSet(
        recoMinEt = cms.double(20),
        l1MinEt = cms.double(30),
        workingPointEt = cms.double(56),
        maxL1DeltaR = cms.double(0.5),
        maxHLTDeltaR = cms.double(0.3)
    )
)
