import FWCore.ParameterSet.Config as cms

patMuonsWithRelIso = cms.EDFilter("PatRelIsoMuonSelector",
    rho = cms.InputTag("kt6PFJets", "rho"),
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(""),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
    coneSize = cms.double(0.3),
)

patElectronsWithRelIso = cms.EDFilter("PatRelIsoElectronSelector",
    rho = cms.InputTag("kt6PFJets", "rho"),
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string(""),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
    coneSize = cms.double(0.3),
)
