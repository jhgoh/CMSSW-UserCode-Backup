import FWCore.ParameterSet.Config as cms

# Lepton producer 
goodPatMuons = cms.EDProducer("PATMuonSelector",
    src = cms.InputTag("cleanPatMuons"),
    cut = cms.string('pt > 5 & abs(eta) < 2.5 & isGood("GlobalMuonPromptTight") & isGood("TrackerMuonArbitrated")')
)

goodPatElectrons = cms.EDProducer("PATElectronSelector",
    src = cms.InputTag("cleanPatElectrons"),
    cut = cms.string('pt > 5 & abs(eta) < 2.5 & electronID("eidRobustLoose") > 0.5')
)

leptonSelectionSeq = cms.Sequence(goodPatElectrons+goodPatMuons)

