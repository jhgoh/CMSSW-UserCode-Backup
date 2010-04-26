import FWCore.ParameterSet.Config as cms

dhCandProducerToEE = cms.EDProducer("DileptonProducer",
    lepton1 = cms.PSet(
        src = cms.InputTag("goodPatElectrons"),
        type = cms.string("electron"),
        charge = cms.int32(1)
    ),
    lepton2 = cms.PSet(
        src = cms.InputTag("goodPatElectrons"),
        type = cms.string("electron"),
        charge = cms.int32(1)
    ),
    chargeConjugation = cms.bool(True)
)


dhCandProducerToMM = cms.EDProducer("DileptonProducer",
    lepton1 = cms.PSet(
        src = cms.InputTag("goodPatMuons"),
        type = cms.string("muon"),
        charge = cms.int32(1)
    ),
    lepton2 = cms.PSet(
        src = cms.InputTag("goodPatMuons"),
        type = cms.string("muon"),
        charge = cms.int32(1)
    ),
    chargeConjugation = cms.bool(True)
)


dhCandProducerToEM = cms.EDProducer("DileptonProducer",
    lepton1 = cms.PSet(
        src = cms.InputTag("goodPatElectrons"),
        type = cms.string("electron"),
        charge = cms.int32(1)
    ),
    lepton2 = cms.PSet(
        src = cms.InputTag("goodPatMuons"),
        type = cms.string("muon"),
        charge = cms.int32(1)
    ),
    chargeConjugation = cms.bool(True)
)


