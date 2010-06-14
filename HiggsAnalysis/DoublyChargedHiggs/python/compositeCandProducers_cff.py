import FWCore.ParameterSet.Config as cms

dhCandProducerToEE = cms.EDProducer("DileptonProducer",
    lepton1 = cms.PSet(
        src = cms.InputTag("goodPatElectrons"),
        type = cms.string("electron"),
        charge = cms.int32(1),
        minPt = cms.double(5),
        maxEta = cms.double(2.5)
    ),
    lepton2 = cms.PSet(
        src = cms.InputTag("goodPatElectrons"),
        type = cms.string("electron"),
        charge = cms.int32(1),
        minPt = cms.double(5),
        maxEta = cms.double(2.5)
    ),
    chargeConjugation = cms.bool(True),
    minMass = cms.double(10),
    minPt = cms.double(5),
    maxEta = cms.double(3)
)


dhCandProducerToMM = cms.EDProducer("DileptonProducer",
    lepton1 = cms.PSet(
        src = cms.InputTag("goodPatMuons"),
        type = cms.string("muon"),
        charge = cms.int32(1),
        minPt = cms.double(5),
        maxEta = cms.double(2.5)
    ),
    lepton2 = cms.PSet(
        src = cms.InputTag("goodPatMuons"),
        type = cms.string("muon"),
        charge = cms.int32(1),
        minPt = cms.double(5),
        maxEta = cms.double(2.5)
    ),
    chargeConjugation = cms.bool(True),
    minMass = cms.double(10),
    minPt = cms.double(5),
    maxEta = cms.double(3)
)


dhCandProducerToEM = cms.EDProducer("DileptonProducer",
    lepton1 = cms.PSet(
        src = cms.InputTag("goodPatElectrons"),
        type = cms.string("electron"),
        charge = cms.int32(1),
        minPt = cms.double(5),
        maxEta = cms.double(2.5)
    ),
    lepton2 = cms.PSet(
        src = cms.InputTag("goodPatMuons"),
        type = cms.string("muon"),
        charge = cms.int32(1),
        minPt = cms.double(5),
        maxEta = cms.double(2.5)
    ),
    chargeConjugation = cms.bool(True),
    minMass = cms.double(10),
    minPt = cms.double(5),
    maxEta = cms.double(3)
)


