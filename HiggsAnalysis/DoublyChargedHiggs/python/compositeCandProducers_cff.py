import FWCore.ParameterSet.Config as cms

dhCandToEE = cms.EDProducer("DileptonProducer",
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


dhCandToMM = cms.EDProducer("DileptonProducer",
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
    minPt = cms.double(10),
    maxEta = cms.double(3)
)


dhCandToEM = cms.EDProducer("DileptonProducer",
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
    minPt = cms.double(10),
    maxEta = cms.double(3)
)

atLeastOneDHCandFilter = cms.EDFilter("MultipleCandCounterFilter",
    cands = cms.VInputTag(
        cms.InputTag("dhCandToEE"),
        cms.InputTag("dhCandToMM"),
        cms.InputTag("dhCandToEM")
    ),
    ptThresholds = cms.vdouble(10)
)
