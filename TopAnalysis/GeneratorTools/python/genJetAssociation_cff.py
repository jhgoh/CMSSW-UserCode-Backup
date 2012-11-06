import FWCore.ParameterSet.Config as cms

recoToGenJet = cms.EDProducer("RecoToGenJetAssociator",
    genJets = cms.InputTag("ak5GenJetsNoNu"),
    recoJets = cms.InputTag("cmgPFJetSelCHS"),
    cuts = cms.PSet(
        maxDR = cms.double(0.5),
    ),
)

genJetToPartons = cms.EDProducer("GenJetPartonAssociator",
    genJets = cms.InputTag("ak5GenJetsNoNu"),
    genParticles = cms.InputTag("genParticlesPruned"),
#    genParticles = cms.InputTag("genParticles"),
    cuts = cms.PSet(
        minNConstituent = cms.int32(1),
        minFracConstituent = cms.double(0.5),
    ),
)

