import FWCore.ParameterSet.Config as cms

recoToGenJet = cms.EDProducer("RecoToGenJetAssociator",
    genJets = cms.InputTag("genJetSel"),
#    genJets = cms.InputTag("ak5GenJetsNoNu"),
    recoJets = cms.InputTag("cmgPFJetSelCHS"),
    cuts = cms.PSet(
        maxDR = cms.double(0.5),
    ),
)

genJetToPartons = cms.EDProducer("GenJetPartonAssociator",
    genJets = cms.InputTag("genJetSel"),
    genParticles = cms.InputTag("genParticlesPruned"),
#    genJets = cms.InputTag("ak5GenJetsNoNu"),
#    genParticles = cms.InputTag("genParticles"),
    matchAlgo = cms.string("deltaR"),
    cuts = cms.PSet(
        maxDR = cms.double(0.5),
        #maxDPt = cms.double(1e9),
    ),
    pdgIdsToMatch = cms.vuint32(1,2,3,4,5,6,21,), # Jets from Quarks and gluons
#    matchAlgo = cms.string("constituent"),
#    cuts = cms.PSet(
#        minNConstituent = cms.int32(1),
#        minFracConstituent = cms.double(0.1),
#    ),
#    pdgIdsToMatch = cms.vuint32(), # Empty list to "match to any partons"
)

