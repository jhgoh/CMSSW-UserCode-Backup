import FWCore.ParameterSet.Config as cms

eventTuple = cms.EDAnalyzer("EventTupleProducer",
    eventCounters = cms.vstring(
        'prePathCounter',
        'postPathCounter',
    ),

    doMCMatch = cms.bool(False),
    gen = cms.InputTag("genParticlesPruned"),
    genJetToPartonsMap = cms.InputTag("genJetToPartons"),
    recoToGenJetMap = cms.InputTag("recoToGenJet"),

    weight = cms.string("lumiWeight"),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    met = cms.InputTag("cmgPFMET"),

    electron = cms.PSet(
        src = cms.InputTag("cmgElectronSel"),
        cut = cms.string(
            "pt>20 && abs(eta) < 2.5"
            " && sourcePtr.get.gsfTrack.isNonnull && sourcePtr.get.gsfTrack.trackerExpectedHitsInner.numberOfLostHits<2"
            " && sourcePtr.get.gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1 "
            " && relIso(0.5, 0, 0.3) < 0.15 && sourcePtr.get.dB < 0.04"
            ' && passConversionVeto && sourcePtr.get.electronID("mvaTrigV0") >= 0'
        ),
        dz = cms.double(999),
        minNumber = cms.uint32(0),
        maxNumber = cms.uint32(999),
    ),
    muon = cms.PSet(
        src = cms.InputTag("cmgMuonSel"),
        cut = cms.string(
            "abs(eta) < 2.5 && pt > 20" 
            " && sourcePtr.get.isPFMuon && (sourcePtr.get.isGlobalMuon || sourcePtr.get.isTrackerMuon)"
            " && relIso(0.5, 0, 0.3) < 0.15"
        ),
        dz = cms.double(999),
        minNumber = cms.uint32(0),
        maxNumber = cms.uint32(999),
    ),
    jet = cms.PSet(
        src = cms.InputTag("cmgPFJetSelCHS"),
        cut = cms.string(
            " abs(eta) < 2.5 && pt > 30 && nConstituents > 1"
            " && component(5).fraction < 0.99 && component(4).fraction < 0.99"
            " && (abs(eta) >= 2.4 || component(2).fraction < 0.99 )"
            " && (abs(eta) >= 2.4 || component(1).fraction > 0 )"
            " && (abs(eta) >= 2.4 || component(1).number > 0 ) "
        ),
        leptonDeltaR = cms.double(0.3),
        bTagType = cms.string("combinedSecondaryVertexBJetTags"),
    ),
)

