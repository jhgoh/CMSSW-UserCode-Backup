import FWCore.ParameterSet.Config as cms

muonHLTAnalyzer = cms.EDAnalyzer("MuonHLTAnalyzer",
  interestedFilterName = cms.string("HLT_Mu9"),
  cut = cms.PSet(
    minEt = cms.double(15),
    maxL1DeltaR = cms.double(0.4),
    maxHLTDeltaR = cms.double(0.3)
  ),
  l1MatcherConfig = cms.PSet(
    maxDeltaR = cms.double(0.3),
    useTrack = cms.string('tracker'),
    useState = cms.string('atVertex'),
    useSimpleGeometry = cms.bool(True),
    cosmicPropagationHypothesis = cms.bool(False)
  )
)

