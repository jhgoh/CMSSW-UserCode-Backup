import FWCore.ParameterSet.Config as cms

muonHLTAnalyzer = cms.EDAnalyzer("MuonHLTAnalyzer",
  trigNames = cms.vstring(
    "HLT_L1MuOpen",
    "HLT_L1Mu",
    "HLT_L1Mu20",
    "HLT_L2Mu0",
    "HLT_L2Mu3",
    "HLT_L2Mu9",
    "HLT_L2Mu11",
    "HLT_Mu3",
    "HLT_IsoMu3",
    "HLT_Mu5",
    "HLT_Mu9",
    "HLT_L1DoubleMuOpen",
    "HLT_L2DoubleMu0",
    "HLT_DoubleMu0",
    "HLT_DoubleMu3"
  ),
  l1Muon = cms.InputTag("hltL1extraParticles"),
  recoMuon = cms.InputTag("muons"),
  minPt = cms.double(15.0),
  maxRelIso = cms.double(0.15),
  l1MatcherConfig = cms.PSet(
    maxDeltaR = cms.double(0.3),
    useTrack = cms.string('tracker'),
    useState = cms.string('atVertex'),
    useSimpleGeometry = cms.bool(True),
    cosmicPropagationHypothesis = cms.bool(False)
  ),
  propagatorName = cms.string('SteppingHelixPropagatorAny')
)

