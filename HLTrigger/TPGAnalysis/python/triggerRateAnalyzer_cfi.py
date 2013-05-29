import FWCore.ParameterSet.Config as cms

muonHLTRateAnalyzer = cms.EDAnalyzer("TriggerRateAnalyzer",
    myTrigNames = cms.vstring(
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
    )
)

