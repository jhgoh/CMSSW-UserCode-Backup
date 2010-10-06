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
        "HLT_Mu7",
        "HLT_Mu9",
        "HLT_Mu11",
        "HLT_Mu15",
        "HLT_L1DoubleMuOpen",
        "HLT_L2DoubleMu0",
        "HLT_DoubleMu0",
        "HLT_DoubleMu3"
    )
)

jetHLTRateAnalyzer = cms.EDAnalyzer("TriggerRateAnalyzer",
    myTrigNames = cms.vstring(
        "HLT_L1Jet6U",
        "HLT_L1Jet10U",
        "HLT_Jet15U",
        "HLT_Jet30U",
        "HLT_Jet50U",
        "HLT_Jet70U",
        "HLT_Jet70U_v1",
        "HLT_Jet70U_v2",
        "HLT_Jet100U",
        "HLT_Jet100U_v1",
        "HLT_Jet100U_v2",
        "HLT_Jet140U",
        "HLT_Jet140U_v1",
        "HLT_Jet140U_v2"
    )
)

