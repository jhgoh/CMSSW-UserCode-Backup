import FWCore.ParameterSet.Config as cms

HLTPaths = {
    #Run2012A : Run2012A : Run range  190456 - 193805
    #Run2012B : Run2012B : Run range  193806 - 198021
    #Run2012C-D : Run2012B and Run2012C : Run range   198022 - YYYYYY

    "Run2012A-DoubleMu":["HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*",],
    "Run2012B-DoubleMu":["HLT_Mu17_Mu8_v*", "LLT_Mu17_TkMu8_v*",],
    "Run2012C-DoubleMu":["HLT_Mu17_Mu8_v*", "LLT_Mu17_TkMu8_v*",],
    "Run2012D-DoubleMu":["HLT_Mu17_Mu8_v*", "LLT_Mu17_TkMu8_v*",],

    "Run2012A-DoubleElectron":["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],
    "Run2012B-DoubleElectron":["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],
    "Run2012C-DoubleElectron":["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],
    "Run2012D-DoubleElectron":["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],

    "Run2012A-MuEG":["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],
    "Run2012B-MuEG":["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],
    "Run2012C-MuEG":["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],
    "Run2012D-MuEG":["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],

    "Run2012A-MuJet":["HLT_IsoMu20_eta2p1_TriCentralPFJet30_v*", "HLT_IsoMu20_eta2p1_TriCentralPFNoPUJet30_v*",],
    "Run2012B-MuJet":["HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_v*", "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v*",],
    "Run2012C-MuJet":["HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v*", "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet45_35_25_v*",],
    "Run2012D-MuJet":["HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v*", "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet45_35_25_v*",],

    "Run2012A-EleJet":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v*", "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_v*"],
    "Run2012B-EleJet":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_v*", "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_30_20_v*",],
    "Run2012C-EleJet":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_30_20_v*", "HLT_Ele25_CaloIdVT_CaloIsoVL_TrkIdVL_TrkIsoT_TriCentralPFNoPUJet45_35_25_v*"],
    "Run2012D-EleJet":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_30_20_v*", "HLT_Ele25_CaloIdVT_CaloIsoVL_TrkIdVL_TrkIsoT_TriCentralPFNoPUJet45_35_25_v*"],

    "Summer12_51X-DoubleMu"     :["HLT_Mu17_Mu8_v12", "HLT_Mu17_TkMu8_v5",],
    "Summer12_52X_GTv5-DoubleMu":["HLT_Mu17_Mu8_v16", "HLT_Mu17_TkMu8_v9",],
    "Summer12_52X_GTv9-DoubleMu":["HLT_Mu17_Mu8_v16", "HLT_Mu17_TkMu8_v9",],
    "Summer12-DoubleMu"         :["HLT_Mu17_Mu8_v17", "HLT_Mu17_TkMu8_v10",],

    "Summer12_51X-DoubleElectron"     :["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v11",],
    "Summer12_52X_GTv7-DoubleElectron":["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15",],
    "Summer12_52X_GTv9-DoubleElectron":["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17",],
    "Summer12-DoubleElectron"         :["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17",],

    
    "Summer12_51X-MuEG"     :["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v9", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v9",],
    "Summer12_52X_GTv7-MuEG":["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4",],
    "Summer12-52X_GTv9-MuEG":["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",],
    "Summer12-MuEG"         :["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",],

    "Summer12_51X-MuJet"     :["HLT_IsoMu17_eta2p1_TriCentralPFJet30_v4",],
    "Summer12_52X_GTv5-MuJet":["HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2",],
    "Summer12_52X_GTv9-MuJet":["HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2", "HLT_IsoMu20_eta2p1_TriCentralPFNoPUJet30_v2",],
    "Summer12-MuJet"         :["HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet50_40_30_v1",],

    "Summer12_51X-EleJet"     :["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v4",],
    "Summer12_52X_GTv5-EleJet":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8",],
    "Summer12_52X_GTv9-EleJet":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_v3",],
    "Summer12-EleJet"         :["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet50_40_30_v5", "HLT_Ele25_CaloIdVT_CaloIsoVL_TrkIdVL_TrkIsoT_TriCentralPFNoPUJet50_40_30_v1",],

}

