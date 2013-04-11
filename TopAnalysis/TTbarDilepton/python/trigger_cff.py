import FWCore.ParameterSet.Config as cms

HLTPaths = {
    #Run2012A : Run2012A : Run range  190456 - 193805
    #Run2012B : Run2012B : Run range  193806 - 198021
    #Run2012C-D : Run2012B and Run2012C : Run range   198022 - YYYYYY

    "DoubleMu-Run2012A":["HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*",],
    "DoubleMu-Run2012B":["HLT_Mu17_Mu8_v*", "LLT_Mu17_TkMu8_v*",],
    "DoubleMu-Run2012C":["HLT_Mu17_Mu8_v*", "LLT_Mu17_TkMu8_v*",],
    "DoubleMu-Run2012D":["HLT_Mu17_Mu8_v*", "LLT_Mu17_TkMu8_v*",],

    "DoubleElectron-Run2012A":["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],
    "DoubleElectron-Run2012B":["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],
    "DoubleElectron-Run2012C":["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],
    "DoubleElectron-Run2012D":["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],

    "MuEG-Run2012A":["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],
    "MuEG-Run2012B":["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],
    "MuEG-Run2012C":["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],
    "MuEG-Run2012D":["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",],

    "MuJet-Run2012A":["HLT_IsoMu20_eta2p1_TriCentralPFJet30_v*", "HLT_IsoMu20_eta2p1_TriCentralPFNoPUJet30_v*",],
    "MuJet-Run2012B":["HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_v*", "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v*",],
    "MuJet-Run2012C":["HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v*", "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet45_35_25_v*",],
    "MuJet-Run2012D":["HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v*", "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet45_35_25_v*",],

    "EleJet-Run2012A":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v*", "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_v*"],
    "EleJet-Run2012B":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_v*", "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_30_20_v*",],
    "EleJet-Run2012C":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_30_20_v*", "HLT_Ele25_CaloIdVT_CaloIsoVL_TrkIdVL_TrkIsoT_TriCentralPFNoPUJet45_35_25_v*"],
    "EleJet-Run2012D":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_30_20_v*", "HLT_Ele25_CaloIdVT_CaloIsoVL_TrkIdVL_TrkIsoT_TriCentralPFNoPUJet45_35_25_v*"],

    "DoubleMu-Summer12_51X"     :["HLT_Mu17_Mu8_v12", "HLT_Mu17_TkMu8_v5",],
    "DoubleMu-Summer12_52X_GTv5":["HLT_Mu17_Mu8_v16", "HLT_Mu17_TkMu8_v9",],
    "DoubleMu-Summer12_52X_GTv9":["HLT_Mu17_Mu8_v16", "HLT_Mu17_TkMu8_v9",],
    "DoubleMu-Summer12"         :["HLT_Mu17_Mu8_v17", "HLT_Mu17_TkMu8_v10",],

    "DoubleElectron-Summer12_51X"     :["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v11",],
    "DoubleElectron-Summer12_52X_GTv7":["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15",],
    "DoubleElectron-Summer12_52X_GTv9":["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17",],
    "DoubleElectron-Summer12"         :["HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17",],


    "MuEG-Summer12_51X"     :["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v9", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v9",],
    "MuEG-Summer12_52X_GTv7":["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4",],
    "MuEG-Summer12_52X_GTv9":["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",],
    "MuEG-Summer12"         :["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",],

    "MuJet-Summer12_51X"     :["HLT_IsoMu17_eta2p1_TriCentralPFJet30_v4",],
    "MuJet-Summer12_52X_GTv5":["HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2",],
    "MuJet-Summer12_52X_GTv9":["HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2", "HLT_IsoMu20_eta2p1_TriCentralPFNoPUJet30_v2",],
    "MuJet-Summer12"         :["HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet50_40_30_v1",],

    "EleJet-Summer12_51X"     :["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v4",],
    "EleJet-Summer12_52X_GTv5":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8",],
    "EleJet-Summer12_52X_GTv9":["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_v3",],
    "EleJet-Summer12"         :["HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet50_40_30_v5", "HLT_Ele25_CaloIdVT_CaloIsoVL_TrkIdVL_TrkIsoT_TriCentralPFNoPUJet50_40_30_v1",],


    ## Run2011
    "DoubleElectron-Run2011A":['HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*', 'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*', 'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*'],
    "DoubleElectron-Run2011B":['HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*', 'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*', 'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*'],
    "DoubleMu-Run2011A":['HLT_DoubleMu7_*', 'HLT_Mu13_Mu8_v*', 'HLT_Mu17_Mu8_v*' ],
    "DoubleMu-Run2011B":['HLT_DoubleMu7_*', 'HLT_Mu13_Mu8_v*', 'HLT_Mu17_Mu8_v*' ],
    "MuEG-Run2011A":['HLT_Mu10_Ele10_CaloIdL_v*', 'HLT_Mu17_Ele8_CaloIdL_v*', 'HLT_Mu8_Ele17_CaloIdL_v*', 'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*', 'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*',],
    "MuEG-Run2011B":['HLT_Mu10_Ele10_CaloIdL_v*', 'HLT_Mu17_Ele8_CaloIdL_v*', 'HLT_Mu8_Ele17_CaloIdL_v*', 'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*', 'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*',],

    ## Fall11
    "DoubleMu-Fall11"      :['HLT_DoubleMu3_v*', 'HLT_DoubleMu6_v*', 'HLT_DoubleMu7_v*',],
    "DoubleElectron-Fall11":['HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*','HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*', ],
    "MuEG-Fall11"          :['HLT_Mu10_Ele10_CaloIdL_v*','HLT_Mu17_Ele8_CaloIdL_v*', 'HLT_Mu8_Ele17_CaloIdL_v*',],
}

