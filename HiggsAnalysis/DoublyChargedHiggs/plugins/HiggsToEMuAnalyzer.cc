#include "HiggsAnalysis/DoublyChargedHiggs/interface/HiggsToEMuAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
#include "PhysicsTools/Utilities/interface/PtComparator.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <cmath>
#include <algorithm>

using namespace std;

// Purpose of this analyzer : 
// analyze pat::Electrons and pat::Muons to be used in the HiggsToEMu combination

HiggsToEMuAnalyzer::HiggsToEMuAnalyzer(const edm::ParameterSet& pset)
{
  edm::Service<TFileService> fs;

  muon1Label_ = pset.getParameter<edm::InputTag>("muon1");
  e1Label_ = pset.getParameter<edm::InputTag>("e1");
  nInterested_ = pset.getUntrackedParameter<unsigned int>("nInterested", 3);

  // Muon ID cut selection
  const string muonSelectionType = pset.getParameter<string>("muonSelection");
  if ( muonSelectionType == "GlobalMuonPromptTight" )
  {
    muonSelectionType_ = reco::Muon::GlobalMuonPromptTight;
  }
  else if ( muonSelectionType == "AllGlobalMuons" )
  {
    muonSelectionType_ = reco::Muon::AllGlobalMuons;
  }
  else
  {
    cerr << "HiggsToEMuAnalyzer::HiggsToEMuAnalyzer() : Disables good muon selection\n";
    muonSelectionType_ = reco::Muon::All;
  }

  // Histograms for muons
  TFileDirectory m1Dir = fs->mkdir("m1", "m1");

  h1_["m1_pt"] = m1Dir.make<TH1F>("pt", "Transverse momentum;Transverse momentum [GeV/c]", 100, 0, 1000);
  h1_["m1_eta"] = m1Dir.make<TH1F>("eta", "Pseudorapidity;Pseudorapidity #eta", 100, -2.5, 2.5);
  h1_["m1_phi"] = m1Dir.make<TH1F>("phi", "Azimuthal angle;Azimuthal angle [Radian]", 100, -TMath::Pi(), TMath::Pi());
  h1_["m1_trackIso"] = m1Dir.make<TH1F>("trackIso", "Track isolation;Track isolation", 100, 0, 20);
  h1_["m1_caloIso"] = m1Dir.make<TH1F>("caloIso", "Calo isolation;Calo isolation", 100, 0, 20);
  h1_["m1_relIso"] = m1Dir.make<TH1F>("relIso", "Relative isolation;Relative isolation", 100, 0, 20);

  // Leading muons
  h1_["m1_pt1_pt"] = m1Dir.make<TH1F>("pt1_pt", "Transverse momentum of leading muon;Transverse momentum [GeV/c]", 100, 0, 1000);
  h1_["m1_pt1_eta"] = m1Dir.make<TH1F>("pt1_eta", "Pseudorapidity of leading muon;Pseudorapidity #eta", 100, -2.5, 2.5);
  h1_["m1_pt1_phi"] = m1Dir.make<TH1F>("pt1_phi", "Azimuthal angle of leading muon;Azimuthal angle [Radian]", 100, -TMath::Pi(), TMath::Pi());
  h1_["m1_pt1_trackIso"] = m1Dir.make<TH1F>("pt1_trackIso", "Track isolation of leading muon;Track isolation", 100, 0, 20);
  h1_["m1_pt1_caloIso"] = m1Dir.make<TH1F>("pt1_caloIso", "Calo isolation of leading muon;Calo isolation", 100, 0, 20);
  h1_["m1_pt1_relIso"] = m1Dir.make<TH1F>("pt1_relIso", "Relative isolation of leading muon;Relative isolation", 100, 0, 20);

  // 2nd leading muons
  h1_["m1_pt2_pt"] = m1Dir.make<TH1F>("pt2_pt", "Transverse momentum of 2^{nd} leading muon;Transverse momentum [GeV/c]", 100, 0, 1000);
  h1_["m1_pt2_eta"] = m1Dir.make<TH1F>("pt2_eta", "Pseudorapidity of 2^{nd} leading muon;Pseudorapidity #eta", 100, -2.5, 2.5);
  h1_["m1_pt2_phi"] = m1Dir.make<TH1F>("pt2_phi", "Azimuthal angle of 2^{nd} leading muon;Azimuthal angle [Radian]", 100, -TMath::Pi(), TMath::Pi());
  h1_["m1_pt2_trackIso"] = m1Dir.make<TH1F>("pt2_trackIso", "Track isolation of 2^{nd} leading muon;Track isolation", 100, 0, 20);
  h1_["m1_pt2_caloIso"] = m1Dir.make<TH1F>("pt2_caloIso", "Calo isolation of 2^{nd} leading muon;Calo isolation", 100, 0, 20);
  h1_["m1_pt2_relIso"] = m1Dir.make<TH1F>("pt2_relIso", "Relative isolation of 2^{nd} leading muon;Relative isolation", 100, 0, 20);

  // 3rd leading muons
  h1_["m1_pt3_pt"] = m1Dir.make<TH1F>("pt3_pt", "Transverse momentum of 3^{rd} leading muon;Transverse momentum [GeV/c]", 100, 0, 1000);
  h1_["m1_pt3_eta"] = m1Dir.make<TH1F>("pt3_eta", "Pseudorapidity of 3^{rd} leading muon;Pseudorapidity #eta", 100, -2.5, 2.5);
  h1_["m1_pt3_phi"] = m1Dir.make<TH1F>("pt3_phi", "Azimuthal angle of 3^{rd} leading muon;Azimuthal angle [Radian]", 100, -TMath::Pi(), TMath::Pi());
  h1_["m1_pt3_trackIso"] = m1Dir.make<TH1F>("pt3_trackIso", "Track isolation of 3^{rd} leading muon;Track isolation", 100, 0, 20);
  h1_["m1_pt3_caloIso"] = m1Dir.make<TH1F>("pt3_caloIso", "Calo isolation of 3^{rd} leading muon;Calo isolation", 100, 0, 20);
  h1_["m1_pt3_relIso"] = m1Dir.make<TH1F>("pt3_relIso", "Relative isolation of 3^{rd} leading muon;Relative isolation", 100, 0, 20);

  // Histograms for electrons
  TFileDirectory e1Dir = fs->mkdir("e1", "e1");

  h1_["e1_pt"] = e1Dir.make<TH1F>("pt", "Transverse momentum;Transverse momentum [GeV/c]", 100, 0, 1000);
  h1_["e1_eta"] = e1Dir.make<TH1F>("eta", "Pseudorapidity;Pseudorapidity #eta", 100, -2.5, 2.5);
  h1_["e1_phi"] = e1Dir.make<TH1F>("phi", "Azimuthal angle;Azimuthal angle [Radian]", 100, -TMath::Pi(), TMath::Pi());
  h1_["e1_trackIso"] = e1Dir.make<TH1F>("trackIso", "Track isolation;Track isolation", 100, 0, 20);
  h1_["e1_caloIso"] = e1Dir.make<TH1F>("caloIso", "Calo isolation;Calo isolation", 100, 0, 20);
  h1_["e1_relIso"] = e1Dir.make<TH1F>("relIso", "Relative isolation;Relative isolation", 100, 0, 20);

  h1_["e1_robustLoose"] = e1Dir.make<TH1F>("robustLoose", "ElectronID with RobustLoose;ElectronID", 11, 0, 1.1);
  h1_["e1_robustTight"] = e1Dir.make<TH1F>("robustTight", "ElectronID with RobustTight;ElectronID", 11, 0, 1.1);

  // Leading electrons
  h1_["e1_pt1_pt"] = e1Dir.make<TH1F>("pt1_pt", "Transverse momentum of leading muon;Transverse momentum [GeV/c]", 100, 0, 1000);
  h1_["e1_pt1_eta"] = e1Dir.make<TH1F>("pt1_eta", "Pseudorapidity of leading muon;Pseudorapidity #eta", 100, -2.5, 2.5);
  h1_["e1_pt1_phi"] = e1Dir.make<TH1F>("pt1_phi", "Azimuthal angle of leading muon;Azimuthal angle [Radian]", 100, -TMath::Pi(), TMath::Pi());
  h1_["e1_pt1_trackIso"] = e1Dir.make<TH1F>("pt1_trackIso", "Track isolation of leading muon;Track isolation", 100, 0, 20);
  h1_["e1_pt1_caloIso"] = e1Dir.make<TH1F>("pt1_caloIso", "Calo isolation of leading muon;Calo isolation", 100, 0, 20);
  h1_["e1_pt1_relIso"] = e1Dir.make<TH1F>("pt1_relIso", "Relative isolation of leading muon;Relative isolation", 100, 0, 20);
  h1_["e1_pt1_robustLoose"] = e1Dir.make<TH1F>("pt1_robustLoose", "ElectronID of leading muon with RobustLoose;ElectronID", 11, 0, 1.1);
  h1_["e1_pt1_robustTight"] = e1Dir.make<TH1F>("pt1_robustTight", "ElectronID of leading muon with RobustTight;ElectronID", 11, 0, 1.1);

  // 2nd leading electrons
  h1_["e1_pt2_pt"] = e1Dir.make<TH1F>("pt2_pt", "Transverse momentum of 2^{nd} leading muon;Transverse momentum [GeV/c]", 100, 0, 1000);
  h1_["e1_pt2_eta"] = e1Dir.make<TH1F>("pt2_eta", "Pseudorapidity of 2^{nd} leading muon;Pseudorapidity #eta", 100, -2.5, 2.5);
  h1_["e1_pt2_phi"] = e1Dir.make<TH1F>("pt2_phi", "Azimuthal angle of 2^{nd} leading muon;Azimuthal angle [Radian]", 100, -TMath::Pi(), TMath::Pi());
  h1_["e1_pt2_trackIso"] = e1Dir.make<TH1F>("pt2_trackIso", "Track isolation of 2^{nd} leading muon;Track isolation", 100, 0, 20);
  h1_["e1_pt2_caloIso"] = e1Dir.make<TH1F>("pt2_caloIso", "Calo isolation of 2^{nd} leading muon;Calo isolation", 100, 0, 20);
  h1_["e1_pt2_relIso"] = e1Dir.make<TH1F>("pt2_relIso", "Relative isolation of 2^{nd} leading muon;Relative isolation", 100, 0, 20);
  h1_["e1_pt2_robustLoose"] = e1Dir.make<TH1F>("pt2_robustLoose", "ElectronID of 2^{nd} leading muon with RobustLoose;ElectronID", 11, 0, 1.1);
  h1_["e1_pt2_robustTight"] = e1Dir.make<TH1F>("pt2_robustTight", "ElectronID of 2^{nd} leading muon with RobustTight;ElectronID", 11, 0, 1.1);

  // 3rd leading electrons
  h1_["e1_pt3_pt"] = e1Dir.make<TH1F>("pt3_pt", "Transverse momentum of 3^{rd} leading muon;Transverse momentum [GeV/c]", 100, 0, 1000);
  h1_["e1_pt3_eta"] = e1Dir.make<TH1F>("pt3_eta", "Pseudorapidity of 3^{rd} leading muon;Pseudorapidity #eta", 100, -2.5, 2.5);
  h1_["e1_pt3_phi"] = e1Dir.make<TH1F>("pt3_phi", "Azimuthal angle of 3^{rd} leading muon;Azimuthal angle [Radian]", 100, -TMath::Pi(), TMath::Pi());
  h1_["e1_pt3_trackIso"] = e1Dir.make<TH1F>("pt3_trackIso", "Track isolation of 3^{rd} leading muon;Track isolation", 100, 0, 20);
  h1_["e1_pt3_caloIso"] = e1Dir.make<TH1F>("pt3_caloIso", "Calo isolation of 3^{rd} leading muon;Calo isolation", 100, 0, 20);
  h1_["e1_pt3_relIso"] = e1Dir.make<TH1F>("pt3_relIso", "Relative isolation of 3^{rd} leading muon;Relative isolation", 100, 0, 20);
  h1_["e1_pt3_robustLoose"] = e1Dir.make<TH1F>("pt3_robustLoose", "ElectronID of 3^{rd} leading muon with RobustLoose;ElectronID", 11, 0, 1.1);
  h1_["e1_pt3_robustTight"] = e1Dir.make<TH1F>("pt3_robustTight", "ElectronID of 3^{rd} leading muon with RobustTight;ElectronID", 11, 0, 1.1);

  // Histograms for EMu combination
  TFileDirectory emuDir = fs->mkdir("emu", "emu");

  h1_["emu_rawMass"] = emuDir.make<TH1F>("rawMass", "Raw mass of Electron-Muon pair;Mass [GeV/c^{2}]", 100, 0, 1000);
  h1_["emu_rawPt"] = emuDir.make<TH1F>("rawPt", "Raw transverse momentum of Electorn-Muon pair;Transverse momentum [GeV/c]", 100, 0, 1000);
  h1_["emu_rawEta"] = emuDir.make<TH1F>("rawEta", "Raw pseudorapidity of Electron-Muon pair;Pseudorapidity", 100, -2.5, 2.5);

  h1_["emu_mass"] = emuDir.make<TH1F>("mass", "Mass of Electron-Muon pair;Mass [GeV/c^{2}]", 100, 0, 1000);
  h1_["emu_pt"] = emuDir.make<TH1F>("pt", "Transverse momentum of Electorn-Muon pair;Transverse momentum [GeV/c]", 100, 0, 1000);
  h1_["emu_eta"] = emuDir.make<TH1F>("eta", "Pseudorapidity of Electron-Muon pair;Pseudorapidity", 100, -2.5, 2.5);

  h1_["emu_cosDeltaPhi"] = emuDir.make<TH1F>("cosDeltaPhi", "Cosine of Electron-Muon pair opening angle;Cos(#Delta#phi)", 100, -1.01, 1.01);
  h1_["emu_resolM"] = emuDir.make<TH1F>("resolM", "Mass resolution of Electron-Muon pair;Mass resolution [GeV/c^{2}]", 100, 0, 100);
  h1_["emu_resolPt"] = emuDir.make<TH1F>("resolPt", "Transverse momentum resolution of Electron-Muon pair;Transverse momentum resolution [GeV/c]", 100, 0, 100);
  h1_["emu_fitChi2"] = emuDir.make<TH1F>("fitChi2", "Fit #Chi^{2} of Electron-Muon pair;#Chi^{2}", 100, 0, 30);
  h1_["emu_trackDZ"] = emuDir.make<TH1F>("trackDZ", "Track #Delta Z of Electron-Muon pair;#Delta Z [cm]", 100, 0, 25);

  // Histograms for 4Leptons
  TFileDirectory emuemuDir = fs->mkdir("emuemu", "emuemu");

  h1_["emuemu_mass"] = emuemuDir.make<TH1F>("mass", "#sqrt{s} of (e#mu)(e#mu);#sqrt{s} [GeV/c^{2}]", 100, 0, 1000);
  h1_["emuemu_massDiff"] = emuemuDir.make<TH1F>("massDiff", "Mass difference between e#mu pairs;Mass difference [GeV/c^{2}]", 100, 0, 200);
  h1_["emuemu_massDiffSigif"] = emuemuDir.make<TH1F>("massDiffSignif", "Significance of mass difference between e#mu pairs", 100, 0, 10);

  h1_["emuemu_eeMass"] = emuemuDir.make<TH1F>("eeMass", "Mass of electron-electorn pair;Mass [GeV/c^{2}", 100, 0, 200);
  h1_["emuemu_mumuMass"] = emuemuDir.make<TH1F>("mumuMass", "Mass of muon-muon pair;Mass [GeV/c^{2}]", 100, 0, 200);

}

HiggsToEMuAnalyzer::~HiggsToEMuAnalyzer()
{
}

void HiggsToEMuAnalyzer::beginJob(const edm::EventSetup& eventSetup)
{
}

void HiggsToEMuAnalyzer::endJob()
{
}

void HiggsToEMuAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  // Retrieve data from the event
  edm::Handle<pat::MuonCollection> muon1Handle;
  edm::Handle<pat::ElectronCollection> e1Handle;
  event.getByLabel(muon1Label_, muon1Handle);
  event.getByLabel(e1Label_, e1Handle);

  // Make sorted list of candidates by decreasing pT
  pat::MuonCollection muon1Cands(muon1Handle->size());
  pat::ElectronCollection e1Cands(e1Handle->size());

  std::copy(muon1Handle->begin(), muon1Handle->end(), muon1Cands.begin());
  std::copy(e1Handle->begin(), e1Handle->end(), e1Cands.begin());

  std::sort(muon1Cands.begin(), muon1Cands.end(), GreaterByPt<pat::Muon>());
  std::sort(e1Cands.begin(), e1Cands.end(), GreaterByPt<pat::Electron>());

  const unsigned int nInterestedMuon1 = std::min(muon1Cands.size(), nInterested_);
  const unsigned int nInterestedElectron1 = std::min(e1Cands.size(), nInterested_);

  // Draw muon plots
  for ( unsigned int i=0; i<nInterestedMuon1; ++i )
  {
    const pat::Muon& muon1 = muon1Cands[i];

    h1_["m1_pt"]->Fill(muon1.pt());
    h1_["m1_eta"]->Fill(muon1.eta());
    h1_["m1_phi"]->Fill(muon1.phi());

    h1_["m1_trackIso"]->Fill(muon1.trackIso());
    h1_["m1_caloIso"]->Fill(muon1.caloIso());
    h1_["m1_relIso"]->Fill((muon1.trackIso()+muon1.caloIso())/muon1.pt());

    switch(i)
    {
      case 0:
        h1_["m1_pt1_pt"]->Fill(muon1.pt());
        h1_["m1_pt1_eta"]->Fill(muon1.eta());
        h1_["m1_pt1_phi"]->Fill(muon1.phi());

        h1_["m1_pt1_trackIso"]->Fill(muon1.trackIso());
        h1_["m1_pt1_caloIso"]->Fill(muon1.caloIso());
        h1_["m1_pt1_relIso"]->Fill((muon1.trackIso()+muon1.caloIso())/muon1.pt());
        break;
      case 1:
        h1_["m1_pt2_pt"]->Fill(muon1.pt());
        h1_["m1_pt2_eta"]->Fill(muon1.eta());
        h1_["m1_pt2_phi"]->Fill(muon1.phi());

        h1_["m1_pt2_trackIso"]->Fill(muon1.trackIso());
        h1_["m1_pt2_caloIso"]->Fill(muon1.caloIso());
        h1_["m1_pt2_relIso"]->Fill((muon1.trackIso()+muon1.caloIso())/muon1.pt());
        break;
      case 2:
        h1_["m1_pt3_pt"]->Fill(muon1.pt());
        h1_["m1_pt3_eta"]->Fill(muon1.eta());
        h1_["m1_pt3_phi"]->Fill(muon1.phi());

        h1_["m1_pt3_trackIso"]->Fill(muon1.trackIso());
        h1_["m1_pt3_caloIso"]->Fill(muon1.caloIso());
        h1_["m1_pt3_relIso"]->Fill((muon1.trackIso()+muon1.caloIso())/muon1.pt());
        break;
    }
  }

  for ( unsigned int i=0; i<nInterestedElectron1; ++i )
  {
    const pat::Electron& e1 = e1Cands[i];

    h1_["e1_pt"]->Fill(e1.pt());
    h1_["e1_eta"]->Fill(e1.eta());
    h1_["e1_phi"]->Fill(e1.phi());

    h1_["e1_trackIso"]->Fill(e1.trackIso());
    h1_["e1_caloIso"]->Fill(e1.caloIso());
    h1_["e1_relIso"]->Fill((e1.trackIso()+e1.caloIso())/e1.pt());

    h1_["e1_robustLoose"]->Fill(e1.electronID("eidRobustLoose"));
    h1_["e1_robustTight"]->Fill(e1.electronID("eidRobustTight"));

    switch(i)
    {
      case 0:
        h1_["e1_pt1_pt"]->Fill(e1.pt());
        h1_["e1_pt1_eta"]->Fill(e1.eta());
        h1_["e1_pt1_phi"]->Fill(e1.phi());

        h1_["e1_pt1_trackIso"]->Fill(e1.trackIso());
        h1_["e1_pt1_caloIso"]->Fill(e1.caloIso());
        h1_["e1_pt1_relIso"]->Fill((e1.trackIso()+e1.caloIso())/e1.pt());

        h1_["e1_pt1_robustLoose"]->Fill(e1.electronID("eidRobustLoose"));
        h1_["e1_pt1_robustTight"]->Fill(e1.electronID("eidRobustTight"));
        break;
      case 1:
        h1_["e1_pt2_pt"]->Fill(e1.pt());
        h1_["e1_pt2_eta"]->Fill(e1.eta());
        h1_["e1_pt2_phi"]->Fill(e1.phi());

        h1_["e1_pt2_trackIso"]->Fill(e1.trackIso());
        h1_["e1_pt2_caloIso"]->Fill(e1.caloIso());
        h1_["e1_pt2_relIso"]->Fill((e1.trackIso()+e1.caloIso())/e1.pt());

        h1_["e1_pt2_robustLoose"]->Fill(e1.electronID("eidRobustLoose"));
        h1_["e1_pt2_robustTight"]->Fill(e1.electronID("eidRobustTight"));
        break;
      case 2:
        h1_["e1_pt3_pt"]->Fill(e1.pt());
        h1_["e1_pt3_eta"]->Fill(e1.eta());
        h1_["e1_pt3_phi"]->Fill(e1.phi());

        h1_["e1_pt3_trackIso"]->Fill(e1.trackIso());
        h1_["e1_pt3_caloIso"]->Fill(e1.caloIso());
        h1_["e1_pt3_relIso"]->Fill((e1.trackIso()+e1.caloIso())/e1.pt());

        h1_["e1_pt3_robustLoose"]->Fill(e1.electronID("eidRobustLoose"));
        h1_["e1_pt3_robustTight"]->Fill(e1.electronID("eidRobustTight"));
        break;
    }
  }

  // Make emu combinations
  pat::CompositeCandidateCollection emuCands;
  OverlapChecker isOverlap;
  AddFourMomenta addP4;

  unsigned int nPosMuon1 = 0, nNegMuon1 = 0;
  for ( pat::MuonCollection::const_iterator muon1Iter = muon1Cands.begin();
        muon1Iter != muon1Cands.end(); ++muon1Iter )
  {
    if ( nPosMuon1 > nInterestedMuon1 || nNegMuon1 > nInterestedMuon1 ) continue;

    const pat::Muon& muon1 = *muon1Iter;

    if ( muon1.charge() > 0 ) ++nPosMuon1;
    else ++nNegMuon1;

    unsigned int nPosElectron1 = 0, nNegElectron1 = 0;
    for ( pat::ElectronCollection::const_iterator e1Iter = e1Cands.begin();
          e1Iter != e1Cands.end(); ++e1Iter )
    {
      if ( nPosElectron1 > nInterestedElectron1 || nNegElectron1 > nInterestedElectron1 ) continue;

      const pat::Electron& e1 = *e1Iter;

      if ( muon1.charge()*e1.charge() < 0 ) continue;
      if ( isOverlap(muon1, e1) ) continue;

      pat::CompositeCandidate emuCand;
      emuCand.addDaughter(muon1, "mu");
      emuCand.addDaughter(e1, "e");

      addP4.set(emuCand);

      h1_["emu_rawMass"]->Fill(emuCand.mass());
      h1_["emu_rawPt"]->Fill(emuCand.pt());
      h1_["emu_rawEta"]->Fill(emuCand.eta());

      // Do vertex fit
//      vtxFitter.set(emuCand);

/*
      h1_["emu_mass"]->Fill(emuCand.mass());
      h1_["emu_pt"]->Fill(emuCand.pt());
      h1_["emu_eta"]->Fill(emuCand.eta());
      
      h1_["emu_resolM"]->Fill(emuCand.resolM());
      h1_["emu_resolPt"]->Fill(emuCand.resolPt());
      h1_["emu_fitChi2"]->Fill(emuCand.normalizedChi2());
*/
      h1_["emu_cosDeltaPhi"]->Fill(cos(muon1.phi()-e1.phi()));
      h1_["emu_trackDZ"]->Fill(fabs(muon1.combinedMuon()->dz()-e1.gsfTrack()->dz()));

      emuCands.push_back(emuCand);
    }
  }

  // Now consider 4-lepton combination
  for ( pat::CompositeCandidateCollection::const_iterator emuCandIter1 = emuCands.begin();
        emuCandIter1 != emuCands.end(); ++emuCandIter1 )
  {
    const pat::CompositeCandidate& emuCand1 = *emuCandIter1;

    for ( pat::CompositeCandidateCollection::const_iterator emuCandIter2 = emuCandIter1+1;
          emuCandIter2 != emuCands.end(); ++emuCandIter2 )
    {
      const pat::CompositeCandidate& emuCand2 = *emuCandIter2;

      if ( emuCand1.charge()+emuCand2.charge() != 0 ) continue;
      if ( isOverlap(emuCand1, emuCand2) ) continue;

      pat::CompositeCandidate fourLeptons;
      fourLeptons.addDaughter(emuCand1, "emu1");
      fourLeptons.addDaughter(emuCand2, "emu2");

      addP4.set(fourLeptons);

      h1_["emuemu_mass"]->Fill(fourLeptons.mass());
      h1_["emuemu_massDiff"]->Fill(fabs(emuCand1.mass()-emuCand2.mass()));

      const pat::Electron* e1 = dynamic_cast<const pat::Electron*>(emuCand1.daughter("e"));
      const pat::Electron* e2 = dynamic_cast<const pat::Electron*>(emuCand2.daughter("e"));
      const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(emuCand1.daughter("mu"));
      const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(emuCand2.daughter("mu"));

      if ( !e1 or !e2 or !muon1 or !muon2 ) continue;

      pat::CompositeCandidate eeCand;
      pat::CompositeCandidate mumuCand;

      eeCand.addDaughter(*e1, "e1");
      eeCand.addDaughter(*e2, "e2");

      mumuCand.addDaughter(*muon1, "muon1");
      mumuCand.addDaughter(*muon2, "muon2");

      addP4.set(eeCand);
      addP4.set(mumuCand);

      h1_["emuemu_eeMass"]->Fill(eeCand.mass());
      h1_["emuemu_mumuMass"]->Fill(mumuCand.mass());
    }
  }
}

