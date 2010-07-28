#include "HLTrigger/TPGAnalysis/plugins/MuonHLTAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <TString.h>
#include <iostream>
#include <memory>

const static int nRecoMuonCutStep = 3;
const static char* recoMuonCutStepNames[nRecoMuonCutStep] = {
  "Prompt+Muon station", "Prompt+Muon station+MisHit", "EWK"
};

//// Basic constants
const static int ptNBin = 11;
const static double ptBins[ptNBin+1] = {
  0.,4.,8.,10.,15.,18.,21.,25.,30.,40.,70.,100.
};

/*
  const int l1EtaNBin = 62;
  const double l1EtaBins[l1EtaNBin+1] = {
    -2.40, -2.35, -2.30, -2.25, -2.20, -2.15, -2.10, -2.05,-2.00, -1.95, -1.90, -1.85, -1.80, -1.75, -1.70, -1.60,-1.50, -1.40, -1.30, -1.20, -1.10, -1.00, -0.90, -0.80,-0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, -0.00,0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60,1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00, 2.05,2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40
  };
*/

const static int etaNBin = 8;
const static double etaBins[etaNBin+1] = {
  -2.40, -1.95, -1.20, -0.90, 0.00, 0.90, 1.20, 1.95, 2.40
};

struct Histograms
{
  typedef TH1F* TH1FP;

  Histograms(const std::string dirName)
  {
    edm::Service<TFileService> fs;

    TFileDirectory dir = fs->mkdir(dirName);

    hNRecoMuon = dir.make<TH1F>("hNRecoMuon", "Number of reco muons per event", 4, 1, 5);
    hNL1Muon = dir.make<TH1F>("hNL1Muon", "Number of L1 matched reco muons per event", 4, 1, 5);
    hNHLTMuon = dir.make<TH1F>("hNHLTMuon", "Number of HLT matched reco muons per event", 4, 1, 5);

    // Basic Kinematic variables
    hPt = dir.make<TH1F>("hPt", "Global muon Transverse momentum;Global muon Transverse momentum [GeV/c]", ptNBin, ptBins);
    hEta = dir.make<TH1F>("hEta", "Global muon Pseudorapidity;Global muon Pseudorapidity #eta", etaNBin, etaBins);
    hPhi = dir.make<TH1F>("hPhi", "Global muon Azimuthal angle;Global muon Azimuthal angle #phi [Radian]", 50, -3.15, 3.15);
    hQ = dir.make<TH1F>("hQ", "Global muon charge;Global muon charge", 3, -1.5, 1.5);

    hPt_B = dir.make<TH1F>("hPt_B", "Global muon Transverse momentum;Global muon Transverse momentum [GeV/c]", ptNBin, ptBins);
    hPhi_B = dir.make<TH1F>("hPhi_B", "Global muon Azimuthal angle;Global muon Azimuthal angle #phi [Radian]", 50, -3.15, 3.15);
    hQ_B = dir.make<TH1F>("hQ_B", "Global muon charge;Global muon charge", 3, -1.5, 1.5);

    hPt_O = dir.make<TH1F>("hPt_O", "Global muon Transverse momentum;Global muon Transverse momentum [GeV/c]", ptNBin, ptBins);
    hPhi_O = dir.make<TH1F>("hPhi_O", "Global muon Azimuthal angle;Global muon Azimuthal angle #phi [Radian]", 50, -3.15, 3.15);
    hQ_O = dir.make<TH1F>("hQ_O", "Global muon charge;Global muon charge", 3, -1.5, 1.5);

    hPt_E = dir.make<TH1F>("hPt_E", "Global muon Transverse momentum;Global muon Transverse momentum [GeV/c]", ptNBin, ptBins);
    hPhi_E = dir.make<TH1F>("hPhi_E", "Global muon Azimuthal angle;Global muon Azimuthal angle #phi [Radian]", 50, -3.15, 3.15);
    hQ_E = dir.make<TH1F>("hQ_E", "Global muon charge;Global muon charge", 3, -1.5, 1.5);

    // Kinematic variables matched to L1 objects
    hL1Pt = dir.make<TH1F>("hL1Pt", "L1 matched global muon Transverse momentum;Global muon Transverse momentum [GeV/c]", ptNBin, ptBins);
    hL1Eta = dir.make<TH1F>("hL1Eta", "L1 matched global muon Pseudorapidity;Global muon Pseudorapidity #eta", etaNBin, etaBins);
    hL1Phi = dir.make<TH1F>("hL1Phi", "L1 matched global muon Azumuthal angle;Global muon Azimuthal angle #phi [Radian]", 50, -3.15, 3.15);
    hL1Q = dir.make<TH1F>("hL1Q", "L1 matched global muon charge;Global muon charge", 3, -1.5, 1.5);

    hL1Pt_B = dir.make<TH1F>("hL1Pt_B", "L1 matched global muon Transverse momentum;Global muon Transverse momentum [GeV/c]", ptNBin, ptBins);
    hL1Phi_B = dir.make<TH1F>("hL1Phi_B", "L1 matched global muon Azumuthal angle;Global muon Azimuthal angle #phi [Radian]", 50, -3.15, 3.15);
    hL1Q_B = dir.make<TH1F>("hL1Q_B", "L1 matched global muon charge;Global muon charge", 3, -1.5, 1.5);

    hL1Pt_O = dir.make<TH1F>("hL1Pt_O", "L1 matched global muon Transverse momentum;Global muon Transverse momentum [GeV/c]", ptNBin, ptBins);
    hL1Phi_O = dir.make<TH1F>("hL1Phi_O", "L1 matched global muon Azumuthal angle;Global muon Azimuthal angle #phi [Radian]", 50, -3.15, 3.15);
    hL1Q_O = dir.make<TH1F>("hL1Q_O", "L1 matched global muon charge;Global muon charge", 3, -1.5, 1.5);

    hL1Pt_E = dir.make<TH1F>("hL1Pt_E", "L1 matched global muon Transverse momentum;Global muon Transverse momentum [GeV/c]", ptNBin, ptBins);
    hL1Phi_E = dir.make<TH1F>("hL1Phi_E", "L1 matched global muon Azumuthal angle;Global muon Azimuthal angle #phi [Radian]", 50, -3.15, 3.15);
    hL1Q_E = dir.make<TH1F>("hL1Q_E", "L1 matched global muon charge;Global muon charge", 3, -1.5, 1.5);

    // Kinematic variables matched to HLT objects
    hHLTPt = dir.make<TH1F>("hHLTPt", "HLT matched global muon Transverse momentum;Global muon Transverse momentum [GeV/c];Global muon Transverse momentum [GeV/c]", ptNBin, ptBins);
    hHLTEta = dir.make<TH1F>("hHLTEta", "HLT matched global muon Pseudorapidity;Global muon Pseudorapidity #eta", etaNBin, etaBins);
    hHLTPhi = dir.make<TH1F>("hHLTPhi", "HLT matched global muon Azumuthal angle;Global muon Azimuthal angle #phi [Radian]", 50, -3.15, 3.15);
    hHLTQ = dir.make<TH1F>("hHLTQ", "HLT matched global muon charge;Global muon charge", 3, -1.5, 1.5);

    hHLTPt_B = dir.make<TH1F>("hHLTPt_B", "HLT matched global muon Transverse momentum;Global muon Transverse momentum [GeV/c]", ptNBin, ptBins);
    hHLTPhi_B = dir.make<TH1F>("hHLTPhi_B", "HLT matched global muon Azumuthal angle;Global muon Azimuthal angle #phi [Radian]", 50, -3.15, 3.15);
    hHLTQ_B = dir.make<TH1F>("hHLTQ_B", "HLT matched global muon charge;Global muon charge", 3, -1.5, 1.5);

    hHLTPt_O = dir.make<TH1F>("hHLTPt_O", "HLT matched global muon Transverse momentum;Global muon Transverse momentum [GeV/c]", ptNBin, ptBins);
    hHLTPhi_O = dir.make<TH1F>("hHLTPhi_O", "HLT matched global muon Azumuthal angle;Global muon Azimuthal angle #phi [Radian]", 50, -3.15, 3.15);
    hHLTQ_O = dir.make<TH1F>("hHLTQ_O", "HLT matched global muon charge;Global muon charge", 3, -1.5, 1.5);

    hHLTPt_E = dir.make<TH1F>("hHLTPt_E", "HLT matched global muon Transverse momentum;Global muon Transverse momentum [GeV/c]", ptNBin, ptBins);
    hHLTPhi_E = dir.make<TH1F>("hHLTPhi_E", "HLT matched global muon Azumuthal angle;Global muon Azimuthal angle #phi [Radian]", 50, -3.15, 3.15);
    hHLTQ_E = dir.make<TH1F>("hHLTQ_E", "HLT matched global muon charge;Global muon charge", 3, -1.5, 1.5);

    // Number of hits and Chi^2
    hGlbNHit = dir.make<TH1F>("hGlbNHit", "Global muon Number of valid muon hits", 100, 0, 100);
    hGlbX2 = dir.make<TH1F>("hGlbX2", "Global muon Normalized #Chi^{2} of global track", 50, 0, 50);
    hTrkNHit = dir.make<TH1F>("hTrkNHit", "Global muon Number of valid tracker hits", 100, 0, 100);
    hTrkX2 = dir.make<TH1F>("hTrkX2", "Global muon Normalized #Chi^{2} of tracker track", 50, 0, 50);

    hGlbNHit_B = dir.make<TH1F>("hGlbNHit_B", "Global muon Number of valid muon hits", 100, 0, 100);
    hGlbX2_B = dir.make<TH1F>("hGlbX2_B", "Global muon Normalized #Chi^{2} of global track", 50, 0, 50);
    hTrkNHit_B = dir.make<TH1F>("hTrkNHit_B", "Global muon Number of valid tracker hits", 100, 0, 100);
    hTrkX2_B = dir.make<TH1F>("hTrkX2_B", "Global muon Normalized #Chi^{2} of tracker track", 50, 0, 50);

    hGlbNHit_O = dir.make<TH1F>("hGlbNHit_O", "Global muon Number of valid muon hits", 100, 0, 100);
    hGlbX2_O = dir.make<TH1F>("hGlbX2_O", "Global muon Normalized #Chi^{2} of global track", 50, 0, 50);
    hTrkNHit_O = dir.make<TH1F>("hTrkNHit_O", "Global muon Number of valid tracker hits", 100, 0, 100);
    hTrkX2_O = dir.make<TH1F>("hTrkX2_O", "Global muon Normalized #Chi^{2} of tracker track", 50, 0, 50);

    hGlbNHit_E = dir.make<TH1F>("hGlbNHit_E", "Global muon Number of valid muon hits", 100, 0, 100);
    hGlbX2_E = dir.make<TH1F>("hGlbX2_E", "Global muon Normalized #Chi^{2} of global track", 50, 0, 50);
    hTrkNHit_E = dir.make<TH1F>("hTrkNHit_E", "Global muon Number of valid tracker hits", 100, 0, 100);
    hTrkX2_E = dir.make<TH1F>("hTrkX2_E", "Global muon Normalized #Chi^{2} of tracker track", 50, 0, 50);

    // L1 Matching information
    hL1DeltaR = dir.make<TH1F>("hL1DeltaR", "Global muon - L1 matching #DeltaR;Position #DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}", 100, 0, 1);
    hL1DeltaPhi = dir.make<TH1F>("hL1DeltaPhi", "Global muon - L1 matching #Delta#phi;Position #Delta#phi", 100, -1, 1);
    hL1DeltaEta = dir.make<TH1F>("hL1DeltaEta", "Global muon - L1 matching #Delta#eta;Position #Delta#eta", 100, -1, 1);

    hL1DeltaR_B = dir.make<TH1F>("hL1DeltaR_B", "Global muon - L1 matching #DeltaR;Position #DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}", 100, 0, 1);
    hL1DeltaPhi_B = dir.make<TH1F>("hL1DeltaPhi_B", "Global muon - L1 matching #Delta#phi;Position #Delta#phi [Radian]", 100, -1, 1);
    hL1DeltaEta_B = dir.make<TH1F>("hL1DeltaEta_B", "Global muon - L1 matching #Delta#eta;Position #Delta#eta", 100, -1, 1);

    hL1DeltaR_O = dir.make<TH1F>("hL1DeltaR_O", "Global muon - L1 matching #DeltaR;Position #DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}", 100, 0, 1);
    hL1DeltaPhi_O = dir.make<TH1F>("hL1DeltaPhi_O", "Global muon - L1 matching #Delta#phi;Position #Delta#phi [Radian]", 100, -1, 1);
    hL1DeltaEta_O = dir.make<TH1F>("hL1DeltaEta_O", "Global muon - L1 matching #Delta#eta;Position #Delta#eta", 100, -1, 1);

    hL1DeltaR_E = dir.make<TH1F>("hL1DeltaR_E", "Global muon - L1 matching #DeltaR;Position #DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}", 100, 0, 1);
    hL1DeltaPhi_E = dir.make<TH1F>("hL1DeltaPhi_E", "Global muon - L1 matching #Delta#phi;Position #Delta#phi [Radian]", 100, -1, 1);
    hL1DeltaEta_E = dir.make<TH1F>("hL1DeltaEta_E", "Global muon - L1 matching #Delta#eta;Position #Delta#eta", 100, -1, 1);

    // HLT Matching information
    hHLTDeltaR = dir.make<TH1F>("hHLTDeltaR", "Global muon - HLT matching #DeltaR;#DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}", 100, 0, 1);
    hHLTDeltaPhi = dir.make<TH1F>("hHLTDeltaPhi", "Global muon - HLT matching #Delta#phi;#Delta#phi [Radian]", 100, -1, 1);
    hHLTDeltaEta = dir.make<TH1F>("hHLTDeltaEta", "Global muon - HLT matching #Delta#eta;#Delta#eta", 100, -1, 1);

    hHLTDeltaR_B = dir.make<TH1F>("hHLTDeltaR_B", "Global muon - HLT matching #DeltaR;#DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}", 100, 0, 1);
    hHLTDeltaPhi_B = dir.make<TH1F>("hHLTDeltaPhi_B", "Global muon - HLT matching #Delta#phi;#Delta#phi [Radian]", 100, -1, 1);
    hHLTDeltaEta_B = dir.make<TH1F>("hHLTDeltaEta_B", "Global muon - HLT matching #Delta#eta;#Delta#eta", 100, -1, 1);

    hHLTDeltaR_O = dir.make<TH1F>("hHLTDeltaR_O", "Global muon - HLT matching #DeltaR;#DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}", 100, 0, 1);
    hHLTDeltaPhi_O = dir.make<TH1F>("hHLTDeltaPhi_O", "Global muon - HLT matching #Delta#phi;#Delta#phi [Radian]", 100, -1, 1);
    hHLTDeltaEta_O = dir.make<TH1F>("hHLTDeltaEta_O", "Global muon - HLT matching #Delta#eta;#Delta#eta", 100, -1, 1);

    hHLTDeltaR_E = dir.make<TH1F>("hHLTDeltaR_E", "Global muon - HLT matching #DeltaR;#DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}", 100, 0, 1);
    hHLTDeltaPhi_E = dir.make<TH1F>("hHLTDeltaPhi_E", "Global muon - HLT matching #Delta#phi;#Delta#phi [Radian]", 100, -1, 1);
    hHLTDeltaEta_E = dir.make<TH1F>("hHLTDeltaEta_E", "Global muon - HLT matching #Delta#eta;#Delta#eta", 100, -1, 1);

    // Other cut variables
    hRelIso   = dir.make<TH1F>("hRelIso", "Global muon relative isolation", 100, 0, 10);
    hRelIso_B = dir.make<TH1F>("hRelIsoBarrel", "Global muon relative isolation", 100, 0, 10);
    hRelIso_O = dir.make<TH1F>("hRelIsoOverlap", "Global muon relative isolation", 100, 0, 10);
    hRelIso_E = dir.make<TH1F>("hRelIsoEndcap", "Global muon relative isolation", 100, 0, 10);
  }

  TH1FP hNRecoMuon, hNL1Muon, hNHLTMuon;

  // Basic kinematic variables
  TH1FP hPt, hEta, hPhi, hQ;
  TH1FP hPt_B, hPhi_B, hQ_B;
  TH1FP hPt_O, hPhi_O, hQ_O;
  TH1FP hPt_E, hPhi_E, hQ_E;

  // Kinematic variables matched to L1 object
  TH1FP hL1Pt, hL1Eta, hL1Phi, hL1Q;
  TH1FP hL1Pt_B, hL1Phi_B, hL1Q_B;
  TH1FP hL1Pt_O, hL1Phi_O, hL1Q_O;
  TH1FP hL1Pt_E, hL1Phi_E, hL1Q_E;

  // Kinematic variables matched to HLT object
  TH1FP hHLTPt, hHLTEta, hHLTPhi, hHLTQ;
  TH1FP hHLTPt_B, hHLTPhi_B, hHLTQ_B;
  TH1FP hHLTPt_O, hHLTPhi_O, hHLTQ_O;
  TH1FP hHLTPt_E, hHLTPhi_E, hHLTQ_E;

  // Number of hits and Chi^2
  TH1FP hGlbNHit, hGlbX2, hTrkNHit, hTrkX2;
  TH1FP hGlbNHit_B, hGlbX2_B, hTrkNHit_B, hTrkX2_B;
  TH1FP hGlbNHit_O, hGlbX2_O, hTrkNHit_O, hTrkX2_O;
  TH1FP hGlbNHit_E, hGlbX2_E, hTrkNHit_E, hTrkX2_E;

  // Matching information
  TH1FP hL1DeltaR, hL1DeltaPhi, hL1DeltaEta, hHLTDeltaR, hHLTDeltaPhi, hHLTDeltaEta;
  TH1FP hL1DeltaR_B, hL1DeltaPhi_B, hL1DeltaEta_B, hHLTDeltaR_B, hHLTDeltaPhi_B, hHLTDeltaEta_B;
  TH1FP hL1DeltaR_O, hL1DeltaPhi_O, hL1DeltaEta_O, hHLTDeltaR_O, hHLTDeltaPhi_O, hHLTDeltaEta_O;
  TH1FP hL1DeltaR_E, hL1DeltaPhi_E, hL1DeltaEta_E, hHLTDeltaR_E, hHLTDeltaPhi_E, hHLTDeltaEta_E;

  // Other cut variables
  TH1FP hRelIso, hRelIso_B, hRelIso_O, hRelIso_E;

};

MuonHLTAnalyzer::MuonHLTAnalyzer(const edm::ParameterSet& pset):
  l1Matcher_(pset.getParameter<edm::ParameterSet>("l1MatcherConfig"))
{
  muonL1TNames_ = pset.getParameter<VString>("muonL1TNames");
  const int nMuonL1TNames = muonL1TNames_.size();

  l1MuonTag_ = pset.getParameter<edm::InputTag>("l1Muon");
  triggerEventTag_ = pset.getParameter<edm::InputTag>("triggerEvent");
  recoMuonTag_ = pset.getParameter<edm::InputTag>("recoMuon");
  minPt_ = pset.getParameter<double>("minPt");
  maxRelIso_ = pset.getParameter<double>("maxRelIso");

  // Book histograms
  edm::Service<TFileService> fs;

  hNEvent_ = fs->make<TH1F>("hNEvent", "Number of events passing trigger paths;Trigger path", nMuonL1TNames+1, 0, nMuonL1TNames+1);
  hNEvent_->GetXaxis()->SetBinLabel(1, "All");
  for ( int muonL1TIdx=0; muonL1TIdx<nMuonL1TNames; ++muonL1TIdx )
  {
    hNEvent_->GetXaxis()->SetBinLabel(muonL1TIdx+2, muonL1TNames_[muonL1TIdx].c_str());
  }

  // Book histograms for each cut steps
  for ( int recoCutStep=0; recoCutStep<nRecoMuonCutStep; ++recoCutStep )
  {
    histograms_.push_back(Histograms(Form("CutStep%d_%s", recoCutStep, recoMuonCutStepNames[recoCutStep])));
  }
}

MuonHLTAnalyzer::~MuonHLTAnalyzer()
{
}

void MuonHLTAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  l1Matcher_.init(eventSetup);
}

void MuonHLTAnalyzer::endRun()
{
}

void MuonHLTAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::TriggerResults> trigResultHandle;
  if ( !event.getByLabel(edm::InputTag("TriggerResults", "", "HLT"), trigResultHandle) )
  {
    edm::LogError("MuonHLTAnalyzer") << "Cannot find TriggerResults\n";
    return;
  }
  const edm::TriggerResults* trigResult = trigResultHandle.product();

  edm::Handle<edm::View<l1extra::L1MuonParticle> > l1MuonHandle;
  if ( !event.getByLabel(l1MuonTag_, l1MuonHandle) )
  {
    edm::LogError("MuonHLTAnalyzer") << "Cannot find L1MuonParticle\n";
    return;
  }

  edm::Handle<trigger::TriggerEvent> triggerEventHandle;
  if ( !event.getByLabel(triggerEventTag_, triggerEventHandle) )
  {
    edm::LogError("MuonHLTAnalyzer") << "Cannot find TriggerEvvent\n";
    return;
  }
  const trigger::TriggerObjectCollection& triggerObjects = triggerEventHandle->getObjects();

  edm::Handle<edm::View<reco::Muon> > recoMuonHandle;
  if ( !event.getByLabel(recoMuonTag_, recoMuonHandle) )
  {
    edm::LogError("MuonHLTAnalyzer") << "Cannot find recoMuon\n";
    return;
  }

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  if ( !event.getByLabel("offlineBeamSpot", beamSpotHandle) )
  {
    edm::LogError("MuonHLTAnalyzer") << "Cannot find offlineBeamSpot\n";
    return;
  }

  // Count how much event passed each trigger paths
  if ( trigResult->wasrun() and trigResult->accept() )
  {
    hNEvent_->Fill(0);

    const int nTrigResult = trigResult->size();
    const edm::TriggerNames& trigNames = event.triggerNames(*trigResult);

    for ( int idxL1T=0; idxL1T<nTrigResult; ++idxL1T )
    {
      if ( !trigResult->accept(idxL1T) ) continue;

      const std::string triggerName = trigNames.triggerName(idxL1T);

      // Find Muon-HLT in the triggerResult
      int bin = -1;
      const int nMuonL1TNames = muonL1TNames_.size();
      for ( int muonL1TIdx=0; muonL1TIdx<nMuonL1TNames; ++muonL1TIdx )
      {
        if ( muonL1TNames_[muonL1TIdx] == triggerName )
        {
          bin = muonL1TIdx;
          break;
        }
      }

      if ( bin != -1 ) hNEvent_->Fill(bin+1);
    }
  } // Loop over trigger results

  // Loop over all global muons
  std::vector<int> nRecoMuon(nRecoMuonCutStep);
  std::vector<int> nL1MatchedRecoMuon(nRecoMuonCutStep);
  std::vector<int> nHLTMatchedRecoMuon(nRecoMuonCutStep);
  for ( edm::View<reco::Muon>::const_iterator recoMuon = recoMuonHandle->begin();
        recoMuon != recoMuonHandle->end(); ++recoMuon )
  {
    if ( !recoMuon->isGlobalMuon() ) continue;

    const reco::TrackRef trkTrack = recoMuon->innerTrack();
    const reco::TrackRef staTrack = recoMuon->outerTrack();
    const reco::TrackRef glbTrack = recoMuon->globalTrack();

    const reco::HitPattern& trkHit = trkTrack->hitPattern();
    //const reco::HitPattern& staHit = staTrack->hitPattern();
    const reco::HitPattern& glbHit = glbTrack->hitPattern();

    const int misHitInner = trkTrack->trackerExpectedHitsInner().numberOfHits();
    const int misHitOuter = trkTrack->trackerExpectedHitsOuter().numberOfHits();

    const double recoPt = recoMuon->pt();
    const double recoEta = recoMuon->eta();
    const double recoPhi = recoMuon->phi();
    const double recoQ = recoMuon->charge();

    const double glbX2 = glbTrack->normalizedChi2();
    const double trkX2 = trkTrack->normalizedChi2();
    const int nMuonHit = glbHit.numberOfValidMuonHits();
    const int nTrkHit = trkHit.numberOfValidTrackerHits();
    const int nPixelHit = trkHit.numberOfValidPixelHits();
    const int nMatches = recoMuon->numberOfMatches();

    const double trackIso = recoMuon->isolationR03().sumPt;
    const double caloIso = recoMuon->isolationR03().emEt + recoMuon->isolationR03().hadEt;
    const double relIso = (trackIso+caloIso)/recoPt;

    const double dxy = glbTrack->dxy(beamSpotHandle->position());

    // Check cut steps
    std::vector<bool> muonQuality(nRecoMuonCutStep);
    if ( recoPt > minPt_ )
    {
      muonQuality[0] = muon::isGoodMuon(*recoMuon, muon::GlobalMuonPromptTight) and nMatches > 1;
      muonQuality[1] = muonQuality[0] and ( misHitInner < 2 and misHitOuter < 2 );
      muonQuality[2] = nMatches > 1 and fabs(recoEta) < 2.1 and fabs(dxy) < 0.2 and glbX2 < 10 and
                       nPixelHit > 0 and nTrkHit > 10 and nMuonHit > 0;
    }

    for ( int recoCutStep=0; recoCutStep<nRecoMuonCutStep; ++recoCutStep )
    {
      if ( !muonQuality[recoCutStep] ) continue;

      ++nRecoMuon[recoCutStep];

      Histograms& h = histograms_[recoCutStep];

      h.hPt->Fill(recoPt);
      h.hEta->Fill(recoEta);
      h.hPhi->Fill(recoPhi);
      h.hQ->Fill(recoQ);

      h.hGlbNHit->Fill(nMuonHit);
      h.hTrkNHit->Fill(nTrkHit);
      h.hGlbX2->Fill(glbX2);
      h.hTrkX2->Fill(trkX2);
      h.hRelIso->Fill(relIso);

      if ( fabs(recoEta) < 0.9 )
      {
        h.hPt_B->Fill(recoPt);
        h.hPhi_B->Fill(recoPhi);

        h.hGlbNHit_B->Fill(nMuonHit);
        h.hTrkNHit_B->Fill(nTrkHit);
        h.hGlbX2_B->Fill(glbX2);
        h.hTrkX2_B->Fill(trkX2);
        h.hRelIso_B->Fill(relIso);
      }
      else if ( fabs(recoEta) < 1.2 )
      {
        h.hPt_O->Fill(recoPt);
        h.hPhi_O->Fill(recoPhi);

        h.hGlbNHit_O->Fill(nMuonHit);
        h.hTrkNHit_O->Fill(nTrkHit);
        h.hGlbX2_O->Fill(glbX2);
        h.hTrkX2_O->Fill(trkX2);
        h.hRelIso_O->Fill(relIso);
      }
      else
      {
        h.hPt_E->Fill(recoPt);
        h.hPhi_E->Fill(recoPhi);

        h.hGlbNHit_E->Fill(nMuonHit);
        h.hTrkNHit_E->Fill(nTrkHit);
        h.hGlbX2_E->Fill(glbX2);
        h.hTrkX2_E->Fill(trkX2);
        h.hRelIso_E->Fill(relIso);
      }
    }

    // Now start Matching
    TrajectoryStateOnSurface tsos = l1Matcher_.extrapolate(*recoMuon);
    if ( !tsos.isValid() ) continue;

    const double recoPosEta = tsos.globalPosition().eta();
    const double recoPosPhi = tsos.globalPosition().phi();

    double matchedL1DeltaR = -999, matchedL1DeltaPhi = -999, matchedL1DeltaEta = -999;
    l1extra::L1MuonParticle bestMatchingL1Muon;

    for ( edm::View<l1extra::L1MuonParticle>::const_iterator l1Muon = l1MuonHandle->begin();
          l1Muon != l1MuonHandle->end(); ++l1Muon )
    {
      const L1MuGMTExtendedCand& gmtCand = l1Muon->gmtMuonCand();
      const unsigned int l1Quality = gmtCand.quality();

      if ( l1Quality < 4 or l1Muon->pt() < 7 ) continue;

      const double l1PosEta = l1Muon->eta();
      const double l1PosPhi = l1Muon->phi();
      const double dR = deltaR(recoPosEta, recoPosPhi, l1PosEta, l1PosPhi);
      if ( matchedL1DeltaR < 0 or dR < matchedL1DeltaR )
      {
        matchedL1DeltaR = dR;
        matchedL1DeltaPhi = recoPosPhi - l1PosPhi;
        matchedL1DeltaEta = recoEta - l1PosEta;
        bestMatchingL1Muon = *l1Muon;
      }
    }

    for ( int recoCutStep=0; recoCutStep<nRecoMuonCutStep; ++recoCutStep )
    {
      if ( !muonQuality[recoCutStep] ) continue;

      Histograms& h = histograms_[recoCutStep];

      h.hL1DeltaR->Fill(matchedL1DeltaR);
      h.hL1DeltaPhi->Fill(matchedL1DeltaPhi);
      h.hL1DeltaEta->Fill(matchedL1DeltaEta);

      if ( fabs(recoEta) < 0.9 )
      {
        h.hL1DeltaR_B->Fill(matchedL1DeltaR);
        h.hL1DeltaPhi_B->Fill(matchedL1DeltaPhi);
        h.hL1DeltaEta_B->Fill(matchedL1DeltaEta);
      }
      else if ( fabs(recoEta) < 1.2 )
      {
        h.hL1DeltaR_O->Fill(matchedL1DeltaR);
        h.hL1DeltaPhi_O->Fill(matchedL1DeltaPhi);
        h.hL1DeltaEta_O->Fill(matchedL1DeltaEta);
      }
      else
      {
        h.hL1DeltaR_E->Fill(matchedL1DeltaR);
        h.hL1DeltaPhi_B->Fill(matchedL1DeltaPhi);
        h.hL1DeltaEta_B->Fill(matchedL1DeltaEta);
      }
    }

    if ( matchedL1DeltaR < 0 or matchedL1DeltaR > 0.3 ) continue;

    // Now we have best matching l1Muon
    for ( int recoCutStep=0; recoCutStep<nRecoMuonCutStep; ++recoCutStep )
    {
      if ( !muonQuality[recoCutStep] ) continue;

      Histograms& h = histograms_[recoCutStep];

      h.hL1Pt->Fill(recoPt);
      h.hL1Eta->Fill(recoEta);
      h.hL1Phi->Fill(recoPhi);

      if ( fabs(recoEta) < 0.9 )
      {
        h.hL1Pt_B->Fill(recoPt);
        h.hL1Phi_B->Fill(recoPhi);
      }
      else if ( fabs(recoEta) < 1.2 )
      {
        h.hL1Pt_O->Fill(recoPt);
        h.hL1Phi_O->Fill(recoPhi);
      }
      else
      {
        h.hL1Pt_E->Fill(recoPt);
        h.hL1Phi_E->Fill(recoPhi);
      }

      ++nL1MatchedRecoMuon[recoCutStep];
    }

    // Now start recoMuon-HLT matching
    for ( unsigned int filterIdx = 0; filterIdx < triggerEventHandle->sizeFilters(); ++filterIdx )
    {
      const std::string filterFullName = triggerEventHandle->filterTag(filterIdx).encode();
      const size_t fsPos = filterFullName.find_first_of(':');
      const std::string filterName = ( fsPos == std::string::npos ) ? filterFullName : filterFullName.substr(0, fsPos);

      if ( filterName != "hltSingleMu9L3Filtered9" ) continue;

      double matchedHLTDeltaR = -999;
      double matchedHLTDeltaPhi = -999;
      double matchedHLTDeltaEta = -999;
      double matchedHLTPtRes = -999;
      const trigger::Keys& trgKeys = triggerEventHandle->filterKeys(filterIdx);
      for ( trigger::Keys::const_iterator trgKey = trgKeys.begin();
            trgKey != trgKeys.end(); ++trgKey )
      {
        const double hltPt = triggerObjects[*trgKey].pt();
        const double hltEta = triggerObjects[*trgKey].eta();
        const double hltPhi = triggerObjects[*trgKey].phi();

        const double dR = deltaR(recoEta, recoPhi, hltEta, hltPhi);

        if ( matchedHLTDeltaR < 0 or dR < matchedHLTDeltaR )
        {
          matchedHLTDeltaR = dR;
          matchedHLTPtRes = hltPt == 0 ? 1e14 : (hltPt-recoPt)/hltPt;
        }
      }

      for ( int recoCutStep=0; recoCutStep<nRecoMuonCutStep; ++recoCutStep )
      {
        if ( !muonQuality[recoCutStep] ) continue;

        Histograms& h = histograms_[recoCutStep];

        h.hHLTDeltaR->Fill(matchedHLTDeltaR);
        h.hHLTDeltaPhi->Fill(matchedHLTDeltaPhi);
        h.hHLTDeltaEta->Fill(matchedHLTDeltaEta);

        if ( fabs(recoEta) < 0.9 )
        {
          h.hHLTDeltaR_B->Fill(matchedHLTDeltaR);
          h.hHLTDeltaPhi_B->Fill(matchedHLTDeltaPhi);
          h.hHLTDeltaEta_B->Fill(matchedHLTDeltaEta);
        }
        else if ( fabs(recoEta) < 1.2 )
        {
          h.hHLTDeltaR_O->Fill(matchedHLTDeltaR);
          h.hHLTDeltaPhi_O->Fill(matchedHLTDeltaPhi);
          h.hHLTDeltaEta_O->Fill(matchedHLTDeltaEta);
        }
        else
        {
          h.hHLTDeltaR_E->Fill(matchedHLTDeltaR);
          h.hHLTDeltaPhi_E->Fill(matchedHLTDeltaPhi);
          h.hHLTDeltaEta_E->Fill(matchedHLTDeltaEta);
        }
      }

      if ( matchedHLTDeltaR < 0 or matchedHLTDeltaR > 0.5 or matchedHLTPtRes > 10 ) continue;

      // Now we have best matching candidate for HLT-global muon pair
      for ( int recoCutStep=0; recoCutStep<nRecoMuonCutStep; ++recoCutStep )
      {
        if ( !muonQuality[recoCutStep] ) continue;

        Histograms& h = histograms_[recoCutStep];

        h.hHLTPt->Fill(recoPt);
        h.hHLTEta->Fill(recoEta);
        h.hHLTPhi->Fill(recoPhi);
        if ( fabs(recoEta) < 0.9 )
        {
          h.hHLTPt_B->Fill(recoPt);
          h.hHLTPhi_B->Fill(recoPhi);
        }
        else if ( fabs(recoEta) < 1.2 )
        {
          h.hHLTPt_O->Fill(recoPt);
          h.hHLTPhi_O->Fill(recoPhi);
        }
        else
        {
          h.hHLTPt_E->Fill(recoPt);
          h.hHLTPhi_E->Fill(recoPhi);
        }

        ++nHLTMatchedRecoMuon[recoCutStep];
      }
    } // Loop over HLTEvent
  } // Loop over global muons

  for ( int recoCutStep=0; recoCutStep<nRecoMuonCutStep; ++recoCutStep )
  {
    Histograms& h = histograms_[recoCutStep];

    h.hNRecoMuon->Fill(nRecoMuon[recoCutStep]);
    h.hNL1Muon->Fill(nL1MatchedRecoMuon[recoCutStep]);
    h.hNHLTMuon->Fill(nHLTMatchedRecoMuon[recoCutStep]);
  }
}

