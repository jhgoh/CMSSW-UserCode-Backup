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
#include "DataFormats/Math/interface/deltaPhi.h"

#include <TString.h>
#include <iostream>
#include <memory>

const static int nRecoMuonCutStep = 6;
const static char* recoMuonCutStepNames[nRecoMuonCutStep] = {
  "All", "Prompt", "MuonSTN", "MisHit", "EWK", "EWK+Iso"
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

  //// Basic constants
  const int l1PtNBin = 11;
  const double l1PtBins[l1PtNBin+1] = {
    0.,4.,8.,10.,15.,18.,21.,25.,30.,40.,70.,100.
  };

  const int l1EtaNBin = 62;
  const double l1EtaBins[l1EtaNBin+1] = {
    -2.40, -2.35, -2.30, -2.25, -2.20, -2.15, -2.10, -2.05,-2.00, -1.95, -1.90, -1.85, -1.80, -1.75, -1.70, -1.60,-1.50, -1.40, -1.30, -1.20, -1.10, -1.00, -0.90, -0.80,-0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, -0.00,0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60,1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00, 2.05,2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40
  };

  hNEvent_ = fs->make<TH1F>("hNEvent", "Number of events passing trigger paths;Trigger path", nMuonL1TNames+1, 0, nMuonL1TNames+1);
  hNEvent_->GetXaxis()->SetBinLabel(1, "All");
  for ( int muonL1TIdx=0; muonL1TIdx<nMuonL1TNames; ++muonL1TIdx )
  {
    hNEvent_->GetXaxis()->SetBinLabel(muonL1TIdx+2, muonL1TNames_[muonL1TIdx].c_str());
  }

  // Book histograms for each cut steps
  for ( int recoCutStep=0; recoCutStep<nRecoMuonCutStep; ++recoCutStep )
  {
    TFileDirectory dir = fs->mkdir(Form("CutStep%d_%s", recoCutStep, recoMuonCutStepNames[recoCutStep]));

    hNRecoMuon_.push_back(dir.make<TH1F>("hNRecoMuon", "Number of recoMuons per event", nRecoMuonCutStep, 0, nRecoMuonCutStep));
    hNL1MatchedRecoMuon_.push_back(dir.make<TH1F>("hNL1MatchedRecoMuon", "Number of L1 matched recoMuons per event", nRecoMuonCutStep, 0, nRecoMuonCutStep));
    hNHLTMatchedRecoMuon_.push_back(dir.make<TH1F>("hNHLTMatchedRecoMuon", "Number of HLT matched recoMuons per event", nRecoMuonCutStep, 0, nRecoMuonCutStep));

    // Histograms for track variables
    hPt_.push_back(dir.make<TH1F>("hPt", "Global muon Transverse momentum", 50, 0, 50));
    hEta_.push_back(dir.make<TH1F>("hEta", "Global muon Pseudorapidity", 100, -2.5, 2.5));
    hPhi_.push_back(dir.make<TH1F>("hPhi", "Global muon Azimuthal angle", 100, -3.15, 3.15));
    hQ_.push_back(dir.make<TH1F>("hQ", "Global muon charge", 3, -1.5, 1.5));
    hPtWithL1Bin_.push_back(dir.make<TH1F>("hPtWithL1Bin", "Global muon Transverse momentum", l1PtNBin, l1PtBins));
    hEtaWithL1Bin_.push_back(dir.make<TH1F>("hEtaWithL1Bin", "Global muon Pseudorapidity", l1EtaNBin, l1EtaBins));

    hNGlbHit_.push_back(dir.make<TH1F>("hNGlbHit", "Global muon Number of valid muon hits", 100, 0, 100));
    hNTrkHit_.push_back(dir.make<TH1F>("hNTrkHit", "Global muon Number of valid tracker hits", 100, 0, 100));
    hGlbX2_.push_back(dir.make<TH1F>("hGlbX2", "Global muon Normalized #Chi^{2} of global track", 50, 0, 50));
    hTrkX2_.push_back(dir.make<TH1F>("hTrkX2", "Global muon Normalized #Chi^{2} of tracker track", 50, 0, 50));
    hRelIso_.push_back(dir.make<TH1F>("hRelIso", "Global muon relative isolation", 100, 0, 20));

    // Histograms for L1 matching
    hL1DeltaR_.push_back(dir.make<TH1F>("hL1DeltaR", "Global muon - L1 muon matching #DeltaR;#DeltaR = #sqrt{#Delta#phi^{2} + #Delta#eta^{2}}", 50, 0, 2));
    hL1DeltaPhi_.push_back(dir.make<TH1F>("hL1DeltaPhi", "Global muon - L1 muon matching #DeltaR;#Delta#phi [Radian]", 50, 0, 2));
    hL1DeltaEta_.push_back(dir.make<TH1F>("hL1DeltaEta", "Global muon - L1 muon matching #Delta#eta;#Delta#eta", 50, 0, 2));

    // Histograms for L1 variables
    hMatchedL1Pt_.push_back(dir.make<TH1F>("hMatchedL1Pt", "Matched L1 p_{T};L1 muon p_{T} [GeV/c]", l1PtNBin, l1PtBins));
    hMatchedL1Eta_.push_back(dir.make<TH1F>("hMatchedL1Eta", "Matched L1 #eta;#eta", l1EtaNBin, l1EtaBins));
    hMatchedL1Phi_.push_back(dir.make<TH1F>("hMatchedL1Phi", "MatchedL1 #phi;#phi [Radian]", 50, -3.15, 3.15));
    hL1MatchedGlbPt_.push_back(dir.make<TH1F>("hL1MatchedGlbPt", "L1 matched global muon p_{T};global muon p_{T} [GeV/c]", l1PtNBin, l1PtBins));
    hL1MatchedGlbEta_.push_back(dir.make<TH1F>("hL1MatchedGlbEta", "L1 matched global muon #eta;global muon #eta", l1EtaNBin, l1EtaBins));
    hL1MatchedGlbPhi_.push_back(dir.make<TH1F>("hL1MatchedGlbPhi", "L1 matched global muon #phi;global muon #phi", 50, -3.15, 3.15));

    // Histograms for HLT matching
    hHLTDeltaR_.push_back(dir.make<TH1F>("hHLTDeltaR", "Global muon - HLT muon matching #DeltaR;#DeltaR = #sqrt{#Delta#phi^{2} + #Delta#eta^{2}}", 50, 0, 2));
    hHLTDeltaPhi_.push_back(dir.make<TH1F>("hHLTDeltaPhi", "Global muon - HLT muon matching #Delta#phi;#Delta#phi [Radian]", 50, 0, 2));
    hHLTDeltaEta_.push_back(dir.make<TH1F>("hHLTDeltaEta", "Global muon - HLT muon matching #Delta#eta;#Delta#eta", 50, 0, 2));

    // Histograms for HLT variables
    hHLTMatchedGlbPt_.push_back(dir.make<TH1F>("hHLTMatchedGlbPt", "HLT matched global muon p_{T};global muon p_{T} [GeV/c]", 50, 0, 50));
    hHLTMatchedGlbEta_.push_back(dir.make<TH1F>("hHLTMatchedGlbEta", "HLT matched global muon #eta;global muon #eta", 50, -2.5, 2.5));
    hHLTMatchedGlbPhi_.push_back(dir.make<TH1F>("hHLTMatchedGlbPhi", "HLT matched global muon #phi;global muon #phi", 50, -3.15, 3.15));

    // Histograms in Barrel region
    hPtBarrel_.push_back(dir.make<TH1F>("hPtBarrel", "Global muon Transverse momentum", 50, 0, 50));
    hEtaBarrel_.push_back(dir.make<TH1F>("hEtaBarrel", "Global muon Pseudorapidity", 100, -2.5, 2.5));
    hPhiBarrel_.push_back(dir.make<TH1F>("hPhiBarrel", "Global muon Azimuthal angle", 100, -3.15, 3.15));
    hPtWithL1BinBarrel_.push_back(dir.make<TH1F>("hPtWithL1BinBarrel", "Global muon Transverse momentum", l1PtNBin, l1PtBins));
    hEtaWithL1BinBarrel_.push_back(dir.make<TH1F>("hEtaWithL1BinBarrel", "Global muon Pseudorapidity", l1EtaNBin, l1EtaBins));

    hNGlbHitBarrel_.push_back(dir.make<TH1F>("hNGlbHitBarrel", "Global muon Number of valid muon hits", 100, 0, 100));
    hNTrkHitBarrel_.push_back(dir.make<TH1F>("hNTrkHitBarrel", "Global muon Number of valid tracker hits", 100, 0, 100));
    hGlbX2Barrel_.push_back(dir.make<TH1F>("hGlbX2Barrel", "Global muon Normalized #Chi^{2} of global track", 50, 0, 50));
    hTrkX2Barrel_.push_back(dir.make<TH1F>("hTrkX2Barrel", "Global muon Normalized #Chi^{2} of tracker track", 50, 0, 50));
    hRelIsoBarrel_.push_back(dir.make<TH1F>("hRelIsoBarrel", "Global muon relative isolation", 100, 0, 20));

    hL1DeltaRBarrel_.push_back(dir.make<TH1F>("hL1DeltaRBarrel", "Global muon - L1 muon matching #DeltaR;#DeltaR = #sqrt{#Delta#phi^{2} + #Delta#eta^{2}}", 50, 0, 2));
    hL1DeltaPhiBarrel_.push_back(dir.make<TH1F>("hL1DeltaPhiBarrel", "Global muon - L1 muon matching #DeltaR;#Delta#phi [Radian]", 50, 0, 2));
    hL1DeltaEtaBarrel_.push_back(dir.make<TH1F>("hL1DeltaEtaBarrel", "Global muon - L1 muon matching #Delta#eta;#Delta#eta", 50, 0, 2));

    hMatchedL1PtBarrel_.push_back(dir.make<TH1F>("hMatchedL1PtBarrel", "Matched L1 p_{T};L1 muon p_{T} [GeV/c]", l1PtNBin, l1PtBins));
    hMatchedL1EtaBarrel_.push_back(dir.make<TH1F>("hMatchedL1EtaBarrel", "Matched L1 #eta;#eta", l1EtaNBin, l1EtaBins));
    hMatchedL1PhiBarrel_.push_back(dir.make<TH1F>("hMatchedL1PhiBarrel", "MatchedL1 #phi;#phi [Radian]", 50, -3.15, 3.15));
    hL1MatchedGlbPtBarrel_.push_back(dir.make<TH1F>("hL1MatchedGlbPtBarrel", "L1 matched global muon p_{T};global muon p_{T} [GeV/c]", l1PtNBin, l1PtBins));
    hL1MatchedGlbEtaBarrel_.push_back(dir.make<TH1F>("hL1MatchedGlbEtaBarrel", "L1 matched global muon #eta;global muon #eta", l1EtaNBin, l1EtaBins));
    hL1MatchedGlbPhiBarrel_.push_back(dir.make<TH1F>("hL1MatchedGlbPhiBarrel", "L1 matched global muon #phi;global muon #phi", 50, -3.15, 3.15));

    hHLTDeltaRBarrel_.push_back(dir.make<TH1F>("hHLTDeltaRBarrel", "Global muon - HLT muon matching #DeltaR;#DeltaR = #sqrt{#Delta#phi^{2} + #Delta#eta^{2}}", 50, 0, 2));
    hHLTDeltaPhiBarrel_.push_back(dir.make<TH1F>("hHLTDeltaPhiBarrel", "Global muon - HLT muon matching #Delta#phi;#Delta#phi [Radian]", 50, 0, 2));
    hHLTDeltaEtaBarrel_.push_back(dir.make<TH1F>("hHLTDeltaEtaBarrel", "Global muon - HLT muon matching #Delta#eta;#Delta#eta", 50, 0, 2));

    hHLTMatchedGlbPtBarrel_.push_back(dir.make<TH1F>("hHLTMatchedGlbPtBarrel", "HLT matched global muon p_{T};global muon p_{T} [GeV/c]", 50, 0, 50));
    hHLTMatchedGlbEtaBarrel_.push_back(dir.make<TH1F>("hHLTMatchedGlbEtaBarrel", "HLT matched global muon #eta;global muon #eta", 50, -2.5, 2.5));
    hHLTMatchedGlbPhiBarrel_.push_back(dir.make<TH1F>("hHLTMatchedGlbPhiBarrel", "HLT matched global muon #phi;global muon #phi", 50, -3.15, 3.15));

    // Histograms in Overlap region
    hPtOverlap_.push_back(dir.make<TH1F>("hPtOverlap", "Global muon Transverse momentum", 50, 0, 50));
    hEtaOverlap_.push_back(dir.make<TH1F>("hEtaOverlap", "Global muon Pseudorapidity", 100, -2.5, 2.5));
    hPhiOverlap_.push_back(dir.make<TH1F>("hPhiOverlap", "Global muon Azimuthal angle", 100, -3.15, 3.15));
    hPtWithL1BinOverlap_.push_back(dir.make<TH1F>("hPtWithL1BinOverlap", "Global muon Transverse momentum", l1PtNBin, l1PtBins));
    hEtaWithL1BinOverlap_.push_back(dir.make<TH1F>("hEtaWithL1BinOverlap", "Global muon Pseudorapidity", l1EtaNBin, l1EtaBins));

    hNGlbHitOverlap_.push_back(dir.make<TH1F>("hNGlbHitOverlap", "Global muon Number of valid muon hits", 100, 0, 100));
    hNTrkHitOverlap_.push_back(dir.make<TH1F>("hNTrkHitOverlap", "Global muon Number of valid tracker hits", 100, 0, 100));
    hGlbX2Overlap_.push_back(dir.make<TH1F>("hGlbX2Overlap", "Global muon Normalized #Chi^{2} of global track", 50, 0, 50));
    hTrkX2Overlap_.push_back(dir.make<TH1F>("hTrkX2Overlap", "Global muon Normalized #Chi^{2} of tracker track", 50, 0, 50));
    hRelIsoOverlap_.push_back(dir.make<TH1F>("hRelIsoOverlap", "Global muon relative isolation", 100, 0, 20));

    hL1DeltaROverlap_.push_back(dir.make<TH1F>("hL1DeltaROverlap", "Global muon - L1 muon matching #DeltaR;#DeltaR = #sqrt{#Delta#phi^{2} + #Delta#eta^{2}}", 50, 0, 2));
    hL1DeltaPhiOverlap_.push_back(dir.make<TH1F>("hL1DeltaPhiOverlap", "Global muon - L1 muon matching #DeltaR;#Delta#phi [Radian]", 50, 0, 2));
    hL1DeltaEtaOverlap_.push_back(dir.make<TH1F>("hL1DeltaEtaOverlap", "Global muon - L1 muon matching #Delta#eta;#Delta#eta", 50, 0, 2));

    hMatchedL1PtOverlap_.push_back(dir.make<TH1F>("hMatchedL1PtOverlap", "Matched L1 p_{T};L1 muon p_{T} [GeV/c]", l1PtNBin, l1PtBins));
    hMatchedL1EtaOverlap_.push_back(dir.make<TH1F>("hMatchedL1EtaOverlap", "Matched L1 #eta;#eta", l1EtaNBin, l1EtaBins));
    hMatchedL1PhiOverlap_.push_back(dir.make<TH1F>("hMatchedL1PhiOverlap", "MatchedL1 #phi;#phi [Radian]", 50, -3.15, 3.15));
    hL1MatchedGlbPtOverlap_.push_back(dir.make<TH1F>("hL1MatchedGlbPtOverlap", "L1 matched global muon p_{T};global muon p_{T} [GeV/c]", l1PtNBin, l1PtBins));
    hL1MatchedGlbEtaOverlap_.push_back(dir.make<TH1F>("hL1MatchedGlbEtaOverlap", "L1 matched global muon #eta;global muon #eta", l1EtaNBin, l1EtaBins));
    hL1MatchedGlbPhiOverlap_.push_back(dir.make<TH1F>("hL1MatchedGlbPhiOverlap", "L1 matched global muon #phi;global muon #phi", 50, -3.15, 3.15));

    hHLTDeltaROverlap_.push_back(dir.make<TH1F>("hHLTDeltaROverlap", "Global muon - HLT muon matching #DeltaR;#DeltaR = #sqrt{#Delta#phi^{2} + #Delta#eta^{2}}", 50, 0, 2));
    hHLTDeltaPhiOverlap_.push_back(dir.make<TH1F>("hHLTDeltaPhiOverlap", "Global muon - HLT muon matching #Delta#phi;#Delta#phi [Radian]", 50, 0, 2));
    hHLTDeltaEtaOverlap_.push_back(dir.make<TH1F>("hHLTDeltaEtaOverlap", "Global muon - HLT muon matching #Delta#eta;#Delta#eta", 50, 0, 2));

    hHLTMatchedGlbPtOverlap_.push_back(dir.make<TH1F>("hHLTMatchedGlbPtOverlap", "HLT matched global muon p_{T};global muon p_{T} [GeV/c]", 50, 0, 50));
    hHLTMatchedGlbEtaOverlap_.push_back(dir.make<TH1F>("hHLTMatchedGlbEtaOverlap", "HLT matched global muon #eta;global muon #eta", 50, -2.5, 2.5));
    hHLTMatchedGlbPhiOverlap_.push_back(dir.make<TH1F>("hHLTMatchedGlbPhiOverlap", "HLT matched global muon #phi;global muon #phi", 50, -3.15, 3.15));

    // Histograms in Endcap region
    hPtEndcap_.push_back(dir.make<TH1F>("hPtEndcap", "Global muon Transverse momentum", 50, 0, 50));
    hEtaEndcap_.push_back(dir.make<TH1F>("hEtaEndcap", "Global muon Pseudorapidity", 100, -2.5, 2.5));
    hPhiEndcap_.push_back(dir.make<TH1F>("hPhiEndcap", "Global muon Azimuthal angle", 100, -3.15, 3.15));
    hPtWithL1BinEndcap_.push_back(dir.make<TH1F>("hPtWithL1BinEndcap", "Global muon Transverse momentum", l1PtNBin, l1PtBins));
    hEtaWithL1BinEndcap_.push_back(dir.make<TH1F>("hEtaWithL1BinEndcap", "Global muon Pseudorapidity", l1EtaNBin, l1EtaBins));

    hNGlbHitEndcap_.push_back(dir.make<TH1F>("hNGlbHitEndcap", "Global muon Number of valid muon hits", 100, 0, 100));
    hNTrkHitEndcap_.push_back(dir.make<TH1F>("hNTrkHitEndcap", "Global muon Number of valid tracker hits", 100, 0, 100));
    hGlbX2Endcap_.push_back(dir.make<TH1F>("hGlbX2Endcap", "Global muon Normalized #Chi^{2} of global track", 50, 0, 50));
    hTrkX2Endcap_.push_back(dir.make<TH1F>("hTrkX2Endcap", "Global muon Normalized #Chi^{2} of tracker track", 50, 0, 50));
    hRelIsoEndcap_.push_back(dir.make<TH1F>("hRelIsoEndcap", "Global muon relative isolation", 100, 0, 20));

    hL1DeltaREndcap_.push_back(dir.make<TH1F>("hL1DeltaREndcap", "Global muon - L1 muon matching #DeltaR;#DeltaR = #sqrt{#Delta#phi^{2} + #Delta#eta^{2}}", 50, 0, 2));
    hL1DeltaPhiEndcap_.push_back(dir.make<TH1F>("hL1DeltaPhiEndcap", "Global muon - L1 muon matching #DeltaR;#Delta#phi [Radian]", 50, 0, 2));
    hEtaWithL1BinBarrel_.push_back(dir.make<TH1F>("hEtaWithL1BinBarrel", "Global muon Pseudorapidity", l1EtaNBin, l1EtaBins));
    hL1DeltaEtaEndcap_.push_back(dir.make<TH1F>("hL1DeltaEtaEndcap", "Global muon - L1 muon matching #Delta#eta;#Delta#eta", 50, 0, 2));

    hMatchedL1PtEndcap_.push_back(dir.make<TH1F>("hMatchedL1PtEndcap", "Matched L1 p_{T};L1 muon p_{T} [GeV/c]", l1PtNBin, l1PtBins));
    hMatchedL1EtaEndcap_.push_back(dir.make<TH1F>("hMatchedL1EtaEndcap", "Matched L1 #eta;#eta", l1EtaNBin, l1EtaBins));
    hMatchedL1PhiEndcap_.push_back(dir.make<TH1F>("hMatchedL1PhiEndcap", "MatchedL1 #phi;#phi [Radian]", 50, -3.15, 3.15));
    hL1MatchedGlbPtEndcap_.push_back(dir.make<TH1F>("hL1MatchedGlbPtEndcap", "L1 matched global muon p_{T};global muon p_{T} [GeV/c]", l1PtNBin, l1PtBins));
    hL1MatchedGlbEtaEndcap_.push_back(dir.make<TH1F>("hL1MatchedGlbEtaEndcap", "L1 matched global muon #eta;global muon #eta", l1EtaNBin, l1EtaBins));
    hL1MatchedGlbPhiEndcap_.push_back(dir.make<TH1F>("hL1MatchedGlbPhiEndcap", "L1 matched global muon #phi;global muon #phi", 50, -3.15, 3.15));

    hHLTDeltaREndcap_.push_back(dir.make<TH1F>("hHLTDeltaREndcap", "Global muon - HLT muon matching #DeltaR;#DeltaR = #sqrt{#Delta#phi^{2} + #Delta#eta^{2}}", 50, 0, 2));
    hHLTDeltaPhiEndcap_.push_back(dir.make<TH1F>("hHLTDeltaPhiEndcap", "Global muon - HLT muon matching #Delta#phi;#Delta#phi [Radian]", 50, 0, 2));
    hHLTDeltaEtaEndcap_.push_back(dir.make<TH1F>("hHLTDeltaEtaEndcap", "Global muon - HLT muon matching #Delta#eta;#Delta#eta", 50, 0, 2));

    hHLTMatchedGlbPtEndcap_.push_back(dir.make<TH1F>("hHLTMatchedGlbPtEndcap", "HLT matched global muon p_{T};global muon p_{T} [GeV/c]", 50, 0, 50));
    hHLTMatchedGlbEtaEndcap_.push_back(dir.make<TH1F>("hHLTMatchedGlbEtaEndcap", "HLT matched global muon #eta;global muon #eta", 50, -2.5, 2.5));
    hHLTMatchedGlbPhiEndcap_.push_back(dir.make<TH1F>("hHLTMatchedGlbPhiEndcap", "HLT matched global muon #phi;global muon #phi", 50, -3.15, 3.15));
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
    const reco::HitPattern& staHit = staTrack->hitPattern();
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
      muonQuality[0] = true;
      muonQuality[1] = muon::isGoodMuon(*recoMuon, muon::GlobalMuonPromptTight);
      muonQuality[2] = nMatches > 1;
      muonQuality[3] = ( misHitInner < 2 and misHitOuter < 2 );
      muonQuality[4] = ( muonQuality[2] and fabs(recoEta) < 2.1 and fabs(dxy) < 0.2 and glbX2 < 10 and 
                         nPixelHit > 0 and nTrkHit > 10 and nMuonHit > 0 );
      muonQuality[5] = ( muonQuality[4] and relIso < maxRelIso_ );
    }

    for ( int recoCutStep=0; recoCutStep<nRecoMuonCutStep; ++recoCutStep )
    {
      if ( !muonQuality[recoCutStep] ) continue;

      ++nRecoMuon[recoCutStep];

      hPt_[recoCutStep]->Fill(recoPt);
      hEta_[recoCutStep]->Fill(recoEta);
      hPhi_[recoCutStep]->Fill(recoPhi);
      hQ_[recoCutStep]->Fill(recoQ);
      hPtWithL1Bin_[recoCutStep]->Fill(recoPt);
      hEtaWithL1Bin_[recoCutStep]->Fill(recoEta);

      hNGlbHit_[recoCutStep]->Fill(nMuonHit);
      hNTrkHit_[recoCutStep]->Fill(nTrkHit);
      hGlbX2_[recoCutStep]->Fill(glbX2);
      hTrkX2_[recoCutStep]->Fill(trkX2);
      hRelIso_[recoCutStep]->Fill(relIso);

      if ( fabs(recoEta) < 0.9 )
      {
        hPtBarrel_[recoCutStep]->Fill(recoPt);
        hEtaBarrel_[recoCutStep]->Fill(recoEta);
        hPhiBarrel_[recoCutStep]->Fill(recoPhi);
        hPtWithL1BinBarrel_[recoCutStep]->Fill(recoPt);
        hEtaWithL1BinBarrel_[recoCutStep]->Fill(recoEta);

        hNGlbHitBarrel_[recoCutStep]->Fill(nMuonHit);
        hNTrkHitBarrel_[recoCutStep]->Fill(nTrkHit);
        hGlbX2Barrel_[recoCutStep]->Fill(glbX2);
        hTrkX2Barrel_[recoCutStep]->Fill(trkX2);
        hRelIsoBarrel_[recoCutStep]->Fill(relIso);
      }
      else if ( fabs(recoEta) < 1.2 )
      {
        hPtOverlap_[recoCutStep]->Fill(recoPt);
        hEtaOverlap_[recoCutStep]->Fill(recoEta);
        hPhiOverlap_[recoCutStep]->Fill(recoPhi);
        hPtWithL1BinOverlap_[recoCutStep]->Fill(recoPt);
        hEtaWithL1BinOverlap_[recoCutStep]->Fill(recoEta);

        hNGlbHitOverlap_[recoCutStep]->Fill(nMuonHit);
        hNTrkHitOverlap_[recoCutStep]->Fill(nTrkHit);
        hGlbX2Overlap_[recoCutStep]->Fill(glbX2);
        hTrkX2Overlap_[recoCutStep]->Fill(trkX2);
        hRelIsoOverlap_[recoCutStep]->Fill(relIso);
      }
      else
      {
        hPtEndcap_[recoCutStep]->Fill(recoPt);
        hEtaEndcap_[recoCutStep]->Fill(recoEta);
        hPhiEndcap_[recoCutStep]->Fill(recoPhi);
        hPtWithL1BinEndcap_[recoCutStep]->Fill(recoPt);
        hEtaWithL1BinEndcap_[recoCutStep]->Fill(recoEta);

        hNGlbHitEndcap_[recoCutStep]->Fill(nMuonHit);
        hNTrkHitEndcap_[recoCutStep]->Fill(nTrkHit);
        hGlbX2Endcap_[recoCutStep]->Fill(glbX2);
        hTrkX2Endcap_[recoCutStep]->Fill(trkX2);
        hRelIsoEndcap_[recoCutStep]->Fill(relIso);
      }
    }

    // Now start Matching
    TrajectoryStateOnSurface tsos = l1Matcher_.extrapolate(*recoMuon);
    if ( !tsos.isValid() ) continue;

    const double recoPosEta = tsos.globalPosition().eta();
    const double recoPosPhi = tsos.globalPosition().phi();

    double matchedDeltaR = -999, matchedDeltaPhi = -999, matchedDeltaEta = -999;
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
      if ( matchedDeltaR < 0 or (dR < matchedDeltaR and dR < 0.3) ) 
      {
        matchedDeltaR = dR;
        matchedDeltaPhi = deltaPhi(recoPosPhi, l1PosPhi);
        matchedDeltaEta = fabs(recoEta - l1PosEta);
        bestMatchingL1Muon = *l1Muon;
      }
    }

    if ( matchedDeltaR < 0 ) continue;

    // Now we have best matching l1Muon
    for ( int recoCutStep=0; recoCutStep<nRecoMuonCutStep; ++recoCutStep )
    {
      if ( !muonQuality[recoCutStep] ) continue;

      hL1DeltaR_[recoCutStep]->Fill(matchedDeltaR);
      hL1DeltaPhi_[recoCutStep]->Fill(matchedDeltaPhi);
      hL1DeltaEta_[recoCutStep]->Fill(matchedDeltaEta);
      hMatchedL1Pt_[recoCutStep]->Fill(recoPt);
      hMatchedL1Eta_[recoCutStep]->Fill(recoEta);
      hMatchedL1Phi_[recoCutStep]->Fill(recoPhi);
      hL1MatchedGlbPt_[recoCutStep]->Fill(recoPt);
      hL1MatchedGlbEta_[recoCutStep]->Fill(recoEta);
      hL1MatchedGlbPhi_[recoCutStep]->Fill(recoPhi);

      if ( fabs(recoEta) < 0.9 )
      {
        hL1DeltaRBarrel_[recoCutStep]->Fill(matchedDeltaR);
        hL1DeltaPhiBarrel_[recoCutStep]->Fill(matchedDeltaPhi);
        hL1DeltaEtaBarrel_[recoCutStep]->Fill(matchedDeltaEta);
        hMatchedL1PtBarrel_[recoCutStep]->Fill(recoPt);
        hMatchedL1EtaBarrel_[recoCutStep]->Fill(recoEta);
        hMatchedL1PhiBarrel_[recoCutStep]->Fill(recoPhi);
        hL1MatchedGlbPtBarrel_[recoCutStep]->Fill(recoPt);
        hL1MatchedGlbEtaBarrel_[recoCutStep]->Fill(recoEta);
        hL1MatchedGlbPhiBarrel_[recoCutStep]->Fill(recoPhi);
      }
      else if ( fabs(recoEta) < 1.2 )
      {
        hL1DeltaROverlap_[recoCutStep]->Fill(matchedDeltaR);
        hL1DeltaPhiOverlap_[recoCutStep]->Fill(matchedDeltaPhi);
        hL1DeltaEtaOverlap_[recoCutStep]->Fill(matchedDeltaEta);
        hMatchedL1PtOverlap_[recoCutStep]->Fill(recoPt);
        hMatchedL1EtaOverlap_[recoCutStep]->Fill(recoEta);
        hMatchedL1PhiOverlap_[recoCutStep]->Fill(recoPhi);
        hL1MatchedGlbPtOverlap_[recoCutStep]->Fill(recoPt);
        hL1MatchedGlbEtaOverlap_[recoCutStep]->Fill(recoEta);
        hL1MatchedGlbPhiOverlap_[recoCutStep]->Fill(recoPhi);
      }
      else
      {
        hL1DeltaREndcap_[recoCutStep]->Fill(matchedDeltaR);
        hL1DeltaPhiEndcap_[recoCutStep]->Fill(matchedDeltaPhi);
        hL1DeltaEtaEndcap_[recoCutStep]->Fill(matchedDeltaEta);
        hMatchedL1PtEndcap_[recoCutStep]->Fill(recoPt);
        hMatchedL1EtaEndcap_[recoCutStep]->Fill(recoEta);
        hMatchedL1PhiEndcap_[recoCutStep]->Fill(recoPhi);
        hL1MatchedGlbPtEndcap_[recoCutStep]->Fill(recoPt);
        hL1MatchedGlbEtaEndcap_[recoCutStep]->Fill(recoEta);
        hL1MatchedGlbPhiEndcap_[recoCutStep]->Fill(recoPhi);
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

      double matchedDeltaR = -999;
      double matchedPtRes = -999;
      const trigger::Keys& trgKeys = triggerEventHandle->filterKeys(filterIdx);
      for ( trigger::Keys::const_iterator trgKey = trgKeys.begin(); 
            trgKey != trgKeys.end(); ++trgKey )
      {
        const double hltPt = triggerObjects[*trgKey].pt();
        const double hltEta = triggerObjects[*trgKey].eta();
        const double hltPhi = triggerObjects[*trgKey].phi();

        const double dR = deltaR(recoEta, hltEta, recoPhi, hltPhi);
        const double ptRes = hltPt == 0 ? 1e14 : (hltPt-recoPt)/hltPt;

        if ( matchedDeltaR < 0 or (matchedDeltaR > dR and dR < 0.5 and ptRes < 10) ) 
        {
          matchedDeltaR = dR;
          matchedPtRes = ptRes;
        }
      }

      if ( matchedDeltaR < 0 ) continue;

      // Now we have best matching candidate for HLT-global muon pair
      for ( int recoCutStep=0; recoCutStep<nRecoMuonCutStep; ++recoCutStep )
      {
        if ( !muonQuality[recoCutStep] ) continue;

        hHLTDeltaR_[recoCutStep]->Fill(matchedDeltaR);
        hHLTDeltaPhi_[recoCutStep]->Fill(matchedDeltaPhi);
        hHLTDeltaEta_[recoCutStep]->Fill(matchedDeltaEta);
        hHLTMatchedGlbPt_[recoCutStep]->Fill(recoPt);
        hHLTMatchedGlbEta_[recoCutStep]->Fill(recoEta);
        hHLTMatchedGlbPhi_[recoCutStep]->Fill(recoPhi);
        if ( fabs(recoEta) < 0.9 )
        {
          hHLTDeltaRBarrel_[recoCutStep]->Fill(matchedDeltaR);
          hHLTDeltaPhiBarrel_[recoCutStep]->Fill(matchedDeltaPhi);
          hHLTDeltaEtaBarrel_[recoCutStep]->Fill(matchedDeltaEta);
          hHLTMatchedGlbPtBarrel_[recoCutStep]->Fill(recoPt);
          hHLTMatchedGlbEtaBarrel_[recoCutStep]->Fill(recoEta);
          hHLTMatchedGlbPhiBarrel_[recoCutStep]->Fill(recoPhi);
        }
        else if ( fabs(recoEta) < 1.2 )
        {
          hHLTDeltaROverlap_[recoCutStep]->Fill(matchedDeltaR);
          hHLTDeltaPhiOverlap_[recoCutStep]->Fill(matchedDeltaPhi);
          hHLTDeltaEtaOverlap_[recoCutStep]->Fill(matchedDeltaEta);
          hHLTMatchedGlbPtOverlap_[recoCutStep]->Fill(recoPt);
          hHLTMatchedGlbEtaOverlap_[recoCutStep]->Fill(recoEta);
          hHLTMatchedGlbPhiOverlap_[recoCutStep]->Fill(recoPhi);
        }
        else
        {
          hHLTDeltaREndcap_[recoCutStep]->Fill(matchedDeltaR);
          hHLTDeltaPhiEndcap_[recoCutStep]->Fill(matchedDeltaPhi);
          hHLTDeltaEtaEndcap_[recoCutStep]->Fill(matchedDeltaEta);
          hHLTMatchedGlbPtEndcap_[recoCutStep]->Fill(recoPt);
          hHLTMatchedGlbEtaEndcap_[recoCutStep]->Fill(recoEta);
          hHLTMatchedGlbPhiEndcap_[recoCutStep]->Fill(recoPhi);
        }

        ++nHLTMatchedRecoMuon[recoCutStep];
      }
    } // Loop over HLTEvent
  } // Loop over global muons

  for ( int recoCutStep=0; recoCutStep<nRecoMuonCutStep; ++recoCutStep )
  {
    hNRecoMuon_[recoCutStep]->Fill(nRecoMuon[recoCutStep]);
    hNL1MatchedRecoMuon_[recoCutStep]->Fill(nL1MatchedRecoMuon[recoCutStep]);
    hNHLTMatchedRecoMuon_[recoCutStep]->Fill(nHLTMatchedRecoMuon[recoCutStep]);
  }
}

