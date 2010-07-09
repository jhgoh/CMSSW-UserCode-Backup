#ifndef HLTrigger_TPGAnalysis_MuonHLTAnalyzer_H
#define HLTrigger_TPGAnalysis_MuonHLTAnalyzer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MuonAnalysis/MuonAssociators/interface/L1MuonMatcherAlgo.h"

#include <TH1F.h>
#include <TH2F.h>

#include <vector>

namespace reco
{
  class Muon;
}

class TrajectoryStateOnSurface;
class MagneticField;
class Propagator;
class GlobalTrackingGeometry;

class MuonHLTAnalyzer : public edm::EDAnalyzer
{
public:
  MuonHLTAnalyzer(const edm::ParameterSet& pset);
  ~MuonHLTAnalyzer();

  void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup);
  void endRun();
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  typedef std::vector<std::string> VString;

  edm::InputTag l1MuonTag_, recoMuonTag_;
  VString muonTrigNames_;
  
  L1MuonMatcherAlgo l1Matcher_;
  edm::ESHandle<MagneticField> bField_;
  edm::ESHandle<GlobalTrackingGeometry> geometry_;
  edm::ESHandle<Propagator> propagator_;

  // Cut variables
  double minPt_, maxRelIso_;

  // List of histograms
  TH1F* hNEvent_;
  TH1F* hNRecoMuon_;

  // List of histograms : track varialbes for each cut steps
  std::vector<TH1F*> hPt_, hEta_, hPhi_, hQ_;
  std::vector<TH1F*> hNGlbHit_, hGlbX2_, hNTrkHit_, hTrkX2_;
  std::vector<TH1F*> hRelIso_;
  std::vector<TH1F*> hDeltaR_, hDeltaPhi_, hDeltaEta_;
  std::vector<TH1F*> hL1Pt_, hL1Eta_, hL1Phi_;
  std::vector<TH1F*> hMatchedL1Pt_, hMatchedL1Eta_, hMatchedL1Phi_;
  std::vector<TH1F*> hL1MatchedGlbPt_, hL1MatchedGlbEta_, hL1MatchedGlbPhi_;

  std::vector<TH1F*> hPtBarrel_, hEtaBarrel_, hPhiBarrel_, hQBarrel_;
  std::vector<TH1F*> hNGlbHitBarrel_, hGlbX2Barrel_, hNTrkHitBarrel_, hTrkX2Barrel_;
  std::vector<TH1F*> hRelIsoBarrel_;

  std::vector<TH1F*> hPtOverlap_, hEtaOverlap_, hPhiOverlap_, hQOverlap_;
  std::vector<TH1F*> hNGlbHitOverlap_, hGlbX2Overlap_, hNTrkHitOverlap_, hTrkX2Overlap_;
  std::vector<TH1F*> hRelIsoOverlap_;

  std::vector<TH1F*> hPtEndcap_, hEtaEndcap_, hPhiEndcap_, hQEndcap_;
  std::vector<TH1F*> hNGlbHitEndcap_, hGlbX2Endcap_, hNTrkHitEndcap_, hTrkX2Endcap_;
  std::vector<TH1F*> hRelIsoEndcap_;
};

#endif

