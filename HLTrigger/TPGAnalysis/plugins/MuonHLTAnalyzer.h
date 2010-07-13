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

  edm::InputTag l1MuonTag_, triggerEventTag_, recoMuonTag_;
  VString muonL1TNames_;
  
  L1MuonMatcherAlgo l1Matcher_;
  edm::ESHandle<MagneticField> bField_;
  edm::ESHandle<GlobalTrackingGeometry> geometry_;
  edm::ESHandle<Propagator> propagator_;

  // Cut variables
  double minPt_, maxRelIso_;

  // List of histograms
  TH1F* hNEvent_;
  std::vector<TH1F*> hNRecoMuon_;
  std::vector<TH1F*> hNL1MatchedRecoMuon_;
  std::vector<TH1F*> hNHLTMatchedRecoMuon_;

  // List of histograms : track varialbes for each cut steps
  std::vector<TH1F*> hPt_, hEta_, hPhi_, hQ_;
  std::vector<TH1F*> hPtWithL1Bin_, hEtaWithL1Bin_;
  std::vector<TH1F*> hNGlbHit_, hGlbX2_, hNTrkHit_, hTrkX2_;
  std::vector<TH1F*> hRelIso_;
  std::vector<TH1F*> hL1DeltaR_, hL1DeltaPhi_, hL1DeltaEta_;
  std::vector<TH1F*> hL1Pt_, hL1Eta_, hL1Phi_;
  std::vector<TH1F*> hMatchedL1Pt_, hMatchedL1Eta_, hMatchedL1Phi_;
  std::vector<TH1F*> hHLTDeltaR_, hHLTDeltaPhi_, hHLTDeltaEta_;
  std::vector<TH1F*> hL1MatchedGlbPt_, hL1MatchedGlbEta_, hL1MatchedGlbPhi_;
  std::vector<TH1F*> hHLTMatchedGlbPt_, hHLTMatchedGlbEta_, hHLTMatchedGlbPhi_;

  std::vector<TH1F*> hPtBarrel_, hEtaBarrel_, hPhiBarrel_, hQBarrel_;
  std::vector<TH1F*> hPtWithL1BinBarrel_, hEtaWithL1BinBarrel_;
  std::vector<TH1F*> hNGlbHitBarrel_, hGlbX2Barrel_, hNTrkHitBarrel_, hTrkX2Barrel_;
  std::vector<TH1F*> hRelIsoBarrel_;
  std::vector<TH1F*> hL1DeltaRBarrel_, hL1DeltaPhiBarrel_, hL1DeltaEtaBarrel_;
  std::vector<TH1F*> hL1PtBarrel_, hL1EtaBarrel_, hL1PhiBarrel_;
  std::vector<TH1F*> hMatchedL1PtBarrel_, hMatchedL1EtaBarrel_, hMatchedL1PhiBarrel_;
  std::vector<TH1F*> hL1MatchedGlbPtBarrel_, hL1MatchedGlbEtaBarrel_, hL1MatchedGlbPhiBarrel_;
  std::vector<TH1F*> hHLTDeltaRBarrel_, hHLTDeltaPhiBarrel_, hHLTDeltaEtaBarrel_;
  std::vector<TH1F*> hHLTMatchedGlbPtBarrel_, hHLTMatchedGlbEtaBarrel_, hHLTMatchedGlbPhiBarrel_;

  std::vector<TH1F*> hPtOverlap_, hEtaOverlap_, hPhiOverlap_, hQOverlap_;
  std::vector<TH1F*> hPtWithL1BinOverlap_, hEtaWithL1BinOverlap_;
  std::vector<TH1F*> hNGlbHitOverlap_, hGlbX2Overlap_, hNTrkHitOverlap_, hTrkX2Overlap_;
  std::vector<TH1F*> hRelIsoOverlap_;
  std::vector<TH1F*> hL1DeltaROverlap_, hL1DeltaPhiOverlap_, hL1DeltaEtaOverlap_;
  std::vector<TH1F*> hL1PtOverlap_, hL1EtaOverlap_, hL1PhiOverlap_;
  std::vector<TH1F*> hMatchedL1PtOverlap_, hMatchedL1EtaOverlap_, hMatchedL1PhiOverlap_;
  std::vector<TH1F*> hHLTDeltaROverlap_, hHLTDeltaPhiOverlap_, hHLTDeltaEtaOverlap_;
  std::vector<TH1F*> hL1MatchedGlbPtOverlap_, hL1MatchedGlbEtaOverlap_, hL1MatchedGlbPhiOverlap_;
  std::vector<TH1F*> hHLTMatchedGlbPtOverlap_, hHLTMatchedGlbEtaOverlap_, hHLTMatchedGlbPhiOverlap_;

  std::vector<TH1F*> hPtEndcap_, hEtaEndcap_, hPhiEndcap_, hQEndcap_;
  std::vector<TH1F*> hPtWithL1BinEndcap_, hEtaWithL1BinEndcap_;
  std::vector<TH1F*> hNGlbHitEndcap_, hGlbX2Endcap_, hNTrkHitEndcap_, hTrkX2Endcap_;
  std::vector<TH1F*> hRelIsoEndcap_;
  std::vector<TH1F*> hL1DeltaREndcap_, hL1DeltaPhiEndcap_, hL1DeltaEtaEndcap_;
  std::vector<TH1F*> hL1PtEndcap_, hL1EtaEndcap_, hL1PhiEndcap_;
  std::vector<TH1F*> hMatchedL1PtEndcap_, hMatchedL1EtaEndcap_, hMatchedL1PhiEndcap_;
  std::vector<TH1F*> hHLTDeltaREndcap_, hHLTDeltaPhiEndcap_, hHLTDeltaEtaEndcap_;
  std::vector<TH1F*> hL1MatchedGlbPtEndcap_, hL1MatchedGlbEtaEndcap_, hL1MatchedGlbPhiEndcap_;
  std::vector<TH1F*> hHLTMatchedGlbPtEndcap_, hHLTMatchedGlbEtaEndcap_, hHLTMatchedGlbPhiEndcap_;
};

#endif

