#ifndef HLTrigger_TPGAnalysis_MuonHLTAnalyzer_H
#define HLTrigger_TPGAnalysis_MuonHLTAnalyzer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "MuonAnalysis/MuonAssociators/interface/L1MuonMatcherAlgo.h"
#include "DataFormats/Math/interface/Point3D.h"

#include <TH1F.h>
#include <TH2F.h>

#include <map>

class Histograms;

class MuonHLTAnalyzer : public edm::EDAnalyzer
{
public:
  MuonHLTAnalyzer(const edm::ParameterSet& pset);
  ~MuonHLTAnalyzer();

  void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup);
  void endRun();
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  typedef l1extra::L1MuonParticleCollection L1Muons;
  typedef std::vector<trigger::TriggerObject> TriggerObjects;
  typedef L1Muons::const_iterator L1Iter;
  typedef TriggerObjects::const_iterator HLTIter;

  const l1extra::L1MuonParticle* getBestMatch(const double recoPosEta, const double recoPosPhi, L1Muons& l1Particles);
  const trigger::TriggerObject* getBestMatch(const reco::Candidate& recoCand, TriggerObjects& hltObjects);
  bool isGoodMuon(const reco::Muon& recoMuon);

  typedef std::vector<std::string> VString;
  typedef TH1F* TH1FP;

  edm::ParameterSet muonCutSet_;
  double recoMinPt_, l1MinPt_, maxL1DeltaR_;
  unsigned int l1MinQuality_;

  std::string interestedFilterName_;

  //edm::InputTag l1MuonTag_;
  //edm::InputTag triggerEventTag_;
  edm::InputTag recoMuonTag_;
  
  // List of histograms
  std::map<int, Histograms*> hAllMuon_ByRun_, hBarrelMuon_ByRun_, hOverlapMuon_ByRun_, hEndcapMuon_ByRun_;
  std::map<int, Histograms*> hAllLeadingMuon_ByRun_, hBarrelLeadingMuon_ByRun_, hOverlapLeadingMuon_ByRun_, hEndcapLeadingMuon_ByRun_;
  Histograms* hAllMuon_, * hBarrelMuon_, * hOverlapMuon_, * hEndcapMuon_;
  Histograms* hAllLeadingMuon_, * hBarrelLeadingMuon_, * hOverlapLeadingMuon_, * hEndcapLeadingMuon_;

  L1MuonMatcherAlgo* l1Matcher_;
//  edm::ESHandle<MagneticField> bField_;
//  edm::ESHandle<GlobalTrackingGeometry> geometry_;
//  edm::ESHandle<Propagator> propagator_;

  const double maxEtaBarrel_, maxEtaOverlap_;

  math::XYZPoint beamSpotPosition_;
};

#endif

