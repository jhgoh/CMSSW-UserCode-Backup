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

#include <TH1F.h>
#include <TH2F.h>

#include <map>

class MuonHistograms;

class MuonHLTAnalyzer : public edm::EDAnalyzer
{
public:
  MuonHLTAnalyzer(const edm::ParameterSet& pset);
  ~MuonHLTAnalyzer();

  void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup);
  void endRun();
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  typedef l1extra::L1MuonParticleCollection::const_iterator L1Iter;
  typedef std::vector<trigger::TriggerObject>::const_iterator HLTIter;

  const l1extra::L1MuonParticle* getBestMatch(const reco::Candidate& recoCand, L1Iter l1Begin, L1Iter l1End);
  const trigger::TriggerObject* getBestMatch(const reco::Candidate& recoCand, HLTIter hltBegin, HLTIter hltEnd);
  bool isGoodMuon(const reco::Muon& recoMuon);

  typedef std::vector<std::string> VString;
  typedef TH1F* TH1FP;

  edm::ParameterSet muonCutSet_;

  std::string interestedFilterName_;

  //edm::InputTag l1MuonTag_;
  //edm::InputTag triggerEventTag_;
  edm::InputTag recoMuonTag_;
  
  // List of histograms
  std::map<int, MuonHistograms*> hMuon_ByRun_, hBarrelMuon_ByRun_, hOverlapMuon_ByRun_, hEndcapMuon_ByRun_;
  MuonHistograms* hMuon_, * hBarrelMuon_, * hOverlapMuon_, * hEndcapMuon_;

  L1MuonMatcherAlgo* l1Matcher_;
//  edm::ESHandle<MagneticField> bField_;
//  edm::ESHandle<GlobalTrackingGeometry> geometry_;
//  edm::ESHandle<Propagator> propagator_;

  const double maxEtaBarrel_, maxEtaOverlap_;
};

#endif

