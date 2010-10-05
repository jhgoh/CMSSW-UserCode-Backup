#ifndef HLTrigger_TPGAnalysis_JetHLTAnalyzer_H
#define HLTrigger_TPGAnalysis_JetHLTAnalyzer_H

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"

#include <TH1F.h>
#include <TH2F.h>

#include <map>

class Histograms;

class JetHLTAnalyzer : public edm::EDAnalyzer
{
public:
  JetHLTAnalyzer(const edm::ParameterSet& pset);
  ~JetHLTAnalyzer();

  void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup);
  void endRun();
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  template<typename T>
  const T* getBestMatch(const reco::Candidate& recoCand, std::vector<T>& targetCands);
  bool isGoodJet(const reco::CaloJet& recoJet, const edm::Event& event);

  typedef std::vector<std::string> VString;
  typedef TH1F* TH1FP;

  edm::ParameterSet jetCutSet_;
  double recoMinEt_, l1MinEt_, maxL1DeltaR_;

  std::string interestedFilterName_;

  //edm::InputTag l1tJetTag_;
  //edm::InputTag triggerEventTag_;
  edm::InputTag recoJetTag_;

  TH1FP hNReco_;
  TH1FP hNCentralL1T_, hNOverlapL1T_, hNForwardL1T_, hNDuplicatedL1T_; 
  TH1FP hNCentralHLT_, hNOverlapHLT_, hNForwardHLT_;

  std::map<int, Histograms*> hAllJet_ByRun_, hAllLeadingJet_ByRun_;
  std::map<int, Histograms*> hCentralJet_ByRun_, hCentralLeadingJet_ByRun_;
  std::map<int, Histograms*> hOverlapJet_ByRun_, hOverlapLeadingJet_ByRun_;
  std::map<int, Histograms*> hForwardJet_ByRun_, hForwardLeadingJet_ByRun_;

  Histograms* hAllJet_, * hAllLeadingJet_;
  Histograms* hCentralJet_, * hCentralLeadingJet_;
  Histograms* hOverlapJet_, * hOverlapLeadingJet_;
  Histograms* hForwardJet_, * hForwardLeadingJet_;

  Histograms* hAllJet_AllRun_, * hAllLeadingJet_AllRun_;
  Histograms* hCentralJet_AllRun_, * hCentralLeadingJet_AllRun_;
  Histograms* hOverlapJet_AllRun_, * hOverlapLeadingJet_AllRun_;
  Histograms* hForwardJet_AllRun_, * hForwardLeadingJet_AllRun_;

  Histograms* hAllJetNoL1_;
  Histograms* hCentralJetNoL1_;
  Histograms* hOverlapJetNoL1_;
  Histograms* hForwardJetNoL1_;
};

#endif

