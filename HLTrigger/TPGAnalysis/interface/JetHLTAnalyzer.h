#ifndef HLTrigger_TPGAnalysis_JetHLTAnalyzer_H
#define HLTrigger_TPGAnalysis_JetHLTAnalyzer_H

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
#include "RecoJets/JetAlgorithms/interface/JetIDHelper.h"

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
  typedef l1extra::L1JetParticleCollection::const_iterator L1Iter;
  typedef std::vector<trigger::TriggerObject>::const_iterator HLTIter;

  const l1extra::L1JetParticle* getBestMatch(const reco::Candidate& recoCand, L1Iter l1Begin, L1Iter l1End);
  const trigger::TriggerObject* getBestMatch(const reco::Candidate& recoCand, HLTIter hltBegin, HLTIter hltEnd);
  bool isGoodJet(const reco::CaloJet& recoJet, const edm::Event& event);

  typedef std::vector<std::string> VString;
  typedef TH1F* TH1FP;

  edm::ParameterSet jetCutSet_;
  double recoMinEt_, l1MinEt_, maxL1DeltaR_;

  std::string interestedFilterName_;

  //edm::InputTag l1tJetTag_;
  //edm::InputTag triggerEventTag_;
  edm::InputTag recoJetTag_;

  std::map<int, TH1FP> hNReco_ByRun_;
  std::map<int, TH1FP> hNCentralL1T_ByRun_, hNForwardL1T_ByRun_, hNDuplicatedL1T_ByRun_;
  std::map<int, TH1FP> hNCentralHLT_ByRun_, hNForwardHLT_ByRun_;

  TH1FP hNReco_;
  TH1FP hNCentralL1T_, hNForwardL1T_, hNDuplicatedL1T_; 
  TH1FP hNCentralHLT_, hNForwardHLT_;

  std::map<int, Histograms*> hAllJet_ByRun_, hCentralJet_ByRun_, hForwardJet_ByRun_;
  std::map<int, Histograms*> hAllLeadingJet_ByRun_, hCentralLeadingJet_ByRun_, hForwardLeadingJet_ByRun_;
  Histograms* hAllJet_, * hCentralJet_, * hForwardJet_;
  Histograms* hAllLeadingJet_, * hCentralLeadingJet_, * hForwardLeadingJet_;

  reco::helper::JetIDHelper* jetIDHelper_;

  TH2F * hCentralL1EtVsRecoEt_, * hCentralHLTEtVsRecoEt_;
};

#endif

