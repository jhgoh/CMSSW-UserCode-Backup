#ifndef HLTrigger_TPGAnalysis_Histograms_H
#define HLTrigger_TPGAnalysis_Histograms_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>

struct Histograms
{
  Histograms(TFileDirectory& dir, TString prefix, edm::ParameterSet& cutSet, int objectType);

  void setRecoCand(const reco::Candidate* recoCand,
                   const double recoPosEta = 0, const double recoPosPhi = 0);
  void setL1Cand(const reco::LeafCandidate* l1Cand);
  void setHLTCand(const trigger::TriggerObject* hltCand);
  void init();
  void fill();

  // L1 Muon association is special, associated by by position deltaR
  void setL1MuonCand(const reco::LeafCandidate* l1Cand); 

  struct ObjectType
  {
    enum
    {
      Muon, Jet, Electron
    };
  };

  typedef TH1F* TH1FP;
  typedef TH2F* TH2FP;

  TH1FP hNReco, hNL1T, hNHLT;

  TH1FP hEtReco, hEtaReco, hPhiReco;
  TH1FP hEtL1T, hEtaL1T, hPhiL1T;
  TH1FP hL1EtL1T, hL1EtaL1T, hL1PhiL1T;
  TH1FP hEtHLT, hEtaHLT, hPhiHLT;
  TH1FP hHLTEtHLT, hHLTEtaHLT, hHLTPhiHLT;

  TH2FP hEtVsL1Et, hEtVsHLTEt;
  TH2FP hEtaVsL1Eta, hEtaVsHLTEta;
  TH2FP hPhiVsL1Phi, hPhiVsHLTPhi;

  TH1FP hDeltaRL1T, hDeltaPhiL1T, hDeltaEtaL1T;
  TH1FP hDeltaRHLT, hDeltaPhiHLT, hDeltaEtaHLT;

  TH2FP hDeltaEtaVsDeltaPhiL1T;
  TH2FP hDeltaEtaVsDeltaPhiHLT;

  // For the Muons
  TH1FP hQ;
  TH1FP hGlbNHit, hTrkNHit, hGlbX2, hTrkX2;
  TH1FP hRelIso;

  const int objectType_;

protected:
  const reco::Candidate* recoCand_;
  const reco::LeafCandidate* l1Cand_;
  const trigger::TriggerObject* hltCand_;

  double workingPointEt_, maxL1DeltaR_, maxHLTDeltaR_;
  double recoPosEta_, recoPosPhi_;
};

#endif
