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

struct Histograms
{
  Histograms(TFileDirectory& dir, TString prefix, edm::ParameterSet& cutSet, int objectType);

  virtual void FillReco(const reco::Candidate& recoCand);
  virtual void FillL1T(const reco::Candidate& recoCand, const reco::LeafCandidate& l1Cand);
  virtual void FillHLT(const reco::Candidate& recoCand, const trigger::TriggerObject& triggerObject);

  struct ObjectType
  {
    enum
    {
      Muon, Jet, Electron
    };
  };

  typedef TH1F* TH1FP;

  TH1FP hNReco, hNL1T, hNHLT;

  TH1FP hEtReco, hEtaReco, hPhiReco;
  TH1FP hEtL1T, hEtaL1T, hPhiL1T;
  TH1FP hL1EtL1T, hL1EtaL1T, hL1PhiL1T;
  TH1FP hEtHLT, hEtaHLT, hPhiHLT;
  TH1FP hHLTEtHLT, hHLTEtaHLT, hHLTPhiHLT;

  TH1FP hDeltaRL1T, hDeltaPhiL1T, hDeltaEtaL1T;
  TH1FP hDeltaRHLT, hDeltaPhiHLT, hDeltaEtaHLT;

  // For the Muons
  TH1FP hQ;
  TH1FP hGlbNHit, hTrkNHit, hGlbX2, hTrkX2;
  TH1FP hRelIso;

  void FillL1T(const reco::Muon& recoCand, const reco::LeafCandidate& l1Cand, const double recoPosEta, const double recoPosPhi);

  const int objectType_;

protected:
  double workingPointEt_, maxL1DeltaR_, maxHLTDeltaR_;
};

#endif
