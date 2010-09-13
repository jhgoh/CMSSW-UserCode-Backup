#ifndef HLTrigger_TPGAnalysis_Histograms_H
#define HLTrigger_TPGAnalysis_Histograms_H

#include <TString.h>
#include <TH1F.h>

class TFileDirectory;

namespace edm
{
  class ParameterSet;
}

namespace reco
{
  class LeafCandidate;
  class Candidate;
  class Muon;
}

namespace trigger
{
  class TriggerObject;
}

struct Histograms
{
  Histograms(TFileDirectory& dir, TString prefix, edm::ParameterSet& cutSet);

  virtual void FillReco(const reco::Candidate& recoCand);
  virtual void FillL1T(const reco::Candidate& recoCand, const reco::LeafCandidate& l1Cand);
  virtual void FillHLT(const reco::Candidate& recoCand, const trigger::TriggerObject& triggerObject);

  typedef TH1F* TH1FP;

  TH1FP hNReco, hNL1T, hNHLT;

  TH1FP hEtReco, hEtaReco, hPhiReco;
  TH1FP hEtL1T, hEtaL1T, hPhiL1T;
  TH1FP hL1EtL1T, hL1EtaL1T, hL1PhiL1T;
  TH1FP hEtHLT, hEtaHLT, hPhiHLT;
  TH1FP hHLTEtHLT, hHLTEtaHLT, hHLTPhiHLT;

  TH1FP hDeltaRL1T, hDeltaPhiL1T, hDeltaEtaL1T;
  TH1FP hDeltaRHLT, hDeltaPhiHLT, hDeltaEtaHLT;

protected:
  virtual void SetHistogramBins();

  std::vector<double> binsEt_;
  std::vector<double> binsEta_;

  double workingPointEt_, maxL1DeltaR_, maxHLTDeltaR_;
};

struct MuonHistograms : public Histograms
{
  MuonHistograms(TFileDirectory& dir, TString prefix, edm::ParameterSet& cutSet);
  void FillReco(const reco::Candidate& recoCand);
  void FillL1T(const reco::Candidate& recoCand, const reco::LeafCandidate& l1Cand);
  void FillL1T(const reco::Muon& recoCand, const reco::LeafCandidate& l1Cand, const double recoPosEta, const double recoPosPhi);

  // Muon specific histograms
  typedef TH1F* TH1FP;

  TH1FP hQ;
  TH1FP hGlbNHit, hTrkNHit, hGlbX2, hTrkX2;
  TH1FP hRelIso;

protected:
  void SetHistogramBins();
};

struct JetHistograms : public Histograms
{
  JetHistograms(TFileDirectory& dir, TString prefix, edm::ParameterSet& cutSet);

protected:
  void SetHistogramBins();
};

#endif
