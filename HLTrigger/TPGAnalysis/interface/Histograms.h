#ifndef HLTrigger_TPGAnalysis_Histograms_H
#define HLTrigger_TPGAnalysis_Histograms_H

#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>

#include <string>

struct Histograms;
class HTrigger;

class HTrigger
{
public:
  HTrigger(TString subDir, const std::string prefix_,
           const double workingPointEt, const double maxL1DeltaR, const double maxHLTDeltaR,
           int objectType);
  void init(const edm::EventID& eventID);
  void fill(const reco::Candidate* recoCand, const reco::LeafCandidate* l1Cand, const trigger::TriggerObject* hltCand,
            const double recoPosEta = 0, const double recoPosPhi = 0);
  void fillNCand(const int nReco, const int nL1, const int nHLT);

private:
  // Histogram collections
  Histograms* hist_; // All run, menu merged histogram
  std::map<int, Histograms*> runToHistMap_; // Run-by-run histograms
  //std::map<int, Histograms*> prescaleToHistMap_; // Prescale-by-prescale histograms

  TString subDir_;
  std::string prefix_;
  const double workingPointEt_, maxL1DeltaR_, maxHLTDeltaR_;
  int objectType_;

  // Cashed histograms
  int runNumber_;
  Histograms* runHist_;
};

struct Histograms
{
  Histograms(TString dirName, TString prefix,
             const double workingPointEt, const double maxL1DeltaR, const double maxHLTDeltaR,
             int objectType);
  void fill(const reco::Candidate* recoCand, const reco::LeafCandidate* l1Cand, const trigger::TriggerObject* hltCand,
            const double recoPosEta = 0, const double recoPosPhi = 0);
  void fillNCand(const int nReco, const int nL1, const int nHLT);

  struct ObjectType
  {
    enum
    {
      Muon, Jet, Electron
    };
  };

  typedef TH1F* TH1FP;
  typedef TH2F* TH2FP;

  TH1FP hNReco, hNL1, hNHLT;

  TH1FP hEtReco, hEtaReco, hPhiReco;
  TH1FP hEtL1, hEtaL1, hPhiL1;
  TH1FP hL1EtL1, hL1EtaL1, hL1PhiL1;
  TH1FP hEtHLT, hEtaHLT, hPhiHLT;
  TH1FP hHLTEtHLT, hHLTEtaHLT, hHLTPhiHLT;

  TH2FP hEtVsL1Et, hEtVsHLTEt;
  TH2FP hEtaVsL1Eta, hEtaVsHLTEta;
  TH2FP hPhiVsL1Phi, hPhiVsHLTPhi;

  TH1FP hDeltaRL1, hDeltaPhiL1, hDeltaEtaL1;
  TH1FP hDeltaRHLT, hDeltaPhiHLT, hDeltaEtaHLT;

  TH2FP hDeltaEtaVsDeltaPhiL1;
  TH2FP hDeltaEtaVsDeltaPhiHLT;

  // For the Muons
  TH1FP hQ;
  TH1FP hGlbNHit, hTrkNHit, hGlbX2, hTrkX2;
  TH1FP hRelIso;

  const int objectType_;

protected:
  const double workingPointEt_, maxL1DeltaR_, maxHLTDeltaR_;

};

#endif
