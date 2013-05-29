#ifndef HiggsAnalysis_DoublyChargedHiggs_FourMuonAnalyzer_H
#define HiggsAnalysis_DoublyChargedHiggs_FourMuonAnalyzer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TH1F;
class TH2F;

class FourMuonAnalyzer : public edm::EDAnalyzer
{
public:
  explicit FourMuonAnalyzer(const edm::ParameterSet& pset);
  ~FourMuonAnalyzer();

protected:
  virtual void beginJob(const edm::EventSetup& eventSetup);
  virtual void endJob();
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  edm::InputTag posDeltaLabel_, negDeltaLabel_;
  bool fsStatus_;

  double delta_minPt_;
  double delta_maxNormalizedChi2_;
  unsigned int nInterested_;

  std::map<int, TH1F*> h1_;
  std::map<int, TH2F*> h2_;
};

/* vim:set ts=2 sts=2 sw=2 expandtab: */

#endif
