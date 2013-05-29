#ifndef HiggsAnalysis_DoublyChargedHiggs_PATMuonCandSelector_H
#define HiggsAnalysis_DoublyChargedHiggs_PATMuonCandSelector_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class PATMuonCandSelector : public edm::EDProducer
{
public:
  explicit PATMuonCandSelector(const edm::ParameterSet& pset);
  ~PATMuonCandSelector() {};

protected:
  virtual void beginJob(const edm::EventSetup& eventSetup) {};
  virtual void produce(edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void endJob() {};

private:
  edm::InputTag muonLabel_;

  int charge_;
  bool useGlobalMuonsOnly_;
  double minPt_;
};

/* vim:set ts=2 sts=2 sw=2 expandtab: */

#endif
