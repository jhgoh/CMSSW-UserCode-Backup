#ifndef HiggsAnalysis_DoublyChargedHiggs_MultipleCandCounterFilter_H
#define HiggsAnalysis_DoublyChargedHiggs_MultipleCandCounterFilter_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class MultipleCandCounterFilter : public edm::EDFilter
{
public:
  explicit MultipleCandCounterFilter(const edm::ParameterSet& pset);
  ~MultipleCandCounterFilter();

protected:
  virtual bool filter(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  std::vector<edm::InputTag> candTags_;
  std::vector<double> ptThresholds_;
};

#endif

