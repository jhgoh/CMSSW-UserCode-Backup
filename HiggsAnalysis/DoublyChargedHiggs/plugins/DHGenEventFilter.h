#ifndef HiggsAnalysis_DoublyChargedHiggs_DHGenEventFilter_H
#define HiggsAnalysis_DoublyChargedHiggs_DHGenEventFilter_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class DHGenEventFilter : public edm::EDFilter
{
public:
  explicit DHGenEventFilter(const edm::ParameterSet& pset);
  ~DHGenEventFilter();

protected:
  virtual void beginJob();
  virtual bool filter(edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void endJob();

private:
  edm::InputTag genLabel_;

  int decay1_, decay2_;
  int hL_pdgId_, hR_pdgId_;
};

#endif
