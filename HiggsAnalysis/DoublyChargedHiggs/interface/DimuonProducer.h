#ifndef HiggsAnalysis_DoublyChargedHiggs_DimuonProducer_H
#define HiggsAnalysis_DoublyChargedHiggs_DimuonProducer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"

class DimuonProducer : public edm::EDProducer
{
public:
  explicit DimuonProducer(const edm::ParameterSet& pset);
  ~DimuonProducer();

protected:
  virtual void beginJob(const edm::EventSetup& eventSetup);
  virtual void produce(edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void endJob();

private:
  edm::InputTag muon1Label_, muon2Label_;
  
  bool isSameCollection_;
  CandCommonVertexFitter<KalmanVertexFitter>* kvFitter_;
};

/* vim:set ts=2 sts=2 sw=2 expandtab: */

#endif
