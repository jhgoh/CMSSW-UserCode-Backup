#ifndef HiggsAnalysis_DoublyChargedHiggs_DileptonProducer_H
#define HiggsAnalysis_DoublyChargedHiggs_DileptonProducer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class DileptonProducer : public edm::EDProducer
{
public:
  explicit DileptonProducer(const edm::ParameterSet& pset);
  ~DileptonProducer();

protected:
  virtual void beginJob();
  virtual void produce(edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void endJob();

  template<typename LeptonIter, typename OutCollection> 
  void combineLeptons(const LeptonIter begin, const LeptonIter end,
                      OutCollection& dileptonCands);
  template <typename LeptonIter1, typename LeptonIter2, typename OutCollection>
  void combineLeptons(const LeptonIter1 begin1, const LeptonIter1 end1,
                      const LeptonIter2 begin2, const LeptonIter2 end2,
                      OutCollection& dileptonCands);
private:
  // Lepton1 source sets
  edm::InputTag lepton1Label_;
  int lepton1Type_, lepton1Charge_;

  // Lepton2 source sets
  edm::InputTag lepton2Label_;
  int lepton2Type_, lepton2Charge_;
  
  // Candidate combination
  bool chargeConj_, isSameCollection_;
  int dileptonCharge_;

  // Cut values
  double lepton1MinPt_, lepton1MaxEta_;
  double lepton2MinPt_, lepton2MaxEta_;
  double dileptonMinMass_, dileptonMinPt_, dileptonMaxEta_;
};

/* vim:set ts=2 sts=2 sw=2 expandtab: */

#endif
