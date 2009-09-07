#ifndef HiggsAnalysis_DoublyChargedHiggs_DileptonProducer_H
#define HiggsAnalysis_DoublyChargedHiggs_DileptonProducer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"

class DileptonProducer : public edm::EDProducer
{
public:
  explicit DileptonProducer(const edm::ParameterSet& pset);
  ~DileptonProducer();

protected:
  virtual void beginJob(const edm::EventSetup& eventSetup);
  virtual void produce(edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void endJob();

  template<typename LeptonIter, typename OutCollection> 
  void combineLeptons(const LeptonIter lepton_begin, const LeptonIter lepton_end,
                      OutCollection& dileptonCands);
  template <typename LeptonIter, typename OutCollection>
  void combineLeptons(const LeptonIter lepton1_begin, const LeptonIter lepton1_end,
                      const LeptonIter lepton2_begin, const LeptonIter lepton2_end,
                      OutCollection& dileptonCands);
  template <typename LeptonIter1, typename LeptonIter2, typename OutCollection>
  void combineLeptons(const LeptonIter1 lepton1_begin, const LeptonIter1 lepton1_end,
                      const LeptonIter2 lepton2_begin, const LeptonIter2 lepton2_end,
                      OutCollection& dileptonCands);

private:
  edm::InputTag lepton1Label_, lepton2Label_;
  
  int lepton1Type_, lepton2Type_;
  
  bool isSameCollection_;
  CandCommonVertexFitter<KalmanVertexFitter>* kvFitter_;
};

/* vim:set ts=2 sts=2 sw=2 expandtab: */

#endif
