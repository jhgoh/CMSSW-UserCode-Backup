#ifndef HiggsAnalysis_DoublyChargedHiggs_HiggsToEMuAnalyzer_H
#define HiggsAnalysis_DoublyChargedHiggs_HiggsToEMuAnalyzer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>

class TH1F;
class TH2F;

class HiggsToEMuAnalyzer : public edm::EDAnalyzer
{
public:
  HiggsToEMuAnalyzer(const edm::ParameterSet& pset);
  ~HiggsToEMuAnalyzer();

protected:
  virtual void beginJob(const edm::EventSetup& eventSetup);
  virtual void endJob();
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  edm::InputTag higgs1Label_;
  edm::InputTag higgs2Label_;
  
  std::map<std::string, TH1F*> h1_;
  std::map<std::string, TH2F*> h2_;
};

#endif
