#ifndef HLTrigger_TPGAnalysis_TriggerRateAnalyzer_H
#define HLTrigger_TPGAnalysis_TriggerRateAnalyser_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include <TH1F.h>

#include <string>
#include <vector>
#include <map>

class TriggerRateAnalyzer : public edm::EDAnalyzer
{
public:
  TriggerRateAnalyzer(const edm::ParameterSet& pset);
  ~TriggerRateAnalyzer();

  void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup);
  void endRun();
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  std::vector<std::string> myTrigNames_;

  int runNumber_;
  
  TH1F* hNEvent_;
  std::map<int, TH1F*> hNEventInRun_;
};

#endif

