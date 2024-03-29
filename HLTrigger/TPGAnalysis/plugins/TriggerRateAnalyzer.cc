#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include <TString.h>
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
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  std::vector<std::string> myTrigNames_;

  int runNumber_;
  
  TH1F* hNEvent_;
  std::map<int, TH1F*> hNEventInRun_;
};

TriggerRateAnalyzer::TriggerRateAnalyzer(const edm::ParameterSet& pset)
{
  myTrigNames_ = pset.getParameter<std::vector<std::string> >("myTrigNames");

  edm::Service<TFileService> fs;

  const int nMyTrigNames = myTrigNames_.size();
  hNEvent_ = fs->make<TH1F>("hNEvent", "Number of events passing trigger paths;Trigger path", nMyTrigNames+1, 0, nMyTrigNames+1);
  hNEvent_->GetXaxis()->SetBinLabel(1, "All");
  for ( int i=0; i<nMyTrigNames; ++i )
  {
    hNEvent_->GetXaxis()->SetBinLabel(i+2, myTrigNames_[i].c_str());
  }
}

TriggerRateAnalyzer::~TriggerRateAnalyzer()
{
}

void TriggerRateAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  runNumber_ = run.run();

  if ( hNEventInRun_.find(runNumber_) == hNEventInRun_.end() )
  {
    edm::Service<TFileService> fs;
    TFileDirectory dir = fs->mkdir(Form("Run %d", runNumber_));

    const int nMyTrigNames = myTrigNames_.size();
    hNEventInRun_[runNumber_] = dir.make<TH1F>("hNEvent", "Number of events passing trigger paths;Trigger path", nMyTrigNames+1, 0, nMyTrigNames+1);
    hNEventInRun_[runNumber_]->GetXaxis()->SetBinLabel(1, "All");
    for ( int i=0; i<nMyTrigNames; ++i )
    {
      hNEventInRun_[runNumber_]->GetXaxis()->SetBinLabel(i+2, myTrigNames_[i].c_str());
    }
  }
}

void TriggerRateAnalyzer::endRun()
{
}

void TriggerRateAnalyzer::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup)
{
}

void TriggerRateAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::TriggerResults> trigResultHandle;
  if ( !event.getByLabel(edm::InputTag("TriggerResults", "", "HLT"), trigResultHandle) )
  {
    edm::LogError("TriggerRateAnalyzer") << "Cnanot find TriggerResults\n";
    return;
  }
  const edm::TriggerResults* trigResult = trigResultHandle.product();

  if ( trigResult->wasrun() and trigResult->accept() )
  {
    hNEvent_->Fill(0);
    hNEventInRun_[runNumber_]->Fill(0);

    const int nTrigResult = trigResult->size();
    const edm::TriggerNames& trigNames = event.triggerNames(*trigResult);

    for ( int trigIndex=0; trigIndex<nTrigResult; ++trigIndex )
    {
      if ( !trigResult->accept(trigIndex) ) continue;

      const std::string triggerName = trigNames.triggerName(trigIndex);

      // Find Trigger name in the triggerResult
      int bin = -1;
      const int nMyTrigNames = myTrigNames_.size();
      for ( int myTrigIndex=0; myTrigIndex<nMyTrigNames; ++myTrigIndex )
      {
        if ( myTrigNames_[myTrigIndex] == triggerName )
        {
          bin = myTrigIndex;
          break;
        }
      }

      if ( bin != -1 )
      {
        hNEvent_->Fill(bin+1);
        hNEventInRun_[runNumber_]->Fill(bin+1);
      }
    }
  } // Loop over trigger results
}

DEFINE_FWK_MODULE(TriggerRateAnalyzer);
