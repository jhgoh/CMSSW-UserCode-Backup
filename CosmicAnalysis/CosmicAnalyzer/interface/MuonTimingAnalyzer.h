#ifndef CosmicAnalysis_CosmicAnalyzer_MuonTimingAnalyzer_H
#define CosmicAnalysis_CosmicAnalyzer_MuonTimingAnalyzer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TString.h"
#include <string>
#include <map>

class TH1F;
class TH2F;
class TProfile;

class MuonTimingAnalyzer : public edm::EDAnalyzer
{
public:
  MuonTimingAnalyzer(const edm::ParameterSet& pset);
  ~MuonTimingAnalyzer();

  virtual void beginJob(const edm::EventSetup& eventSetup);
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void endJob();

private:
  edm::Service<TFileService> fs_;

  edm::InputTag digiLabel_;
  edm::InputTag recHitLabel_;

  std::map<std::string, TH1F*> h1_;
  std::map<std::string, TH2F*> h2_;
  std::map<std::string, TProfile*> prf_;

  int minNDigiCut_;
  int minBxNumber_, maxBxNumber_;

  // Store detector cell names and its average of IOV values
  //typedef std::map<std::string, std::pair<unsigned int, double> > DetIOVMap;
  //DetIOVMap rpcIValues_, rpcVValues_, rpcTValues_;

  //uint64_t rpcIMinTime_, rpcVMinTime_, rpcTMinTime_;
  //uint64_t rpcIMaxTime_, rpcVMaxTime_, rpcTMaxTime_;

  //uint64_t timeOffset_;
};

#endif
