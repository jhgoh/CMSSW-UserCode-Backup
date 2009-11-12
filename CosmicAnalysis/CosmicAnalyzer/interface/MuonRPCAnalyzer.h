#ifndef CosmicAnalysis_CosmicAnalyzer_MuonRPCAnalyzer_H
#define CosmicAnalysis_CosmicAnalyzer_MuonRPCAnalyzer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/RPCObjects/interface/RPCObCond.h"

#include <map>

class TH1F;

class MuonRPCAnalyzer : public edm::EDAnalyzer
{
public:
  MuonRPCAnalyzer(const edm::ParameterSet& pset);
  ~MuonRPCAnalyzer();

  virtual void beginJob(const edm::EventSetup& eventSetup);
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void endJob();

private:
  edm::Service<TFileService> fs_;

  edm::InputTag digiLabel_;

  std::map<int, TH1F*> h1_;

  std::vector<float> eventNumbers_;
  std::vector<float> rpcAvgT_, rpcAvgV_, rpcAvgI_;
  std::vector<float> rpcErrT_, rpcErrV_, rpcErrI_;
};

#endif
