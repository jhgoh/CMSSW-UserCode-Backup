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
class TString;

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

  std::map<std::string, TH1F*> h1_;

  // Store detector cell names and its average of IOV values
  typedef std::map<std::string, std::pair<unsigned int, double> > DetIOVMap;
  DetIOVMap rpcIValues_, rpcVValues_, rpcTValues_;
};

#endif
