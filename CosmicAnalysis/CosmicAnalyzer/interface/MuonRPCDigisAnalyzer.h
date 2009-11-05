#ifndef CosmicAnalysis_CosmicAnalyzer_MuonRPCDigisAnalyzer_H
#define CosmicAnalysis_CosmicAnalyzer_MuonRPCDigisAnalyzer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"

#include <map>

class TH1F;
class TProfile2D;

class MuonRPCDigisAnalyzer : public edm::EDAnalyzer
{
public:
  MuonRPCDigisAnalyzer(const edm::ParameterSet& pset);
  ~MuonRPCDigisAnalyzer();

  virtual void beginJob(const edm::EventSetup& eventSetup);
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void endJob();

private:
  edm::InputTag digiLabel_;

  std::map<int, TH1F*> h1_;
  std::map<int, TProfile2D*> hPrf2_;
};

#endif
