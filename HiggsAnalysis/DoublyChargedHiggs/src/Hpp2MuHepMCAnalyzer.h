#ifndef Hpp2MuHepMCAnalyzer_H
#define Hpp2MuHepMCAnalyzer_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

class TH1D;

class Hpp2MuHepMCAnalyzer : public edm::EDAnalyzer
{
 public:
  Hpp2MuHepMCAnalyzer(const edm::ParameterSet& pset);
  ~Hpp2MuHepMCAnalyzer();

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);
  void beginJob(const edm::EventSetup& eventSetup);
  void endJob();

 private:
  typedef TH1D* H1P;
  H1P hMuPt_, hMuEta_;
  H1P hMuMuM_;
};

#endif

/* vim:set ts=2 sts=2 sw=2 expandtab: */
