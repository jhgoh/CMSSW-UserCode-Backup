#ifndef Hpp2MuHepMCAnalyzer_H
#define Hpp2MuHepMCAnalyzer_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"

#include <string>

class HPtcl;
class HTT;

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
  typedef HPtcl* HPtclP;
  typedef HTT* HTTP;
  HPtclP hTrk_, hMu_, hGoodMu_;
  HPtclP hHpp_;
  HPtclP hHppMu_;
  HPtclP hHppGoodMuP_, hHppGoodMuM_;

  HTTP hDimuonPP_, hDimuonMM_;
  HTTP hGoodDimuonPP_, hGoodDimuonMM_;

  H1P hNMuP_, hNMuM_, hNGoodMuP_, hNGoodMuM_;
  H1P hNHiggsMu_, hNHiggsGoodMu_;
};


#endif

/* vim:set ts=2 sts=2 sw=2 expandtab: */
