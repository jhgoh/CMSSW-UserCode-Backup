#ifndef Hpp2MuHepMCAnalyzer_H
#define Hpp2MuHepMCAnalyzer_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"

#include <string>

struct HPtcl
{
public:
  HPtcl(const std::string namePrefix, const std::string titleSuffix)
  {
    edm::Service<TFileService> fs;
    hPt_ = hEta_ = 0;

    hPt_ = fs->make<TH1D>(("h"+namePrefix+"Pt").c_str(), ("p_{T} of "+titleSuffix).c_str(), 50, 0, 200);
    hEta_ = fs->make<TH1D>(("h"+namePrefix+"Eta").c_str(), ("#eta of "+titleSuffix).c_str(), 50, -2.5, 2.5);
  };

  void operator()(HepMC::GenParticle* ptcl)
  {
    if ( !hPt_ || !hEta_ ) return;
    hPt_->Fill(ptcl->momentum().perp());
    hEta_->Fill(ptcl->momentum().eta());
  };

private:
  TH1D * hPt_, * hEta_;

};

struct HTT
{
public:
  HTT(const std::string namePrefix, const std::string titleSuffix) {
    edm::Service<TFileService> fs;
    hM_ = hPt_ = hEta_ = 0;    

    hM_ = fs->make<TH1D>(("h"+namePrefix+"M").c_str(), ("Mass of "+titleSuffix).c_str(), 50, 50, 200);
    hPt_ = fs->make<TH1D>(("h"+namePrefix+"Pt").c_str(), ("p_{T} of "+titleSuffix).c_str(), 50, 0, 200);
    hEta_ = fs->make<TH1D>(("h"+namePrefix+"Eta").c_str(), ("#eta of "+titleSuffix).c_str(), 50, -2.5, 2.5);
  };

  typedef std::pair<HepMC::GenParticle*, HepMC::GenParticle*> PtclPair;
  void operator()(PtclPair pPair)
  {
    if ( !hM_ || !hPt_ || !hEta_ ) return;

    const HepMC::FourVector p1 = pPair.first->momentum();
    const HepMC::FourVector p2 = pPair.second->momentum();
    
    HepMC::FourVector p(p1.px()+p2.px(), p1.py()+p2.py(), p1.pz()+p2.pz(), p1.e()+p2.e());

    hM_->Fill(p.m());
    hPt_->Fill(p.perp());
    hEta_->Fill(p.eta());
  };

private:
  TH1D * hM_, * hPt_, * hEta_;
};

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
  HPtcl hTrk_, hMu_, hGoodMu_;
  HPtcl hHpp_;
  HPtcl hHppMu_;
  HPtcl hHppGoodMuP_, hHppGoodMuM_;

  HTT hDimuonPP_, hDimuonMM_;
  HTT hGoodDimuonPP_, hGoodDimuonMM_;

  H1P hNMuP_, hNMuM_, hNGoodMuP_, hNGoodMuM_;
  H1P hNHiggsMu_, hNHiggsGoodMu_;
};


#endif

/* vim:set ts=2 sts=2 sw=2 expandtab: */
