#include "HiggsAnalysis/DoublyChargedHiggs/src/Hpp2MuHepMCAnalyzer.h"

#include "DataFormats/Common/interface/Handle.h"

#include "TMath.h"

#include <iostream>

using namespace std;
using namespace edm;

class MuonFilter
{
public:
  const static int Minus = -1;
  const static int Plus = 1;

  MuonFilter(const double minPt, const double maxEta, const int charge=0):
    minPt_(minPt), maxEta_(maxEta), charge_(charge)
  {
  };
  // charge = 0 for charge conjugation
  bool operator()(const HepMC::GenParticle* ptcl)
  {
    if ( ptcl->status() != 1 ) return false;

    const int pdg_id = ptcl->pdg_id();
    if ( abs(pdg_id) != 13 ) return false;

    const int charge = pdg_id/-13;
    if ( charge*charge_ == -1 ) return false;

    const double pt = ptcl->momentum().perp();
    const double eta = ptcl->momentum().eta();

    if ( pt < minPt_ || fabs(eta) > maxEta_ ) return false;

    return true;
  };

  const double minPt_, maxEta_;
  const int charge_;
};



Hpp2MuHepMCAnalyzer::Hpp2MuHepMCAnalyzer(const ParameterSet& pset):
  hTrk_("Trk", "all tracks"),
  hMu_("Mu", "all muons"),
  hGoodMu_("GoodMu", "good muons"),
  hHpp_("Hpp", "Higgs"),
  hHppMu_("HppMu", "muons from Higgs decay"),
  hHppGoodMuP_("HppGoodMuP", "good #mu^{+} from Higgs decay"),
  hHppGoodMuM_("HppGoodMuM", "good #mu^{-} from Higgs decay"),
  hDimuonPP_("DimuonPP", "#mu^{+}#mu^{+}"),
  hDimuonMM_("DimuonMM", "#mu^{-}#mu^{-}"),
  hGoodDimuonPP_("GoodDimuonPP", "Good #mu^{+}#mu^{+}"),
  hGoodDimuonMM_("GoodDimuonMM", "Good #mu^{-}#mu^{-}")
{
  edm::Service<TFileService> fs;
  
  hNMuP_ = fs->make<TH1D>("hNMuP", "# of #mu^{+}", 10, 0, 10);
  hNMuM_ = fs->make<TH1D>("hNMuM", "# of #mu^{-}", 10, 0, 10);
  hNGoodMuP_ = fs->make<TH1D>("hNGoodMuP", "# of good #mu^{+}", 10, 0, 10);
  hNGoodMuM_ = fs->make<TH1D>("hNGoodMuM", "# of good #mu^{-}", 10, 0, 10);
  hNHiggsMu_ = fs->make<TH1D>("hNHiggsMu", "# of #mu from Higgs", 10, 0, 10);
  hNHiggsGoodMu_ = fs->make<TH1D>("hNHiggsGoodMu", "# of Good #mu from Higgs", 10, 0, 10);
}

Hpp2MuHepMCAnalyzer::~Hpp2MuHepMCAnalyzer()
{
}

void Hpp2MuHepMCAnalyzer::beginJob(const EventSetup& eventSetup)
{
}

void Hpp2MuHepMCAnalyzer::analyze(const Event& event, const EventSetup& eventSetup)
{
  Handle<HepMCProduct> genEvtHandle;
  event.getByLabel("source", genEvtHandle);
  const HepMC::GenEvent* genEvt = genEvtHandle->GetEvent();

  vector<HepMC::GenParticle*> muPs;
  vector<HepMC::GenParticle*> muMs;

  MuonFilter isMuon(0, 999), isMuM(0, 999, MuonFilter::Minus), isMuP(0, 999, MuonFilter::Plus);
  MuonFilter isGoodMuon(20, 2.0), isGoodMuP(20, 2.0, MuonFilter::Plus), isGoodMuM(20, 2.0, MuonFilter::Minus);

  // Fill distributions for all tracks
  int nGoodMuP = 0, nGoodMuM = 0;
  for(HepMC::GenEvent::particle_const_iterator iGenPtcl = genEvt->particles_begin();
      iGenPtcl != genEvt->particles_end(); ++iGenPtcl) {
    HepMC::GenParticle* genPtcl = *iGenPtcl;

    if ( genPtcl->status() != 1 ) continue;

    hTrk_(genPtcl);

    if ( isMuP(genPtcl) ) muPs.push_back(genPtcl);
    else if ( isMuM(genPtcl) ) muMs.push_back(genPtcl);

    if ( isMuon(genPtcl) ) hMu_(genPtcl);
    if ( isGoodMuon(genPtcl) ) hGoodMu_(genPtcl);

    if ( isGoodMuM(genPtcl) ) ++nGoodMuM;
    if ( isGoodMuP(genPtcl) ) ++nGoodMuP;
  }
  hNMuP_->Fill(muPs.size());
  hNMuM_->Fill(muMs.size());
  hNGoodMuP_->Fill(nGoodMuP);
  hNGoodMuM_->Fill(nGoodMuM);

  vector<HepMC::GenVertex*> higgsVtxs;

  for(HepMC::GenEvent::particle_const_iterator iGenPtcl = genEvt->particles_begin();
      iGenPtcl != genEvt->particles_end(); ++iGenPtcl) {
    HepMC::GenParticle* genPtcl = *iGenPtcl;

    if ( abs(genPtcl->pdg_id()) == 9900041 ||
         abs(genPtcl->pdg_id()) == 9900042 ) {
      HepMC::GenVertex* higgsDecVtx = genPtcl->end_vertex();

      higgsVtxs.push_back(higgsDecVtx);

      hHpp_(genPtcl);
    }
  }

  if ( higgsVtxs.empty() ) {
    LogError("Hpp2MuHepMCAnalyzer") << "No Higgs in this event\n";
    return;
  }

  for(vector<HepMC::GenVertex*>::const_iterator iGenVtx = higgsVtxs.begin();
      iGenVtx != higgsVtxs.end(); ++iGenVtx) {
    HepMC::GenVertex* genVtx = *iGenVtx;
    if ( genVtx == 0 ) continue;

    vector<HepMC::GenParticle*> higgsMuPs, higgsMuMs;
    vector<HepMC::GenParticle*> higgsGoodMuPs, higgsGoodMuMs;

    for(HepMC::GenVertex::particle_iterator iDesc = genVtx->particles_begin(HepMC::descendants);
        iDesc != genVtx->particles_end(HepMC::descendants); ++iDesc) {
      HepMC::GenParticle* desc = *iDesc;

      if ( desc->status() != 1 ) continue;

      if ( isMuM(desc) ) higgsMuMs.push_back(desc);
      else if ( isMuP(desc) ) higgsMuPs.push_back(desc);
      
      if ( isGoodMuM(desc) ) {
        higgsGoodMuMs.push_back(desc);
        hHppGoodMuM_(desc);
      }
      else if ( isGoodMuP(desc) ) {
        higgsGoodMuPs.push_back(desc);
        hHppGoodMuP_(desc);
      }
    }

    if ( higgsMuPs.size() == 2 ) hDimuonPP_(make_pair(higgsMuPs[0], higgsMuPs[1]));
    if ( higgsMuMs.size() == 2 ) hDimuonMM_(make_pair(higgsMuMs[0], higgsMuMs[1]));

    if ( higgsGoodMuPs.size() == 2 ) hGoodDimuonPP_(make_pair(higgsGoodMuPs[0], higgsMuPs[1]));
    if ( higgsGoodMuMs.size() == 2 ) hGoodDimuonMM_(make_pair(higgsGoodMuMs[0], higgsMuMs[1]));

    hNHiggsMu_->Fill(higgsMuPs.size()+higgsMuMs.size());
    hNHiggsGoodMu_->Fill(higgsGoodMuPs.size()+higgsGoodMuMs.size());
  }

}

void Hpp2MuHepMCAnalyzer::endJob()
{
}

/* vim:set ts=2 sts=2 sw=2 expandtab: */
