#include "HiggsAnalysis/DoublyChargedHiggs/src/Hpp2MuHepMCAnalyzer.h"

#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TH1D.h"
#include "TMath.h"

#include <iostream>

using namespace std;
using namespace edm;

Hpp2MuHepMCAnalyzer::Hpp2MuHepMCAnalyzer(const ParameterSet& pset)
{
  string outputFileName = pset.getUntrackedParameter<string>("outputFileName");

  edm::Service<TFileService> fs;

  hMuPt_ = fs->make<TH1D>("hMuPt", "Muon p_{T} from H++", 50, 0, 200);
  hMuEta_ = fs->make<TH1D>("hMuEta", "Muon #eta from H++", 50, -5.0, 5.0);
  hMuMuM_ = fs->make<TH1D>("hMuMuM", "Dimuon mass", 100, 50, 300);
  
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

  vector<HepMC::GenVertex*> higgsVtxs;

  for(HepMC::GenEvent::particle_const_iterator iGenPtcl = genEvt->particles_begin();
      iGenPtcl != genEvt->particles_end(); ++iGenPtcl) {
    HepMC::GenParticle* genPtcl = *iGenPtcl;

    if ( abs(genPtcl->pdg_id()) == 9900041 ||
         abs(genPtcl->pdg_id()) == 9900042 ) {
      HepMC::GenVertex* higgsDecVtx = genPtcl->end_vertex();
      higgsVtxs.push_back(higgsDecVtx);
    }
  }

  if ( higgsVtxs.empty() ) {
    LogError("Hpp2MuHepMCAnalyzer") << "No Higgs in this event\n";
    return;
  }

  vector<HepMC::GenParticle*> stableHiggsDescs;
  for(vector<HepMC::GenVertex*>::const_iterator iGenVtx = higgsVtxs.begin();
      iGenVtx != higgsVtxs.end(); ++iGenVtx) {
    HepMC::GenVertex* genVtx = *iGenVtx;
    if ( genVtx == 0 ) continue;

    for(HepMC::GenVertex::particle_iterator iDesc = genVtx->particles_begin(HepMC::descendants);
        iDesc != genVtx->particles_end(HepMC::descendants); ++iDesc) {
      HepMC::GenParticle* desc = *iDesc;

      if ( desc->status() == 1 ) {
        stableHiggsDescs.push_back(desc);
      }
    }
  }
}

void Hpp2MuHepMCAnalyzer::endJob()
{
}

/* vim:set ts=2 sts=2 sw=2 expandtab: */
