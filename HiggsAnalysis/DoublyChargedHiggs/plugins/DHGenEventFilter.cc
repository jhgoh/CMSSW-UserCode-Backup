#include "HiggsAnalysis/DoublyChargedHiggs/plugins/DHGenEventFilter.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "HiggsAnalysis/DoublyChargedHiggs/interface/LeptonTypes.h"

using namespace std;

DHGenEventFilter::DHGenEventFilter(const edm::ParameterSet& pset)
{
  genLabel_ = pset.getParameter<edm::InputTag>("genLabel");

  // decay modes input : ee, mm, tt, em, et, mt
  const std::string decayExpr1 = pset.getParameter<std::string>("decay1");
  const std::string decayExpr2 = pset.getParameter<std::string>("decay2");

  decay1_ = LeptonTypes::getType(&decayExpr1[0]) | LeptonTypes::getType(&decayExpr1[1]);
  decay2_ = LeptonTypes::getType(&decayExpr2[0]) | LeptonTypes::getType(&decayExpr2[1]);

  hL_pdgId_ = abs(pset.getUntrackedParameter<int>("hL_pdgId", 9900041));
  hR_pdgId_ = abs(pset.getUntrackedParameter<int>("hR_pdgId", 9900042));
}

DHGenEventFilter::~DHGenEventFilter()
{
}

void DHGenEventFilter::beginJob()
{
}

void DHGenEventFilter::endJob()
{
}

bool DHGenEventFilter::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::GenParticleCollection> genHandle;
  // Pass event when no genParticle exists - make ineffective in case of real data
  if ( !event.getByLabel(genLabel_, genHandle) ) return true;

  int nHpp = 0, nHmm = 0;
  int nDecayMode1 = 0, nDecayMode2 = 0;
  for ( reco::GenParticleCollection::const_iterator genIter = genHandle->begin();
        genIter != genHandle->end(); ++genIter )
  {
    const reco::GenParticle& ptcl = *genIter;

    if ( ptcl.status() != 3 ) continue;
    if ( abs(ptcl.pdgId()) != hL_pdgId_ &&
         abs(ptcl.pdgId()) != hR_pdgId_ ) continue;

    int nDau = 0, dau1Idx = 0, dau2Idx = 0;
    const int nDauCand = ptcl.numberOfDaughters();
    if ( nDauCand < 2 ) continue;
    for ( int i=0; i<nDauCand; ++i )
    {
      const reco::GenParticle* dauCand = dynamic_cast<const reco::GenParticle*>(ptcl.daughter(i));
      if ( !dauCand ) continue;
      if ( dauCand->status() != 3 ) continue;

      ++nDau;
      if ( nDau == 1 ) dau1Idx = i;
      else if ( nDau == 2 ) dau2Idx = i;
    }
    if ( nDau != 2 ) continue;

    const reco::GenParticle* dau1 = dynamic_cast<const reco::GenParticle*>(ptcl.daughter(dau1Idx));
    const reco::GenParticle* dau2 = dynamic_cast<const reco::GenParticle*>(ptcl.daughter(dau2Idx));

    if ( !dau1 || !dau2 ) continue;

    const int decay = (LeptonTypes::getType(dau1->pdgId()) | LeptonTypes::getType(dau2->pdgId()));
    if ( decay & LeptonTypes::None ) continue;
  
    if ( decay == decay1_ ) ++nDecayMode1;
    if ( decay == decay2_ ) ++nDecayMode2;
  
    ptcl.pdgId() > 0 ? ++nHpp : ++nHmm;
  }

  if ( nHpp < 1 or nHmm < 1 ) return false;
  if ( nDecayMode1 < 1 or nDecayMode2 < 1 ) return false;
  else return true;
}

