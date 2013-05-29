#include "HiggsAnalysis/DoublyChargedHiggs/interface/DimuonProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
//#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

using namespace edm;
using namespace std;

DimuonProducer::DimuonProducer(const edm::ParameterSet& pset):
  kvFitter_(0)
{
  muon1Label_ = pset.getParameter<edm::InputTag>("muon1");
  muon2Label_ = pset.getParameter<edm::InputTag>("muon2");
  isSameCollection_ = muon1Label_ == muon2Label_ ? true : false;

  produces<pat::CompositeCandidateCollection>();

  const string fitterType = pset.getParameter<string>("fitterType");
  kvFitter_ = 0;
  if ( fitterType == "KalmanVertexFitter" )
  {
    kvFitter_ = new CandCommonVertexFitter<KalmanVertexFitter>(pset.getParameter<edm::ParameterSet>("vertexFitSet"));
  }
}

DimuonProducer::~DimuonProducer()
{
  if ( kvFitter_ ) delete kvFitter_;
}

void DimuonProducer::beginJob(const edm::EventSetup& eventSetup)
{
}

void DimuonProducer::endJob()
{
}

void DimuonProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<pat::MuonCollection> muon1Handle, muon2Handle;
  event.getByLabel(muon1Label_, muon1Handle);
  if ( !isSameCollection_ ) event.getByLabel(muon2Label_, muon2Handle);

  // B-field information for vertex fit
  edm::ESHandle<MagneticField> magneticField;
  eventSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  kvFitter_->set(magneticField.product());

  // Manually combine two muons
  auto_ptr<pat::CompositeCandidateCollection> dimuonCands(new pat::CompositeCandidateCollection);

  // Start candidate combination
  typedef pat::MuonCollection::const_iterator MuonIter;
  const MuonIter muon1_begin = muon1Handle->begin();
  const MuonIter muon1_end = muon1Handle->end();

  for ( MuonIter iMuon1 = muon1_begin; iMuon1 != muon1_end; ++iMuon1 )
  {
    const pat::Muon& muon1 = *iMuon1;

    const MuonIter muon2_begin = isSameCollection_ ? iMuon1+1 : muon2Handle->begin();
    const MuonIter muon2_end = isSameCollection_ ? muon1_end : muon2Handle->end();
    for ( MuonIter iMuon2 = muon2_begin; iMuon2 != muon2_end; ++iMuon2 )
    {
      const pat::Muon& muon2 = *iMuon2;
      //if ( muon1 == muon2 ) continue;

      pat::CompositeCandidate dimuonCand;
      dimuonCand.addDaughter(muon1);
      dimuonCand.addDaughter(muon2);

      AddFourMomenta addP4;
      addP4.set(dimuonCand);

      // FIXME : Add vertex fit routine here
      dimuonCands->push_back(dimuonCand);
    }
  }

  event.put(dimuonCands);
}

/* vim:set ts=2 sts=2 sw=2 expandtab: */

