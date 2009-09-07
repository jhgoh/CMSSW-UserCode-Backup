#include "HiggsAnalysis/DoublyChargedHiggs/interface/DileptonProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
//#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

using namespace edm;
using namespace std;

namespace LeptonTypes
{
  enum 
  {
    None = 0, Muon, Electron
  };
}

DileptonProducer::DileptonProducer(const edm::ParameterSet& pset):
  kvFitter_(0)
{
  // Lepton1 Input
  edm::ParameterSet lepton1Set = pset.getParameter<edm::ParameterSet>("lepton1");

  lepton1Label_ = lepton1Set.getParameter<edm::InputTag>("src");

  const string lepton1Type = lepton1Set.getParameter<string>("type");
  if ( lepton1Type == "muon" ) lepton1Type_ = LeptonTypes::Muon;
  else if ( lepton1Type == "electron" ) lepton1Type_ = LeptonTypes::Electron;
  else lepton1Type_ = LeptonTypes::None;

  // Lepton2 Input
  edm::ParameterSet lepton2Set = pset.getParameter<edm::ParameterSet>("lepton2");
  
  lepton2Label_ = lepton2Set.getParameter<edm::InputTag>("src");
  
  const string lepton2Type = lepton2Set.getParameter<string>("type");
  if ( lepton2Type == "muon" ) lepton2Type_ = LeptonTypes::Muon;
  else if ( lepton2Type == "electron" ) lepton2Type_ = LeptonTypes::Electron;
  else lepton2Type_ = LeptonTypes::None;

  // Compositing options
  isSameCollection_ = lepton1Label_ == lepton2Label_ ? true : false;

  // Register output product
  produces<pat::CompositeCandidateCollection>();

  // Vertex fitter setup
  const string fitterType = pset.getParameter<string>("fitterType");
  kvFitter_ = 0;
  if ( fitterType == "KalmanVertexFitter" )
  {
    kvFitter_ = new CandCommonVertexFitter<KalmanVertexFitter>(pset.getParameter<edm::ParameterSet>("vertexFitSet"));
  }
}

DileptonProducer::~DileptonProducer()
{
  if ( kvFitter_ ) delete kvFitter_;
}

void DileptonProducer::beginJob(const edm::EventSetup& eventSetup)
{
}

void DileptonProducer::endJob()
{
}

// Candidate combiner : Same types, same sources
template <typename LeptonIter, typename OutCollection>
void DileptonProducer::combineLeptons(const LeptonIter lepton_begin, const LeptonIter lepton_end,
                                      OutCollection& dileptonCands)
{
  for ( LeptonIter lepton1 = lepton_begin; lepton1 != lepton_end; ++lepton1 )
  {
    for ( LeptonIter lepton2 = lepton1+1; lepton2 != lepton_end; ++lepton2 )
    {
      pat::CompositeCandidate dileptonCand;
      dileptonCand.addDaughter(*lepton1);
      dileptonCand.addDaughter(*lepton2);

      AddFourMomenta addP4;
      addP4.set(dileptonCand);

      // FIXME : Add vertex fit routine here
      dileptonCands->push_back(dileptonCand);
    }
  }
}

// Candidate combiner : same types, different sources
template <typename LeptonIter, typename OutCollection>
void DileptonProducer::combineLeptons(const LeptonIter lepton1_begin, const LeptonIter lepton1_end,
                                      const LeptonIter lepton2_begin, const LeptonIter lepton2_end,
                                      OutCollection& dileptonCands)
{
  for ( LeptonIter lepton1 = lepton1_begin; lepton1 != lepton1_end; ++lepton1 )
  {
    for ( LeptonIter lepton2 = lepton2_begin; lepton2 != lepton2_end; ++lepton2 )
    {
      pat::CompositeCandidate dileptonCand;
      dileptonCand.addDaughter(*lepton1);
      dileptonCand.addDaughter(*lepton2);

      AddFourMomenta addP4;
      addP4.set(dileptonCand);

      // FIXME : Add vertex fit routine here
      dileptonCands->push_back(dileptonCand);
    }
  }
}

// Candidate combiner : different types and different sources
template <typename LeptonIter1, typename LeptonIter2, typename OutCollection>
void DileptonProducer::combineLeptons(const LeptonIter1 lepton1_begin, const LeptonIter1 lepton1_end,
                                      const LeptonIter2 lepton2_begin, const LeptonIter2 lepton2_end,
                                      OutCollection& dileptonCands)
{
  for ( LeptonIter1 lepton1 = lepton1_begin; lepton1 != lepton1_end; ++lepton1 )
  {
    for ( LeptonIter2 lepton2 = lepton2_begin; lepton2 != lepton2_end; ++lepton2 )
    {
      pat::CompositeCandidate dileptonCand;
      dileptonCand.addDaughter(*lepton1);
      dileptonCand.addDaughter(*lepton2);

      AddFourMomenta addP4;
      addP4.set(dileptonCand);

      // FIXME : Add vertex fit routine here
      dileptonCands->push_back(dileptonCand);
    }
  }
}

void DileptonProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  typedef pat::MuonCollection::const_iterator MuonIter;
  typedef pat::ElectronCollection::const_iterator ElectronIter;

  // B-field information for vertex fit
  edm::ESHandle<MagneticField> magneticField;
  eventSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  if ( kvFitter_ ) kvFitter_->set(magneticField.product());

  // Manually combine two leptons
  typedef auto_ptr<pat::CompositeCandidateCollection> CompositeCandCollPtr;
  CompositeCandCollPtr dileptonCands(new pat::CompositeCandidateCollection);

  if ( lepton1Type_ == LeptonTypes::Muon && lepton2Type_ == LeptonTypes::Muon )
  {
    if ( isSameCollection_ ) 
    {
      edm::Handle<pat::MuonCollection> lepton1Handle;
      event.getByLabel(lepton1Label_, lepton1Handle);

      combineLeptons<MuonIter, CompositeCandCollPtr>(lepton1Handle->begin(), lepton1Handle->end(), dileptonCands);
    }
    else
    {
      edm::Handle<pat::MuonCollection> lepton1Handle, lepton2Handle;
      event.getByLabel(lepton1Label_, lepton1Handle);
      event.getByLabel(lepton2Label_, lepton2Handle);

      combineLeptons<MuonIter, CompositeCandCollPtr>(lepton1Handle->begin(), lepton1Handle->end(),
                               lepton2Handle->begin(), lepton2Handle->end(), dileptonCands);
    }
  }
  else if ( lepton1Type_ == LeptonTypes::Electron && lepton2Type_ == LeptonTypes::Electron )
  {
    if ( isSameCollection_ )
    {
      edm::Handle<pat::ElectronCollection> lepton1Handle;
      event.getByLabel(lepton1Label_, lepton1Handle);

      combineLeptons<ElectronIter, CompositeCandCollPtr>(lepton1Handle->begin(), lepton1Handle->end(), dileptonCands);
    }
    else
    {
      edm::Handle<pat::ElectronCollection> lepton1Handle, lepton2Handle;
      event.getByLabel(lepton1Label_, lepton1Handle);
      event.getByLabel(lepton2Label_, lepton2Handle);

      combineLeptons<ElectronIter, CompositeCandCollPtr>(lepton1Handle->begin(), lepton1Handle->end(),
                                   lepton2Handle->begin(), lepton2Handle->end(), dileptonCands);
    }
  }
  else if ( lepton1Type_ == LeptonTypes::Electron && lepton2Type_ == LeptonTypes::Muon )
  {
    edm::Handle<pat::ElectronCollection> lepton1Handle;
    edm::Handle<pat::MuonCollection> lepton2Handle;
    event.getByLabel(lepton1Label_, lepton1Handle);
    event.getByLabel(lepton2Label_, lepton2Handle);

    combineLeptons<ElectronIter, MuonIter, CompositeCandCollPtr>(lepton1Handle->begin(), lepton1Handle->end(),
                                           lepton2Handle->begin(), lepton2Handle->end(), dileptonCands);
  }
  else if ( lepton1Type_ == LeptonTypes::Muon && lepton2Type_ == LeptonTypes::Electron )
  {
    edm::Handle<pat::MuonCollection> lepton1Handle;
    edm::Handle<pat::ElectronCollection> lepton2Handle;
    event.getByLabel(lepton1Label_, lepton1Handle);
    event.getByLabel(lepton2Label_, lepton2Handle);

    combineLeptons<ElectronIter, MuonIter, CompositeCandCollPtr>(lepton2Handle->begin(), lepton2Handle->end(),
                                           lepton1Handle->begin(), lepton1Handle->end(), dileptonCands);
  }

  event.put(dileptonCands);
}

/* vim:set ts=2 sts=2 sw=2 expandtab: */

