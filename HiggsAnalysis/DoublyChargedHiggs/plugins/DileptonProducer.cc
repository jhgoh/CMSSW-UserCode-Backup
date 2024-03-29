#include "HiggsAnalysis/DoublyChargedHiggs/plugins/DileptonProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "HiggsAnalysis/DoublyChargedHiggs/interface/LeptonTypes.h"

using namespace std;

DileptonProducer::DileptonProducer(const edm::ParameterSet& pset)
{
  // Set lepton1 input
  edm::ParameterSet lepton1Set = pset.getParameter<edm::ParameterSet>("lepton1");
  
  lepton1Label_ = lepton1Set.getParameter<edm::InputTag>("src");
  lepton1Type_ = LeptonTypes::getType(lepton1Set.getParameter<string>("type"));
  lepton1Charge_ = lepton1Set.getParameter<int>("charge");

  // Set lepton2 input
  edm::ParameterSet lepton2Set = pset.getParameter<edm::ParameterSet>("lepton2");
  
  lepton2Label_ = lepton2Set.getParameter<edm::InputTag>("src");
  lepton2Type_ = LeptonTypes::getType(lepton2Set.getParameter<string>("type"));
  lepton2Charge_ = lepton2Set.getParameter<int>("charge");


  if ( lepton1Type_ == LeptonTypes::None or lepton2Type_ == LeptonTypes::None ) 
  {
    edm::LogError("DileptonProducer") << "Source lepton type is wrong or not supported\n"
                                      << "-- lepton1 = " << lepton1Set.getParameter<string>("type") << '\n'
                                      << "-- lepton2 = " << lepton2Set.getParameter<string>("type") << endl;
    return;
  }

  // Compositing options
  chargeConj_ = pset.getParameter<bool>("chargeConjugation"); // Do charge conjugation?
  isSameCollection_ = (lepton1Label_ == lepton2Label_);

  dileptonCharge_ = lepton1Charge_ + lepton2Charge_;
  if ( chargeConj_ ) dileptonCharge_ = abs(dileptonCharge_);

  // Raw cuts
  lepton1MinPt_ = lepton1Set.getParameter<double>("minPt");
  lepton1MaxEta_ = lepton1Set.getParameter<double>("maxEta");
  
  lepton2MinPt_ = lepton2Set.getParameter<double>("minPt");
  lepton2MaxEta_ = lepton2Set.getParameter<double>("maxEta");

  dileptonMinMass_ = pset.getParameter<double>("minMass");
  dileptonMinPt_ = pset.getParameter<double>("minPt");
  dileptonMaxEta_ = pset.getParameter<double>("maxEta");

  // Register output product
  produces<pat::CompositeCandidateCollection>();
}

DileptonProducer::~DileptonProducer()
{
}

void DileptonProducer::beginJob()
{
}

void DileptonProducer::endJob()
{
}

void DileptonProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  typedef pat::MuonCollection::const_iterator MuonIter;
  typedef pat::ElectronCollection::const_iterator ElectronIter;

  // Create dilepton candidate collection
  typedef auto_ptr<pat::CompositeCandidateCollection> CompositeCandCollPtr;
  CompositeCandCollPtr dileptonCands(new pat::CompositeCandidateCollection);

  // Combine leptons
  // -- EE channel
  if ( (lepton1Type_ | lepton2Type_ ) == LeptonTypes::Electron )
  {
    edm::Handle<pat::ElectronCollection> e1Handle;
    event.getByLabel(lepton1Label_, e1Handle);

    if ( isSameCollection_ )
    {
      combineLeptons(e1Handle->begin(), e1Handle->end(), dileptonCands);
    }
    else
    {
      edm::Handle<pat::ElectronCollection> e2Handle;
      event.getByLabel(lepton2Label_, e2Handle);
    
      combineLeptons(e1Handle->begin(), e1Handle->end(), e2Handle->begin(), e2Handle->end(), dileptonCands);
    }
  }
  // -- MuMu channel
  else if ( (lepton1Type_ | lepton2Type_) == LeptonTypes::Muon )
  {
    edm::Handle<pat::MuonCollection> mu1Handle;
    event.getByLabel(lepton1Label_, mu1Handle);

    if ( isSameCollection_ )
    {
      combineLeptons(mu1Handle->begin(), mu1Handle->end(), dileptonCands);
    }
    else
    {
      edm::Handle<pat::MuonCollection> mu2Handle;
      event.getByLabel(lepton2Label_, mu2Handle);

      combineLeptons(mu1Handle->begin(), mu1Handle->end(), mu2Handle->begin(), mu2Handle->end(), dileptonCands);
    }
  }
  // -- MuE or EMu channel
  else if ( (lepton1Type_ | lepton2Type_) == (LeptonTypes::Electron | LeptonTypes::Muon) )
  {
    edm::Handle<pat::ElectronCollection> eHandle;
    edm::Handle<pat::MuonCollection> muHandle;

    if ( lepton1Type_ == LeptonTypes::Electron )
    {
      event.getByLabel(lepton1Label_, eHandle); 
      event.getByLabel(lepton2Label_, muHandle);
    }
    else
    {
      event.getByLabel(lepton1Label_, muHandle);
      event.getByLabel(lepton2Label_, eHandle);
    }
    
    combineLeptons(eHandle->begin(), eHandle->end(), muHandle->begin(), muHandle->end(), dileptonCands);
  }
  // -- Otherwise, do nothing
  else
  {
    edm::LogError("DileptonProducer::produce()") << "Input collection must be electron or muons. "
                                                    "Thus this should not happen\n";
  }

  // Now put the candidate collection to the event
  event.put(dileptonCands);
}

// Candidate combiner : Same sources 
template<typename LeptonIter, typename OutCollection>
void DileptonProducer::combineLeptons(const LeptonIter begin, const LeptonIter end,
                                      OutCollection& dileptonCands)
{
  for ( LeptonIter lepton1 = begin; lepton1 != end; ++lepton1 )
  {
    // Precuts for lepton1
    if ( lepton1->pt() < lepton1MinPt_ ) continue;
    if ( fabs(lepton1->eta()) > lepton1MaxEta_ ) continue;

    // Check charges
    const int lepton1Charge = lepton1->charge();
    if (  chargeConj_ && abs(lepton1Charge) != abs(lepton1Charge_) ) continue;
    if ( !chargeConj_ && lepton1Charge != lepton1Charge_ ) continue;

    for ( LeptonIter lepton2 = lepton1+1; lepton2 != end; ++lepton2 )
    {
      // Precuts for lepton2
      if ( lepton2->pt() < lepton2MinPt_ ) continue;
      if ( fabs(lepton2->eta()) > lepton2MaxEta_ ) continue;

      const int lepton2Charge = lepton2->charge();
      if (  chargeConj_ && abs(lepton1Charge+lepton2Charge) != dileptonCharge_ ) continue;
      if ( !chargeConj_ && lepton2Charge != lepton2Charge_ ) continue;

      pat::CompositeCandidate dileptonCand;
      dileptonCand.addDaughter(*lepton1);
      dileptonCand.addDaughter(*lepton2);

      AddFourMomenta addP4;
      addP4.set(dileptonCand);

      // Precuts for dilepton candidate
      if ( dileptonCand.mass() < dileptonMinMass_ ) continue;
      if ( dileptonCand.pt() < dileptonMinPt_ ) continue;
      if ( fabs(dileptonCand.eta()) > dileptonMaxEta_ ) continue;

      // FIXME :: Add vertex fit routine here

      dileptonCands->push_back(dileptonCand);
    }
  }
}

// Candidate combiner :: Different sources
template<typename LeptonIter1, typename LeptonIter2, typename OutCollection>
void DileptonProducer::combineLeptons(const LeptonIter1 begin1, const LeptonIter1 end1, 
                                      const LeptonIter2 begin2, const LeptonIter2 end2, 
                                      OutCollection& dileptonCands)
{
  for ( LeptonIter1 lepton1 = begin1; lepton1 != end1; ++lepton1 )
  {
    // Precuts for lepton1
    if ( lepton1->pt() < lepton1MinPt_ ) continue;
    if ( fabs(lepton1->eta()) > lepton1MaxEta_ ) continue;

    // Check charges
    const int lepton1Charge = lepton1->charge();
    if (  chargeConj_ && abs(lepton1Charge) != abs(lepton1Charge_) ) continue;
    if ( !chargeConj_ && lepton1Charge != lepton1Charge_ ) continue;

    for ( LeptonIter2 lepton2 = begin2; lepton2 != end2; ++lepton2 )
    {
      // Precuts for lepton2
      if ( lepton2->pt() < lepton2MinPt_ ) continue;
      if ( fabs(lepton2->eta()) > lepton2MaxEta_ ) continue;

      const int lepton2Charge = lepton2->charge();
      if (  chargeConj_ && abs(lepton1Charge+lepton2Charge) != dileptonCharge_ ) continue;
      if ( !chargeConj_ && lepton2Charge != lepton2Charge_ ) continue;

      pat::CompositeCandidate dileptonCand;
      dileptonCand.addDaughter(*lepton1); 
      dileptonCand.addDaughter(*lepton2);

      AddFourMomenta addP4;
      addP4.set(dileptonCand);

      // Precuts for dilepton candidate
      if ( dileptonCand.mass() < dileptonMinMass_ ) continue;
      if ( dileptonCand.pt() < dileptonMinPt_ ) continue;
      if ( fabs(dileptonCand.eta()) > dileptonMaxEta_ ) continue;

      // FIXME :: Add vertex fit routine here

      dileptonCands->push_back(dileptonCand);
    }
  }
}
