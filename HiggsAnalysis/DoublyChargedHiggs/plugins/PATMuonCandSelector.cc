#include "HiggsAnalysis/DoublyChargedHiggs/interface/PATMuonCandSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

using namespace std; 

PATMuonCandSelector::PATMuonCandSelector(const edm::ParameterSet& pset)
{
  muonLabel_ = pset.getParameter<edm::InputTag>("muon");

  edm::ParameterSet muonCutSet = pset.getParameter<edm::ParameterSet>("muonCutSet");
  charge_ = muonCutSet.getUntrackedParameter<int>("charge");
  useGlobalMuonsOnly_ = muonCutSet.getUntrackedParameter<bool>("useGlobalMuonsOnly");
  minPt_ = muonCutSet.getUntrackedParameter<double>("minPt");

  produces<pat::MuonCollection>();
}

void PATMuonCandSelector::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<pat::MuonCollection> muonHandle;
  event.getByLabel(muonLabel_, muonHandle);

  auto_ptr<pat::MuonCollection> muonCands(new pat::MuonCollection);

  for ( pat::MuonCollection::const_iterator iMuon = muonHandle->begin();
        iMuon != muonHandle->end(); ++iMuon )
  {
    const pat::Muon& muon = *iMuon;

    if ( charge_ != 0 && muon.charge() != charge_ ) continue;
    if ( useGlobalMuonsOnly_ && !muon.isGlobalMuon() ) continue;
    if ( muon.pt() < minPt_ ) continue;

    muonCands->push_back(muon);
  }

  event.put(muonCands);
}

/* vim:set ts=2 sts=2 sw=2 expandtab: */
