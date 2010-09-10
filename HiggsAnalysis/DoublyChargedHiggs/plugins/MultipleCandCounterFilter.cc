#include "HiggsAnalysis/DoublyChargedHiggs/plugins/MultipleCandCounterFilter.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Candidate/interface/Candidate.h"

using namespace std;

MultipleCandCounterFilter::MultipleCandCounterFilter(const edm::ParameterSet& pset)
{
  candTags_ = pset.getParameter<std::vector<edm::InputTag> >("cands");
  ptThresholds_ = pset.getParameter<std::vector<double> >("ptThresholds");
  sort(ptThresholds_.rbegin(), ptThresholds_.rend());
}

MultipleCandCounterFilter::~MultipleCandCounterFilter()
{
}

bool MultipleCandCounterFilter::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  const unsigned int nCandMin = ptThresholds_.size();
  unsigned int nCand = 0;

  std::vector<double> ptValues;
  for ( std::vector<edm::InputTag>::const_iterator candTag = candTags_.begin();
        candTag != candTags_.end(); ++candTag )
  {
    edm::Handle<edm::View<reco::Candidate> > candHandle;
    if ( !event.getByLabel(*candTag, candHandle) )
    {
      continue;
    }

    nCand += candHandle->size();

    for ( edm::View<reco::Candidate>::const_iterator cand = candHandle->begin();
          cand != candHandle->end(); ++cand )
    {
      ptValues.push_back(cand->pt());
    }
  }

  if ( nCand < nCandMin ) return false;

  // Check pt thresholds
  if ( ptValues.size() < nCandMin ) return false;
  sort(ptValues.rbegin(), ptValues.rend());

  for ( unsigned int i=0; i<nCandMin; ++i )
  {
    if ( ptValues[i] < ptThresholds_[i] ) return false;
  }

  return true;
}

