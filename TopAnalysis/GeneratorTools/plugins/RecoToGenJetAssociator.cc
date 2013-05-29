#ifndef TopAnalysis_GeneratorTools_RecoToGenJetAssociator_H
#define TopAnalysis_GeneratorTools_RecoToGenJetAssociator_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "AnalysisDataFormats/CMGTools/interface/GenericTypes.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToOne.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <memory>
#include <vector>
#include <string>

class RecoToGenJetAssociator : public edm::EDProducer
{
public:
  RecoToGenJetAssociator(const edm::ParameterSet& pset);
  ~RecoToGenJetAssociator() {};

  void produce(edm::Event& event, const edm::EventSetup& eventSetup);

  //typedef edm::AssociationMap<edm::OneToOne<std::vector<cmg::PFJet>, std::vector<reco::GenJet> > > RecoToGenJetMap;
  typedef edm::AssociationMap<edm::OneToOne<std::vector<cmg::PFJet>, std::vector<cmg::GenJet> > > CMGRecoToGenJetMap;

private:
  edm::InputTag recoJetLabel_;
  edm::InputTag genJetLabel_;
  double cut_maxDR_, cut_maxDPt_;

};

RecoToGenJetAssociator::RecoToGenJetAssociator(const edm::ParameterSet& pset)
{
  recoJetLabel_ = pset.getParameter<edm::InputTag>("recoJets");
  genJetLabel_ = pset.getParameter<edm::InputTag>("genJets");

  cut_maxDR_ = cut_maxDPt_ = 1e9;
  edm::ParameterSet cuts = pset.getParameter<edm::ParameterSet>("cuts");
  if ( cuts.exists("maxDR") ) cut_maxDR_ = cuts.getParameter<double>("maxDR");
  if ( cuts.exists("maxDPt") ) cut_maxDPt_ = cuts.getParameter<double>("maxDPt");

  produces<CMGRecoToGenJetMap>();
}

void RecoToGenJetAssociator::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  std::auto_ptr<CMGRecoToGenJetMap> recoToGenJetMap(new CMGRecoToGenJetMap);

  edm::Handle<std::vector<cmg::PFJet> > recoJetHandle;
  event.getByLabel(recoJetLabel_, recoJetHandle);

  //edm::Handle<std::vector<reco::GenJet> > genJetHandle;
  edm::Handle<std::vector<cmg::GenJet> > genJetHandle;
  event.getByLabel(genJetLabel_, genJetHandle);

  for ( int i=0, n=recoJetHandle->size(); i<n; ++i )
  {
    const cmg::PFJet& recoJet = recoJetHandle->at(i);
    int matchedJetIndex = -1;
    double maxDR = cut_maxDR_;
    double maxDPt = cut_maxDPt_;

    // Find best pair with dR or dPt cut
    for ( int j=0, m=genJetHandle->size(); j<m; ++j )
    {
      //const reco::GenJet& genJet = genJetHandle->at(j);
      const cmg::GenJet& genJet = genJetHandle->at(j);

      const double dR = deltaR(recoJet, genJet);
      const double dPt = abs(recoJet.pt() - genJet.pt());
      if ( dR > maxDR and cut_maxDR_ < 1e9 ) continue;
      if ( dPt > maxDPt and cut_maxDPt_ < 1e9 ) continue;

      maxDR = dR;
      maxDPt = dPt;
      matchedJetIndex = j;
    }

    if ( matchedJetIndex == -1 ) continue;
      
    // Now we have best matching reco->gen jet pair
    edm::Ref<std::vector<cmg::PFJet> > recoJetRef(recoJetHandle, i);
    //edm::Ref<std::vector<reco::GenJet> > genJetRef(genJetHandle, matchedJetIndex);
    edm::Ref<std::vector<cmg::GenJet> > genJetRef(genJetHandle, matchedJetIndex);
    recoToGenJetMap->insert(recoJetRef, genJetRef);
  }

  event.put(recoToGenJetMap);
}

DEFINE_FWK_MODULE(RecoToGenJetAssociator);

#endif

