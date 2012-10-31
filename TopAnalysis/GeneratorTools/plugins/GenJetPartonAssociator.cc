#ifndef TopAnalysis_GeneratorTools_GenJetPartonAssociator_H
#define TopAnalysis_GeneratorTools_GenJetPartonAssociator_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToMany.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <memory>
#include <vector>
#include <string>

using namespace std;

class GenJetPartonAssociator : public edm::EDProducer
{
public:
  GenJetPartonAssociator(const edm::ParameterSet& pset);
  ~GenJetPartonAssociator() {};

  void produce(edm::Event& event, const edm::EventSetup& eventSetup);

  bool hasMother(const reco::Candidate* p, const reco::Candidate* mother);
  typedef edm::AssociationMap<edm::OneToMany<std::vector<reco::GenJet>, reco::GenParticleCollection> > GenJetToGenParticlesMap;

private:
  edm::InputTag genParticleLabel_;
  edm::InputTag genJetLabel_;
  //double maxDR_, maxDPt_;
  int cut_minNConstituent_;
  double cut_minFracConstituent_;

};

GenJetPartonAssociator::GenJetPartonAssociator(const edm::ParameterSet& pset)
{
  genParticleLabel_ = pset.getParameter<edm::InputTag>("genParticles");
  genJetLabel_ = pset.getParameter<edm::InputTag>("genJets");

  cut_minNConstituent_ = 1; // Minimum requirement to say "matched"
  cut_minFracConstituent_ = 0.; 

  //minDR_ = minDPt_ = 1e9;
  edm::ParameterSet cuts = pset.getParameter<edm::ParameterSet>("cuts");
  //if ( cuts.exists("minDR") ) minDR_ = cuts.getParameter<double>("minDR");
  //if ( cuts.exists("minDPt") ) minDPt_ = cuts.getParameter<double>("minDPt");
  if ( cuts.exists("minNConstituent") ) cut_minNConstituent_ = cuts.getParameter<int>("minNConstituent");
  if ( cuts.exists("minFracConstituent") ) cut_minFracConstituent_ = cuts.getParameter<double>("minFracConstituent");

  produces<GenJetToGenParticlesMap>();
}

void GenJetPartonAssociator::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  std::auto_ptr<GenJetToGenParticlesMap> genJetToGenParticlesMap(new GenJetToGenParticlesMap);

  edm::Handle<std::vector<reco::GenJet> > genJetHandle;
  event.getByLabel(genJetLabel_, genJetHandle);

  edm::Handle<std::vector<reco::GenParticle> > genParticleHandle;
  event.getByLabel(genParticleLabel_, genParticleHandle);

  // Collect list of gen particles in hard process
  std::vector<unsigned int> genPartonIndicies;
  for ( int i=0, n=genParticleHandle->size(); i<n; ++i )
  {
    const reco::GenParticle& p = genParticleHandle->at(i);
    if ( p.status() != 3 ) continue;
    if ( abs(p.pdgId()) == 2212 and p.pt() == 0 and abs(p.pz()) > 1000 ) continue; // ignore initial incident protons

    genPartonIndicies.push_back(i);
  }
  const int nGenParton = genPartonIndicies.size();

  // Find matching between genJet to genParticle
  for ( int i=0, n=genJetHandle->size(); i<n; ++i )
  {
    std::vector<unsigned int> matchedPartons;

    const reco::GenJet& genJet = genJetHandle->at(i);
    std::vector<const reco::GenParticle*> genConstituents = genJet.getGenConstituents();
    const int nConstituent = genConstituents.size();
    if ( nConstituent < cut_minNConstituent_ ) continue;

    for ( int k=0; k<nGenParton; ++k )
    {
      const reco::GenParticle* genParton = &genParticleHandle->at(genPartonIndicies[k]);

      int nMatched = 0;

      for ( int j=0; j<nConstituent; ++j )
      {
        const reco::GenParticle* p = genConstituents[j];
        if ( !hasMother(p, genParton) ) continue;

        ++nMatched;
      }

      if ( nMatched < cut_minNConstituent_ ) continue;
      if ( 1.*nMatched/nConstituent < cut_minFracConstituent_ ) continue;

      cout << genJet.p4() << " : " << genParton->pdgId() << ':' << genParton->p4();
      cout << ' ' << nMatched << '/' << genConstituents.size() << '=' << 1.*nMatched/nConstituent << endl;

      matchedPartons.push_back(k);
    }

    if ( matchedPartons.empty() ) continue;

    edm::Ref<std::vector<reco::GenJet> > genJetRef(genJetHandle, i);
    for ( int j=0, m=matchedPartons.size(); j<m; ++j )
    {
      const unsigned int k = matchedPartons[j];
      edm::Ref<std::vector<reco::GenParticle> > genParticleRef(genParticleHandle, k);
      genJetToGenParticlesMap->insert(genJetRef, genParticleRef);
    }
  }

  event.put(genJetToGenParticlesMap);
}

bool GenJetPartonAssociator::hasMother(const reco::Candidate* p, const reco::Candidate* mother)
{
  if ( !p or !mother ) return false;
  if ( p == mother ) return false;

  const reco::Candidate* m = p->mother();
  if ( m == mother ) return true;

  return hasMother(m, mother);
}

DEFINE_FWK_MODULE(GenJetPartonAssociator);

#endif

