#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToOne.h"

namespace cmg {
  edm::Wrapper<edm::AssociationMap<edm::OneToOne<std::vector<cmg::PFJet>,std::vector<reco::GenJet>,unsigned int> > > dummy00;
  edm::AssociationMap<edm::OneToOne<std::vector<cmg::PFJet>, std::vector<reco::GenJet> > > dummy01;
  edm::RefProd<std::vector<cmg::PFJet> > dummy02;

  edm::Wrapper<edm::AssociationMap<edm::OneToMany<std::vector<reco::GenJet>,std::vector<reco::GenParticle>,unsigned int> > > dummy03;

}
