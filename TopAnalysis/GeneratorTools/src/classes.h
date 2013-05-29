#include "AnalysisDataFormats/CMGTools/interface/GenericTypes.h"
#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToOne.h"

namespace cmg {
  edm::RefProd<std::vector<cmg::PFJet> > dummy00;
  edm::RefProd<std::vector<cmg::GenJet> > dummy01;

  // cmg::PFJet -> reco::GenJet mapping
  edm::Wrapper<edm::AssociationMap<edm::OneToOne<std::vector<cmg::PFJet>,std::vector<reco::GenJet>,unsigned int> > > dummy10;
  edm::AssociationMap<edm::OneToOne<std::vector<cmg::PFJet>, std::vector<reco::GenJet> > > dummy11;

  // reco::GenJet -> reco::GenParticle mapping
  edm::Wrapper<edm::AssociationMap<edm::OneToMany<std::vector<reco::GenJet>,std::vector<reco::GenParticle>,unsigned int> > > dummy12;

  // cmg::PFJet -> cmg::GenJet mapping
  edm::Wrapper<edm::AssociationMap<edm::OneToOne<std::vector<cmg::PFJet>,std::vector<cmg::GenJet>,unsigned int> > > dummy20;
  edm::AssociationMap<edm::OneToOne<std::vector<cmg::PFJet>, std::vector<cmg::GenJet> > > dummy21;

  // cmg::GenJet -> reco::GenParticle mapping
  edm::Wrapper<edm::AssociationMap<edm::OneToMany<std::vector<cmg::GenJet>,std::vector<reco::GenParticle>,unsigned int> > > dummy22;

}
