#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToMany.h"
#include "DataFormats/Common/interface/OneToOne.h"
#include "AnalysisDataFormats/CMGTools/interface/Electron.h"
#include "AnalysisDataFormats/CMGTools/interface/Muon.h"
#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "AnalysisDataFormats/CMGTools/interface/BaseMET.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"
#include "TH1F.h"

#include <memory>
#include <vector>
#include <string>

class EventTupleProducer : public edm::EDAnalyzer
{
public:
  EventTupleProducer(const edm::ParameterSet& pset);
  ~EventTupleProducer();

  //void beginJob();
  //void beginLuminosityBlock(edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  bool doMCMatch_;

  // Input objects
  edm::InputTag genLabel_;
  edm::InputTag genJetLabel_;
  edm::InputTag recoToGenJetMapLabel_;
  edm::InputTag genJetToPartonMapLabel_;
  std::string weightStr_;
  edm::InputTag vertexLabel_;

  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag metLabel_;
  std::string bTagType_;

  // Cuts
  StringCutObjectSelector<cmg::Muon, true>* isGoodMuon_;
  StringCutObjectSelector<cmg::Electron, true>* isGoodElectron_;
  StringCutObjectSelector<cmg::PFJet, true>* isGoodJet_;

  // Output tree
  TTree* eventTree_;
  int run_, lumi_, event_;
  double weight_, weightUp_, weightDn_;
  int nVertex_;
  math::XYZTLorentzVector met_;
  std::vector<math::XYZTLorentzVector> muons_, electrons_;
  std::vector<int> muons_Q_, electrons_Q_;
  std::vector<double> muons_Iso_, electrons_Iso_;
  std::vector<math::XYZTLorentzVector> jets_;
  std::vector<double> jets_bTag_;

  // Generator level information
  TTree* genTree_;
  // jet MC matching
  std::vector<int> jets_motherId_;
  //std::vector<math::XYZTLorentzVector> genJet_;
  //std::vector<int> genJetMotherId_;
  std::vector<math::XYZTLorentzVector> genMuons_, genElectrons_;
  std::vector<math::XYZTLorentzVector> genMuonNus_, genElectronNus_;

};

EventTupleProducer::EventTupleProducer(const edm::ParameterSet& pset)
{
  doMCMatch_ = pset.getParameter<bool>("doMCMatch");

  // Input labels
  genLabel_ = pset.getParameter<edm::InputTag>("gen");
  recoToGenJetMapLabel_ = pset.getParameter<edm::InputTag>("recoToGenJetMap");
  genJetToPartonMapLabel_ = pset.getParameter<edm::InputTag>("genJetToPartonsMap");

  weightStr_ = pset.getParameter<std::string>("weight");
  vertexLabel_ = pset.getParameter<edm::InputTag>("vertex");
  metLabel_ = pset.getParameter<edm::InputTag>("met");

  edm::ParameterSet electronPSet = pset.getParameter<edm::ParameterSet>("electron");
  std::string electronCut = electronPSet.getParameter<std::string>("cut");
  isGoodElectron_ = new StringCutObjectSelector<cmg::Electron, true>(electronCut);
  electronLabel_ = electronPSet.getParameter<edm::InputTag>("src");

  edm::ParameterSet muonPSet = pset.getParameter<edm::ParameterSet>("muon");
  std::string muonCut = muonPSet.getParameter<std::string>("cut");
  isGoodMuon_ = new StringCutObjectSelector<cmg::Muon, true>(muonCut);
  muonLabel_ = muonPSet.getParameter<edm::InputTag>("src");

  edm::ParameterSet jetPSet = pset.getParameter<edm::ParameterSet>("jet");
  std::string jetCut = jetPSet.getParameter<std::string>("cut");
  isGoodJet_ = new StringCutObjectSelector<cmg::PFJet, true>(jetCut);
  jetLabel_ = jetPSet.getParameter<edm::InputTag>("src");
  bTagType_ = jetPSet.getParameter<std::string>("bTagType");

  // Output histograms and tree
  edm::Service<TFileService> fs;
  eventTree_ = fs->make<TTree>("event", "Mixed event tree");
  eventTree_->Branch("run", &run_, "run/I");
  eventTree_->Branch("lumi", &lumi_, "lumi/I");
  eventTree_->Branch("event", &event_, "event/I");

  eventTree_->Branch("weight", &weight_, "weight/D");
  eventTree_->Branch("weightUp", &weightUp_, "weightUp/D");
  eventTree_->Branch("weightDn", &weightDn_, "weightDn/D");
  eventTree_->Branch("nVertex", &nVertex_, "nVertex/I"); 

  eventTree_->Branch("electrons", &electrons_);
  eventTree_->Branch("electrons_Q", &electrons_Q_);
  eventTree_->Branch("electrons_Iso", &electrons_Iso_);

  eventTree_->Branch("muons", &muons_);
  eventTree_->Branch("muons_Q", &muons_Q_);
  eventTree_->Branch("muons_Iso", &muons_Iso_);

  eventTree_->Branch("met", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>", &met_);

  eventTree_->Branch("jets", &jets_);
  eventTree_->Branch("jets_bTag", &jets_bTag_);

  genTree_ = fs->make<TTree>("gen", "Gen level event tree");
  genTree_->Branch("muons", &genMuons_);
  genTree_->Branch("electrons", &genElectrons_);
  genTree_->Branch("muonNus", &genMuonNus_);
  genTree_->Branch("electronNus", &genElectronNus_);

  genTree_->Branch("jets_motherId", &jets_motherId_);

}

EventTupleProducer::~EventTupleProducer()
{
}

void EventTupleProducer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace std;

  // Clear up
  electrons_.clear();
  electrons_Q_.clear();
  electrons_Iso_.clear();
  muons_.clear();
  muons_Q_.clear();
  muons_Iso_.clear();

  jets_.clear();
  jets_bTag_.clear();
  jets_motherId_.clear();

  genMuons_.clear();
  genElectrons_.clear();
  genMuonNus_.clear();
  genElectronNus_.clear();

  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByLabel(vertexLabel_, vertexHandle);
  nVertex_ = vertexHandle->size();

  if ( event.isRealData() )
  {
    weight_ = weightUp_ = weightDn_ = 1.0;
  }
  else
  {
    edm::Handle<double> weightHandle, weightUpHandle, weightDnHandle;
    event.getByLabel(edm::InputTag(weightStr_), weightHandle);
    event.getByLabel(edm::InputTag(weightStr_, "up"), weightUpHandle);
    event.getByLabel(edm::InputTag(weightStr_, "dn"), weightDnHandle);
    weight_ = *(weightHandle.product());
    weightUp_ = *(weightUpHandle.product());
    weightDn_ = *(weightDnHandle.product());
  }

  edm::Handle<std::vector<cmg::Electron> > electronHandle;
  event.getByLabel(electronLabel_, electronHandle);
  for ( int i=0, n=electronHandle->size(); i<n; ++i )
  {
    const cmg::Electron& e = electronHandle->at(i);
    if ( !(*isGoodElectron_)(e) ) continue;

    electrons_.push_back(e.p4());
    electrons_Q_.push_back(e.charge());
    electrons_Iso_.push_back(e.relIso(0.5, 0, 0.3)); // Default isolation cone size to be checked
  }

  edm::Handle<std::vector<cmg::Muon> > muonHandle;
  event.getByLabel(muonLabel_, muonHandle);
  for ( int i=0, n=muonHandle->size(); i<n; ++i )
  {
    const cmg::Muon& mu = muonHandle->at(i);
    if ( !(*isGoodMuon_)(mu) ) continue;

    muons_.push_back(mu.p4());
    muons_Q_.push_back(mu.charge());
    muons_Iso_.push_back(mu.relIso(0.5, 0, 0.3)); // Default isolation cone size to be checked
  }

  if ( muons_.size() + electrons_.size() < 2 ) return; // Dilepton channel only

  edm::Handle<std::vector<cmg::BaseMET> > metHandle;
  event.getByLabel(metLabel_, metHandle);
  met_ = metHandle->at(0).p4();

  typedef edm::AssociationMap<edm::OneToMany<std::vector<reco::GenJet>, reco::GenParticleCollection> > GenJetToGenParticlesMap;
  typedef edm::AssociationMap<edm::OneToOne<std::vector<cmg::PFJet>, std::vector<reco::GenJet> > > RecoToGenJetMap;
  edm::Handle<GenJetToGenParticlesMap> genJetToPartonMapHandle;
  edm::Handle<RecoToGenJetMap> recoToGenJetMapHandle;

  // This while loop runs just for once, a "break" statement must be kept in the end of loop
  // It reduces nested loop
  while ( doMCMatch_ ) 
  {
    if ( event.isRealData() ) doMCMatch_ = false;

    edm::Handle<reco::GenParticleCollection> genHandle;
    event.getByLabel(genLabel_, genHandle);
    event.getByLabel(genJetToPartonMapLabel_, genJetToPartonMapHandle);
    event.getByLabel(recoToGenJetMapLabel_, recoToGenJetMapHandle);
    if ( !genHandle.isValid() or
         !genJetToPartonMapHandle.isValid() or 
         !recoToGenJetMapHandle.isValid() )
    {
      doMCMatch_ = false;
      break;
    }

    // Find top quark from the genParticles
    for ( int i=0, n=genHandle->size(); i<n; ++i )
    {
      const reco::GenParticle& p = genHandle->at(i);
      if ( p.status() != 3 ) continue;

      switch(abs(p.pdgId()))
      {
        case 11: genElectrons_.push_back(p.p4())  ; break;
        case 13: genMuons_.push_back(p.p4())      ; break;
        case 12: genElectronNus_.push_back(p.p4()); break;
        case 14: genMuonNus_.push_back(p.p4())    ; break;
        default: break;
      }
    }

    break;
  }

  edm::Handle<std::vector<cmg::PFJet> > jetHandle;
  event.getByLabel(jetLabel_, jetHandle);
  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
    const cmg::PFJet& jet = jetHandle->at(i);

    if ( !(*isGoodJet_)(jet) ) continue;

    jets_.push_back(jet.p4());
    jets_bTag_.push_back(jet.btag(bTagType_.c_str()));

    int jetMotherId = 0;

    while ( doMCMatch_ )
    {
      edm::Ref<std::vector<cmg::PFJet> > jetRef(jetHandle, i);
      RecoToGenJetMap::const_iterator recoToGenJet = recoToGenJetMapHandle->find(jetRef);
      if ( recoToGenJet == recoToGenJetMapHandle->end() ) break;

      const edm::Ref<std::vector<reco::GenJet> >& genJet = recoToGenJet->val;
      GenJetToGenParticlesMap::const_iterator genJetToParton = genJetToPartonMapHandle->find(genJet);
      if ( genJetToParton == genJetToPartonMapHandle->end() ) break;

      const edm::RefVector<reco::GenParticleCollection>& genPartons = genJetToParton->val;
      for ( int j=0, m=genPartons.size(); j<m; ++j )
      {
        const int partonId = genPartons.at(j)->pdgId();
        // NOTE : Maybe there can be better way to set mother parton's id
        if ( abs(partonId) == 24 or partonId == 23 ) jetMotherId = partonId;
        else if ( jetMotherId == 0 and abs(partonId) == 6 ) jetMotherId = partonId;
        else if ( jetMotherId == 0 and partonId == 25 ) jetMotherId = partonId;
      }

      break;
    }

    jets_motherId_.push_back(jetMotherId);

  }

  // Now put jets in current event to the event cache
  run_ = event.run();
  lumi_ = event.luminosityBlock();
  event_ = event.id().event();

  eventTree_->Fill();
  genTree_->Fill();
}

DEFINE_FWK_MODULE(EventTupleProducer);

