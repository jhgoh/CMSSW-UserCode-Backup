#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisDataFormats/CMGTools/interface/BaseJet.h"
#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "AnalysisDataFormats/CMGTools/interface/BaseMET.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <memory>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "AnalysisDataFormats/CMGTools/interface/GenericTypes.h"
#include <TH1F.h>
#include <TH2F.h>

using namespace edm;
using namespace std;

class TopJetSelector : public edm::EDFilter
{
public:
  TopJetSelector(const edm::ParameterSet& pset);
  ~TopJetSelector() {};

private:
  void beginJob() {};
  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);
  void endJob() {};

  bool doFilter_;
  double overlapDeltaR_;
  unsigned int minNumber_;
  unsigned int maxNumber_;

  edm::InputTag jetLabel_, metLabel_;
  std::vector<edm::InputTag> overlapCandLabels_;

  JetCorrectionUncertainty *jecUncCalculator_;
  StringCutObjectSelector<cmg::PFJet, true>* isGoodJet_;
  double minPt_, maxEta_;

  int cleanMethod_;

private:
  bool debug_;
  edm::InputTag genJetLabel_, genLeptonLabel_;

  TH1F* hNPFJet_;
  TH1F* hNGenJet_;
  TH2F* hNGenJetVsNPFJet_;

};

TopJetSelector::TopJetSelector(const edm::ParameterSet& pset)
{
  doFilter_ = pset.getParameter<bool>("doFilter");
  debug_ = pset.getUntrackedParameter<bool>("debug",false);

  jetLabel_ = pset.getParameter<edm::InputTag>("jet");
  metLabel_ = pset.getParameter<edm::InputTag>("met");

  // Selection cuts
  edm::ParameterSet selectionPSet = pset.getParameter<edm::ParameterSet>("selection");
  std::string jetCut = selectionPSet.getParameter<std::string>("cut");
  isGoodJet_ = new StringCutObjectSelector<cmg::PFJet, true>(jetCut);
  minPt_ = selectionPSet.getParameter<double>("minPt");
  maxEta_ = selectionPSet.getParameter<double>("maxEta");

  // Cleaning
  edm::ParameterSet cleanPSet = pset.getParameter<edm::ParameterSet>("cleaning");
  overlapDeltaR_ = cleanPSet.getParameter<double>("overlapDeltaR");
  overlapCandLabels_ = cleanPSet.getParameter<std::vector<edm::InputTag> >("overlapCands");
  const std::string cleanMethodName = cleanPSet.getParameter<std::string>("cleanMethod");
  // Cleaning methods:
  //  subtract    =  1: Check acceptance cut after lepton p4 subtraction and store subtacted jet
  //  subtractAll =  2: Check acceptance cut after lepton p4 subtraction but keep original 4 momentum
  //  cleanAll    =  0: No lepton subtraction
  //  Default     = -1: No jet cleaning
  if ( cleanMethodName == "subtract" ) cleanMethod_ = 2; 
  else if ( cleanMethodName == "subtractAll" ) cleanMethod_ = 1; 
  else if ( cleanMethodName == "cleanAll" ) cleanMethod_ = 0; 
  else cleanMethod_ = -1; 

  minNumber_ = pset.getParameter<unsigned int>("minNumber");
  maxNumber_ = pset.getParameter<unsigned int>("maxNumber");

  produces<std::vector<cmg::PFJet> >("jet");
  produces<std::vector<cmg::PFJet> >("jetUp");
  produces<std::vector<cmg::PFJet> >("jetDn");
  produces<std::vector<cmg::BaseMET> >("met");
  produces<std::vector<cmg::BaseMET> >("metUp");
  produces<std::vector<cmg::BaseMET> >("metDn");

  if ( debug_ )
  {
    genJetLabel_ = pset.getParameter<edm::InputTag>("genJet");
    genLeptonLabel_ = pset.getParameter<edm::InputTag>("genLepton");

    edm::Service<TFileService> fs;
    hNPFJet_ = fs->make<TH1F>("hNPFJet", "nPFJet;Number of PF Jet;Events", 10, 0, 10);
    hNGenJet_ = fs->make<TH1F>("hNGenJet", "nGenJet;Number of Gen Jet;Events", 10, 0, 10);
    hNGenJetVsNPFJet_ = fs->make<TH2F>("hNGenJetVsNPFJet", "nGenJet vs nPFJet;Number of Gen Jet;Number of PF Jet", 10, 0, 10, 10, 0, 10);
  }
}

bool TopJetSelector::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<std::vector<cmg::PFJet> > jetHandle;
  event.getByLabel(jetLabel_, jetHandle);

  edm::Handle<std::vector<cmg::BaseMET> > metHandle;
  event.getByLabel(metLabel_, metHandle);

  std::auto_ptr<std::vector<cmg::PFJet> > corrJets(new std::vector<cmg::PFJet>());
  std::auto_ptr<std::vector<cmg::PFJet> > corrJetsUp(new std::vector<cmg::PFJet>());
  std::auto_ptr<std::vector<cmg::PFJet> > corrJetsDn(new std::vector<cmg::PFJet>());
  std::auto_ptr<std::vector<cmg::BaseMET> > corrMets(new std::vector<cmg::BaseMET>());
  std::auto_ptr<std::vector<cmg::BaseMET> > corrMetsUp(new std::vector<cmg::BaseMET>());
  std::auto_ptr<std::vector<cmg::BaseMET> > corrMetsDn(new std::vector<cmg::BaseMET>());

  double metUpX = 0, metUpY = 0;
  double metDnX = 0, metDnY = 0;

  std::vector<const reco::Candidate*> overlapCands;
  if ( cleanMethod_ != -1 )
  {
    for ( int iLabel=0, nLabel=overlapCandLabels_.size(); iLabel<nLabel; ++iLabel )
    {
      edm::Handle<edm::View<reco::Candidate> > overlapCandHandle;
      event.getByLabel(overlapCandLabels_.at(iLabel), overlapCandHandle);

      //if ( !overlapCandHandle.isValid() ) continue;

      for ( int i=0, n=overlapCandHandle->size(); i<n; ++i )
      {
        overlapCands.push_back(&(overlapCandHandle->at(i)));
      }
    }
  }

  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
    cmg::PFJet jet = jetHandle->at(i);

    if ( !(*isGoodJet_)(jet) ) continue;

    reco::Candidate::LorentzVector jetP4 = jet.p4();
    jecUncCalculator_->setJetPt(jetP4.pt());
    jecUncCalculator_->setJetEta(jetP4.eta());
    const double jecUncUp = jecUncCalculator_->getUncertainty(true);
    const double jecUncDn = jecUncCalculator_->getUncertainty(false);
    reco::Candidate::LorentzVector jetUpP4, jetDnP4;

    jetUpP4 = jetP4*(1+jecUncUp);
    jetDnP4 = jetP4*(1-jecUncDn);
    metUpX += jetP4.px() - jetUpP4.px();
    metUpY += jetP4.py() - jetUpP4.py();
    metDnX += jetP4.px() - jetDnP4.px();
    metDnY += jetP4.py() - jetDnP4.py();

    cmg::PFJet jetUp = jet, jetDn = jet;
    jetUp.setP4(jetUpP4);
    jetDn.setP4(jetDnP4);

    bool isOverlap = false;
    for ( int j=0, m=overlapCands.size(); j<m; ++j )
    {
      if ( deltaR(jet.p4(), overlapCands.at(j)->p4()) < overlapDeltaR_ )
      {
        isOverlap = true;
        if ( cleanMethod_ == 0 ) break;
        else if ( cleanMethod_ == 1 or cleanMethod_ == 2 )
        {
          jetP4 -= overlapCands.at(j)->p4();
          jetUpP4 -= overlapCands.at(j)->p4();
          jetDnP4 -= overlapCands.at(j)->p4();
        }
      }
    }

    if ( isOverlap )
    {
      if ( cleanMethod_ == 0 ) continue;
      else
      {
        if ( cleanMethod_ == 1 )
        {
          jet.setP4(jetP4);
          jetUp.setP4(jetUpP4);
          jetDn.setP4(jetDnP4);
        }
        if ( jetP4.pt() >= minPt_ and abs(jetP4.eta()) <= maxEta_ ) corrJets->push_back(jet);
        if ( jetUpP4.pt() >= minPt_ and abs(jetUpP4.eta()) <= maxEta_ ) corrJetsUp->push_back(jetUp);
        if ( jetDnP4.pt() >= minPt_ and abs(jetDnP4.eta()) <= maxEta_ ) corrJetsDn->push_back(jetDn);
      }
    }
    else 
    {
      corrJets->push_back(jet);
      corrJetsUp->push_back(jetUp);
      corrJetsDn->push_back(jetDn);
    }
  }

  const unsigned int nCleanJet = corrJets->size();
  std::sort(corrJets->begin(), corrJets->end());
  std::sort(corrJetsUp->begin(), corrJetsDn->end());
  std::sort(corrJetsDn->begin(), corrJetsDn->end());

  cmg::BaseMET met = metHandle->at(0);
  cmg::BaseMET metUp, metDn;
  metUp.setP4(reco::Candidate::LorentzVector(metUpX, metUpY, 0, hypot(metUpX, metUpY)));
  metDn.setP4(reco::Candidate::LorentzVector(metDnX, metDnY, 0, hypot(metDnX, metDnY)));
  corrMets->push_back(met);
  corrMetsUp->push_back(metUp);
  corrMetsDn->push_back(metDn);

  event.put(corrJets);
  event.put(corrJetsUp, "jetsUp");
  event.put(corrJetsDn, "jetsDn");
  event.put(corrMets);
  event.put(corrMetsUp, "metsUp");
  event.put(corrMetsDn, "metsDn");

  if ( debug_ and !event.isRealData() )
  {
    edm::Handle<std::vector<cmg::GenJet> > genJetHandle;
    event.getByLabel(genJetLabel_, genJetHandle);

    edm::Handle<std::vector<reco::GenParticle> > genLeptonHandle;
    event.getByLabel(genLeptonLabel_, genLeptonHandle);

    int nGenJet = 0;
    for ( int i=0, n=genJetHandle->size(); i<n; ++i )
    {
      const cmg::GenJet& genJet = genJetHandle->at(i);
      if ( genJet.pt() < minPt_ or abs(genJet.eta()) > maxEta_ ) continue;
      const std::vector<const reco::GenParticle*> genConstituents = genJet.sourcePtr()->get()->getGenConstituents();
      if ( genConstituents.size() <= 1 ) continue;
      bool isOverlap = false;
      for ( int j=0, m=genConstituents.size(); j<m; ++j )
      {
        const reco::Candidate* p = genConstituents.at(j);
        if ( abs(p->pdgId()) != 11 and abs(p->pdgId()) != 13 ) continue;

        while ( (p = p->mother()) != 0 )
        {
          if ( p->status() != 3 ) continue;
          if ( abs(p->pdgId()) == 11 or abs(p->pdgId()) == 13 )
          {
            isOverlap = true;
            break;
          }
        }
      }
      if ( isOverlap ) continue;
      
      for ( int j=0, m=genLeptonHandle->size(); j<m; ++j )
      {
        const reco::GenParticle& genLepton = genLeptonHandle->at(j);
        if ( deltaR(genJet, genLepton) > 0.05 ) continue;
        //if ( abs(genJet.pt()-genLepton.pt())/genJet.pt() < 0.1 ) continue;
        //cout << genJet.pt() << "   " << genLepton.pt() << " " << deltaR(genJet, genLepton) << endl;
        isOverlap = true;
        break;
      }
      ++nGenJet;
    }
    hNPFJet_->Fill(nCleanJet);
    hNGenJet_->Fill(nGenJet);
    hNGenJetVsNPFJet_->Fill(nGenJet, nCleanJet);
  }

  if ( !doFilter_ ) return true;
  else if ( nCleanJet >= minNumber_ and nCleanJet <= maxNumber_ ) return true;

  return false;
}

DEFINE_FWK_MODULE(TopJetSelector);

