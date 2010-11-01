#include "HLTrigger/TPGAnalysis/interface/JetHLTAnalyzer.h"
#include "HLTrigger/TPGAnalysis/interface/Histograms.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
    
#include "DataFormats/Math/interface/deltaR.h"

#include <TString.h>
#include <memory>
#include <iostream>

using namespace std;

template<typename T>
const T* JetHLTAnalyzer::getBestMatch(const reco::Candidate& recoCand, std::vector<T>& candColl)
{
  const T* matchedCand = 0;
  double matchedDeltaR = 1e14;

  for ( typename std::vector<T>::const_iterator cand = candColl.begin();
        cand != candColl.end(); ++cand )
  {
    const double candDeltaR = deltaR(recoCand, *cand);

    if ( candDeltaR < matchedDeltaR )
    {
      matchedCand = &(*cand);
      matchedDeltaR = candDeltaR;
    }
  }

  return matchedCand;
}

JetHLTAnalyzer::JetHLTAnalyzer(const edm::ParameterSet& pset)
{
  interestedFilterName_ = pset.getParameter<std::string>("interestedFilterName");
  recoJetTag_ = pset.getParameter<edm::InputTag>("recoJet");

  edm::ParameterSet jetCutSet = pset.getParameter<edm::ParameterSet>("cut");

  recoMinEt_ = jetCutSet.getParameter<double>("recoMinEt");
  l1MinEt_ = jetCutSet.getParameter<double>("l1MinEt");
  maxL1DeltaR_ = jetCutSet.getParameter<double>("maxL1DeltaR");
  const double workingPointEt = jetCutSet.getParameter<double>("workingPointEt");
  const double maxHLTDeltaR = jetCutSet.getParameter<double>("maxHLTDeltaR");
  const double maxL1DeltaR = jetCutSet.getParameter<double>("maxL1DeltaR");

  const int objectType = Histograms::ObjectType::Jet;

  hAll_ = new HTrigger("All", "All", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);
  hCentral_ = new HTrigger("Central", "Central", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);
  hForward_ = new HTrigger("Forward", "Forward", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);

  hLeadingAll_ = new HTrigger("All1", "All leading", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);
  hLeadingCentral_ = new HTrigger("Central1", "Central leading", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);
  hLeadingForward_ = new HTrigger("Forward1", "Forward leading", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);
}

JetHLTAnalyzer::~JetHLTAnalyzer()
{
}

void JetHLTAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
}

void JetHLTAnalyzer::endRun()
{
}

void JetHLTAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<l1extra::L1JetParticleCollection> centralL1JetHandle;
  edm::Handle<l1extra::L1JetParticleCollection> forwardL1JetHandle;

  if ( !event.getByLabel(edm::InputTag("l1extraParticles", "Central"), centralL1JetHandle) )
  {
    edm::LogError("JetHLTAnalyzer") << "Cannot find central l1extra Jets\n";
    return;
  }
  if ( !event.getByLabel(edm::InputTag("l1extraParticles", "Forward"), forwardL1JetHandle) )
  {
    edm::LogError("JetHLTAnalyzer") << "Cannot find forward l1extra Jets\n";
    return;
  }

  edm::Handle<trigger::TriggerEvent> triggerEventHandle;
  if ( !event.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", "HLT"), triggerEventHandle) )
  {
    edm::LogError("JetHLTAnalyzer") << "Cannot find TriggerEvent\n";
    return;
  }
  const trigger::TriggerObjectCollection& allTriggerObjects = triggerEventHandle->getObjects();

  edm::Handle<edm::View<reco::CaloJet> > recoJetHandle;
  if ( !event.getByLabel(recoJetTag_, recoJetHandle) )
  {
    edm::LogError("JetHLTAnalyzer") << "Cannot find reco Jets\n";
    return;
  }

  edm::Handle<reco::JetIDValueMap> jetIDMapHandle;
  if ( !event.getByLabel("ak5JetID", jetIDMapHandle) )
  {
    edm::LogError("JetHLTAnalyzer") << "Cannot find Jet ID map\n";
    return; 
  }
  JetIDSelectionFunctor jetIDSelector(JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::LOOSE);

  // Initialize histogram objects
  hAll_->init(event.id());
  hCentral_->init(event.id());
  hForward_->init(event.id());

  hLeadingAll_->init(event.id());
  hLeadingCentral_->init(event.id());
  hLeadingForward_->init(event.id());

  // Collect L1 objects
  l1extra::L1JetParticleCollection centralL1Jets;
  l1extra::L1JetParticleCollection forwardL1Jets;

  centralL1Jets.reserve(centralL1JetHandle->size());
  forwardL1Jets.reserve(forwardL1JetHandle->size());

  for ( l1extra::L1JetParticleCollection::const_iterator l1Jet = centralL1JetHandle->begin();
        l1Jet != centralL1JetHandle->end(); ++l1Jet )
  {
    if ( l1Jet->et() < l1MinEt_ ) continue;
    centralL1Jets.push_back(*l1Jet);
  }

  for ( l1extra::L1JetParticleCollection::const_iterator l1Jet = forwardL1JetHandle->begin();
        l1Jet != forwardL1JetHandle->end(); ++l1Jet )
  {
    if ( l1Jet->et() < l1MinEt_ ) continue;
    forwardL1Jets.push_back(*l1Jet);
  }

  // Collect HLT objects
  trigger::TriggerObjectCollection triggerObjects;
  for ( unsigned int filterIdx = 0; filterIdx < triggerEventHandle->sizeFilters(); ++filterIdx )
  {
    const std::string filterFullName = triggerEventHandle->filterTag(filterIdx).encode();
    const size_t fsPos = filterFullName.find_first_of(':');
    const std::string filterName = ( fsPos == std::string::npos ) ? filterFullName : filterFullName.substr(0, fsPos);

    if ( filterName != interestedFilterName_ ) continue;

    const trigger::Keys& trgKeys = triggerEventHandle->filterKeys(filterIdx);
    for ( trigger::Keys::const_iterator trgKey = trgKeys.begin();
          trgKey != trgKeys.end(); ++trgKey )
    {
      triggerObjects.push_back(allTriggerObjects[*trgKey]);
    }
  }

  const reco::CaloJet* leadingJet = 0;
  const reco::CaloJet* leadingCentralJet = 0;
  const reco::CaloJet* leadingForwardJet = 0;
  for ( edm::View<reco::CaloJet>::const_iterator recoJet = recoJetHandle->begin();
        recoJet != recoJetHandle->end(); ++recoJet )
  {
    // Do some basic cuts
    //if ( !isGoodJet(*recoJet, event) ) continue;
    edm::RefToBase<reco::CaloJet> recoJetRef = recoJetHandle->refAt(recoJet-recoJetHandle->begin());
    reco::JetID const& jetID = jetIDMapHandle->operator[](recoJetRef);
    if ( !jetIDSelector(*recoJet, jetID) ) continue;

    const double recoJetAbsEta = fabs(recoJet->eta());
    const double recoJetEt = recoJet->et();

    // All jets
    const reco::CaloJet* recoJetP = &*recoJet;
    const l1extra::L1JetParticle* matchedL1 = getBestMatch(*recoJet, recoJetAbsEta < 3.0 ? centralL1Jets : forwardL1Jets); 
    const trigger::TriggerObject* matchedHLT = getBestMatch(*recoJet, triggerObjects);

    const double l1DeltaR = deltaR(*recoJet, *matchedL1);
    const double hltDeltaR = deltaR(*recoJet, *matchedHLT);

    if ( !leadingJet or leadingJet->et() < recoJetEt ) leadingJet = recoJetP;
    if ( recoJetAbsEta < 3.0 and (!leadingCentralJet or leadingCentralJet->et() < recoJetEt ) ) leadingCentralJet = recoJetP;
    if ( recoJetAbsEta >= 3.0 and (!leadingForwardJet or leadingForwardJet->et() < recoJetEt ) ) leadingForwardJet = recoJetP;

    hAll_->fill(recoJetP, matchedL1, matchedHLT);
    if ( recoJetAbsEta < 3.0 ) hCentral_->fill(recoJetP, matchedL1, matchedHLT);
    else hForward_->fill(recoJetP, matchedL1, matchedHLT);
  }

  if ( leadingJet )
  {
    const double recoJetAbsEta = fabs(leadingJet->eta());
    const l1extra::L1JetParticle* matchedL1 = getBestMatch(*leadingJet, recoJetAbsEta < 3.0 ? centralL1Jets : forwardL1Jets);
    const trigger::TriggerObject* matchedHLT = getBestMatch(*leadingJet, triggerObjects);
    hLeadingAll_->fill(leadingJet, matchedL1, matchedHLT);
  }

  if ( leadingCentralJet )
  {
    const l1extra::L1JetParticle* matchedL1 = getBestMatch(*leadingCentralJet, centralL1Jets);
    const trigger::TriggerObject* matchedHLT = getBestMatch(*leadingCentralJet, triggerObjects);
    hLeadingCentral_->fill(leadingCentralJet, matchedL1, matchedHLT);
  }

  if ( leadingForwardJet )
  {
    const l1extra::L1JetParticle* matchedL1 = getBestMatch(*leadingForwardJet, forwardL1Jets);
    const trigger::TriggerObject* matchedHLT = getBestMatch(*leadingForwardJet, triggerObjects);
    hForward_->fill(leadingForwardJet, matchedL1, matchedHLT);
  }

}

bool JetHLTAnalyzer::isGoodJet(const reco::CaloJet& recoJet, const edm::Event& event)
{
  return true;
}
