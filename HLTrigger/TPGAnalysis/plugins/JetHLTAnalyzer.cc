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

  jetCutSet_ = pset.getParameter<edm::ParameterSet>("cut");

  recoMinEt_ = jetCutSet_.getParameter<double>("recoMinEt");
  l1MinEt_ = jetCutSet_.getParameter<double>("l1MinEt");
  maxL1DeltaR_ = jetCutSet_.getParameter<double>("maxL1DeltaR");

  // Book run independent histograms
  edm::Service<TFileService> fs;
  const int objectType = Histograms::ObjectType::Jet;

  TFileDirectory allJetDir = fs->mkdir("All");
  TFileDirectory centralJetDir = fs->mkdir("Central");
  TFileDirectory forwardJetDir = fs->mkdir("Forward");

  TFileDirectory allLeadingJetDir = fs->mkdir("AllLeading");
  TFileDirectory centralLeadingJetDir = fs->mkdir("CentralLeading");
  TFileDirectory forwardLeadingJetDir = fs->mkdir("ForwardLeading");

  hAllJet_AllRun_ = new Histograms(allJetDir, "All", jetCutSet_, objectType);
  hCentralJet_AllRun_ = new Histograms(centralJetDir, "Central", jetCutSet_, objectType);
  hForwardJet_AllRun_ = new Histograms(forwardJetDir, "Forward", jetCutSet_, objectType);

  hAllLeadingJet_AllRun_ = new Histograms(allLeadingJetDir, "All Leading", jetCutSet_, objectType);
  hCentralLeadingJet_AllRun_ = new Histograms(centralLeadingJetDir, "Central Leading", jetCutSet_, objectType);
  hForwardLeadingJet_AllRun_ = new Histograms(forwardLeadingJetDir, "Forward Leading", jetCutSet_, objectType);  

  // Special histograms
  hNReco_ = fs->make<TH1F>("hNReco", "Number of reco object", 6, -0.5, 5.5);
  hNCentralL1T_ = fs->make<TH1F>("hNCentralL1T", "Number of matched Central L1T object", 6, -0.5, 5.5);
  hNForwardL1T_ = fs->make<TH1F>("hNForwardL1T", "Number of matched Forward L1T object", 6, -0.5, 5.5);
  hNCentralHLT_ = fs->make<TH1F>("hNCentralHLT", "Number of matched Central HLT object", 6, -0.5, 5.5);
  hNForwardHLT_ = fs->make<TH1F>("hNForwardHLT", "Number of matched Forward HLT object", 6, -0.5, 5.5);
  hNDuplicatedL1T_ = fs->make<TH1F>("hNDuplicatedL1T", "Number of duplicated L1T-reco matching in Central L1T - Forward L1T", 6, -0.5, 5.5);

  TFileDirectory allJetNoL1Dir = fs->mkdir("AllJetNoL1");
  TFileDirectory centralJetNoL1Dir = fs->mkdir("CentralJetNoL1");
  TFileDirectory forwardJetNoL1Dir = fs->mkdir("ForwardJetNoL1");

  hAllJetNoL1_ = new Histograms(allJetNoL1Dir, "All Jets no L1 matching", jetCutSet_, objectType);
  hCentralJetNoL1_ = new Histograms(centralJetNoL1Dir, "Central Jets no L1 matching", jetCutSet_, objectType);
  hForwardJetNoL1_ = new Histograms(forwardJetNoL1Dir, "Forward Jets no L1 matching", jetCutSet_, objectType);
}

JetHLTAnalyzer::~JetHLTAnalyzer()
{
}

void JetHLTAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  const int runNumber = run.run();

  if ( hCentralJet_ByRun_.find(runNumber) == hCentralJet_ByRun_.end() )
  {
    edm::Service<TFileService> fs;

    TFileDirectory runDir = fs->mkdir(Form("Run %d", runNumber));

    TFileDirectory allJetDir = runDir.mkdir("All");
    TFileDirectory centralJetDir = runDir.mkdir("Central");
    TFileDirectory forwardJetDir = runDir.mkdir("Forward");

    TFileDirectory allLeadingJetDir = runDir.mkdir("AllLeading");
    TFileDirectory centralLeadingJetDir = runDir.mkdir("CentralLeading");
    TFileDirectory forwardLeadingJetDir = runDir.mkdir("ForwardLeading");

    const int objectType = Histograms::ObjectType::Jet;

    hAllJet_ByRun_[runNumber] = new Histograms(allJetDir, "All", jetCutSet_, objectType);
    hCentralJet_ByRun_[runNumber] = new Histograms(centralJetDir, "Central", jetCutSet_, objectType);
    hForwardJet_ByRun_[runNumber] = new Histograms(forwardJetDir, "Forward", jetCutSet_, objectType);

    hAllLeadingJet_ByRun_[runNumber] = new Histograms(allLeadingJetDir, "All Leading", jetCutSet_, objectType);
    hCentralLeadingJet_ByRun_[runNumber] = new Histograms(centralLeadingJetDir, "Central Leading", jetCutSet_, objectType);
    hForwardLeadingJet_ByRun_[runNumber] = new Histograms(forwardLeadingJetDir, "Forward Leading", jetCutSet_, objectType);   
  }

  hAllJet_ = hAllJet_ByRun_[runNumber];
  hCentralJet_ = hCentralJet_ByRun_[runNumber];
  hForwardJet_ = hForwardJet_ByRun_[runNumber];

  hAllLeadingJet_ = hAllLeadingJet_ByRun_[runNumber];
  hCentralLeadingJet_ = hCentralLeadingJet_ByRun_[runNumber];
  hForwardLeadingJet_ = hForwardLeadingJet_ByRun_[runNumber];
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
  hAllJet_->init();
  hCentralJet_->init();
  hForwardJet_->init();

  hAllLeadingJet_->init();
  hCentralLeadingJet_->init();
  hForwardLeadingJet_->init();  

  hAllJet_AllRun_->init();
  hCentralJet_AllRun_->init();
  hForwardJet_AllRun_->init();

  hAllLeadingJet_AllRun_->init();
  hCentralLeadingJet_AllRun_->init();
  hForwardLeadingJet_AllRun_->init();

  hAllJetNoL1_->init();
  hCentralJetNoL1_->init();
  hForwardJetNoL1_->init();

  // Loop over all reco jets
  int nReco = 0;
  int nDuplicatedL1T = 0;
  int nCentralL1T = 0, nForwardL1T = 0;
  int nCentralHLT = 0, nForwardHLT = 0;

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
  const reco::CaloJet* centralLeadingJet = 0;
  const reco::CaloJet* forwardLeadingJet = 0;
  for ( edm::View<reco::CaloJet>::const_iterator recoJet = recoJetHandle->begin();
        recoJet != recoJetHandle->end(); ++recoJet )
  {
    // Do some basic cuts
    //if ( !isGoodJet(*recoJet, event) ) continue;
    edm::RefToBase<reco::CaloJet> recoJetRef = recoJetHandle->refAt(recoJet-recoJetHandle->begin());
    reco::JetID const& jetID = jetIDMapHandle->operator[](recoJetRef);
    if ( !jetIDSelector(*recoJet, jetID) ) continue;

    ++nReco;
    const double recoJetAbsEta = fabs(recoJet->eta());
    const double recoJetEt = recoJet->et();

    // Fill basic information of recoJet
    // and find the leading Jet
    hAllJet_->setRecoCand(&*recoJet);
    hAllJetNoL1_->setRecoCand(&*recoJet);
    hAllJet_AllRun_->setRecoCand(&*recoJet);
    if ( !leadingJet or leadingJet->et() < recoJetEt ) leadingJet = &(*recoJet);
    if ( recoJetAbsEta < 3.0 )
    {
      hCentralJet_->setRecoCand(&*recoJet);
      hCentralJetNoL1_->setRecoCand(&*recoJet);
      hCentralJet_AllRun_->setRecoCand(&*recoJet);
      if ( !centralLeadingJet or centralLeadingJet->et() < recoJetEt ) centralLeadingJet = &(*recoJet);
    }
    else
    {
      hForwardJet_->setRecoCand(&*recoJet);
      hForwardJetNoL1_->setRecoCand(&*recoJet);
      hForwardJet_AllRun_->setRecoCand(&*recoJet);
      if ( !forwardLeadingJet or forwardLeadingJet->et() < recoJetEt ) forwardLeadingJet = &(*recoJet);
    }

    // Try matching
    const trigger::TriggerObject* matchedHLTJet = getBestMatch(*recoJet, triggerObjects);

    if ( recoJetAbsEta < 3.0 ) 
    {
      // Barrel region
      const l1extra::L1JetParticle* matchedL1Jet = getBestMatch(*recoJet, centralL1Jets);
      if ( matchedL1Jet )
      {
        ++nCentralL1T;

        hAllJet_->setL1Cand(matchedL1Jet);
        hCentralJet_->setL1Cand(matchedL1Jet);

        hAllJet_AllRun_->setL1Cand(matchedL1Jet);
        hCentralJet_AllRun_->setL1Cand(matchedL1Jet);

        if ( matchedHLTJet and deltaR(*recoJet, *matchedL1Jet) < maxL1DeltaR_ )
        {
          ++nCentralHLT;

          hAllJet_->setHLTCand(matchedHLTJet);
          hCentralJet_->setHLTCand(matchedHLTJet);

          hAllJet_AllRun_->setHLTCand(matchedHLTJet);
          hCentralJet_AllRun_->setHLTCand(matchedHLTJet);
        }
      }
    }
    else
    {
      // Forward region. They are comming from HF
      const l1extra::L1JetParticle* matchedL1Jet = getBestMatch(*recoJet, forwardL1Jets);
      if ( matchedL1Jet )
      {
        ++nForwardL1T;

        hAllJet_->setL1Cand(matchedL1Jet);
        hForwardJet_->setL1Cand(matchedL1Jet);

        hAllJet_AllRun_->setL1Cand(matchedL1Jet);
        hForwardJet_AllRun_->setL1Cand(matchedL1Jet);

        if ( matchedHLTJet and deltaR(*recoJet, *matchedL1Jet) < maxL1DeltaR_ )
        {
          ++nForwardHLT;

          hAllJet_->setHLTCand(matchedHLTJet);
          hForwardJet_->setHLTCand(matchedHLTJet);

          hAllJet_AllRun_->setHLTCand(matchedHLTJet);
          hForwardJet_AllRun_->setHLTCand(matchedHLTJet);
        }
      }
    }

    if ( matchedHLTJet )
    {
      hAllJetNoL1_->setHLTCand(matchedHLTJet);
      if ( recoJetAbsEta < 3 ) hCentralJetNoL1_->setHLTCand(matchedHLTJet);
      else hForwardJetNoL1_->setHLTCand(matchedHLTJet);
    }
  }

  // We found leading jets. do the trigger object matching
  if ( leadingJet )
  {
    hAllLeadingJet_->setRecoCand(&*leadingJet);
    hAllLeadingJet_AllRun_->setRecoCand(&*leadingJet);

    // Retry matching
    const l1extra::L1JetParticle* matchedCentralL1Jet = getBestMatch(*leadingJet, centralL1Jets);
    const l1extra::L1JetParticle* matchedForwardL1Jet = getBestMatch(*leadingJet, forwardL1Jets);

    if ( matchedCentralL1Jet )
    {
      hAllLeadingJet_->setL1Cand(matchedCentralL1Jet);
      hAllLeadingJet_AllRun_->setL1Cand(matchedCentralL1Jet);

      const trigger::TriggerObject* matchedHLTJet = getBestMatch(*leadingJet, triggerObjects);
      const double l1DeltaR = deltaR(*leadingJet, *matchedCentralL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hAllLeadingJet_->setHLTCand(matchedHLTJet);
        hAllLeadingJet_AllRun_->setHLTCand(matchedHLTJet);
      }
    }
    if ( matchedForwardL1Jet )
    {
      hAllLeadingJet_->setL1Cand(matchedForwardL1Jet);
      hAllLeadingJet_AllRun_->setL1Cand(matchedForwardL1Jet);

      const trigger::TriggerObject* matchedHLTJet = getBestMatch(*leadingJet, triggerObjects);
      const double l1DeltaR = deltaR(*leadingJet, *matchedForwardL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hAllLeadingJet_->setHLTCand(matchedHLTJet);
        hAllLeadingJet_AllRun_->setHLTCand(matchedHLTJet);
      }
    }
  }

  if ( centralLeadingJet )
  {
    hCentralLeadingJet_->setRecoCand(&*centralLeadingJet);
    hCentralLeadingJet_AllRun_->setRecoCand(&*centralLeadingJet);

    // Retry matching
    const l1extra::L1JetParticle* matchedCentralL1Jet = getBestMatch(*centralLeadingJet, centralL1Jets);

    if ( matchedCentralL1Jet )
    {
      hCentralLeadingJet_->setL1Cand(matchedCentralL1Jet);
      hCentralLeadingJet_AllRun_->setL1Cand(matchedCentralL1Jet);

      const trigger::TriggerObject* matchedHLTJet = getBestMatch(*centralLeadingJet, triggerObjects);
      const double l1DeltaR = deltaR(*centralLeadingJet, *matchedCentralL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hCentralLeadingJet_->setHLTCand(matchedHLTJet);
        hCentralLeadingJet_AllRun_->setHLTCand(matchedHLTJet);
      }
    }
  }

  if ( forwardLeadingJet )
  {
    hForwardLeadingJet_->setRecoCand(&*forwardLeadingJet);
    hForwardLeadingJet_AllRun_->setRecoCand(&*forwardLeadingJet);

    // Retry matching
    const l1extra::L1JetParticle* matchedForwardL1Jet = getBestMatch(*forwardLeadingJet, forwardL1Jets);

    if ( matchedForwardL1Jet )
    {
      hForwardLeadingJet_->setL1Cand(matchedForwardL1Jet);
      hForwardLeadingJet_AllRun_->setL1Cand(matchedForwardL1Jet);

      const trigger::TriggerObject* matchedHLTJet = getBestMatch(*forwardLeadingJet, triggerObjects);
      const double l1DeltaR = deltaR(*forwardLeadingJet, *matchedForwardL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hForwardLeadingJet_->setHLTCand(matchedHLTJet);
        hForwardLeadingJet_AllRun_->setHLTCand(matchedHLTJet);
      }
    }
  }

  hAllJet_->fill();
  hCentralJet_->fill();
  hForwardJet_->fill();

  hAllLeadingJet_->fill();
  hCentralLeadingJet_->fill();
  hForwardLeadingJet_->fill();  

  hAllJet_AllRun_->fill();
  hCentralJet_AllRun_->fill();
  hForwardJet_AllRun_->fill();

  hAllLeadingJet_AllRun_->fill();
  hCentralLeadingJet_AllRun_->fill();
  hForwardLeadingJet_AllRun_->fill();

  hAllJetNoL1_->fill();
  hCentralJetNoL1_->fill();
  hForwardJetNoL1_->fill();

  hNReco_->Fill(nReco);

  hNDuplicatedL1T_->Fill(nDuplicatedL1T);

  hNCentralL1T_->Fill(nCentralL1T);
  hNForwardL1T_->Fill(nForwardL1T);

  hNCentralHLT_->Fill(nCentralHLT);
  hNForwardHLT_->Fill(nForwardHLT);
}

bool JetHLTAnalyzer::isGoodJet(const reco::CaloJet& recoJet, const edm::Event& event)
{
  return true;
}

