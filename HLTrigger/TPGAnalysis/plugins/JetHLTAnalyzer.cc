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
  TFileDirectory overlapJetDir = fs->mkdir("Overlap");
  TFileDirectory forwardJetDir = fs->mkdir("Forward");

  TFileDirectory allLeadingJetDir = fs->mkdir("AllLeading");
  TFileDirectory centralLeadingJetDir = fs->mkdir("CentralLeading");
  TFileDirectory overlapLeadingJetDir = fs->mkdir("OverlapLeading");
  TFileDirectory forwardLeadingJetDir = fs->mkdir("ForwardLeading");

  hAllJet_AllRun_ = new Histograms(allJetDir, "All", jetCutSet_, objectType);
  hCentralJet_AllRun_ = new Histograms(centralJetDir, "Central", jetCutSet_, objectType);
  hOverlapJet_AllRun_ = new Histograms(overlapJetDir, "Overlap", jetCutSet_, objectType);
  hForwardJet_AllRun_ = new Histograms(forwardJetDir, "Forward", jetCutSet_, objectType);

  hAllLeadingJet_AllRun_ = new Histograms(allLeadingJetDir, "All Leading", jetCutSet_, objectType);
  hCentralLeadingJet_AllRun_ = new Histograms(centralLeadingJetDir, "Central Leading", jetCutSet_, objectType);
  hOverlapLeadingJet_AllRun_ = new Histograms(overlapLeadingJetDir, "Overlap Leading", jetCutSet_, objectType);
  hForwardLeadingJet_AllRun_ = new Histograms(forwardLeadingJetDir, "Forward Leading", jetCutSet_, objectType);  

  // Special histograms
  hNReco_ = fs->make<TH1F>("hNReco", "Number of reco object", 6, -0.5, 5.5);
  hNCentralL1T_ = fs->make<TH1F>("hNCentralL1T", "Number of matched Central L1T object", 6, -0.5, 5.5);
  hNOverlapL1T_ = fs->make<TH1F>("hNOverlapL1T", "Number of matched Overlap L1T object", 6, -0.5, 5.5);
  hNForwardL1T_ = fs->make<TH1F>("hNForwardL1T", "Number of matched Forward L1T object", 6, -0.5, 5.5);
  hNCentralHLT_ = fs->make<TH1F>("hNCentralHLT", "Number of matched Central HLT object", 6, -0.5, 5.5);
  hNOverlapHLT_ = fs->make<TH1F>("hNOverlapHLT", "Number of matched Overlap HLT object", 6, -0.5, 5.5);
  hNForwardHLT_ = fs->make<TH1F>("hNForwardHLT", "Number of matched Forward HLT object", 6, -0.5, 5.5);
  hNDuplicatedL1T_ = fs->make<TH1F>("hNDuplicatedL1T", "Number of duplicated L1T-reco matching in Central L1T - Forward L1T", 6, -0.5, 5.5);

  TFileDirectory allJetNoL1Dir = fs->mkdir("AllJetNoL1");
  TFileDirectory centralJetNoL1Dir = fs->mkdir("CentralJetNoL1");
  TFileDirectory overlapJetNoL1Dir = fs->mkdir("OverlapJetNoL1");
  TFileDirectory forwardJetNoL1Dir = fs->mkdir("ForwardJetNoL1");

  hAllJetNoL1_ = new Histograms(allJetNoL1Dir, "All Jets no L1 matching", jetCutSet_, objectType);
  hCentralJetNoL1_ = new Histograms(centralJetNoL1Dir, "Central Jets no L1 matching", jetCutSet_, objectType);
  hOverlapJetNoL1_ = new Histograms(overlapJetNoL1Dir, "Overlap Jets no L1 matching", jetCutSet_, objectType);
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
    TFileDirectory overlapJetDir = runDir.mkdir("Overlap");
    TFileDirectory forwardJetDir = runDir.mkdir("Forward");

    TFileDirectory allLeadingJetDir = runDir.mkdir("AllLeading");
    TFileDirectory centralLeadingJetDir = runDir.mkdir("CentralLeading");
    TFileDirectory overlapLeadingJetDir = runDir.mkdir("OverlapLeading");
    TFileDirectory forwardLeadingJetDir = runDir.mkdir("ForwardLeading");

    const int objectType = Histograms::ObjectType::Jet;

    hAllJet_ByRun_[runNumber] = new Histograms(allJetDir, "All", jetCutSet_, objectType);
    hCentralJet_ByRun_[runNumber] = new Histograms(centralJetDir, "Central", jetCutSet_, objectType);
    hOverlapJet_ByRun_[runNumber] = new Histograms(overlapJetDir, "Overlap", jetCutSet_, objectType);
    hForwardJet_ByRun_[runNumber] = new Histograms(forwardJetDir, "Forward", jetCutSet_, objectType);

    hAllLeadingJet_ByRun_[runNumber] = new Histograms(allLeadingJetDir, "All Leading", jetCutSet_, objectType);
    hCentralLeadingJet_ByRun_[runNumber] = new Histograms(centralLeadingJetDir, "Central Leading", jetCutSet_, objectType);
    hOverlapLeadingJet_ByRun_[runNumber] = new Histograms(overlapLeadingJetDir, "Overlap Leading", jetCutSet_, objectType);
    hForwardLeadingJet_ByRun_[runNumber] = new Histograms(forwardLeadingJetDir, "Forward Leading", jetCutSet_, objectType);   
  }

  hAllJet_ = hAllJet_ByRun_[runNumber];
  hCentralJet_ = hCentralJet_ByRun_[runNumber];
  hOverlapJet_ = hOverlapJet_ByRun_[runNumber];
  hForwardJet_ = hForwardJet_ByRun_[runNumber];

  hAllLeadingJet_ = hAllLeadingJet_ByRun_[runNumber];
  hCentralLeadingJet_ = hCentralLeadingJet_ByRun_[runNumber];
  hOverlapLeadingJet_ = hOverlapLeadingJet_ByRun_[runNumber];
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

  // Loop over all reco jets
  int nReco = 0;
  int nDuplicatedL1T = 0;
  int nCentralL1T = 0, nOverlapL1T = 0, nForwardL1T = 0;
  int nCentralHLT = 0, nOverlapHLT = 0, nForwardHLT = 0;

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
  const reco::CaloJet* overlapLeadingJet = 0;
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
    hAllJet_->FillReco(*recoJet);
    hAllJetNoL1_->FillReco(*recoJet);
    hAllJet_AllRun_->FillReco(*recoJet);
    if ( !leadingJet or leadingJet->et() < recoJetEt ) leadingJet = &(*recoJet);
    if ( recoJetAbsEta < 2.5 )
    {
      hCentralJet_->FillReco(*recoJet);
      hCentralJetNoL1_->FillReco(*recoJet);
      hCentralJet_AllRun_->FillReco(*recoJet);
      if ( !centralLeadingJet or centralLeadingJet->et() < recoJetEt ) centralLeadingJet = &(*recoJet);
    }
    else if ( recoJetAbsEta < 3.0 )
    {
      hOverlapJet_->FillReco(*recoJet);
      hOverlapJetNoL1_->FillReco(*recoJet);
      hOverlapJet_AllRun_->FillReco(*recoJet);
      if ( !overlapLeadingJet or overlapLeadingJet->et() < recoJetEt ) overlapLeadingJet = &(*recoJet);
    }
    else
    {
      hForwardJet_->FillReco(*recoJet);
      hForwardJetNoL1_->FillReco(*recoJet);
      hForwardJet_AllRun_->FillReco(*recoJet);
      if ( !forwardLeadingJet or forwardLeadingJet->et() < recoJetEt ) forwardLeadingJet = &(*recoJet);
    }

    // Try matching
    const trigger::TriggerObject* matchedHLTJet = getBestMatch(*recoJet, triggerObjects);

    if ( recoJetAbsEta < 2.5 ) 
    {
      // Barrel region
      const l1extra::L1JetParticle* matchedL1Jet = getBestMatch(*recoJet, centralL1Jets);
      if ( matchedL1Jet )
      {
        ++nCentralL1T;

        hAllJet_->FillL1T(*recoJet, *matchedL1Jet);
        hCentralJet_->FillL1T(*recoJet, *matchedL1Jet);

        hAllJet_AllRun_->FillL1T(*recoJet, *matchedL1Jet);
        hCentralJet_AllRun_->FillL1T(*recoJet, *matchedL1Jet);

        if ( matchedHLTJet and deltaR(*recoJet, *matchedL1Jet) < maxL1DeltaR_ )
        {
          ++nCentralHLT;

          hAllJet_->FillHLT(*recoJet, *matchedHLTJet);
          hCentralJet_->FillHLT(*recoJet, *matchedHLTJet);

          hAllJet_AllRun_->FillHLT(*recoJet, *matchedHLTJet);
          hCentralJet_AllRun_->FillHLT(*recoJet, *matchedHLTJet);
        }
      }
    }
    else if ( recoJetAbsEta < 3 )
    {
      // Overlap region. The L1 objects are still comming from "Central" collection
      const l1extra::L1JetParticle* matchedL1Jet = getBestMatch(*recoJet, centralL1Jets);
      if ( matchedL1Jet )
      {
        ++nOverlapL1T;

        hAllJet_->FillL1T(*recoJet, *matchedL1Jet);
        hOverlapJet_->FillL1T(*recoJet, *matchedL1Jet);

        hAllJet_AllRun_->FillL1T(*recoJet, *matchedL1Jet);
        hOverlapJet_AllRun_->FillL1T(*recoJet, *matchedL1Jet);

        if ( matchedHLTJet and deltaR(*recoJet, *matchedL1Jet) < maxL1DeltaR_ )
        {
          ++nOverlapHLT;

          hAllJet_->FillHLT(*recoJet, *matchedHLTJet);
          hOverlapJet_->FillHLT(*recoJet, *matchedHLTJet);

          hAllJet_AllRun_->FillHLT(*recoJet, *matchedHLTJet);
          hOverlapJet_AllRun_->FillHLT(*recoJet, *matchedHLTJet);
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

        hAllJet_->FillL1T(*recoJet, *matchedL1Jet);
        hForwardJet_->FillL1T(*recoJet, *matchedL1Jet);

        hAllJet_AllRun_->FillL1T(*recoJet, *matchedL1Jet);
        hForwardJet_AllRun_->FillL1T(*recoJet, *matchedL1Jet);

        if ( matchedHLTJet and deltaR(*recoJet, *matchedL1Jet) < maxL1DeltaR_ )
        {
          ++nForwardHLT;

          hAllJet_->FillHLT(*recoJet, *matchedHLTJet);
          hForwardJet_->FillHLT(*recoJet, *matchedHLTJet);

          hAllJet_AllRun_->FillHLT(*recoJet, *matchedHLTJet);
          hForwardJet_AllRun_->FillHLT(*recoJet, *matchedHLTJet);
        }
      }
    }

    if ( matchedHLTJet )
    {
      hAllJetNoL1_->FillHLT(*recoJet, *matchedHLTJet);
      if ( recoJetAbsEta < 2.5 ) hCentralJetNoL1_->FillHLT(*recoJet, *matchedHLTJet);
      else if ( recoJetAbsEta < 3 ) hOverlapJetNoL1_->FillHLT(*recoJet, *matchedHLTJet);
      else hForwardJetNoL1_->FillHLT(*recoJet, *matchedHLTJet);
    }
  }

  // We found leading jets. do the trigger object matching
  if ( leadingJet )
  {
    hAllLeadingJet_->FillReco(*leadingJet);
    hAllLeadingJet_AllRun_->FillReco(*leadingJet);

    // Retry matching
    const l1extra::L1JetParticle* matchedCentralL1Jet = getBestMatch(*leadingJet, centralL1Jets);
    const l1extra::L1JetParticle* matchedForwardL1Jet = getBestMatch(*leadingJet, forwardL1Jets);

    if ( matchedCentralL1Jet )
    {
      hAllLeadingJet_->FillL1T(*leadingJet, *matchedCentralL1Jet);
      hAllLeadingJet_AllRun_->FillL1T(*leadingJet, *matchedCentralL1Jet);

      const trigger::TriggerObject* matchedHLTJet = getBestMatch(*leadingJet, triggerObjects);
      const double l1DeltaR = deltaR(*leadingJet, *matchedCentralL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hAllLeadingJet_->FillHLT(*leadingJet, *matchedHLTJet);
        hAllLeadingJet_AllRun_->FillHLT(*leadingJet, *matchedHLTJet);
      }
    }
    if ( matchedForwardL1Jet )
    {
      hAllLeadingJet_->FillL1T(*leadingJet, *matchedForwardL1Jet);
      hAllLeadingJet_AllRun_->FillL1T(*leadingJet, *matchedForwardL1Jet);

      const trigger::TriggerObject* matchedHLTJet = getBestMatch(*leadingJet, triggerObjects);
      const double l1DeltaR = deltaR(*leadingJet, *matchedForwardL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hAllLeadingJet_->FillHLT(*leadingJet, *matchedHLTJet);
        hAllLeadingJet_AllRun_->FillHLT(*leadingJet, *matchedHLTJet);
      }
    }
  }

  if ( centralLeadingJet )
  {
    hCentralLeadingJet_->FillReco(*centralLeadingJet);
    hCentralLeadingJet_AllRun_->FillReco(*centralLeadingJet);

    // Retry matching
    const l1extra::L1JetParticle* matchedCentralL1Jet = getBestMatch(*centralLeadingJet, centralL1Jets);

    if ( matchedCentralL1Jet )
    {
      hCentralLeadingJet_->FillL1T(*centralLeadingJet, *matchedCentralL1Jet);
      hCentralLeadingJet_AllRun_->FillL1T(*centralLeadingJet, *matchedCentralL1Jet);

      const trigger::TriggerObject* matchedHLTJet = getBestMatch(*centralLeadingJet, triggerObjects);
      const double l1DeltaR = deltaR(*centralLeadingJet, *matchedCentralL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hCentralLeadingJet_->FillHLT(*centralLeadingJet, *matchedHLTJet);
        hCentralLeadingJet_AllRun_->FillHLT(*centralLeadingJet, *matchedHLTJet);
      }
    }
  }

  if ( overlapLeadingJet )
  {
    hOverlapLeadingJet_->FillReco(*overlapLeadingJet);
    hOverlapLeadingJet_AllRun_->FillReco(*overlapLeadingJet);

    // Retry matching
    const l1extra::L1JetParticle* matchedL1Jet = getBestMatch(*overlapLeadingJet, centralL1Jets);

    if ( matchedL1Jet )
    {
      hOverlapLeadingJet_->FillL1T(*overlapLeadingJet, *matchedL1Jet);
      hOverlapLeadingJet_AllRun_->FillL1T(*overlapLeadingJet, *matchedL1Jet);

      const trigger::TriggerObject* matchedHLTJet = getBestMatch(*overlapLeadingJet, triggerObjects);
      const double l1DeltaR = deltaR(*overlapLeadingJet, *matchedL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hOverlapLeadingJet_->FillHLT(*overlapLeadingJet, *matchedHLTJet);
        hOverlapLeadingJet_AllRun_->FillHLT(*overlapLeadingJet, *matchedHLTJet);
      }
    }
  }

  if ( forwardLeadingJet )
  {
    hForwardLeadingJet_->FillReco(*forwardLeadingJet);
    hForwardLeadingJet_AllRun_->FillReco(*forwardLeadingJet);

    // Retry matching
    const l1extra::L1JetParticle* matchedForwardL1Jet = getBestMatch(*forwardLeadingJet, forwardL1Jets);

    if ( matchedForwardL1Jet )
    {
      hForwardLeadingJet_->FillL1T(*forwardLeadingJet, *matchedForwardL1Jet);
      hForwardLeadingJet_AllRun_->FillL1T(*forwardLeadingJet, *matchedForwardL1Jet);

      const trigger::TriggerObject* matchedHLTJet = getBestMatch(*forwardLeadingJet, triggerObjects);
      const double l1DeltaR = deltaR(*forwardLeadingJet, *matchedForwardL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hForwardLeadingJet_->FillHLT(*forwardLeadingJet, *matchedHLTJet);
        hForwardLeadingJet_AllRun_->FillHLT(*forwardLeadingJet, *matchedHLTJet);
      }
    }
  }

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

