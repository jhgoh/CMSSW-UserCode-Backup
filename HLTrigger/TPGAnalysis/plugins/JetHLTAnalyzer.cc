#include "HLTrigger/TPGAnalysis/interface/JetHLTAnalyzer.h"
#include "HLTrigger/TPGAnalysis/interface/Histograms.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
    
#include "DataFormats/Math/interface/deltaR.h"

#include <TString.h>
#include <memory>

const l1extra::L1JetParticle* JetHLTAnalyzer::getBestMatch(const reco::Candidate& recoCand, L1Iter l1Begin, L1Iter l1End)
{
  const l1extra::L1JetParticle* matchedL1Cand = 0;
  double matchedDeltaR = 1e14;

  for ( L1Iter l1Cand = l1Begin; l1Cand != l1End; ++l1Cand )
  {
    const double l1DeltaR = deltaR(recoCand, *l1Cand);
    
    if ( l1DeltaR < matchedDeltaR )
    {
      matchedL1Cand = &(*l1Cand);
      matchedDeltaR = l1DeltaR;
    }
  }

  return matchedL1Cand;
}

const trigger::TriggerObject* JetHLTAnalyzer::getBestMatch(const reco::Candidate& recoCand, HLTIter hltBegin, HLTIter hltEnd)
{
  const trigger::TriggerObject* matchedHLTCand = 0;
  double matchedDeltaR = 1e14;

  for ( HLTIter hltCand = hltBegin; hltCand != hltEnd; ++hltCand )
  {
    const double hltDeltaR = deltaR(recoCand, *hltCand);

    if ( hltDeltaR < matchedDeltaR )
    {
      matchedHLTCand = &(*hltCand);
      matchedDeltaR = hltDeltaR;
    }
  }

  return matchedHLTCand;
}

JetHLTAnalyzer::JetHLTAnalyzer(const edm::ParameterSet& pset)
{
  interestedFilterName_ = pset.getParameter<std::string>("interestedFilterName");
  recoJetTag_ = pset.getParameter<edm::InputTag>("recoJet");

  jetCutSet_ = pset.getParameter<edm::ParameterSet>("cut");

  recoMinEt_ = jetCutSet_.getParameter<double>("recoMinEt");
  l1MinEt_ = jetCutSet_.getParameter<double>("l1MinEt");
  maxL1DeltaR_ = jetCutSet_.getParameter<double>("maxL1DeltaR");

  jetIDHelper_ = new reco::helper::JetIDHelper(pset.getParameter<edm::ParameterSet>("JetIDParams"));

//  edm::Service<TFileService> fs;
}

JetHLTAnalyzer::~JetHLTAnalyzer()
{
/*
  for ( std::map<int, Histograms*>::iterator h = h.begin();
        h != h.end(); ++h )
  {
    delete(h->second);
    h->second = 0;
  }
*/
}

void JetHLTAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  const int runNumber = run.run();

  if ( hCentralJet_ByRun_.find(runNumber) == hCentralJet_ByRun_.end() )
  {
    edm::Service<TFileService> fs;

    TFileDirectory runDir = fs->mkdir(Form("Run %d", runNumber));

    hNReco_ByRun_[runNumber] = runDir.make<TH1F>("hNReco", "Number of reco object", 6, -0.5, 5.5);
    hNCentralL1T_ByRun_[runNumber] = runDir.make<TH1F>("hNCentralL1T", "Number of matched Central L1T object", 6, -0.5, 5.5);
    hNForwardL1T_ByRun_[runNumber] = runDir.make<TH1F>("hNForwardL1T", "Number of matched Forward L1T object", 6, -0.5, 5.5);
    hNCentralHLT_ByRun_[runNumber] = runDir.make<TH1F>("hNCentralHLT", "Number of matched Central HLT object", 6, -0.5, 5.5);
    hNForwardHLT_ByRun_[runNumber] = runDir.make<TH1F>("hNForwardHLT", "Number of matched Forward HLT object", 6, -0.5, 5.5);
    hNDuplicatedL1T_ByRun_[runNumber] = runDir.make<TH1F>("hNDuplicatedL1T", "Number of duplicated L1T-reco matching in Central L1T - Forward L1T", 6, -0.5, 5.5);

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

    hAllLeadingJet_ByRun_[runNumber] = new Histograms(allLeadingJetDir, "AllLeading", jetCutSet_, objectType);
    hCentralLeadingJet_ByRun_[runNumber] = new Histograms(centralLeadingJetDir, "CentralLeading", jetCutSet_, objectType);
    hForwardLeadingJet_ByRun_[runNumber] = new Histograms(forwardLeadingJetDir, "ForwardLeading", jetCutSet_, objectType);   
  }

  hNReco_ = hNReco_ByRun_[runNumber];
  hNCentralL1T_ = hNCentralL1T_ByRun_[runNumber];
  hNForwardL1T_ = hNForwardL1T_ByRun_[runNumber];
  hNCentralHLT_ = hNCentralHLT_ByRun_[runNumber];
  hNForwardHLT_ = hNForwardHLT_ByRun_[runNumber];
  hNDuplicatedL1T_ = hNDuplicatedL1T_ByRun_[runNumber];

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

/*
  edm::Handle<edm::View<reco::CaloMET> > offMETHandle;
  if ( !event.getByLabel(offMETTag_, offMETHandle) )
  {
    edm::LogError("JetHLTAnalyzer") << "Cannot find reco METs\n";
    return;
  }
*/

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
    if ( !isGoodJet(*recoJet, event) ) continue;

    ++nReco;
    const double recoJetAbsEta = fabs(recoJet->eta());
    const double recoJetEt = recoJet->et();

    // Fill basic information of recoJet
    // and find the leading Jet
    hAllJet_->FillReco(*recoJet);
    if ( !leadingJet or leadingJet->et() < recoJetEt ) leadingJet = &(*recoJet);
    if ( recoJetAbsEta < 2.5 )
    {
      hCentralJet_->FillReco(*recoJet);
      if ( !centralLeadingJet or centralLeadingJet->et() < recoJetEt ) centralLeadingJet = &(*recoJet);
    }
    else
    {
      hForwardJet_->FillReco(*recoJet);
      if ( !forwardLeadingJet or forwardLeadingJet->et() < recoJetEt ) forwardLeadingJet = &(*recoJet);
    }

    // Try matching
    const l1extra::L1JetParticle* matchedCentralL1Jet = getBestMatch(*recoJet, centralL1Jets.begin(), centralL1Jets.end());
    const l1extra::L1JetParticle* matchedForwardL1Jet = getBestMatch(*recoJet, forwardL1Jets.begin(), forwardL1Jets.end());
    const trigger::TriggerObject* matchedHLTJet = getBestMatch(*recoJet, triggerObjects.begin(), triggerObjects.end());

    if ( matchedCentralL1Jet and matchedForwardL1Jet ) ++nDuplicatedL1T;
    if ( recoJetAbsEta < 2.5 and matchedCentralL1Jet ) 
    {
      ++nCentralL1T;
      hAllJet_->FillL1T(*recoJet, *matchedCentralL1Jet);
      hCentralJet_->FillL1T(*recoJet, *matchedCentralL1Jet);

      if ( matchedHLTJet and deltaR(*recoJet, *matchedCentralL1Jet) < maxL1DeltaR_ )
      {
        ++nCentralHLT;
        hCentralJet_->FillHLT(*recoJet, *matchedHLTJet);
      }
    }
    else if ( recoJetAbsEta >= 2.5 and matchedForwardL1Jet )
    {
      ++nForwardL1T;
      hAllJet_->FillL1T(*recoJet, *matchedForwardL1Jet);
      hForwardJet_->FillL1T(*recoJet, *matchedForwardL1Jet);

      if ( matchedHLTJet and deltaR(*recoJet, *matchedForwardL1Jet) < maxL1DeltaR_ )
      {
        ++nForwardHLT;
        hForwardJet_->FillHLT(*recoJet, *matchedHLTJet);
      }
    }

    if ( matchedHLTJet )
    {
      hAllJet_->FillHLT(*recoJet, *matchedHLTJet);
    }

  }

  // We found leading jets. do the trigger object matching
  if ( leadingJet )
  {
    hAllLeadingJet_->FillReco(*leadingJet);

    // Retry matching
    const l1extra::L1JetParticle* matchedCentralL1Jet = getBestMatch(*leadingJet, centralL1Jets.begin(), centralL1Jets.end());
    const l1extra::L1JetParticle* matchedForwardL1Jet = getBestMatch(*leadingJet, forwardL1Jets.begin(), forwardL1Jets.end());
    const trigger::TriggerObject* matchedHLTJet = getBestMatch(*leadingJet, triggerObjects.begin(), triggerObjects.end());

    if ( matchedCentralL1Jet )
    {
      hAllLeadingJet_->FillL1T(*leadingJet, *matchedCentralL1Jet);

      const double l1DeltaR = deltaR(*leadingJet, *matchedCentralL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hAllLeadingJet_->FillHLT(*leadingJet, *matchedHLTJet);
      }
    }
    if ( matchedForwardL1Jet )
    {
      hAllLeadingJet_->FillL1T(*leadingJet, *matchedForwardL1Jet);

      const double l1DeltaR = deltaR(*leadingJet, *matchedCentralL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hAllLeadingJet_->FillHLT(*leadingJet, *matchedHLTJet);
      }
    }
  }

  if ( centralLeadingJet )
  {
    hCentralLeadingJet_->FillReco(*centralLeadingJet);

    // Retry matching
    const l1extra::L1JetParticle* matchedCentralL1Jet = getBestMatch(*centralLeadingJet, centralL1Jets.begin(), centralL1Jets.end());
    const trigger::TriggerObject* matchedHLTJet = getBestMatch(*centralLeadingJet, triggerObjects.begin(), triggerObjects.end());

    if ( matchedCentralL1Jet )
    {
      hCentralLeadingJet_->FillL1T(*centralLeadingJet, *matchedCentralL1Jet);

      const double l1DeltaR = deltaR(*centralLeadingJet, *matchedCentralL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hCentralLeadingJet_->FillHLT(*centralLeadingJet, *matchedHLTJet);
      }
    }
  }

  if ( forwardLeadingJet )
  {
    hForwardLeadingJet_->FillReco(*forwardLeadingJet);

    // Retry matching
    const l1extra::L1JetParticle* matchedForwardL1Jet = getBestMatch(*forwardLeadingJet, forwardL1Jets.begin(), forwardL1Jets.end());
    const trigger::TriggerObject* matchedHLTJet = getBestMatch(*forwardLeadingJet, triggerObjects.begin(), triggerObjects.end());

    if ( matchedForwardL1Jet )
    {
      hForwardLeadingJet_->FillL1T(*forwardLeadingJet, *matchedForwardL1Jet);

      const double l1DeltaR = deltaR(*forwardLeadingJet, *matchedForwardL1Jet);
      if ( matchedHLTJet and l1DeltaR < maxL1DeltaR_ )
      {
        hForwardLeadingJet_->FillHLT(*forwardLeadingJet, *matchedHLTJet);
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
  if ( recoJet.emEnergyFraction() < 0.01 || recoJet.n90() < 2 ) return false;
  if ( recoJet.et() < recoMinEt_ ) return false;

  jetIDHelper_->calculate(event, recoJet);
  if ( jetIDHelper_->fHPD() > 0.98 ) return false;

  return true;
}

