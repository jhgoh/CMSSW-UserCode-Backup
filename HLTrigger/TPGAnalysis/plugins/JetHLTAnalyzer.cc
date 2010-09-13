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

  l1MinEt_ = jetCutSet_.getParameter<double>("l1MinEt");

  jetIDHelper_ = new reco::helper::JetIDHelper(pset.getParameter<edm::ParameterSet>("JetIDParams"));
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

    TFileDirectory jetDir = runDir.mkdir("All");
    TFileDirectory centralJetDir = runDir.mkdir("Central");
    TFileDirectory forwardJetDir = runDir.mkdir("Forward");

    hJet_ByRun_[runNumber] = new JetHistograms(jetDir, "All", jetCutSet_);
    hCentralJet_ByRun_[runNumber] = new JetHistograms(centralJetDir, "Central", jetCutSet_);
    hForwardJet_ByRun_[runNumber] = new JetHistograms(forwardJetDir, "Forward", jetCutSet_);
  }

  hNReco_ = hNReco_ByRun_[runNumber];
  hNCentralL1T_ = hNCentralL1T_ByRun_[runNumber];
  hNForwardL1T_ = hNForwardL1T_ByRun_[runNumber];
  hNCentralHLT_ = hNCentralHLT_ByRun_[runNumber];
  hNForwardHLT_ = hNForwardHLT_ByRun_[runNumber];
  hNDuplicatedL1T_ = hNDuplicatedL1T_ByRun_[runNumber];

  hJet_ = hJet_ByRun_[runNumber];
  hCentralJet_ = hCentralJet_ByRun_[runNumber];
  hForwardJet_ = hForwardJet_ByRun_[runNumber];
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
  const trigger::TriggerObjectCollection& triggerObjects = triggerEventHandle->getObjects();

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
  trigger::TriggerObjectCollection selectedTriggerObjects;
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
      selectedTriggerObjects.push_back(triggerObjects[*trgKey]);
    }
  }

  for ( edm::View<reco::CaloJet>::const_iterator recoJet = recoJetHandle->begin();
        recoJet != recoJetHandle->end(); ++recoJet )
  {
    // Do some basic cuts
    if ( !isGoodJet(*recoJet, event) ) continue;

    ++nReco;
    const double recoJetAbsEta = fabs(recoJet->eta());

    hJet_->FillReco(*recoJet);
    if ( recoJetAbsEta < 2.5 ) hCentralJet_->FillReco(*recoJet);
    if ( recoJetAbsEta >= 2.5 ) hForwardJet_->FillReco(*recoJet);

    // Try matching
    const l1extra::L1JetParticle* matchedCentralL1Jet = getBestMatch(*recoJet, centralL1JetHandle->begin(), centralL1JetHandle->end());
    const l1extra::L1JetParticle* matchedForwardL1Jet = getBestMatch(*recoJet, forwardL1JetHandle->begin(), forwardL1JetHandle->end());
    const trigger::TriggerObject* matchedHLTJet = getBestMatch(*recoJet, selectedTriggerObjects.begin(), selectedTriggerObjects.end());

    if ( matchedCentralL1Jet and matchedForwardL1Jet ) ++nDuplicatedL1T;
    if ( recoJetAbsEta < 2.5 and matchedCentralL1Jet ) 
    {
      ++nCentralL1T;
      hJet_->FillL1T(*recoJet, *matchedCentralL1Jet);
      hCentralJet_->FillL1T(*recoJet, *matchedCentralL1Jet);

      if ( matchedHLTJet )
      {
        ++nCentralHLT;
        hCentralJet_->FillHLT(*recoJet, *matchedHLTJet);
      }
    }
    if ( recoJetAbsEta >= 2.5 and matchedForwardL1Jet )
    {
      ++nForwardL1T;
      hJet_->FillL1T(*recoJet, *matchedForwardL1Jet);
      hForwardJet_->FillL1T(*recoJet, *matchedForwardL1Jet);

      if ( matchedHLTJet )
      {
        ++nForwardHLT;
        hForwardJet_->FillHLT(*recoJet, *matchedHLTJet);
      }
    }

    if ( matchedHLTJet )
    {
      hJet_->FillHLT(*recoJet, *matchedHLTJet);
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

  jetIDHelper_->calculate(event, recoJet);
  if ( jetIDHelper_->fHPD() > 0.98 ) return false;

  return true;
}


