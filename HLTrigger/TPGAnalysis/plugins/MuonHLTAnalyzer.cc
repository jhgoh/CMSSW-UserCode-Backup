#include "HLTrigger/TPGAnalysis/interface/MuonHLTAnalyzer.h"
#include "HLTrigger/TPGAnalysis/interface/Histograms.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
    
// Muon track extrapolation
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <TString.h>
#include <memory>

const l1extra::L1MuonParticle* MuonHLTAnalyzer::getBestMatch(const reco::Candidate& recoCand, L1Iter l1Begin, L1Iter l1End)
{
  const l1extra::L1MuonParticle* matchedL1Cand = 0;
  double matchedDeltaR = 1e14;

  //const double recoEta = recoCand.eta();
  //const double recoPhi = recoCand.phi();

  // Make track extrapolation to muon station
  TrajectoryStateOnSurface tsos = l1Matcher_->extrapolate(recoCand);
  if ( !tsos.isValid() ) return 0;

  const double recoPosEta = tsos.globalPosition().eta();
  const double recoPosPhi = tsos.globalPosition().phi();

  for ( L1Iter l1Cand = l1Begin; l1Cand != l1End; ++l1Cand )
  {
    const double l1Eta = l1Cand->eta();
    const double l1Phi = l1Cand->phi();

    const double l1DeltaR = deltaR(l1Eta, l1Phi, recoPosEta, recoPosPhi);

    if ( l1DeltaR < matchedDeltaR )
    {
      matchedL1Cand = &(*l1Cand);
      matchedDeltaR = l1DeltaR;
    }
  }

  return matchedL1Cand;
}

const trigger::TriggerObject* MuonHLTAnalyzer::getBestMatch(const reco::Candidate& recoCand, HLTIter hltBegin, HLTIter hltEnd)
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

MuonHLTAnalyzer::MuonHLTAnalyzer(const edm::ParameterSet& pset):
  maxEtaBarrel_(0.9), maxEtaOverlap_(1.2)
{
  interestedFilterName_ = pset.getParameter<std::string>("interestedFilterName");

  muonCutSet_ = pset.getParameter<edm::ParameterSet>("cut");

  l1Matcher_ = new L1MuonMatcherAlgo(pset.getParameter<edm::ParameterSet>("l1MatcherConfig"));
}

MuonHLTAnalyzer::~MuonHLTAnalyzer()
{
}

void MuonHLTAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  l1Matcher_->init(eventSetup);

  const int runNumber = run.run();

  if ( hMuon_ByRun_.find(runNumber) == hMuon_ByRun_.end() )
  {
    edm::Service<TFileService> fs;

    TFileDirectory runDir = fs->mkdir(Form("Run %d", runNumber));

    TFileDirectory muonDir = runDir.mkdir("All");
    TFileDirectory barrelMuonDir = runDir.mkdir("Barrel");
    TFileDirectory overlapMuonDir = runDir.mkdir("Overlap");
    TFileDirectory endcapMuonDir = runDir.mkdir("Endcap");

    hMuon_ByRun_[runNumber] = new Histograms(muonDir, "All", muonCutSet_, Histograms::ObjectType::Muon);
    hBarrelMuon_ByRun_[runNumber] = new Histograms(barrelMuonDir, "Barrel", muonCutSet_, Histograms::ObjectType::Muon);
    hOverlapMuon_ByRun_[runNumber] = new Histograms(overlapMuonDir, "Overlap", muonCutSet_, Histograms::ObjectType::Muon);
    hEndcapMuon_ByRun_[runNumber] = new Histograms(endcapMuonDir, "Endcap", muonCutSet_, Histograms::ObjectType::Muon);
  }

  hMuon_ = hMuon_ByRun_[runNumber];
  hBarrelMuon_ = hBarrelMuon_ByRun_[runNumber];
  hOverlapMuon_ = hOverlapMuon_ByRun_[runNumber];
  hEndcapMuon_ = hEndcapMuon_ByRun_[runNumber];
}

void MuonHLTAnalyzer::endRun()
{
}

void MuonHLTAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<l1extra::L1MuonParticleCollection> l1MuonHandle;
  if ( !event.getByLabel(edm::InputTag("l1extraParticles"), l1MuonHandle) )
  {
    edm::LogError("MuonHLTAnalyzer") << "Cannot find l1extra muons\n";
    return;
  }

  edm::Handle<trigger::TriggerEvent> triggerEventHandle;
  if ( !event.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", "HLT"), triggerEventHandle) )
  {
    edm::LogError("MuonHLTAnalyzer") << "Cannot find TriggerEvent\n";
    return;
  }
  const trigger::TriggerObjectCollection& triggerObjects = triggerEventHandle->getObjects();

  edm::Handle<edm::View<reco::Muon> > recoMuonHandle;
  if ( !event.getByLabel(edm::InputTag("muons"), recoMuonHandle) )
  {
    edm::LogError("MuonHLTAnalyzer") << "Cannot find reco muons\n";
    return;
  }

  // Loop over all reco muons
  int nReco = 0, nRecoBarrel = 0, nRecoOverlap = 0, nRecoEndcap = 0;

  // Collect HLT objects
  std::vector<trigger::TriggerObject> selectedTriggerObjects;
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

  for ( edm::View<reco::Muon>::const_iterator recoMuon = recoMuonHandle->begin();
        recoMuon != recoMuonHandle->end(); ++recoMuon )
  {
    // Do some basic cuts
    if ( !isGoodMuon(*recoMuon) ) continue;

    ++nReco;
    const double recoMuonAbsEta = fabs(recoMuon->eta());

    hMuon_->FillReco(*recoMuon);
    if ( recoMuonAbsEta < maxEtaBarrel_ ) hBarrelMuon_->FillReco(*recoMuon);
    else if ( recoMuonAbsEta < maxEtaOverlap_ ) hOverlapMuon_->FillReco(*recoMuon);
    else hEndcapMuon_->FillReco(*recoMuon);

    const l1extra::L1MuonParticle* matchedL1Muon = getBestMatch(*recoMuon, l1MuonHandle->begin(), l1MuonHandle->end());
    const trigger::TriggerObject* matchedHLTMuon = getBestMatch(*recoMuon, selectedTriggerObjects.begin(), selectedTriggerObjects.end());

    if ( matchedL1Muon )
    {
      hMuon_->FillL1T(*recoMuon, *matchedL1Muon);
      if ( recoMuonAbsEta < maxEtaBarrel_ ) hBarrelMuon_->FillL1T(*recoMuon, *matchedL1Muon);
      else if ( recoMuonAbsEta < maxEtaOverlap_ ) hOverlapMuon_->FillL1T(*recoMuon, *matchedL1Muon);
      else hEndcapMuon_->FillL1T(*recoMuon, *matchedL1Muon);

      if ( matchedHLTMuon )
      {
        hMuon_->FillHLT(*recoMuon, *matchedHLTMuon);
        if ( recoMuonAbsEta < maxEtaBarrel_ ) hBarrelMuon_->FillHLT(*recoMuon, *matchedHLTMuon);
        else if ( recoMuonAbsEta < maxEtaOverlap_ ) hOverlapMuon_->FillHLT(*recoMuon, *matchedHLTMuon);
        else hEndcapMuon_->FillHLT(*recoMuon, *matchedHLTMuon);
      }
    }
  }

  hMuon_->hNReco->Fill(nReco);
  hBarrelMuon_->hNReco->Fill(nRecoBarrel);
  hOverlapMuon_->hNReco->Fill(nRecoOverlap);
  hEndcapMuon_->hNReco->Fill(nRecoEndcap);
}

bool MuonHLTAnalyzer::isGoodMuon(const reco::Muon& recoMuon)
{
  if ( !recoMuon.isGlobalMuon() or !recoMuon.isTrackerMuon() ) return false;
  if ( fabs(recoMuon.eta()) < 2.5 ) return false;

  if ( !muon::isGoodMuon(recoMuon, muon::GlobalMuonPromptTight) or
       !muon::isGoodMuon(recoMuon, muon::TrackerMuonArbitrated) ) return false;

  return true;
}
