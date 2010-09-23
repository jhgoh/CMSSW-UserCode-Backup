#include "HLTrigger/TPGAnalysis/interface/MuonHLTAnalyzer.h"
#include "HLTrigger/TPGAnalysis/interface/Histograms.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
    
// Muon track extrapolation
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <TString.h>
#include <memory>

const l1extra::L1MuonParticle* MuonHLTAnalyzer::getBestMatch(const double recoPosEta, const double recoPosPhi, L1Iter l1Begin, L1Iter l1End)
{
  const l1extra::L1MuonParticle* matchedL1Cand = 0;
  double matchedDeltaR = 1e14;

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
    //const double hltDeltaR = deltaR(recoCand, *hltCand);
    const double hltDeltaR = deltaR(recoCand.eta(), recoCand.phi(), hltCand->eta(), hltCand->phi());

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
  recoMinPt_ = muonCutSet_.getParameter<double>("recoMinPt");
  l1MinPt_ = muonCutSet_.getParameter<double>("l1MinPt");
  l1MinQuality_ = muonCutSet_.getParameter<unsigned int>("l1MinQuality");
  maxL1DeltaR_ = muonCutSet_.getParameter<double>("maxL1DeltaR");

  l1Matcher_ = new L1MuonMatcherAlgo(pset.getParameter<edm::ParameterSet>("l1MatcherConfig"));
}

MuonHLTAnalyzer::~MuonHLTAnalyzer()
{
}

void MuonHLTAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  if ( !l1Matcher_ )
  {
    edm::LogError("MuonHLTAnalyzer") << "L1 matcher not initialized!!!\n";
    return;
  }

  l1Matcher_->init(eventSetup);

  const int runNumber = run.run();

  if ( hAllMuon_ByRun_.find(runNumber) == hAllMuon_ByRun_.end() )
  {
    edm::Service<TFileService> fs;

    TFileDirectory runDir = fs->mkdir(Form("Run %d", runNumber));

    TFileDirectory allMuonDir = runDir.mkdir("All");
    TFileDirectory barrelMuonDir = runDir.mkdir("Barrel");
    TFileDirectory overlapMuonDir = runDir.mkdir("Overlap");
    TFileDirectory endcapMuonDir = runDir.mkdir("Endcap");

    const int objectType = Histograms::ObjectType::Muon;

    hAllMuon_ByRun_[runNumber] = new Histograms(allMuonDir, "All", muonCutSet_, objectType);
    hBarrelMuon_ByRun_[runNumber] = new Histograms(barrelMuonDir, "Barrel", muonCutSet_, objectType);
    hOverlapMuon_ByRun_[runNumber] = new Histograms(overlapMuonDir, "Overlap", muonCutSet_, objectType);
    hEndcapMuon_ByRun_[runNumber] = new Histograms(endcapMuonDir, "Endcap", muonCutSet_, objectType);
  }

  hAllMuon_ = hAllMuon_ByRun_[runNumber];
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
  const trigger::TriggerObjectCollection& allTriggerObjects = triggerEventHandle->getObjects();

  edm::Handle<edm::View<reco::Muon> > recoMuonHandle;
  if ( !event.getByLabel(edm::InputTag("muons"), recoMuonHandle) )
  {
    edm::LogError("MuonHLTAnalyzer") << "Cannot find reco muons\n";
    return;
  }

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  if ( !event.getByLabel("offlineBeamSpot", beamSpotHandle) )
  {
    edm::LogError("MuonHLTAnalyzer") << "Cannot find offlineBeamSpot\n";
    return;
  }
  beamSpotPosition_ = beamSpotHandle->position();

  // Loop over all reco muons
  int nReco = 0, nRecoBarrel = 0, nRecoOverlap = 0, nRecoEndcap = 0;

  // Collect L1 objects
  l1extra::L1MuonParticleCollection l1Muons;
  l1Muons.reserve(l1MuonHandle->size());

  for ( l1extra::L1MuonParticleCollection::const_iterator l1Muon = l1MuonHandle->begin();
        l1Muon != l1MuonHandle->end(); ++l1Muon )
  {
    if ( l1Muon->pt() < l1MinPt_ or 
         l1Muon->gmtMuonCand().quality() < l1MinQuality_ ) continue;

    l1Muons.push_back(*l1Muon);
  }

  // Collect HLT objects
  std::vector<trigger::TriggerObject> triggerObjects;
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

  const reco::Muon* leadingMuon = 0;
  const reco::Muon* barrelLeadingMuon = 0;
  const reco::Muon* overlapLeadingMuon = 0;
  const reco::Muon* endcapLeadingMuon = 0;
  for ( edm::View<reco::Muon>::const_iterator recoMuon = recoMuonHandle->begin();
        recoMuon != recoMuonHandle->end(); ++recoMuon )
  {
    // Do some basic cuts
    if ( !isGoodMuon(*recoMuon) ) continue;

    ++nReco;
    const double recoMuonAbsEta = fabs(recoMuon->eta());
    const double recoMuonPt = recoMuon->pt();

    // Fill basic information of recoMuon
    // and find the leading muon
    hAllMuon_->FillReco(*recoMuon);
    if ( !leadingMuon or leadingMuon->pt() < recoMuonPt ) leadingMuon = &(*recoMuon);
    if ( recoMuonAbsEta < maxEtaBarrel_ )
    {
      hBarrelMuon_->FillReco(*recoMuon);
      if ( !barrelLeadingMuon or barrelLeadingMuon->pt() < recoMuonPt ) barrelLeadingMuon = &(*recoMuon);
    }
    else if ( recoMuonAbsEta < maxEtaOverlap_ )
    {
      hOverlapMuon_->FillReco(*recoMuon);
      if ( !overlapLeadingMuon or overlapLeadingMuon->pt() < recoMuonPt ) overlapLeadingMuon = &(*recoMuon);
    }
    else
    {
      hEndcapMuon_->FillReco(*recoMuon);
      if ( !endcapLeadingMuon or endcapLeadingMuon->pt() < recoMuonPt ) endcapLeadingMuon = &(*recoMuon);
    }

    // Get the recoMuon's eta and phi at Muon station 2 to do the extrapolation
    TrajectoryStateOnSurface tsos = l1Matcher_->extrapolate(*recoMuon);
    if ( !tsos.isValid() ) continue;

    const double recoPosEta = tsos.globalPosition().eta();
    const double recoPosPhi = tsos.globalPosition().phi();

    // Then try matching
    const l1extra::L1MuonParticle* matchedL1Muon = getBestMatch(recoPosEta, recoPosPhi, l1Muons.begin(), l1Muons.end());
    const trigger::TriggerObject* matchedHLTMuon = getBestMatch(*recoMuon, triggerObjects.begin(), triggerObjects.end());

    if ( matchedL1Muon )
    {
      hAllMuon_->FillL1T(*recoMuon, *matchedL1Muon, recoPosEta, recoPosPhi);
      if ( recoMuonAbsEta < maxEtaBarrel_ ) hBarrelMuon_->FillL1T(*recoMuon, *matchedL1Muon, recoPosEta, recoPosPhi);
      else if ( recoMuonAbsEta < maxEtaOverlap_ ) hOverlapMuon_->FillL1T(*recoMuon, *matchedL1Muon, recoPosEta, recoPosPhi);
      else hEndcapMuon_->FillL1T(*recoMuon, *matchedL1Muon, recoPosEta, recoPosPhi);

      const double l1DeltaR = deltaR(recoPosEta, recoPosPhi, matchedL1Muon->eta(), matchedL1Muon->phi());
      if ( matchedHLTMuon and l1DeltaR < maxL1DeltaR_ )
      {
        hAllMuon_->FillHLT(*recoMuon, *matchedHLTMuon);
        if ( recoMuonAbsEta < maxEtaBarrel_ ) hBarrelMuon_->FillHLT(*recoMuon, *matchedHLTMuon);
        else if ( recoMuonAbsEta < maxEtaOverlap_ ) hOverlapMuon_->FillHLT(*recoMuon, *matchedHLTMuon);
        else hEndcapMuon_->FillHLT(*recoMuon, *matchedHLTMuon);
      }
    }
  }

//FIXME//
/*
  if ( leadingMuon )
  {
  }
*/

  hAllMuon_->hNReco->Fill(nReco);
  hBarrelMuon_->hNReco->Fill(nRecoBarrel);
  hOverlapMuon_->hNReco->Fill(nRecoOverlap);
  hEndcapMuon_->hNReco->Fill(nRecoEndcap);
}

bool MuonHLTAnalyzer::isGoodMuon(const reco::Muon& recoMuon)
{
  if ( !recoMuon.isGlobalMuon() or !recoMuon.isTrackerMuon() ) return false;

  if ( recoMuon.pt() < recoMinPt_ ) return false;

  const double absEta = fabs(recoMuon.eta());
  if ( absEta > 2.5 ) return false;

//  if ( !muon::isGoodMuon(recoMuon, muon::GlobalMuonPromptTight) or
//       !muon::isGoodMuon(recoMuon, muon::TrackerMuonArbitrated) ) return false;

  // Apply the VBTF selection
  const int nMatches = recoMuon.numberOfMatches();

  const reco::TrackRef trkTrack = recoMuon.innerTrack();
  const reco::TrackRef glbTrack = recoMuon.globalTrack();

  const reco::HitPattern& trkHit = trkTrack->hitPattern();
  const reco::HitPattern& glbHit = glbTrack->hitPattern();

  const double glbX2 = glbTrack->normalizedChi2();

  const int nMuonHit = glbHit.numberOfValidMuonHits();
  const int nTrkHit = trkHit.numberOfValidTrackerHits();
  const int nPixelHit = trkHit.numberOfValidPixelHits();

  const double dxy = glbTrack->dxy(beamSpotPosition_);

  if ( !(nMatches > 1 and absEta < 2.1 and fabs(dxy) < 0.2 and 
         glbX2 < 10 and nPixelHit > 0 and nTrkHit > 10 and nMuonHit > 0 ) ) return false;

  return true;
}
