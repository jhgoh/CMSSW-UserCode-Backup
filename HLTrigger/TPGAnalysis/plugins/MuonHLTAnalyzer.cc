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

using namespace std;

const l1extra::L1MuonParticle* MuonHLTAnalyzer::getBestMatch(const double recoPosEta, const double recoPosPhi, l1extra::L1MuonParticleCollection& l1Particles)
{
  const l1extra::L1MuonParticle* matchedL1Cand = 0;
  double matchedDeltaR = 1e14;

  for ( L1Iter l1Cand = l1Particles.begin(); l1Cand != l1Particles.end(); ++l1Cand )
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

const trigger::TriggerObject* MuonHLTAnalyzer::getBestMatch(const reco::Candidate& recoCand, std::vector<trigger::TriggerObject>& triggerObjects)
{
  const trigger::TriggerObject* matchedHLTCand = 0;
  double matchedDeltaR = 1e14;

  for ( HLTIter hltCand = triggerObjects.begin(); hltCand != triggerObjects.end(); ++hltCand )
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

  const edm::ParameterSet muonCutSet = pset.getParameter<edm::ParameterSet>("cut");
  recoMinPt_ = muonCutSet.getParameter<double>("recoMinPt");
  l1MinPt_ = muonCutSet.getParameter<double>("l1MinPt");
  l1MinQuality_ = muonCutSet.getParameter<unsigned int>("l1MinQuality");
  maxL1DeltaR_ = muonCutSet.getParameter<double>("maxL1DeltaR");
  const double workingPointEt = muonCutSet.getParameter<double>("workingPointEt");
  const double maxHLTDeltaR = muonCutSet.getParameter<double>("maxHLTDeltaR");
  const double maxL1DeltaR = muonCutSet.getParameter<double>("maxL1DeltaR");

  l1Matcher_ = new L1MuonMatcherAlgo(pset.getParameter<edm::ParameterSet>("l1MatcherConfig"));

  const int objectType = Histograms::ObjectType::Muon;

  hAll_ = new HTrigger("All", "All", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);
  hBarrel_ = new HTrigger("Barrel", "Barrel", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);
  hOverlap_ = new HTrigger("Overlap", "Overlap", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);
  hEndcap_ = new HTrigger("Endcap", "Endcap", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);

  hLeadingAll_ = new HTrigger("LeadingAll", "Leading all", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);
  hLeadingBarrel_ = new HTrigger("LeadingBarrel", "Leading barrel", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);
  hLeadingOverlap_ = new HTrigger("LeadingOverlap", "Leading overlap", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);
  hLeadingEndcap_ = new HTrigger("LeadingEndcap", "Leading endcap", workingPointEt, maxL1DeltaR, maxHLTDeltaR, objectType);

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

  // Initialize histogram objects
  hAll_->init(event.id());
  hBarrel_->init(event.id()); 
  hOverlap_->init(event.id());
  hEndcap_->init(event.id());

  hLeadingAll_->init(event.id());
  hLeadingBarrel_->init(event.id());
  hLeadingOverlap_->init(event.id());
  hLeadingEndcap_->init(event.id());

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
  const reco::Muon* leadingBarrelMuon = 0;
  const reco::Muon* leadingOverlapMuon = 0;
  const reco::Muon* leadingEndcapMuon = 0;

  double leadingMuonPosEta = 1e14, leadingMuonPosPhi = 1e14;
  double leadingBarrelMuonPosEta = 1e14, leadingBarrelMuonPosPhi = 1e14;
  double leadingOverlapMuonPosEta = 1e14, leadingOverlapMuonPosPhi = 1e14;
  double leadingEndcapMuonPosEta = 1e14, leadingEndcapMuonPosPhi = 1e14;

  for ( edm::View<reco::Muon>::const_iterator recoMuon = recoMuonHandle->begin();
        recoMuon != recoMuonHandle->end(); ++recoMuon )
  {
    // Do some basic cuts
    if ( !isGoodMuon(*recoMuon) ) continue;

    const double recoMuonAbsEta = fabs(recoMuon->eta());
    const double recoMuonPt = recoMuon->pt();

    // Get the recoMuon's eta and phi at Muon station 2 to do the extrapolation
    TrajectoryStateOnSurface tsos = l1Matcher_->extrapolate(*recoMuon);
    if ( !tsos.isValid() ) continue;

    const double recoPosEta = tsos.globalPosition().eta();
    const double recoPosPhi = tsos.globalPosition().phi();

    const reco::Muon* recoMuonP = &*recoMuon;
    const l1extra::L1MuonParticle* matchedL1 = getBestMatch(recoPosEta, recoPosPhi, l1Muons);
    const trigger::TriggerObject* matchedHLT = getBestMatch(*recoMuon, triggerObjects);

    if ( !leadingMuon or leadingMuon->pt() < recoMuonPt ) 
    {
      leadingMuon = &(*recoMuon);
      leadingMuonPosEta = recoPosEta;
      leadingMuonPosPhi = recoPosPhi;
    }

    if ( recoMuonAbsEta < maxEtaBarrel_ )
    {
      if ( !leadingBarrelMuon or leadingBarrelMuon->pt() < recoMuonPt ) 
      {
        leadingBarrelMuon = &(*recoMuon);
        leadingBarrelMuonPosEta = recoPosEta;
        leadingBarrelMuonPosPhi = recoPosPhi;
      }
    }
    else if ( recoMuonAbsEta < maxEtaOverlap_ )
    {
      if ( !leadingOverlapMuon or leadingOverlapMuon->pt() < recoMuonPt ) 
      {
        leadingOverlapMuon = &(*recoMuon);
        leadingOverlapMuonPosEta = recoPosEta;
        leadingOverlapMuonPosPhi = recoPosPhi;
      }
    }
    else
    {
      if ( !leadingEndcapMuon or leadingEndcapMuon->pt() < recoMuonPt ) 
      {
        leadingEndcapMuon = &(*recoMuon);
        leadingEndcapMuonPosEta = recoPosEta;
        leadingEndcapMuonPosPhi = recoPosPhi;
      }
    }

    hAll_->fill(recoMuonP, matchedL1, matchedHLT);
    if ( recoMuonAbsEta < maxEtaBarrel_ ) hBarrel_->fill(recoMuonP, matchedL1, matchedHLT);
    else if ( recoMuonAbsEta < maxEtaOverlap_ ) hOverlap_->fill(recoMuonP, matchedL1, matchedHLT);
    else hEndcap_->fill(recoMuonP, matchedL1, matchedHLT);
  }

  // We found leading muons. do the trigger object matching
  if ( leadingMuon )
  {
    // Retry matching
    const l1extra::L1MuonParticle* matchedL1 = getBestMatch(leadingMuonPosEta, leadingMuonPosPhi, l1Muons);
    const trigger::TriggerObject* matchedHLT = getBestMatch(*leadingMuon, triggerObjects);
    hLeadingAll_->fill(leadingMuon, matchedL1, matchedHLT, leadingMuonPosEta, leadingMuonPosPhi);
  }

  if ( leadingBarrelMuon )
  {
    // Retry matching
    const l1extra::L1MuonParticle* matchedL1 = getBestMatch(leadingBarrelMuonPosEta, leadingBarrelMuonPosPhi, l1Muons);
    const trigger::TriggerObject* matchedHLT = getBestMatch(*leadingBarrelMuon, triggerObjects);
    hLeadingBarrel_->fill(leadingBarrelMuon, matchedL1, matchedHLT, leadingBarrelMuonPosEta, leadingBarrelMuonPosPhi);
  }

  if ( leadingOverlapMuon )
  {
    // Retry matching
    const l1extra::L1MuonParticle* matchedL1 = getBestMatch(leadingOverlapMuonPosEta, leadingOverlapMuonPosPhi, l1Muons);
    const trigger::TriggerObject* matchedHLT = getBestMatch(*leadingOverlapMuon, triggerObjects);
    hLeadingOverlap_->fill(leadingOverlapMuon, matchedL1, matchedHLT, leadingOverlapMuonPosEta, leadingOverlapMuonPosPhi);
  }

  if ( leadingEndcapMuon )
  {
    // Retry matching
    const l1extra::L1MuonParticle* matchedL1 = getBestMatch(leadingEndcapMuonPosEta, leadingEndcapMuonPosPhi, l1Muons);
    const trigger::TriggerObject* matchedHLT = getBestMatch(*leadingEndcapMuon, triggerObjects);
    hLeadingEndcap_->fill(leadingEndcapMuon, matchedL1, matchedHLT, leadingEndcapMuonPosEta, leadingEndcapMuonPosPhi);
  }

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
