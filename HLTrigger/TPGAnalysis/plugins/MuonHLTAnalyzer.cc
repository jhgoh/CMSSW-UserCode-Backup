#include "HLTrigger/TPGAnalysis/plugins/MuonHLTAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
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
#include "DataFormats/Math/interface/deltaPhi.h"

#include <TString.h>

const static int nRecoMuonCutStep = 6;
const static char* recoMuonCutStepNames[nRecoMuonCutStep] = {
  "All", "Prompt", "MuonSTN", "MisHit", "EWK", "EWK+Iso"
};

MuonHLTAnalyzer::MuonHLTAnalyzer(const edm::ParameterSet& pset):
  l1Matcher_(pset.getParameter<edm::ParameterSet>("l1MatcherConfig"))
{
  muonTrigNames_ = pset.getParameter<VString>("trigNames");
  const int nTrigNames = muonTrigNames_.size();

  l1MuonTag_ = pset.getParameter<edm::InputTag>("l1Muon");
  recoMuonTag_ = pset.getParameter<edm::InputTag>("recoMuon");
  minPt_ = pset.getParameter<double>("minPt");
  maxRelIso_ = pset.getParameter<double>("maxRelIso");

  // Book histograms
  edm::Service<TFileService> fs;

  //// Basic constants
  const int l1PtNBin = 11;
  const double l1PtBins[l1PtNBin+1] = {
    0.,4.,8.,10.,15.,18.,21.,25.,30.,40.,70.,100.
  };

  const int l1EtaNBin = 62;
  const double l1EtaBins[l1EtaNBin+1] = {
    -2.40, -2.35, -2.30, -2.25, -2.20, -2.15, -2.10, -2.05,-2.00, -1.95, -1.90, -1.85, -1.80, -1.75, -1.70, -1.60,-1.50, -1.40, -1.30, -1.20, -1.10, -1.00, -0.90, -0.80,-0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, -0.00,0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60,1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00, 2.05,2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40
  };

  hNEvent_ = fs->make<TH1F>("hNEvent", "Number of events passing trigger paths;Trigger path", nTrigNames+1, 0, nTrigNames+1);
  hNEvent_->GetXaxis()->SetBinLabel(1, "All");
  for ( int i=0; i<nTrigNames; ++i )
  {
    hNEvent_->GetXaxis()->SetBinLabel(i+2, muonTrigNames_[i].c_str());
  }

  hNRecoMuon_ = fs->make<TH1F>("hNRecoMuon", "Number of recoMuons per event", nRecoMuonCutStep, 0, nRecoMuonCutStep);

  for ( int i=0; i<nRecoMuonCutStep; ++i )
  {
    hNRecoMuon_->GetXaxis()->SetBinLabel(i+1, recoMuonCutStepNames[i]);
  }

  // Book histograms for each cut steps
  for ( int i=0; i<nRecoMuonCutStep; ++i )
  {
    TFileDirectory dir = fs->mkdir(Form("CutStep%d_%s", i, recoMuonCutStepNames[i]));

    // Histograms for track variables
    hPt_.push_back(dir.make<TH1F>("hPt", "Global muon Transverse momentum", 100, 0, 100));
    hEta_.push_back(dir.make<TH1F>("hEta", "Global muon Pseudorapidity", 100, -2.5, 2.5));
    hPhi_.push_back(dir.make<TH1F>("hPhi", "Global muon Azimuthal angle", 100, -3.15, 3.15));
    hQ_.push_back(dir.make<TH1F>("hQ", "Global muon charge", 3, -1.5, 1.5));

    hNGlbHit_.push_back(dir.make<TH1F>("hNMuonHit", "Global muon Number of valid muon hits", 100, 0, 100));
    hNTrkHit_.push_back(dir.make<TH1F>("hNTrkHit", "Global muon Number of valid tracker hits", 100, 0, 100));
    hGlbX2_.push_back(dir.make<TH1F>("hMuonX2", "Global muon Normalized #Chi^{2} of global track", 50, 0, 50));
    hTrkX2_.push_back(dir.make<TH1F>("hTrkX2", "Global muon Normalized #Chi^{2} of tracker track", 50, 0, 50));

    hRelIso_.push_back(dir.make<TH1F>("hRelIso", "Global muon relative isolation", 100, 0, 20));

    // Histograms for Matching
    hDeltaR_.push_back(dir.make<TH1F>("hDeltaR", "Global muon - L1 muon matching #DeltaR;#DeltaR = #sqrt{#Delta#phi^{2} + #Delta#eta^{2}}", 50, 0, 2));
    hDeltaPhi_.push_back(dir.make<TH1F>("hDeltaPhi", "Global muon - L1 muon matching #DeltaR;#Delta#phi [Radian]", 50, 0, 2));
    hDeltaEta_.push_back(dir.make<TH1F>("hDeltaEta", "Global muon - L1 muon matching #Delta#eta;#Delta#eta", 50, 0, 2));

    // Histograms for L1 variables
    hL1Pt_.push_back(dir.make<TH1F>("hL1Pt", "L1 muon p_{T}", l1PtNBin, l1PtBins));
    hL1Eta_.push_back(dir.make<TH1F>("hL1Eta", "L1 muon #eta", l1EtaNBin, l1EtaBins));
    hL1Phi_.push_back(dir.make<TH1F>("hL1Phi", "L1 muon #phi", 50, -3.15, 3.15));

    hMatchedL1Pt_.push_back(dir.make<TH1F>("hMatchedL1Pt", "Matched L1 muon p_{T}", l1PtNBin, l1PtBins));
    hMatchedL1Eta_.push_back(dir.make<TH1F>("hMatchedL1Eta", "Matched L1 muon #eta", l1EtaNBin, l1EtaBins));
    hMatchedL1Phi_.push_back(dir.make<TH1F>("hMatchedL1Phi", "Matched L1 muon #phi", 50, -3.15, 3.15));

    hL1MatchedGlbPt_.push_back(dir.make<TH1F>("hL1MatchedGlbPt", "L1 matched global muon p_{T};global muon p_{T} [GeV/c]", l1PtNBin, l1PtBins));
    hL1MatchedGlbEta_.push_back(dir.make<TH1F>("hL1MatchedGlbEta", "L1 matched global muon #eta;global muon #eta", l1EtaNBin, l1EtaBins));
    hL1MatchedGlbPhi_.push_back(dir.make<TH1F>("hL1MatchedGlbPhi", "L1 matched global muon #phi;global muon #phi", 50, -3.15, 3.15));

    // Histograms in Barrel region
    hPtBarrel_.push_back(dir.make<TH1F>("hPtBarrel", "Global muon Transverse momentum", 100, 0, 100));
    hEtaBarrel_.push_back(dir.make<TH1F>("hEtaBarrel", "Global muon Pseudorapidity", 100, -2.5, 2.5));
    hPhiBarrel_.push_back(dir.make<TH1F>("hPhiBarrel", "Global muon Azimuthal angle", 100, -3.15, 3.15));

    hNGlbHitBarrel_.push_back(dir.make<TH1F>("hNGlbHitBarrel", "Global muon Number of valid muon hits", 100, 0, 100));
    hNTrkHitBarrel_.push_back(dir.make<TH1F>("hNTrkHitBarrel", "Global muon Number of valid tracker hits", 100, 0, 100));
    hGlbX2Barrel_.push_back(dir.make<TH1F>("hGlbX2Barrel", "Global muon Normalized #Chi^{2} of global track", 50, 0, 50));
    hTrkX2Barrel_.push_back(dir.make<TH1F>("hTrkX2Barrel", "Global muon Normalized #Chi^{2} of tracker track", 50, 0, 50));
    hRelIsoBarrel_.push_back(dir.make<TH1F>("hRelIsoBarrel", "Global muon relative isolation", 100, 0, 20));

    // Histograms in Overlap region
    hPtOverlap_.push_back(dir.make<TH1F>("hPtOverlap", "Global muon Transverse momentum", 100, 0, 100));
    hEtaOverlap_.push_back(dir.make<TH1F>("hEtaOverlap", "Global muon Pseudorapidity", 100, -2.5, 2.5));
    hPhiOverlap_.push_back(dir.make<TH1F>("hPhiOverlap", "Global muon Azimuthal angle", 100, -3.15, 3.15));

    hNGlbHitOverlap_.push_back(dir.make<TH1F>("hNGlbHitOverlap", "Global muon Number of valid muon hits", 100, 0, 100));
    hNTrkHitOverlap_.push_back(dir.make<TH1F>("hNTrkHitOverlap", "Global muon Number of valid tracker hits", 100, 0, 100));
    hGlbX2Overlap_.push_back(dir.make<TH1F>("hGlbX2Overlap", "Global muon Normalized #Chi^{2} of global track", 50, 0, 50));
    hTrkX2Overlap_.push_back(dir.make<TH1F>("hTrkX2Overlap", "Global muon Normalized #Chi^{2} of tracker track", 50, 0, 50));
    hRelIsoOverlap_.push_back(dir.make<TH1F>("hRelIsoOverlap", "Global muon relative isolation", 100, 0, 20));

    // Histograms in Endcap region
    hPtEndcap_.push_back(dir.make<TH1F>("hPt_Endcap", "Global muon Transverse momentum", 100, 0, 100));
    hEtaEndcap_.push_back(dir.make<TH1F>("hEta_Endcap", "Global muon Pseudorapidity", 100, -2.5, 2.5));
    hPhiEndcap_.push_back(dir.make<TH1F>("hPhi_Endcap", "Global muon Azimuthal angle", 100, -3.15, 3.15));

    hNGlbHitEndcap_.push_back(dir.make<TH1F>("hNGlbHitEndcap", "Global muon Number of valid muon hits", 100, 0, 100));
    hNTrkHitEndcap_.push_back(dir.make<TH1F>("hNTrkHitEndcap", "Global muon Number of valid tracker hits", 100, 0, 100));
    hGlbX2Endcap_.push_back(dir.make<TH1F>("hGlbX2Endcap", "Global muon Normalized #Chi^{2} of global track", 50, 0, 50));
    hTrkX2Endcap_.push_back(dir.make<TH1F>("hTrkX2Endcap", "Global muon Normalized #Chi^{2} of tracker track", 50, 0, 50));
    hRelIsoEndcap_.push_back(dir.make<TH1F>("hRelIsoEndcap", "Global muon relative isolation", 100, 0, 20));
  }
}

MuonHLTAnalyzer::~MuonHLTAnalyzer()
{
}

void MuonHLTAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
}

void MuonHLTAnalyzer::endRun()
{
}

void MuonHLTAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::TriggerResults> trigResultHandle;
  if ( !event.getByLabel("TriggerResults", trigResultHandle) ) return;
  const edm::TriggerResults* trigResult = trigResultHandle.product();

  edm::Handle<edm::View<l1extra::L1MuonParticle> > l1MuonHandle;
  if ( !event.getByLabel(l1MuonTag_, l1MuonHandle) ) return;

  edm::Handle<edm::View<reco::Muon> > recoMuonHandle;
  if ( !event.getByLabel(recoMuonTag_, recoMuonHandle) ) return;

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  if ( !event.getByLabel("offlineBeamSpot", beamSpotHandle) ) return;

  // Count how much event passed each trigger paths
  if ( trigResult->wasrun() and trigResult->accept() )
  {
    hNEvent_->Fill(0);

    const int nTrigResult = trigResult->size();
    const edm::TriggerNames& trigNames = event.triggerNames(*trigResult);

    for ( int i=0; i<nTrigResult; ++i )
    {
      if ( !trigResult->accept(i) ) continue;

      const std::string triggerName = trigNames.triggerName(i);

      // Find Muon-HLT in the triggerResult
      int bin = -1;
      const int nMuonTrigNames = muonTrigNames_.size();
      for ( int j=0; j<nMuonTrigNames; ++j )
      {
        if ( muonTrigNames_[j] == triggerName )
        {
          bin = j;
          break;
        }
      }

      if ( bin != -1 ) hNEvent_->Fill(bin+1);
    }
  }

  // Start GlobalMuon-L1Muon matching
  for ( edm::View<reco::Muon>::const_iterator recoMuon = recoMuonHandle->begin();
        recoMuon != recoMuonHandle->end(); ++recoMuon )
  {
    if ( !recoMuon->isGlobalMuon() ) continue;

    const reco::TrackRef trkTrack = recoMuon->innerTrack();
    const reco::TrackRef staTrack = recoMuon->outerTrack();
    const reco::TrackRef glbTrack = recoMuon->globalTrack();

    const reco::HitPattern& trkHit = trkTrack->hitPattern();
    const reco::HitPattern& staHit = staTrack->hitPattern();
    const reco::HitPattern& glbHit = glbTrack->hitPattern();

    const int misHitInner = trkTrack->trackerExpectedHitsInner().numberOfHits();
    const int misHitOuter = trkTrack->trackerExpectedHitsOuter().numberOfHits();

    const double recoPt = recoMuon->pt();
    const double recoEta = recoMuon->eta();
    const double recoPhi = recoMuon->phi();
    const double recoQ = recoMuon->charge();

    const double glbX2 = glbTrack->normalizedChi2();
    const double trkX2 = trkTrack->normalizedChi2();
    const int nMuonHit = glbHit.numberOfValidMuonHits();
    const int nTrkHit = trkHit.numberOfValidTrackerHits();    
    const int nPixelHit = trkHit.numberOfValidPixelHits();
    const int nMatches = recoMuon->numberOfMatches();

    const double trackIso = recoMuon->isolationR03().sumPt;
    const double caloIso = recoMuon->isolationR03().emEt + recoMuon->isolationR03().hadEt;
    const double relIso = (trackIso+caloIso)/recoPt;

    const double dxy = glbTrack->dxy(beamSpotHandle->position());

    // Check cut steps
    std::vector<bool> muonQuality(nRecoMuonCutStep);
    if ( recoPt > minPt_ )
    {
      muonQuality[0] = true;
      muonQuality[1] = muon::isGoodMuon(*recoMuon, muon::GlobalMuonPromptTight);
      muonQuality[2] = nMatches > 1;
      muonQuality[3] = ( misHitInner < 2 and misHitOuter < 2 );
      muonQuality[4] = ( muonQuality[2] and fabs(recoEta) < 2.1 and fabs(dxy) < 0.2 and glbX2 < 10 and 
                         nPixelHit > 0 and nTrkHit > 10 and nMuonHit > 0 );
      muonQuality[5] = ( muonQuality[4] and relIso < maxRelIso_ );
    }

    for ( int i=0; i<nRecoMuonCutStep; ++i )
    {
      if ( !muonQuality[i] ) continue;

      hPt_[i]->Fill(recoPt);
      hEta_[i]->Fill(recoEta);
      hPhi_[i]->Fill(recoPhi);
      hQ_[i]->Fill(recoQ);

      hNGlbHit_[i]->Fill(nMuonHit);
      hNTrkHit_[i]->Fill(nTrkHit);
      hGlbX2_[i]->Fill(glbX2);
      hTrkX2_[i]->Fill(trkX2);

      if ( fabs(recoEta) < 0.9 )
      {
        hPtBarrel_[i]->Fill(recoPt);
        hEtaBarrel_[i]->Fill(recoEta);
        hPhiBarrel_[i]->Fill(recoPhi);
        hQBarrel_[i]->Fill(recoQ);
      }
      else if ( fabs(recoEta) < 1.2 )
      {
        hPtOverlap_[i]->Fill(recoPt);
        hEtaOverlap_[i]->Fill(recoEta);
        hPhiOverlap_[i]->Fill(recoPhi);
        hQOverlap_[i]->Fill(recoQ);
      }
      else
      {
        hPtEndcap_[i]->Fill(recoPt);
        hEtaEndcap_[i]->Fill(recoEta);
        hPhiEndcap_[i]->Fill(recoPhi);
        hQEndcap_[i]->Fill(recoQ);
      }
    }

    // Now start Matching
    TrajectoryStateOnSurface tsos = l1Matcher_.extrapolate(*recoMuon);
    if ( !tsos.isValid() ) continue;

    const double recoPosEta = tsos.globalPosition().eta();
    const double recoPosPhi = tsos.globalPosition().phi();

    double matchedDeltaR = -999;
    l1extra::L1MuonParticle bestMatchingL1Muon;

    for ( edm::View<l1extra::L1MuonParticle>::const_iterator l1Muon = l1MuonHandle->begin();
          l1Muon != l1MuonHandle->end(); ++l1Muon )
    {
      const L1MuGMTExtendedCand& gmtCand = l1Muon->gmtMuonCand();
      const unsigned int l1Quality = gmtCand.quality();

      if ( l1Quality < 4 or l1Muon->pt() < 7 ) continue;

      const double l1PosEta = l1Muon->eta();
      const double l1PosPhi = l1Muon->phi();
      const double dR = deltaR(recoPosEta, recoPosPhi, l1PosEta, l1PosPhi);
      if ( matchedDeltaR < 0 or dR < matchedDeltaR ) 
      {
        matchedDeltaR = dR;
        bestMatchingL1Muon = *l1Muon;
      }
    }

    if ( matchedDeltaR < 0 or matchedDeltaR > 0.3 ) continue;

    // Now we have best matching l1Muon
    for ( int i=0; i<nRecoMuonCutStep; ++i )
    {
      if ( !muonQuality[i] ) continue;

      hL1MatchedGlbPt_[i]->Fill(recoPt);
      hL1MatchedGlbEta_[i]->Fill(recoEta);
      hL1MatchedGlbPhi_[i]->Fill(recoPhi);
    }
  }
}

