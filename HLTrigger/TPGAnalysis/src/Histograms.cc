#include "HLTrigger/TPGAnalysis/interface/Histograms.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include <TString.h>
#include <TH1F.h>

Histograms::Histograms(TFileDirectory& dir, TString prefix, edm::ParameterSet& cutSet)
{
  minEt_ = cutSet.getParameter<double>("minEt");
  maxL1DeltaR_ = cutSet.getParameter<double>("maxL1DeltaR");
  maxHLTDeltaR_ = cutSet.getParameter<double>("maxHLTDeltaR");

  if ( prefix.Length() != 0 ) prefix += " ";

  SetHistogramBins();

  const unsigned int nBinEt = binsEt_.size()-1;
  const double* binsEt = &binsEt_[0];

  const unsigned int nBinEta = binsEta_.size()-1;
  const double* binsEta = &binsEta_[0];

  hNReco = dir.make<TH1F>("hNReco", prefix+"Number of reco object per event;Number of reco object", 4, 1, 5);
  hNL1T = dir.make<TH1F>("hNL1T", prefix+"Number of L1T matched reco objects per event;Number of reco object", 4, 1, 5);
  hNHLT = dir.make<TH1F>("hNHLT", prefix+"Number of HLT matched reco objects per event;Number of reco object", 4, 1, 5);

  // Kinematic variables of Reco objects
  hEtReco = dir.make<TH1F>("hEtReco", prefix+"Transverse momentum of reco object;Reco p_{T} [GeV/c]", nBinEt, binsEt);
  hEtaReco = dir.make<TH1F>("hEtaReco", prefix+"Pseudorapidity of reco object;Reco #eta", nBinEta, binsEta);
  hPhiReco = dir.make<TH1F>("hPhiReco", prefix+"Azimuthal angle of reco object;Reco #phi [Radian]", 50, -3.15, 3.15);

  // Kinematic variables of L1T objects
  hEtL1T = dir.make<TH1F>("hEtL1T", prefix+"Transverse momentum of reco object matched to L1T object;Reco p_{T} [GeV/c]", nBinEt, binsEt);
  hEtaL1T = dir.make<TH1F>("hEtaL1T", prefix+"Pseudorapidity of reco object matched to L1T object;Reco #eta", nBinEta, binsEta);
  hPhiL1T = dir.make<TH1F>("hPhiL1T", prefix+"Azimuthal angle of reco object matched to L1T object;Reco #phi [Radian]", 50, -3.15, 3.15);

  hL1EtL1T = dir.make<TH1F>("hL1EtL1T", prefix+"Transverse momentum of L1 object matched to reco object;L1 p_{T} [GeV/c]", nBinEt, binsEt);
  hL1EtaL1T = dir.make<TH1F>("hL1EtaL1T", prefix+"Pseudorapidity of L1 object matched to reco object;L1 #eta", nBinEta, binsEta);
  hL1PhiL1T = dir.make<TH1F>("hL1PhiL1T", prefix+"Azimuthal angle of L1 object matched to reco object;L1 #phi [Radian]", 50, -3.15, 3.15);

  // Kinematic variables of HLT objects
  hEtHLT = dir.make<TH1F>("hEtHLT", prefix+"Transverse momentum of reco object matched to HLT object;Reco p_{T} [GeV/c]", nBinEt, binsEt);
  hEtaHLT = dir.make<TH1F>("hEtaHLT", prefix+"Pseudorapidity of reco object matched to HLT object;Reco #eta", nBinEta, binsEta);
  hPhiHLT = dir.make<TH1F>("hPhiHLT", prefix+"Azimuthal angle of reco object matched to HLT object;Reco #phi [Radian]", 50, -3.15, 3.15);

  hHLTEtHLT = dir.make<TH1F>("hHLTEtHLT", prefix+"Transverse momentum of HLT object matched to reco object;HLT p_{T} [GeV/c]", nBinEt, binsEt);
  hHLTEtaHLT = dir.make<TH1F>("hHLTEtaHLT", prefix+"Pseudorapidity of HLT object matched to reco object;HLT #eta", nBinEta, binsEta);
  hHLTPhiHLT = dir.make<TH1F>("hHLTPhiHLT", prefix+"Azimuthal angle of HLT object matched to reco object;HLT #phi [Radian]", 50, -3.15, 3.15);

  // Matching variables
  hDeltaRL1T = dir.make<TH1F>("hDeltaRL1T", prefix+"#DeltaR between reco object - L1T object;#DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}", 100, 0, 1);
  hDeltaPhiL1T = dir.make<TH1F>("hDeltaPhiL1T", prefix+"#Delta#phi between reco object - L1T object;#Delta#phi [Radian]", 100, -1, 1);
  hDeltaEtaL1T = dir.make<TH1F>("hDeltaEtaL1T", prefix+"#Delta#eta between reco object - L1T object;#Delta#eta", 100, -1, 1);

  hDeltaRHLT = dir.make<TH1F>("hDeltaRHLT", prefix+"#DeltaR between reco object - HLT object;#DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}", 100, 0, 1);
  hDeltaPhiHLT = dir.make<TH1F>("hDeltaPhiHLT", prefix+"#Delta#phi between reco object - HLT object;#Delta#phi [Radian]", 100, -1, 1);
  hDeltaEtaHLT = dir.make<TH1F>("hDeltaEtaHLT", prefix+"#Delta#eta between reco object - HLT object;#Delta#eta", 100, -1, 1);
}

MuonHistograms::MuonHistograms(TFileDirectory& dir, TString prefix, edm::ParameterSet& cutSet):
  Histograms(dir, prefix, cutSet)
{
  if ( prefix.Length() != 0 ) prefix += " ";

  hQ = dir.make<TH1F>("hQ", prefix+"Electric charge;Electric charge", 3, -1.5, 1.5);

  hGlbNHit = dir.make<TH1F>("hGlbNHit", prefix+"Number of valid muon hits", 100, 0, 100);
  hTrkNHit = dir.make<TH1F>("hTrkNHit", prefix+"Number of valid tracker hits", 50, 0, 50);

  hGlbX2 = dir.make<TH1F>("hGlbX2", prefix+"Normalized #Chi^{2} of global track", 50, 0, 50);
  hTrkX2 = dir.make<TH1F>("hTrkX2", prefix+"Normalized #Chi^{2} of tracker track", 50, 0, 50);

  hRelIso = dir.make<TH1F>("hRelIso", prefix+"Relative isolation;Relative isolation", 100, 0, 10);
}

JetHistograms::JetHistograms(TFileDirectory& dir, TString prefix, edm::ParameterSet& cutSet):
  Histograms(dir, prefix, cutSet)
{
  if ( prefix.Length() != 0 ) prefix += " ";
}

void Histograms::FillReco(const reco::Candidate& recoCand)
{
  const double recoEt = recoCand.et();
  const double recoEta = recoCand.eta();
  const double recoPhi = recoCand.phi();

  hEtReco->Fill(recoEt);

  if ( recoEt < minEt_ ) return;

  hEtaReco->Fill(recoEta);
  hPhiReco->Fill(recoPhi);
}

void MuonHistograms::FillReco(const reco::Candidate& recoCand)
{
  Histograms::FillReco(recoCand);

  const reco::Muon* recoMuonP = dynamic_cast<const reco::Muon*>(&recoCand);
  if ( !recoMuonP ) return;

  const reco::Muon& recoMuon = *recoMuonP;

  if ( !recoMuon.isGlobalMuon() or !recoMuon.isTrackerMuon() ) return;
  const double recoPt = recoMuon.pt();
  if ( recoPt < minEt_ ) return;

  const reco::TrackRef trkTrack = recoMuon.innerTrack();
  //const reco::TrackRef staTrack = recoMuon->outerTrack();
  const reco::TrackRef glbTrack = recoMuon.globalTrack();

  const reco::HitPattern& trkHit = trkTrack->hitPattern();
  //const reco::HitPattern& staHit = staTrack->hitPattern();
  const reco::HitPattern& glbHit = glbTrack->hitPattern();

  const double glbX2 = glbTrack->normalizedChi2();
  const double trkX2 = trkTrack->normalizedChi2();
  const int nMuonHit = glbHit.numberOfValidMuonHits();
  const int nTrkHit = trkHit.numberOfValidTrackerHits();
  //const int nPixelHit = trkHit.numberOfValidPixelHits();
  //const int nMatches = recoMuon->numberOfMatches();

  //const int misHitInner = trkTrack->trackerExpectedHitsInner().numberOfHits();
  //const int misHitOuter = trkTrack->trackerExpectedHitsOuter().numberOfHits();

  const double trackIso = recoMuon.isolationR03().sumPt;
  const double caloIso = recoMuon.isolationR03().emEt + recoMuon.isolationR03().hadEt;
  const double relIso = (trackIso+caloIso)/recoPt;

  hGlbNHit->Fill(nMuonHit);
  hGlbX2->Fill(glbX2);

  hTrkNHit->Fill(nTrkHit);
  hTrkX2->Fill(trkX2);

  hRelIso->Fill(relIso);
}

void Histograms::FillL1T(const reco::Candidate& recoCand, const reco::LeafCandidate& l1Cand)
{
  const double recoEt = recoCand.et();
  const double recoEta = recoCand.eta();
  const double recoPhi = recoCand.phi();
  //const double recoCharge = recoCand.charge();

  const double l1Et = l1Cand.et();
  const double l1Eta = l1Cand.eta();
  const double l1Phi = l1Cand.phi();
  //const double l1Charge = l1Cand.charge();

  const double dR = deltaR(recoCand, l1Cand);
  const double dEta = l1Eta - recoEta;
  const double dPhi = deltaPhi(recoPhi, l1Phi);

  hEtL1T->Fill(recoEt);
  hL1EtL1T->Fill(l1Et);

  if ( recoEt < minEt_ ) return;
  
  hDeltaRL1T->Fill(dR);
  hDeltaEtaL1T->Fill(dEta);
  hDeltaPhiL1T->Fill(dPhi);

  if ( maxL1DeltaR_ < dR ) return;

  hEtaL1T->Fill(recoEta);
  hPhiL1T->Fill(recoPhi);

  hL1EtaL1T->Fill(l1Eta);
  hL1PhiL1T->Fill(l1Phi);
}

void MuonHistograms::FillL1T(const reco::Candidate& recoCand, const reco::LeafCandidate& l1Cand)
{
  edm::LogError("MuonHistograms") << "Wrong function. Please replace to FillL1T(const reco::Candidate& recoCand, const reco::LeafCandidate& l1Cand, const double recoPosEta, const double recoPosPhi)\n";
}

void MuonHistograms::FillL1T(const reco::Muon& recoMuon, const reco::LeafCandidate& l1Cand, const double recoPosEta, const double recoPosPhi)
{
  const double recoEt = recoMuon.et();
  const double recoEta = recoMuon.eta();
  const double recoPhi = recoMuon.phi();
  //const double recoCharge = recoMuon.charge();

  const double l1Et = l1Cand.et();
  const double l1Eta = l1Cand.eta();
  const double l1Phi = l1Cand.phi();
  //const double l1Charge = l1Cand.charge();

  const double dR = deltaR(recoPosEta, recoPosPhi, l1Eta, l1Phi);
  const double dEta = l1Eta - recoPosEta;
  const double dPhi = deltaPhi(recoPosPhi, l1Phi);

  hEtL1T->Fill(recoEt);
  hL1EtL1T->Fill(l1Et);

  if ( recoEt < minEt_ ) return;
  
  hDeltaRL1T->Fill(dR);
  hDeltaEtaL1T->Fill(dEta);
  hDeltaPhiL1T->Fill(dPhi);

  if ( maxL1DeltaR_ < dR ) return;

  hEtaL1T->Fill(recoEta);
  hPhiL1T->Fill(recoPhi);

  hL1EtaL1T->Fill(l1Eta);
  hL1PhiL1T->Fill(l1Phi);
}

void Histograms::FillHLT(const reco::Candidate& recoCand, const trigger::TriggerObject& triggerObject)
{
  const double recoEt = recoCand.et();
  const double recoEta = recoCand.eta();
  const double recoPhi = recoCand.phi();
  //const double recoCharge = recoCand.charge();

  const double hltEt = triggerObject.et();
  const double hltEta = triggerObject.eta();
  const double hltPhi = triggerObject.phi();
  //const double hltCharge = triggerObject.charge();

  const double dR = deltaR(recoCand, triggerObject);
  const double dEta = hltEta - recoEta;
  const double dPhi = deltaPhi(recoPhi, hltPhi);

  hEtHLT->Fill(recoEt);
  hHLTEtHLT->Fill(hltEt);

  if ( recoEt < minEt_ ) return;
  
  hDeltaRHLT->Fill(dR);
  hDeltaEtaHLT->Fill(dEta);
  hDeltaPhiHLT->Fill(dPhi);

  if ( maxHLTDeltaR_ < dR ) return;

  hEtaHLT->Fill(recoEta);
  hPhiHLT->Fill(recoPhi);

  hHLTEtaHLT->Fill(hltEta);
  hHLTPhiHLT->Fill(hltPhi);
}

void Histograms::SetHistogramBins()
{
  const unsigned int nEt = 50;
  const double minEt = 0, maxEt = 100;

  const unsigned int nEta = 50;
  const double minEta = -2.5, maxEta = 2.5;

  binsEt_.reserve(nEt);
  binsEta_.reserve(nEta);

  const double dEt = (maxEt-minEt)/nEt;
  for ( unsigned int i=0; i<=nEt; ++i )
  {
    binsEt_.push_back(minEt+dEt*i);
  }

  const double dEta = (maxEta-minEta)/nEta;
  for ( unsigned int i=0; i<=nEta; ++i )
  {
    binsEta_.push_back(minEta+dEta*i);
  }
}

void MuonHistograms::SetHistogramBins()
{
  const unsigned int nEt = 12;
  const double binsEt[nEt] = {
    0.,4.,8.,10.,15.,18.,21.,25.,30.,40.,70.,100.
  };

  const unsigned int nEta = 9;
  const double binsEta[nEta] = {
    -2.40, -1.95, -1.20, -0.90, 0.00, 0.90, 1.20, 1.95, 2.40
  };

  binsEt_.reserve(nEt);
  binsEta_.reserve(nEta);

  for ( unsigned int i=0; i<nEt; ++i )
  {
    binsEt_.push_back(binsEt[i]);
  }

  for ( unsigned int i=0; i<nEta; ++i )
  {
    binsEta_.push_back(binsEta[i]);
  }
}
