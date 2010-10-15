#include "HLTrigger/TPGAnalysis/interface/Histograms.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>

Histograms::Histograms(TFileDirectory& dir, TString prefix, edm::ParameterSet& cutSet, int objectType):
  objectType_(objectType)
{
  workingPointEt_ = cutSet.getParameter<double>("workingPointEt");
  maxL1DeltaR_ = cutSet.getParameter<double>("maxL1DeltaR");
  maxHLTDeltaR_ = cutSet.getParameter<double>("maxHLTDeltaR");

  if ( prefix.Length() != 0 ) prefix += " ";

  std::vector<double> binsEt;
  std::vector<double> binsEta;

  // Set bins
  if ( objectType_ == ObjectType::Muon )
  {
    const unsigned int nEt = 12;
    const double binsEtInput[nEt] = {
      0.,4.,8.,10.,15.,18.,21.,25.,30.,40.,70.,100.
    };

    const unsigned int nEta = 9;
    const double binsEtaInput[nEta] = {
      -2.40, -1.95, -1.20, -0.90, 0.00, 0.90, 1.20, 1.95, 2.40
    };

    std::copy(binsEtInput, binsEtInput+nEt, std::back_inserter(binsEt));
    std::copy(binsEtaInput, binsEtaInput+nEta, std::back_inserter(binsEta));
  }
  else if ( objectType_ == ObjectType::Jet )
  {
    const unsigned int nEt = 200;
    const double minEt = 0, maxEt = 200;

    const unsigned int nEta = 50;
    const double minEta = -5, maxEta = 5;

    binsEt.reserve(nEt);
    binsEta.reserve(nEta);

    const double dEt = (maxEt-minEt)/nEt;
    for ( unsigned int i=0; i<=nEt; ++i )
    {
      binsEt.push_back(minEt+dEt*i);
    }

    const double dEta = (maxEta-minEta)/nEta;
    for ( unsigned int i=0; i<=nEta; ++i )
    {
      binsEta.push_back(minEta+dEta*i);
    }
  }
  else
  {
    const unsigned int nEt = 100;
    const double minEt = 0, maxEt = 100;

    const unsigned int nEta = 50;
    const double minEta = -2.5, maxEta = 2.5;

    binsEt.reserve(nEt);
    binsEta.reserve(nEta);

    const double dEt = (maxEt-minEt)/nEt;
    for ( unsigned int i=0; i<=nEt; ++i )
    {
      binsEt.push_back(minEt+dEt*i);
    }

    const double dEta = (maxEta-minEta)/nEta;
    for ( unsigned int i=0; i<=nEta; ++i )
    {
      binsEta.push_back(minEta+dEta*i);
    }
  }

  const unsigned int nBinEt = binsEt.size()-1;
  const double* binsEtPtr = &binsEt[0];

  const unsigned int nBinEta = binsEta.size()-1;
  const double* binsEtaPtr = &binsEta[0];

  hNReco = dir.make<TH1F>("hNReco", prefix+"Number of reco object per event;Number of reco object", 4, 1, 5);
  hNL1T = dir.make<TH1F>("hNL1T", prefix+"Number of L1T matched reco objects per event;Number of reco object", 4, 1, 5);
  hNHLT = dir.make<TH1F>("hNHLT", prefix+"Number of HLT matched reco objects per event;Number of reco object", 4, 1, 5);

  // Kinematic variables of Reco objects
  hEtReco = dir.make<TH1F>("hEtReco", prefix+"Transverse energy of reco object;Reco p_{T} [GeV/c]", nBinEt, binsEtPtr);
  hEtaReco = dir.make<TH1F>("hEtaReco", prefix+"Pseudorapidity of reco object;Reco #eta", nBinEta, binsEtaPtr);
  hPhiReco = dir.make<TH1F>("hPhiReco", prefix+"Azimuthal angle of reco object;Reco #phi [Radian]", 50, -3.15, 3.15);

  // Kinematic variables of L1T objects
  hEtL1T = dir.make<TH1F>("hEtL1T", prefix+"Transverse energy of reco object matched to L1T object;Reco p_{T} [GeV/c]", nBinEt, binsEtPtr);
  hEtaL1T = dir.make<TH1F>("hEtaL1T", prefix+"Pseudorapidity of reco object matched to L1T object;Reco #eta", nBinEta, binsEtaPtr);
  hPhiL1T = dir.make<TH1F>("hPhiL1T", prefix+"Azimuthal angle of reco object matched to L1T object;Reco #phi [Radian]", 50, -3.15, 3.15);

  hL1EtL1T = dir.make<TH1F>("hL1EtL1T", prefix+"Transverse energy of L1 object matched to reco object;L1 p_{T} [GeV/c]", nBinEt, binsEtPtr);
  hL1EtaL1T = dir.make<TH1F>("hL1EtaL1T", prefix+"Pseudorapidity of L1 object matched to reco object;L1 #eta", nBinEta, binsEtaPtr);
  hL1PhiL1T = dir.make<TH1F>("hL1PhiL1T", prefix+"Azimuthal angle of L1 object matched to reco object;L1 #phi [Radian]", 50, -3.15, 3.15);

  hEtVsL1Et = dir.make<TH2F>("hEtVsL1Et", prefix+"Transverse energy of reco object vs L1 object;Reco p_{T} [GeV/c];L1 p_{T} [GeV/c]", nBinEt, binsEtPtr, nBinEt, binsEtPtr);
  hEtaVsL1Eta = dir.make<TH2F>("hEtaVsL1Eta", prefix+"Pseudorapidity of reco object vs L1 object;Reco #eta;L1 #eta", nBinEta, binsEtaPtr, nBinEta, binsEtaPtr);
  hPhiVsL1Phi = dir.make<TH2F>("hPhiVsL1Phi", prefix+"Azimuthal angle of reco object vs L1 object;Reco #phi [Radian];L1 #phi [Radian]", 50, -3.15, 3.15, 50, -3.15, 3.15);

  // Kinematic variables of HLT objects
  hEtHLT = dir.make<TH1F>("hEtHLT", prefix+"Transverse energy of reco object matched to HLT object;Reco p_{T} [GeV/c]", nBinEt, binsEtPtr);
  hEtaHLT = dir.make<TH1F>("hEtaHLT", prefix+"Pseudorapidity of reco object matched to HLT object;Reco #eta", nBinEta, binsEtaPtr);
  hPhiHLT = dir.make<TH1F>("hPhiHLT", prefix+"Azimuthal angle of reco object matched to HLT object;Reco #phi [Radian]", 50, -3.15, 3.15);

  hHLTEtHLT = dir.make<TH1F>("hHLTEtHLT", prefix+"Transverse energy of HLT object matched to reco object;HLT p_{T} [GeV/c]", nBinEt, binsEtPtr);
  hHLTEtaHLT = dir.make<TH1F>("hHLTEtaHLT", prefix+"Pseudorapidity of HLT object matched to reco object;HLT #eta", nBinEta, binsEtaPtr);
  hHLTPhiHLT = dir.make<TH1F>("hHLTPhiHLT", prefix+"Azimuthal angle of HLT object matched to reco object;HLT #phi [Radian]", 50, -3.15, 3.15);

  hEtVsHLTEt = dir.make<TH2F>("hEtVsHLTEt", prefix+"Transverse energy of reco object vs HLT object;Reco p_{T} [GeV/c];HLT p_{T} [GeV/c]", nBinEt, binsEtPtr, nBinEt, binsEtPtr);
  hEtaVsHLTEta = dir.make<TH2F>("hEtaVsHLTEta", prefix+"Pseudorapidity of reco object vs HLT object;Reco #eta;HLT #eta", nBinEta, binsEtaPtr, nBinEta, binsEtaPtr);
  hPhiVsHLTPhi = dir.make<TH2F>("hPhiVsHLTPhi", prefix+"Azimuthal angle of reco object vs HLT object;Reco #phi [Radian];HLT #phi [Radian]", 50, -3.15, 3.15, 50, -3.15, 3.15);

  // Matching variables
  hDeltaRL1T = dir.make<TH1F>("hDeltaRL1T", prefix+"#DeltaR between reco object - L1T object;#DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}", 100, 0, 1);
  hDeltaPhiL1T = dir.make<TH1F>("hDeltaPhiL1T", prefix+"#Delta#phi between reco object - L1T object;#Delta#phi [Radian]", 100, -1, 1);
  hDeltaEtaL1T = dir.make<TH1F>("hDeltaEtaL1T", prefix+"#Delta#eta between reco object - L1T object;#Delta#eta", 100, -1, 1);
  hDeltaEtaVsDeltaPhiL1T = dir.make<TH2F>("hDeltaEtaVsDeltaPhiL1T", prefix+"#Delta#eta vs #Delta#phi between reco object - L1T object;#Delta#eta;#Delta;#phi [Radian]", 100, -1, 1, 100, -1, 1);

  hDeltaRHLT = dir.make<TH1F>("hDeltaRHLT", prefix+"#DeltaR between reco object - HLT object;#DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}", 100, 0, 1);
  hDeltaPhiHLT = dir.make<TH1F>("hDeltaPhiHLT", prefix+"#Delta#phi between reco object - HLT object;#Delta#phi [Radian]", 100, -1, 1);
  hDeltaEtaHLT = dir.make<TH1F>("hDeltaEtaHLT", prefix+"#Delta#eta between reco object - HLT object;#Delta#eta", 100, -1, 1);
  hDeltaEtaVsDeltaPhiHLT = dir.make<TH2F>("hDeltaEtaVsDeltaPhiHLT", prefix+"#Delta#eta vs #Delta#phi between reco object - HLT object;#Delta#eta;#Delta;#phi [Radian]", 100, -1, 1, 100, -1, 1);

  // Object specific histograms
  if ( objectType == ObjectType::Muon )
  {
    hQ = dir.make<TH1F>("hQ", prefix+"Electric charge;Electric charge", 3, -1.5, 1.5);

    hGlbNHit = dir.make<TH1F>("hGlbNHit", prefix+"Number of valid muon hits", 100, 0, 100);
    hTrkNHit = dir.make<TH1F>("hTrkNHit", prefix+"Number of valid tracker hits", 50, 0, 50);

    hGlbX2 = dir.make<TH1F>("hGlbX2", prefix+"Normalized #Chi^{2} of global track", 50, 0, 50);
    hTrkX2 = dir.make<TH1F>("hTrkX2", prefix+"Normalized #Chi^{2} of tracker track", 50, 0, 50);

    hRelIso = dir.make<TH1F>("hRelIso", prefix+"Relative isolation;Relative isolation", 100, 0, 10);
  }
}

void Histograms::FillReco(const reco::Candidate& recoCand)
{
  const double recoEt = recoCand.et();
  const double recoEta = recoCand.eta();
  const double recoPhi = recoCand.phi();

  hEtReco->Fill(recoEt);

  if ( recoEt < workingPointEt_ ) return;

  hEtaReco->Fill(recoEta);
  hPhiReco->Fill(recoPhi);

  if ( objectType_ == ObjectType::Muon )
  {
    const reco::Muon* recoMuonP = dynamic_cast<const reco::Muon*>(&recoCand);
    if ( !recoMuonP ) return;

    const reco::Muon& recoMuon = *recoMuonP;

    if ( !recoMuon.isGlobalMuon() or !recoMuon.isTrackerMuon() ) return;

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
    const double relIso = (trackIso+caloIso)/recoEt;

    hGlbNHit->Fill(nMuonHit);
    hGlbX2->Fill(glbX2);

    hTrkNHit->Fill(nTrkHit);
    hTrkX2->Fill(trkX2);

    hRelIso->Fill(relIso);
  }
}

void Histograms::FillL1T(const reco::Candidate& recoCand, const reco::LeafCandidate& l1Cand)
{
  if ( objectType_ == ObjectType::Muon )
  {
    edm::LogError("Histograms") << "Wrong function. Please replace to FillL1T(const reco::Candidate& recoCand, const reco::LeafCandidate& l1Cand, const double recoPosEta, const double recoPosPhi)\n";
    return;
  }

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

  hDeltaRL1T->Fill(dR);
  hDeltaEtaL1T->Fill(dEta);
  hDeltaPhiL1T->Fill(dPhi);
  hDeltaEtaVsDeltaPhiL1T->Fill(dEta, dPhi);

  if ( maxL1DeltaR_ < dR ) return;

  hEtL1T->Fill(recoEt);
  hL1EtL1T->Fill(l1Et);

  hEtVsL1Et->Fill(recoEt, l1Et);
  hEtaVsL1Eta->Fill(recoEta, l1Eta);
  hPhiVsL1Phi->Fill(recoPhi, l1Phi);

  if ( recoEt < workingPointEt_ ) return;
  
  hEtaL1T->Fill(recoEta);
  hPhiL1T->Fill(recoPhi);

  hL1EtaL1T->Fill(l1Eta);
  hL1PhiL1T->Fill(l1Phi);
}

void Histograms::FillL1T(const reco::Muon& recoMuon, const reco::LeafCandidate& l1Cand, const double recoPosEta, const double recoPosPhi)
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

  hDeltaRL1T->Fill(dR);
  hDeltaEtaL1T->Fill(dEta);
  hDeltaPhiL1T->Fill(dPhi);
  hDeltaPhiL1T->Fill(dPhi);
  hDeltaEtaVsDeltaPhiL1T->Fill(dEta, dPhi);

  if ( maxL1DeltaR_ < dR ) return;

  hEtL1T->Fill(recoEt);
  hL1EtL1T->Fill(l1Et);

  hEtVsL1Et->Fill(recoEt, l1Et);
  hEtaVsL1Eta->Fill(recoEta, l1Eta);
  hPhiVsL1Phi->Fill(recoPhi, l1Phi);

  if ( recoEt < workingPointEt_ ) return;
  
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

  hDeltaRHLT->Fill(dR);
  hDeltaEtaHLT->Fill(dEta);
  hDeltaPhiHLT->Fill(dPhi);
  hDeltaEtaVsDeltaPhiHLT->Fill(dEta, dPhi);

  if ( maxHLTDeltaR_ < dR ) return;

  hEtHLT->Fill(recoEt);
  hHLTEtHLT->Fill(hltEt);

  hEtVsHLTEt->Fill(recoEt, hltEt);
  hEtaVsHLTEta->Fill(recoEta, hltEta);
  hPhiVsHLTPhi->Fill(recoPhi, hltPhi);

  if ( recoEt < workingPointEt_ ) return;
  
  hEtaHLT->Fill(recoEta);
  hPhiHLT->Fill(recoPhi);

  hHLTEtaHLT->Fill(hltEta);
  hHLTPhiHLT->Fill(hltPhi);
}

void Histograms::reset()
{
}

void Histograms::setRecoCand(const edm::Ref<std::vector<reco::Candidate> > recoCand, 
                             const double recoPosEta, const double recoPosPhi)
{
  recoCand_ = recoCand;

  // reco position eta,phi are needed for recoMuon-L1Muon association
  recoPosEta_ = recoPosEta;
  recoPosPhi_ = recoPosPhi;
}

void Histograms::setL1Cand(const edm::Ref<std::vector<reco::LeafCandidate> > l1Cand)
{
  if ( objectType_ == ObjectType::Muon )
  {
    edm::LogError("Histograms") << "Wrong function. Please replace to FillL1T(const reco::Candidate& recoCand, const reco::LeafCandidate& l1Cand, const double recoPosEta, const double recoPosPhi)\n";
    return;
  }

  l1Cand_ = l1Cand;
}

void Histograms::setHLTCand(const edm::Ref<std::vector<trigger::TriggerObject> > hltCand)
{
  hltCand_ = hltCand;
}

void Histograms::fill()
{
  // Histograms are based on associations of reco->trigger objects
  // So nothing can be done without reco object
  if ( recoCand_.isNull() or !recoCand_.isAvailable() ) return;

  const double recoEt = recoCand_->et();
  const double recoEta = recoCand_->eta();
  const double recoPhi = recoCand_->phi();

  // Fill basic reco histograms
  // This do-while statement is dummy, to make histogram filling code
  // to be independent of other routine
  do
  {
    hEtReco->Fill(recoEt);

    if ( recoEt < workingPointEt_ ) continue;
    hEtaReco->Fill(recoEta);
    hPhiReco->Fill(recoPhi);

    // Muon specific histograms
    if ( objectType_ == ObjectType::Muon and recoCand_->isMuon() ) 
    {
      const reco::Muon* recoMuonP = dynamic_cast<const reco::Muon*>(&*recoCand_);
      if ( !recoMuonP ) continue;
      const reco::Muon& recoMuon = *recoMuonP;

      if ( !recoMuon.isGlobalMuon() or !recoMuon.isTrackerMuon() ) continue;

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
      const double relIso = (trackIso+caloIso)/recoEt;

      hGlbNHit->Fill(nMuonHit);
      hGlbX2->Fill(glbX2);

      hTrkNHit->Fill(nTrkHit);
      hTrkX2->Fill(trkX2);

      hRelIso->Fill(relIso);
    }
  } while ( false );

  // Fill L1 histograms
  do
  {
    if ( l1Cand_.isNull() or !l1Cand_.isAvailable() ) continue;

    const bool isMuon = objectType_ == ObjectType::Muon and recoCand_->isMuon();

    const double l1Et = l1Cand_->et();
    const double l1Eta = l1Cand_->eta();
    const double l1Phi = l1Cand_->phi();

    const double l1DeltaR = isMuon ? deltaR(recoPosEta_, recoPosPhi_, l1Eta, l1Phi) : deltaR(*recoCand_, *l1Cand_);
    const double l1DeltaEta = isMuon ? l1Eta - recoPosEta_ : l1Eta - recoEta;
    const double l1DeltaPhi = isMuon ? deltaPhi(l1Phi, recoPosPhi_) : deltaPhi(l1Phi, recoPhi);

    hDeltaRL1T->Fill(l1DeltaR);
    hDeltaEtaL1T->Fill(l1DeltaEta);
    hDeltaPhiL1T->Fill(l1DeltaPhi);
    hDeltaEtaVsDeltaPhiL1T->Fill(l1DeltaEta, l1DeltaPhi);

    if ( maxL1DeltaR_ < l1DeltaR ) continue;

    hEtL1T->Fill(recoEt);
    hL1EtL1T->Fill(l1Et);

    hEtVsL1Et->Fill(recoEt, l1Et);
    hEtaVsL1Eta->Fill(recoEta, l1Eta);
    hPhiVsL1Phi->Fill(recoPhi, l1Phi);

    if ( recoEt < workingPointEt_ ) continue;

    hEtaL1T->Fill(recoEta);
    hPhiL1T->Fill(recoPhi);

    hL1EtaL1T->Fill(l1Eta);
    hL1PhiL1T->Fill(l1Phi);

  } while ( false );

  // Fill HLT histograms
  do
  {
    if ( hltCand_.isNull() or !hltCand_.isAvailable() ) continue;

    const double hltEt = hltCand_->et();
    const double hltEta = hltCand_->eta();
    const double hltPhi = hltCand_->phi();
    //const double hltCharge = hltCand_->charge();

    const double hltDeltaR = deltaR(*recoCand_, *hltCand_);
    const double hltDeltaEta = hltEta - recoEta;
    const double hltDeltaPhi = deltaPhi(recoPhi, hltPhi);

    hDeltaRHLT->Fill(hltDeltaR);
    hDeltaEtaHLT->Fill(hltDeltaEta);
    hDeltaPhiHLT->Fill(hltDeltaPhi);
    hDeltaEtaVsDeltaPhiHLT->Fill(hltDeltaEta, hltDeltaPhi);

    if ( maxHLTDeltaR_ < hltDeltaR ) return;

    hEtHLT->Fill(recoEt);
    hHLTEtHLT->Fill(hltEt);

    hEtVsHLTEt->Fill(recoEt, hltEt);
    hEtaVsHLTEta->Fill(recoEta, hltEta);
    hPhiVsHLTPhi->Fill(recoPhi, hltPhi);

    if ( recoEt < workingPointEt_ ) return;

    hEtaHLT->Fill(recoEta);
    hPhiHLT->Fill(recoPhi);

    hHLTEtaHLT->Fill(hltEta);
    hHLTPhiHLT->Fill(hltPhi);

  } while ( false );
}

