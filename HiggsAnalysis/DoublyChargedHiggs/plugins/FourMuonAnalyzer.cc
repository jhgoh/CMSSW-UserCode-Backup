#include "HiggsAnalysis/DoublyChargedHiggs/interface/FourMuonAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"

//#include "DataFormats/Candidate/interface/CompositeCandidate.h"
//#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "PhysicsTools/Utilities/interface/PtComparator.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <cmath>

using namespace std;

namespace H1
{
  enum ID
  {
    posDelta_pt, posDelta_eta, posDelta_mass,
    negDelta_pt, negDelta_eta, negDelta_mass,
    
    posDelta_dZ, posDelta_dPhi,
    negDelta_dZ, negDelta_dPhi,
    
    posDelta_trackIso, posDelta_caloIso, posDelta_ecalIso, posDelta_hcalIso,
    negDelta_trackIso, negDelta_caloIso, negDelta_ecalIso, negDelta_hcalIso,

    posDelta_relIso,
    negDelta_relIso,

    zStar_mass,
    zStar_dMass,
    zStar_dPhi,
    zStar_minDZMass
  };
}

FourMuonAnalyzer::FourMuonAnalyzer(const edm::ParameterSet& pset)
{
  // Set composite particle candidate labels
  posDeltaLabel_ = pset.getParameter<edm::InputTag>("posDelta"); 
  negDeltaLabel_ = pset.getParameter<edm::InputTag>("negDelta");

  // Preselection variables
  edm::ParameterSet deltaCutSet = pset.getParameter<edm::ParameterSet>("deltaCutSet");
  delta_maxNormalizedChi2_ = deltaCutSet.getUntrackedParameter<double>("maxNormalizedChi2");
  delta_minPt_ = deltaCutSet.getUntrackedParameter<double>("minPt");

  nInterested_ = pset.getUntrackedParameter<unsigned int>("nInterested");

  // Invoke TFileService
  edm::Service<TFileService> fs;
  fsStatus_ = fs.isAvailable();
  if ( fsStatus_ )
  {
    // Book histograms
    h1_[H1::posDelta_pt] = fs->make<TH1F>("posDelta_pt", "#Delta^{++} p_{T};Transverse momentum [GeV/c]", 100, 0, 500);
    h1_[H1::posDelta_eta] = fs->make<TH1F>("posDelta_eta", "#Delta^{++} #eta;Pseudorapidity [Radian]", 100, -2.4, 2.4);
    h1_[H1::posDelta_mass] = fs->make<TH1F>("posDelta_mass", "#Delta^{++} mass;Mass [GeV/c^{2}]", 100, 0, 300);

    h1_[H1::negDelta_pt] = fs->make<TH1F>("negDelta_pt", "#Delta^{--} p_{T};Transverse momentum [GeV/c]", 100, 0, 500);
    h1_[H1::negDelta_eta] = fs->make<TH1F>("negDelta_eta", "#Delta^{--} #eta;Pseudorapidity [Radian]", 100, -2.4, 2.4);
    h1_[H1::negDelta_mass] = fs->make<TH1F>("negDelta_mass", "#Delta^{--} mass;Mass [GeV/c^{2}]", 100, 0, 300);

    h1_[H1::posDelta_dZ] = fs->make<TH1F>("posDelta_dZ", "#Delta^{++} dimuon #Delta z;#Delta;z [cm]", 100, 0, 1);
    h1_[H1::posDelta_dPhi] = fs->make<TH1F>("posDelta_dPhi", "#Delta^{++} dimuon #Delta#phi;#Delta#phi [Radian]", 100, 0, TMath::Pi());

    h1_[H1::negDelta_dZ] = fs->make<TH1F>("negDelta_dZ", "#Delta^{--} dimuon #Delta z;dz [cm]", 100, 0, 1);
    h1_[H1::negDelta_dPhi] = fs->make<TH1F>("negDelta_dPhi", "#Delta^{--} dimuon #Delta#phi;#Delta#phi [Radian]", 100, 0, TMath::Pi());

    h1_[H1::posDelta_trackIso] = fs->make<TH1F>("posDelta_trackIso", "Muon track isolation from #Delta^{++};Track Isolation [GeV]", 100, -1, 10);
    h1_[H1::posDelta_caloIso] = fs->make<TH1F>("posDelta_caloIso", "Muon calorimeter (ecal+hcal) isolation from #Delta^{++};Calorimeter Isolation [GeV]", 100, -1, 10);
    h1_[H1::posDelta_ecalIso] = fs->make<TH1F>("posDelta_ecalIso", "Muon electromagnetic calorimeter isolation from #Delta^{++};ECal Isolation [GeV]", 100, -1, 10);
    h1_[H1::posDelta_hcalIso] = fs->make<TH1F>("posDelta_hcalIso", "Muon Hadronic calorimeter isolation from #Delta^{++};HCal Isolation [GeV]", 100, -1, 10);
    h1_[H1::posDelta_relIso] = fs->make<TH1F>("posDelta_relIso", "Muon relative isolation (trackIso+caloIso)/pt;Relative isolation", 100, 0, 5);

    h1_[H1::negDelta_trackIso] = fs->make<TH1F>("negDelta_trackIso", "Muon track isolation from #Delta^{--};Track Isolation [GeV]", 100, -1, 10);
    h1_[H1::negDelta_caloIso] = fs->make<TH1F>("negDelta_caloIso", "Muon calorimeter (ecal+hcal) isolation from #Delta^{--};Calorimeter Isolation [GeV]", 100, -1, 10);
    h1_[H1::negDelta_ecalIso] = fs->make<TH1F>("negDelta_ecalIso", "Muon electromagnetic calorimeter isolation from #Delta^{--};ECal Isolation [GeV]", 100, -1, 10);
    h1_[H1::negDelta_hcalIso] = fs->make<TH1F>("negDelta_hcalIso", "Muon Hadronic calorimeter isolation from #Delta^{--};HCal Isolation [GeV]", 100, -1, 10);
    h1_[H1::negDelta_relIso] = fs->make<TH1F>("negDelta_relIso", "Muon relative isolation (trackIso+caloIso)/p_{T};Relative isolation", 100, 0, 5);

    h1_[H1::zStar_mass] = fs->make<TH1F>("zStar_mass", "Z* mass;Mass [GeV/c^{2}]", 100, 0, 1000);
    h1_[H1::zStar_dMass] = fs->make<TH1F>("zStar_dMass", "Mass difference between #Delta^{++}, #Delta^{--}", 100, 0, 300);
    h1_[H1::zStar_dPhi] = fs->make<TH1F>("zStar_dPhi", "#Delta#phi between #Delta^{++}, #Delta^{--}", 100, 0, TMath::Pi());
    h1_[H1::zStar_minDZMass] = fs->make<TH1F>("zStar_minDZMass", "min(Mass(#mu^{+}#mu^{-})-Mass(Z_{PDG}));Mass difference [GeV/c^{2}]", 100, -20, 20);
  }
}


FourMuonAnalyzer::~FourMuonAnalyzer()
{
}

void FourMuonAnalyzer::beginJob(const edm::EventSetup& eventSetup)
{
}

void FourMuonAnalyzer::endJob()
{
}

void FourMuonAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<pat::CompositeCandidateCollection> posDeltaHandle, negDeltaHandle;
  event.getByLabel(posDeltaLabel_, posDeltaHandle);
  event.getByLabel(negDeltaLabel_, negDeltaHandle);

  // Make sorted list by pT
  vector<pat::CompositeCandidate> posDeltaCands(posDeltaHandle->size());
  vector<pat::CompositeCandidate> negDeltaCands(negDeltaHandle->size());
  copy(posDeltaHandle->begin(), posDeltaHandle->end(), posDeltaCands.begin());
  copy(negDeltaHandle->begin(), negDeltaHandle->end(), negDeltaCands.begin());
  sort(posDeltaCands.begin(), posDeltaCands.end(), GreaterByPt<reco::Candidate>());
  sort(negDeltaCands.begin(), negDeltaCands.end(), GreaterByPt<reco::Candidate>());

  // Start 4muon combination
  const unsigned int nPosDelta = min(posDeltaCands.size(), nInterested_);
  const unsigned int nNegDelta = min(negDeltaCands.size(), nInterested_);

  for(unsigned int i=0; i<nPosDelta; ++i)
  {
    // Fill dimuon histograms
    const pat::CompositeCandidate& posDelta = posDeltaCands[i];

    h1_[H1::posDelta_pt]->Fill(posDelta.pt());
    h1_[H1::posDelta_eta]->Fill(posDelta.eta());
    h1_[H1::posDelta_mass]->Fill(posDelta.mass());

    // Fill muon track histograms
    const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(posDelta.daughter(0));
    const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(posDelta.daughter(1));

    // Fill muon isolation information
    h1_[H1::posDelta_trackIso]->Fill(muon1->trackIso());
    h1_[H1::posDelta_trackIso]->Fill(muon2->trackIso());
    h1_[H1::posDelta_caloIso]->Fill(muon1->caloIso());
    h1_[H1::posDelta_caloIso]->Fill(muon2->caloIso());
    h1_[H1::posDelta_ecalIso]->Fill(muon1->ecalIso());
    h1_[H1::posDelta_ecalIso]->Fill(muon2->ecalIso());
    h1_[H1::posDelta_hcalIso]->Fill(muon1->hcalIso());
    h1_[H1::posDelta_hcalIso]->Fill(muon2->hcalIso());

    // Relative isolation suggested in SWGuideMuonIsolation
    h1_[H1::posDelta_relIso]->Fill((muon1->trackIso()+muon1->caloIso())/muon1->pt());
    h1_[H1::posDelta_relIso]->Fill((muon2->trackIso()+muon2->caloIso())/muon2->pt());

    // Correlation between two muons
    h1_[H1::posDelta_dZ]->Fill(fabs(muon1->vz() - muon2->vz()));
    h1_[H1::posDelta_dPhi]->Fill(fabs(deltaPhi(muon1->phi(), muon2->phi())));
  }

  for(unsigned int i=0; i<nNegDelta; ++i)
  {
    // Fill dimuon histograms
    const pat::CompositeCandidate& negDelta = negDeltaCands[i];

    h1_[H1::negDelta_pt]->Fill(negDelta.pt());
    h1_[H1::negDelta_eta]->Fill(negDelta.eta());
    h1_[H1::negDelta_mass]->Fill(negDelta.mass());

    // Fill muon track histograms
    const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(negDelta.daughter(0));
    const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(negDelta.daughter(1));

    // Correlation between two muons
    h1_[H1::negDelta_dZ]->Fill(fabs(muon1->vz()-muon2->vz()));
    h1_[H1::negDelta_dPhi]->Fill(fabs(deltaPhi(muon1->phi(), muon2->phi())));

    // Fill muon isolation information
    h1_[H1::negDelta_trackIso]->Fill(muon1->trackIso());
    h1_[H1::negDelta_trackIso]->Fill(muon2->trackIso());
    h1_[H1::negDelta_caloIso]->Fill(muon1->caloIso());
    h1_[H1::negDelta_caloIso]->Fill(muon2->caloIso());
    h1_[H1::negDelta_ecalIso]->Fill(muon1->ecalIso());
    h1_[H1::negDelta_ecalIso]->Fill(muon2->ecalIso());
    h1_[H1::negDelta_hcalIso]->Fill(muon1->hcalIso());
    h1_[H1::negDelta_hcalIso]->Fill(muon2->hcalIso());

    // Relative isolation suggested in SWGuideMuonIsolation
    h1_[H1::negDelta_relIso]->Fill((muon1->trackIso()+muon1->caloIso())/muon1->pt());
    h1_[H1::negDelta_relIso]->Fill((muon2->trackIso()+muon2->caloIso())/muon2->pt());
  }

  // Z* combination
  for(unsigned int i=0; i<nPosDelta; ++i)
  {
    const pat::CompositeCandidate& posDelta = posDeltaCands[i];

    for(unsigned int j=0; j<nNegDelta; ++j)
    {
      const pat::CompositeCandidate& negDelta = negDeltaCands[j];

      pat::CompositeCandidate zStarCand;
      zStarCand.addDaughter(posDelta);
      zStarCand.addDaughter(negDelta);

      AddFourMomenta addP4;
      addP4.set(zStarCand);

      h1_[H1::zStar_dMass]->Fill(fabs(posDelta.mass() - negDelta.mass()));
      h1_[H1::zStar_dPhi]->Fill(fabs(deltaPhi(posDelta.phi(), negDelta.phi())));
      h1_[H1::zStar_mass]->Fill(zStarCand.mass());

      // Try neutral sign combinations
      double minDZMass = 1e+20;
      const static double zMassPDG = 91.1876;
      for(unsigned int zi=0; zi<2; ++zi)
      {
        for(unsigned int zj=0; zj<2; ++zj)
        {
          pat::CompositeCandidate zCand;
          zCand.addDaughter(*posDelta.daughter(zi));
          zCand.addDaughter(*negDelta.daughter(zj));

          AddFourMomenta addZP4;
          addZP4.set(zCand);

          const double dzMass = zCand.mass() - zMassPDG;
          if ( dzMass < minDZMass ) minDZMass = dzMass;
        }
      }

      h1_[H1::zStar_minDZMass]->Fill(minDZMass);
    }
  }
}

/* vim:set ts=2 sts=2 sw=2 expandtab: */
