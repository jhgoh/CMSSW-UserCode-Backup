#include "HiggsAnalysis/DoublyChargedHiggs/interface/HiggsToEMuAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
#include "PhysicsTools/Utilities/interface/PtComparator.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <cmath>

using namespace std;

// Purpose of this analyzer : 
// analyze pat::Electrons and pat::Muons to be used in the HiggsToEMu combination

HiggsToEMuAnalyzer::HiggsToEMuAnalyzer(const edm::ParameterSet& pset)
{
  higgs1Label_ = pset.getParameter<edm::InputTag>("higgs1Label");
  higgs2Label_ = pset.getParameter<edm::InputTag>("higgs2Label");

  edm::Service<TFileService> fs;
}

HiggsToEMuAnalyzer::~HiggsToEMuAnalyzer()
{
}

void HiggsToEMuAnalyzer::beginJob(const edm::EventSetup& eventSetup)
{
}

void HiggsToEMuAnalyzer::endJob()
{
}

void HiggsToEMuAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  // Retrieve from the event
  edm::Handle<pat::CompositeCandidateCollection> higgs1Handle;
  edm::Handle<pat::CompositeCandidateCollection> higgs2Handle;
  event.getByLabel(higgs1Label_, higgs1Handle);
  event.getByLabel(higgs2Label_, higgs2Handle);

  // Make sorted list by decreating pT
  pat::CompositeCandidateCollection higgs1Cands(higgs1Handle->size());
  pat::CompositeCandidateCollection higgs2Cands(higgs2Handle->size());
  copy(higgs1Handle->begin(), higgs1Handle->end(), higgs1Cands.begin());
  copy(higgs2Handle->begin(), higgs2Handle->end(), higgs2Cands.begin());
  std::sort(higgs1Cands.begin(), higgs1Cands.end(), GreaterByPt<pat::CompositeCandidate>());
  std::sort(higgs2Cands.begin(), higgs2Cands.end(), GreaterByPt<pat::CompositeCandidate>());

  const unsigned int nHiggs1 = higgs1Cands.size();
  const unsigned int nHiggs2 = higgs2Cands.size();

  OverlapChecker isOverlap;

  // Plot histograms
  for ( unsigned int i=0; i<nHiggs1; ++i )
  {
    const pat::CompositeCandidate& higgs1Cand = higgs1Cands[i];

    const double higgs1Pt = higgs1Cand.pt();
    const double higgs1Eta = higgs1Cand.eta();
    const double higgs1Phi = higgs1Cand.phi();
    const double higgs1Mass = higgs1Cand.mass();

    for ( unsigned int j=0; j<nHiggs2; ++j )
    {
      const pat::CompositeCandidate& higgs2Cand = higgs2Cands[j];
      if ( isOverlap(higgs1Cand, higgs2Cand) ) continue;

      const double higgs2Pt = higgs2Cand.pt();
      const double higgs2Eta = higgs2Cand.eta();
      const double higgs2Phi = higgs2Cand.phi();
      const double higgs2Mass = higgs2Cand.mass();

      const double deltaDMass = higgs1Mass - higgs2Mass;

      h1_["higgs12_deltaPhi"]->Fill(deltaPhi(higgs1Phi, higgs2Phi));
      h1_["higgs12_deltaMass"]->Fill(deltaDMass);
    }

    continue;

    h1_["higgs1_pt"]->Fill(higgs1Cand.pt());
    h1_["higgs1_eta"]->Fill(higgs1Cand.eta());
    h1_["higgs1_mass"]->Fill(higgs1Cand.mass());

    // Get daughters
    const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(higgs1Cand.daughter(0));
    const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(higgs1Cand.daughter(0));
    const pat::Electron* electron1 = dynamic_cast<const pat::Electron*>(higgs1Cand.daughter(1));
    const pat::Electron* electron2 = dynamic_cast<const pat::Electron*>(higgs1Cand.daughter(1));

    // Make first daughter to be valid before proceed to daughter combination
    // thus EMu/MuE combinations can be analyzed with single routine
    if ( muon1 == 0 and muon2 != 0 ) muon1 = muon2;
    if ( electron1 == 0 and electron2 != 0 ) electron1 = electron2;

    // Consider daughter combinations
    // E-Mu channel
    if ( muon1 != 0 and electron1 != 0 )
    {
      // Fill muon variables
      h1_["muon1_pt"]->Fill(muon1->pt());
      h1_["muon1_eta"]->Fill(muon1->eta());
      h1_["muon1_relIso"]->Fill((muon1->trackIso()+muon1->caloIso())/muon1->pt());

      h1_["electron1_pt"]->Fill(electron1->pt());
      h1_["electron1_eta"]->Fill(electron1->eta());

      h1_["higgs1_deltaPhi"]->Fill(deltaPhi(muon1->phi(), electron1->phi()));
    }
  }
}
