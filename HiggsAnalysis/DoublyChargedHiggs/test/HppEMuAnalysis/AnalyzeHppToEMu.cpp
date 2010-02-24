#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"

#endif

#include "TROOT.h"
#include "TSystem.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "TMath.h"
#include "TH1F.h"
#include "THStack.h"

#include <iostream>
#include <vector>

using namespace std;

const double zMassPDG = 91.1876;
const double zVetoCut = 5;
const double higgsMassMin = 50;

void AnalyzeHppToEMu(TString sampleName, int verbose=0)
{
  // Output file
  TString fileName = TString("res/")+TString(sampleName)+TString(".root");
  TFile* hFile = TFile::Open(fileName, "RECREATE");

  // Book histograms
  TH1F* hElectronPt = new TH1F("ElectronPt", "Electron transverse momentum;p_{T} [GeV/c]", 200, 0, 1000);
  TH1F* hElectronEta = new TH1F("ElectronEta", "Electron pseudorapidity;#eta [Radian]", 200, -3, 3);
  TH1F* hElectronIso = new TH1F("ElectronIso", "Electron relative isolation;Isolation", 200, 0, 20);

  TH1F* hMuonPt = new TH1F("MuonPt", "Muon transverse momentum;p_{T} [GeV/c]", 200, 0, 1000);
  TH1F* hMuonEta = new TH1F("MuonEta", "Muon pseudorapidity;#eta [Radian]", 200, -3, 3);
  TH1F* hMuonIso = new TH1F("MuonIso", "Muon relative isolation;Isolation", 200, 0, 20);

  TH1F* hHppPt = new TH1F("HppPt", "H^{++} transverse momentum;p_{T} [GeV/c]", 200, 0, 1000);
  TH1F* hHmmPt = new TH1F("HmmPt", "H^{--} transverse mementum;p_{T} [GeV/c]", 200, 0, 1000);

  TH1F* hHppMass = new TH1F("HppMass", "H^{++} mass;Mass [GeV/c^{2}]", 200, 0, 1000);
  TH1F* hHmmMass = new TH1F("HmmMass", "H^{--} mass;Mass [GeV/c^{2}]", 200, 0, 1000);
  TH1F* hMassDiff = new TH1F("MassDiff", "Higgs mass difference;Mass difference [GeV/c^{2}]", 200, 0, 1000);

  TH1F* hZToMuMuMass = new TH1F("ZToMuMuMass", "M(#mu^{+}#mu^{-}));Mass difference [GeV/c^{2}]", 200, 0, 200);
  TH1F* hZToEEMass = new TH1F("ZToEEMass", "M(e^{+}e^{-}));Mass difference [GeV/c^{2}]", 200, 0, 200);

  TH1F* hN4lCand = new TH1F("N4lCand", "Number of 4 lepton candidates per event", 5, -0.5, 4.5);
  TH1F* hHppMassNearest = new TH1F("HppMassNearest", "H^{++} mass with nearest selection;Mass [GeV/c^{2}]", 200, 0, 1000);
  TH1F* hHmmMassNearest = new TH1F("HmmMassNearest", "H^{--} mass with nearest selection;Mass [GeV/c^{2}]", 200, 0, 1000);

  // Input samples
  std::vector<std::string> inputFiles;

  {
    TString datasetPath("HiggsAnalysis/DoublyChargedHiggs/");
    datasetPath += getenv("CMSSW_VERSION");
    datasetPath += "/"+sampleName+"/HppEMu/";

    FILE* rfdirCmd = gSystem->OpenPipe("nsls "+datasetPath, "r");
    const size_t bufferSize = 1000;
    char dirContent[bufferSize];
    while ( fgets(dirContent, bufferSize, rfdirCmd) )
    {
      std::string fileName(dirContent);
      if ( fileName.rfind(".root\n") != fileName.npos )
      {
        fileName.erase(fileName.size()-1);
        fileName = "rfio:/castor/cern.ch/user/j/jhgoh/"+datasetPath+fileName;

        inputFiles.push_back(fileName);
      }
    }
  }

  fwlite::ChainEvent event(inputFiles);

  const int eventSize = event.size();
  int eventNumber = 0;

  // Loop over all events
  for(event.toBegin(); !event.atEnd(); ++event, ++eventNumber)
  {
    if ( verbose > 0 )
    {
      const int eventFrac = TMath::Nint(1000.*eventNumber/eventSize);
      if ( eventFrac%10 == 0 ) cout << "@@ Analyzing " << sampleName << " : " << eventFrac/10 << "\% done\r";
    }

    // Retrieve objects from the event
    fwlite::Handle<pat::ElectronCollection> electronColl;
    fwlite::Handle<pat::MuonCollection> muonColl;
    fwlite::Handle<pat::CompositeCandidateCollection> posDeltaColl;
    fwlite::Handle<pat::CompositeCandidateCollection> negDeltaColl;

    electronColl.getByLabel(event, "goodPatElectrons");
    muonColl.getByLabel(event, "goodPatMuons");
    posDeltaColl.getByLabel(event, "posDeltaToEMu");
    negDeltaColl.getByLabel(event, "negDeltaToEMu");

    // Draw basic kinematic distributions
    typedef reco::CandidateCollection::const_iterator CandIter;
    typedef pat::ElectronCollection::const_iterator ElectronIter;
    typedef pat::MuonCollection::const_iterator MuonIter;
    typedef pat::CompositeCandidateCollection::const_iterator CompCandIter;

    OverlapChecker isOverlap;

    for(ElectronIter electronPtr = electronColl->begin();
        electronPtr != electronColl->end(); ++electronPtr)
    {
      hElectronPt->Fill(electronPtr->pt());
      hElectronEta->Fill(electronPtr->eta());
      hElectronIso->Fill((electronPtr->trackIso()+electronPtr->caloIso())/electronPtr->pt());
    }

    for(MuonIter muonPtr = muonColl->begin();
        muonPtr != muonColl->end(); ++muonPtr)
    {
      hMuonPt->Fill(muonPtr->pt());
      hMuonEta->Fill(muonPtr->eta());
      hMuonIso->Fill((muonPtr->trackIso()+muonPtr->caloIso())/muonPtr->pt());
    }

    for(CompCandIter hppCand = posDeltaColl->begin(); 
        hppCand != posDeltaColl->end(); ++hppCand)
    {
      hHppPt->Fill(hppCand->pt());
      hHppMass->Fill(hppCand->mass());
    }

    for(CompCandIter hmmCand = negDeltaColl->begin(); 
        hmmCand != negDeltaColl->end(); ++hmmCand)
    {
      hHmmPt->Fill(hmmCand->pt());
      hHmmMass->Fill(hmmCand->mass());
    }

    // Consider H++/H-- combination
    int n4lCand = 0;
    double massDiffMin = 1e5;
    CompCandIter hppCandNearest = posDeltaColl->end();
    CompCandIter hmmCandNearest = negDeltaColl->end();

    for(CompCandIter hppCand = posDeltaColl->begin(); 
        hppCand != posDeltaColl->end(); ++hppCand)
    {
      if ( hppCand->numberOfDaughters() != 2 ) continue; 
      if ( hppCand->mass() < higgsMassMin ) continue;

      const int muonIdx = hppCand->daughter(0)->isMuon() ? 0 : 1;
      const int electronIdx = hppCand->daughter(1)->isElectron() ? 1 : 0;

      const pat::Muon* posMuon = dynamic_cast<const pat::Muon*>(hppCand->daughter(muonIdx));
      const pat::Electron* posElectron = dynamic_cast<const pat::Electron*>(hppCand->daughter(electronIdx));

      // Skip if neither E-Mu nor Mu-E pair
      if ( !posMuon or !posElectron ) continue;

      // Confirm Skimming cuts
      if ( posMuon->pt() < 10 or posElectron->pt() < 10 or
           fabs(posMuon->eta()) > 2.5 or fabs(posElectron->eta()) > 2.5 ) continue;

      for(CompCandIter hmmCand = negDeltaColl->begin(); 
          hmmCand != negDeltaColl->end(); ++hmmCand)
      {
        if ( hmmCand->numberOfDaughters() != 2 ) continue;
        if ( hmmCand->mass() < higgsMassMin ) continue;

        const int muonIdx = hmmCand->daughter(0)->isMuon() ? 0 : 1;
        const int electronIdx = hmmCand->daughter(1)->isElectron() ? 1 : 0;

        const pat::Muon* negMuon = dynamic_cast<const pat::Muon*>(hmmCand->daughter(muonIdx));
        const pat::Electron* negElectron = dynamic_cast<const pat::Electron*>(hmmCand->daughter(electronIdx));

        // Skip if neigher E-Mu nor Mu-E pair
        if ( !negMuon or !negElectron ) continue;

        // Confirm Skimming cuts
        if ( negMuon->pt() < 10 or negElectron->pt() < 10 or
             fabs(negMuon->eta()) > 2.5 or fabs(negElectron->eta()) > 2.5 ) continue;

        // Check overlap
        if ( isOverlap(*hppCand, *hmmCand) ) continue;

        const math::XYZTLorentzVector dimuonLVec = negMuon->p4() + posMuon->p4();
        const math::XYZTLorentzVector dielectronLVec = negElectron->p4() + posElectron->p4();
        
        const double zToMuMuMass = dimuonLVec.M();
        const double zToEEMass = dielectronLVec.M();

        hZToMuMuMass->Fill(zToMuMuMass);
        hZToEEMass->Fill(zToEEMass);

        // Z veto
        if ( zToMuMuMass > zMassPDG-zVetoCut and zToMuMuMass < zMassPDG+zVetoCut ) continue;
        if ( zToEEMass > zMassPDG-zVetoCut and zToEEMass < zMassPDG+zVetoCut ) continue;

        const double massDiff = fabs(hppCand->mass()-hmmCand->mass());
        hMassDiff->Fill(massDiff);

        ++n4lCand;
        if ( massDiffMin > massDiff )
        {
          massDiffMin = massDiff;
          hppCandNearest = hppCand;
          hmmCandNearest = hmmCand;
        }
      }
    }

    // Fill best candidate pair
    hN4lCand->Fill(n4lCand);
    if ( n4lCand > 0 )
    {
      hHppMassNearest->Fill(hppCandNearest->mass());
      hHmmMassNearest->Fill(hmmCandNearest->mass());
    }

  }

  hFile->Write();

  if ( verbose > 1 )
  {
    TCanvas* c;

    hHppPt->SetLineColor(kRed)  ; hHmmPt->SetLineColor(kBlue)  ;
    hHppMass->SetLineColor(kRed); hHmmMass->SetLineColor(kBlue);

    c = new TCanvas; hHppPt->Draw()  ; hHmmPt->Draw("same")  ;
    c = new TCanvas; hHppMass->Draw(); hHmmMass->Draw("same");
    c = new TCanvas; hMassDiff->Draw();
  }
  else
  {
    hFile->Close();
  }
}
