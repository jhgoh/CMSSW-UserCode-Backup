#include "TH1F.h"
#include "THStack.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TROOT.h"

#include <iostream>
#include <vector>

#include "HppXSec.cpp"

using namespace std;

void view(const double hppMass = 200, const char* datasetNameForm = "Hpp%d_EMu_10TeV_GEN_HLT")
{
  cout << "=== View cuts ===" << endl;

  const double totalLumi = 1; // Unit in fb^{-1}
  
  // List of samples
  const size_t nBkg = 3;
  TString bkgNames[nBkg] = {"TT_4l", "LLBB_4l", "ZZ_4l"};
  const int bkgColors[nBkg] = {kGreen, kBlue, kMagenta};
  // XSec parameters from https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsZZMCsamples
  const double bkgXSec[nBkg] = {280900*1.46*0.01091, 56200*1.66*0.007, 189*(1+0.35+0.2)*0.3165};
  const double bkgNEvents[nBkg] = {1007062, 1063204, 898940};
  //const double bkgNEvents[nBkg] = {10000, 10000, 10000};
  TFile* bkgFiles[nBkg];

  for ( size_t bkgIdx = 0; bkgIdx < nBkg; ++bkgIdx )
  {
    bkgFiles[bkgIdx] = new TFile("res/"+bkgNames[bkgIdx]+"_10TeV_GEN.root");
  }

  TString datasetName = Form(datasetNameForm, int(hppMass));
  TString sigFileName = TString("res/")+datasetName+".root";
  int sigColor = kRed;
  const double sigNEvents = 10000;
  const double sigXSec = getHppXSec(hppMass);
  TFile* sigFile;
  sigFile = new TFile(sigFileName);
  
  // List of histograms
  //TString histPathPrefix = "hToEMuPromptMu";
  const size_t nHist = 3+3+5+2+3;
  const char* histPaths[nHist] = {
    "MuonPt", "MuonEta", "MuonIso",
    "ElectronPt", "ElectronEta", "ElectronIso",
    "HppMass", "HmmMass", "HppPt", "HmmPt", "MassDiff",
    "ZToEEMass", "ZToMuMuMass",
    "N4lCand", "HppMassNearest", "HmmMassNearest"
  };

  // Draw histograms 
  for ( size_t histIdx = 0; histIdx < nHist; ++histIdx )
  {
    TString histPath = histPaths[histIdx];

    gROOT->cd();
    TCanvas* c = new TCanvas(histPaths[histIdx], histPaths[histIdx]);
    THStack* hStack = new THStack(histPaths[histIdx], histPaths[histIdx]);

    double nEntries = 0;

    for ( size_t bkgIdx = 0; bkgIdx < nBkg; ++bkgIdx )
    {
      TFile* f = bkgFiles[bkgIdx];
      f->cd();
      TH1F* h = (TH1F*)f->Get(histPath);

      if ( !h )
      {
        cerr << "Cannot find histogram " << histPath << endl;
        continue;
      }
      else
      {
        h->Scale(bkgXSec[bkgIdx]/bkgNEvents[bkgIdx]);
        nEntries += h->GetEntries();
        h->SetFillColor(bkgColors[bkgIdx]);
        hStack->Add(h);
      }
    }

    sigFile->cd();
    TH1F* hSig = (TH1F*)sigFile->Get(histPath);

    if ( !hSig )
    {
      cerr << "Cannot find histogram " << histPath << endl;
    }
    else
    {
      hSig->Scale(sigXSec/sigNEvents);
      //hSig->SetFillColor(sigColor);
      hSig->SetLineColor(sigColor);
      hSig->SetLineWidth(2);
    }

    TString histName = datasetName+"_"+hStack->GetName();
    histName.ReplaceAll("/", "_");
    //histName.ToLower();
    if ( nEntries > 0 )
    {
      if ( histName.Contains("Iso") || 
           histName.Contains("Pt") ) 
      {
        cout << "Plot with log scales" << endl;
        c->SetLogy();
      }
    }

    hStack->Draw();
    hSig->Draw("SAME");
    c->Print("gif/"+histName+".gif");
  }
}
