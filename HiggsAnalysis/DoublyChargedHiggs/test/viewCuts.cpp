#include "TH1F.h"
#include "THStack.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TROOT.h"

#include <iostream>
#include <vector>

using namespace std;

void viewCuts()
{
  cout << "=== View cuts ===" << endl;
  
  // List of samples
  TString sigName = "HppEMu/PYTHIA6_HppEMu_M200_10TeV_IDEAL_V12";
  int sigColor = kBlue;
  TFile* sigFile;
  sigFile = new TFile(sigName+"/hist.root");

  const size_t nBkg = 3;
  TString bkgNames[nBkg] = {"TT_4l", "LLBB_4l", "ZZ_4l"};
  int bkgColors[nBkg] = {kGreen, kBlue, kMagenta};
  TFile* bkgFiles[nBkg];

  for ( size_t bkgIdx = 0; bkgIdx < nBkg; ++bkgIdx )
  {
    bkgFiles[bkgIdx] = new TFile("result/hist_"+bkgNames[bkgIdx]+".root");
  }
  
  // List of histograms
  //TString histPathPrefix = "hToEMuPromptMu";
  TString histPathPrefix = "hToEMuAllMu";
  const size_t nHist = 51;
  const char* histPaths[nHist] = {
    "m1/pt", "m1/pt1_pt", "m1/pt2_pt", "m1/pt3_pt",
    "m1/trackIso", "m1/pt1_trackIso", "m1/pt2_trackIso", "m1/pt3_trackIso", 
    "m1/caloIso", "m1/pt1_caloIso", "m1/pt2_caloIso", "m1/pt3_caloIso", 
    "m1/relIso", "m1/pt1_relIso", "m1/pt2_relIso", "m1/pt3_relIso", 
 
    "e1/pt", "e1/pt1_pt", "e1/pt2_pt", "e1/pt3_pt",
    "e1/trackIso", "e1/pt1_trackIso", "e1/pt2_trackIso", "e1/pt3_trackIso", 
    "e1/caloIso", "e1/pt1_caloIso", "e1/pt2_caloIso", "e1/pt3_caloIso", 
    "e1/relIso", "e1/pt1_relIso", "e1/pt2_relIso", "e1/pt3_relIso", 

    "e1/robustLoose", "e1/pt1_robustLoose", "e1/pt2_robustLoose", "e1/pt3_robustLoose",
    "e1/robustTight", "e1/pt1_robustTight", "e1/pt2_robustTight", "e1/pt3_robustTight",

    "emu/rawMass", "emu/rawPt", "emu/mass", "emu/pt", "emu/cosDeltaPhi",
//    "emu/resolM", "emu/resolPt",
    "emu/fitChi2", "emu/trackDZ",

    "emuemu/mass", "emuemu/massDiff",
//    "emuemu/massDiffSignif",
    "emuemu/eeMass", "emuemu/mumuMass"
  };

  // Draw histograms 
  for ( size_t histIdx = 0; histIdx < nHist; ++histIdx )
  {
    TString histPath = histPathPrefix+"/"+histPaths[histIdx];

    gROOT->cd();
    TCanvas* c = new TCanvas(histPaths[histIdx], histPaths[histIdx]);
    THStack* hStack = new THStack(histPaths[histIdx], histPaths[histIdx]);

    double nEntries = 0;

    sigFile->cd();
    TH1F* hSig = (TH1F*)sigFile->Get(histPath);

    if ( !hSig )
    {
      cerr << "Cannot find histogram " << histPath << endl;
    }
    else
    {
      nEntries += hSig->GetEntries();
      hSig->SetFillColor(sigColor);
      hStack->Add(hSig);
    }

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
        nEntries += h->GetEntries();
        h->SetFillColor(bkgColors[bkgIdx]);
        hStack->Add(h);
      }
    }

    TString histName(hStack->GetName());
    histName.ReplaceAll("/", "_");
    histName.ToLower();
    if ( nEntries > 0 )
    {
      if ( histName.Contains("_reliso") || 
           histName.Contains("_caloiso") || histName.Contains("_trackiso") ||
           histName.Contains("_rawpt") || histName.Contains("_pt") ) 
      {
        cout << "Plot with log scales" << endl;
        c->SetLogy();
      }
    }

    hStack->Draw();
    c->Print(histName+".gif");
  }
}
