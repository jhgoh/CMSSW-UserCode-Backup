#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"
#include "THStack.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <vector>

void compareReleases()
{
  using namespace std;

  const TString dataset = "RelValSingleMuPt1000";

  const int cWidth = 1024;
  const int cHeight = 768;

  // Open files
  vector<TString> releases;
  vector<TFile*> files;
  void* dataDir = gSystem->OpenDirectory(".");
  if ( !dataDir ) return;
  while ( true )
  {
    TString dir = gSystem->GetDirEntry(dataDir);

    if ( dir.IsNull() ) break; 
    if ( !dir.Contains("CMSSW_") ) continue;

    TString fileName = dir+"/src/Validation/RPCRecHits/test/dqm_"+dataset+".root";
    TFile* file = TFile::Open(fileName);

    if ( !file || !file->IsOpen() ) continue;

    releases.push_back(dir);
    files.push_back(file);
  }
  const int nReleases = releases.size();
  const int nFiles = files.size();

  // List of histograms to draw
  vector<TString> h1Names;
  h1Names.push_back("globalEfficiencies");
  h1Names.push_back("Effic_Wheel");
  h1Names.push_back("Effic_Disk");
  h1Names.push_back("NoiseRate_Wheel");
  h1Names.push_back("NoiseRate_Disk");
  h1Names.push_back("LostRate_Wheel");
  h1Names.push_back("LostRate_Disk");

  h1Names.push_back("Res");

  h1Names.push_back("Res_WM2");
  h1Names.push_back("Res_WM1");
  h1Names.push_back("Res_W00");
  h1Names.push_back("Res_WP1");
  h1Names.push_back("Res_WP2");

  h1Names.push_back("Res_DM3");
  h1Names.push_back("Res_DM2");
  h1Names.push_back("Res_DM1");
  h1Names.push_back("Res_DP1");
  h1Names.push_back("Res_DP2");
  h1Names.push_back("Res_DP3");

  h1Names.push_back("Pull");

  h1Names.push_back("Pull_WM2");
  h1Names.push_back("Pull_WM1");
  h1Names.push_back("Pull_W00");
  h1Names.push_back("Pull_WP1");
  h1Names.push_back("Pull_WP2");

  h1Names.push_back("Pull_DM3");
  h1Names.push_back("Pull_DM2");
  h1Names.push_back("Pull_DM1");
  h1Names.push_back("Pull_DP1");
  h1Names.push_back("Pull_DP2");
  h1Names.push_back("Pull_DP3");

  // Draw histograms
  for ( vector<TString>::const_iterator h1Name = h1Names.begin();
        h1Name != h1Names.end(); ++h1Name )
  {
    TString canvasName = "c"+*h1Name;
    TCanvas* c = new TCanvas(canvasName, canvasName, cWidth, cHeight);
    c->cd();

    THStack* hStack = new THStack;
    TLegend* leg = new TLegend(0.65, 0.2, 0.85, 0.5);

    for ( int i=0; i<nFiles; ++i )
    {
      TFile* file = files[i];
      if ( !file || !file->IsOpen() )
      {
        cout << "File " << file->GetName() << " not opened" << endl;
        continue;
      }

      TH1F* hist = dynamic_cast<TH1F*>(file->Get("DQMData/RPCRecHitsV/"+*h1Name));
      if ( !hist )
      {
        cout << "Cannot open histogram " << *h1Name << endl;
        continue;
      }

      TString histName = hist->GetName();
      if ( histName.Contains("Effic") or
           histName.Contains("effic") or
           histName.Contains("Rate") or
           histName.Contains("rate") )
      {
        hist->SetMinimum(0);
        hist->SetMaximum(1.1);
      }

      hist->SetLineColor(i+2);

      hStack->Add(hist);
      leg->AddEntry(hist, releases[i]);

      if ( i == 0 )
      {
        hStack->SetTitle(TString(hist->GetTitle())+";"+hist->GetXaxis()->GetTitle());
      }
    }

    hStack->Draw("nostack");
    leg->Draw();

    c->Print(TString("gif/")+c->GetName()+".gif");
  }

  // Resolutions, evolution in SW releases
  vector<TString> fitNames;
  fitNames.push_back("Res");

  fitNames.push_back("Res_WM2");
  fitNames.push_back("Res_WM1");
  fitNames.push_back("Res_W00");
  fitNames.push_back("Res_WP1");
  fitNames.push_back("Res_WP2");

  fitNames.push_back("Res_DM3");
  fitNames.push_back("Res_DM2");
  fitNames.push_back("Res_DM1");
  fitNames.push_back("Res_DP1");
  fitNames.push_back("Res_DP2");
  fitNames.push_back("Res_DP3");

  fitNames.push_back("Pull");

  fitNames.push_back("Pull_WM2");
  fitNames.push_back("Pull_WM1");
  fitNames.push_back("Pull_W00");
  fitNames.push_back("Pull_WP1");
  fitNames.push_back("Pull_WP2");

  fitNames.push_back("Pull_DM3");
  fitNames.push_back("Pull_DM2");
  fitNames.push_back("Pull_DM1");
  fitNames.push_back("Pull_DP1");
  fitNames.push_back("Pull_DP2");
  fitNames.push_back("Pull_DP3");

  // Perform gaussian fit
  for ( std::vector<TString>::const_iterator fitName = fitNames.begin();
        fitName != fitNames.end(); ++fitName )
  {
    gROOT->cd();
  
    TH1F* hFitMean = new TH1F("hFitMean_"+*fitName, "Gaussian mean of "+*fitName, nReleases, 0, nReleases);
    TH1F* hFitSigma = new TH1F("hFitSigma_"+*fitName, "Gaussian sigma of "+*fitName, nReleases, 0, nReleases);

//    hFitMean->SetMinimum(-0.5);
//    hFitMean->SetMaximum(0.5);

//    hFitSigma->SetMinimum(0);
//    hFitSigma->SetMaximum(2);

    for ( int i=0; i<nFiles; ++i )
    {
      TFile* file = files[i];
      if ( !file || !file->IsOpen() )
      {
        cout << "File " << file->GetName() << " not opened" << endl;
        continue;
      }

      TH1F* hist = dynamic_cast<TH1F*>(file->Get("DQMData/RPCRecHitsV/"+*fitName));
      if ( !hist )
      {
        cout << "Cannot open histogram " << *fitName << endl;
        continue;
      }

      TF1 gaus("gaus", "gaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
      hist->Fit("gaus", "QER0");
      TF1* fitResult = hist->GetFunction("gaus");

      hFitMean->SetBinContent(i+1, fitResult->GetParameter(1));
      hFitMean->SetBinError(i+1, fitResult->GetParError(1));
      hFitSigma->SetBinContent(i+1, fitResult->GetParameter(2));
      hFitSigma->SetBinError(i+1, fitResult->GetParError(2));

      hFitMean->GetXaxis()->SetBinLabel(i+1, releases[i]);
      hFitSigma->GetXaxis()->SetBinLabel(i+1, releases[i]);
    }

    TString canvasName = "cFit+"+*fitName;
    TCanvas* c = new TCanvas(canvasName, canvasName, cWidth, cHeight);
    c->Divide(1,2);

    c->cd(1);
    hFitMean->Draw("E");

    c->cd(2);
    hFitSigma->Draw("E");

    c->Print(TString("gif/")+c->GetName()+".gif");
  }
}

