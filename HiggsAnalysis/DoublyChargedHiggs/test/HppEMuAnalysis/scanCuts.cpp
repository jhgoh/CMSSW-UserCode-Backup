#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TString.h"

#include "HppXSec.cpp"

#include <iostream>
#include <vector>

using namespace std;

class CutProcessor
{
public:
  void setSignal(TFile* file, const double xSec, const double nEvent)
  {
    if ( !file || !file->IsOpen() ) return;

    sigFile_ = file;
    sigScale_ = xSec/nEvent;
  }
  void addBackground(TFile* file, const double xSec, const double nEvent)
  {
    if ( !file || !file->IsOpen() ) return;

    bkgFiles_.push_back(file);
    bkgScale_.push_back(xSec/nEvent);
  }
  void scanCut(const char* name, bool greaterThan=true)
  {
    if ( !sigFile_ || bkgFiles_.empty() ) return;

    // Read up histograms
    TH1F* hSig = dynamic_cast<TH1F*>(sigFile_->Get(name));
    if ( !hSig )
    {
      cout << "No signal histogram found\n";
      return;
    }
    const double nSigTotal = hSig->GetEntries()*sigScale_;

    const size_t nBkgFiles = bkgFiles_.size();
    std::vector<TH1F*> hBkgs;
    double nBkgTotal = 0;
    for ( size_t bkgIdx=0; bkgIdx<nBkgFiles; ++bkgIdx )
    {
      TH1F* hBkg = dynamic_cast<TH1F*>(bkgFiles_[bkgIdx]->Get(name));
      if ( !hBkg )
      {
        cout << "No histogram found " << endl;
        return;
      }
      hBkgs.push_back(hBkg);
      nBkgTotal += hBkg->GetEntries()*bkgScale_[bkgIdx];
    }

    if ( hBkgs.size() != nBkgFiles ) return;

    // Build up graph
    const Int_t nBinsX = hSig->GetNbinsX();
    TGraph* grpSignif = new TGraph(nBinsX);
    TGraph* grpEffSig = new TGraph(nBinsX);
    TGraph* grpEffBkg = new TGraph(nBinsX);

    // Count # of candidates
    double nSig = 0, nBkg = 0;
    int bin = greaterThan ? nBinsX : 1;

    for ( int i=0; i<nBinsX; ++i)
    {
      greaterThan ? --bin : ++bin;

      nSig += hSig->GetBinContent(bin)*sigScale_;

      for ( size_t bkgIdx=0; bkgIdx<nBkgFiles; ++bkgIdx )
      {
        TH1F* hBkg = hBkgs[bkgIdx];
        if ( !hBkg ) continue;
        nBkg += hBkg->GetBinContent(bin)*bkgScale_[bkgIdx];
      }

      const double x = hSig->GetBinLowEdge(bin);
      const double signif = nSig+nBkg == 0 ? 0 : nSig/sqrt(double(nSig+nBkg));
      grpSignif->SetPoint(i, x, signif);
      grpEffSig->SetPoint(i, x, double(nSig)/nSigTotal);
      grpEffBkg->SetPoint(i, x, double(nBkg)/nBkgTotal);
    }

    grpSignif->SetTitle(TString(hSig->GetTitle())+";"+hSig->GetXaxis()->GetTitle()+";Significance");
    grpEffSig->SetTitle(TString(hSig->GetTitle())+";"+hSig->GetXaxis()->GetTitle()+";Efficiency");
    grpEffBkg->SetTitle(TString(hSig->GetTitle())+";"+hSig->GetXaxis()->GetTitle()+";Efficiency");

    std::string canvasName = string("can_")+name;
    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str());
    c->Divide(1,2);
    grpSignif->SetLineWidth(2);
    grpEffSig->SetLineWidth(2);
    grpEffBkg->SetLineWidth(2);

    grpEffSig->SetLineColor(kRed);
    grpEffBkg->SetLineColor(kGreen);

    TVirtualPad* pad;

    pad = c->cd(1);
    pad->SetGridx(); pad->SetGridy();
    grpSignif->Draw("ALP");
    grpSignif->SetMinimum(0);

    pad = c->cd(2);
    pad->SetGridx(); pad->SetGridy();
    grpEffSig->SetMinimum(0);
    grpEffSig->SetMaximum(1.1);
    grpEffSig->Draw("ALP");
    grpEffBkg->Draw("LP");
  }
private:
  TFile* sigFile_;
  double sigScale_;
  std::vector<TFile*> bkgFiles_;
  std::vector<double> bkgScale_;
};

void scanCuts(const double mass = 140, const char* sigNameForm = "Hpp%.0f_EMu_10TeV_GEN_HLT")
{
  CutProcessor cutProcessor;

  // Set signal sample
  string sigFileName(string("res/")+Form(sigNameForm, mass)+".root");
  cutProcessor.setSignal(TFile::Open(sigFileName.c_str()), getHppXSec(mass), 10000);

  // Set background sample
  cutProcessor.addBackground(TFile::Open("res/LLBB_4l_10TeV_GEN.root"), 56200*1.66*0.007, 1063204);
  cutProcessor.addBackground(TFile::Open("res/TT_4l_10TeV_GEN.root"), 280900*1.46*0.01091, 1007062);
  cutProcessor.addBackground(TFile::Open("res/ZZ_4l_10TeV_GEN.root"), 189*(1+0.35+0.2)*0.3165, 898940);

  // Start cut scan
  cutProcessor.scanCut("ElectronPt");
  cutProcessor.scanCut("MuonPt");
  cutProcessor.scanCut("HppPt");
  cutProcessor.scanCut("HmmPt");

  cutProcessor.scanCut("ZToEEMass");
  cutProcessor.scanCut("ZToMuMuMass");

  cutProcessor.scanCut("ElectronIso", false);
  cutProcessor.scanCut("MuonIso", false);
}

