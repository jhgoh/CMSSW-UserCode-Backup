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
    sigFile_ = file;
    sigScale_ = xSec/nEvent;
  }
  void addBackground(TFile* file, const double xSec, const double nEvent)
  {
    bkgFiles_.push_back(file);
    bkgScale_.push_back(xSec/nEvent);
  }
  void scanCut(TString name, bool greaterThan=true)
  {
    // Read up histograms
    TH1F* hSig = (TH1F*)(sigFile_->Get(name));
    const int nSigTotal = hSig->GetEntries()*sigScale_;

    std::vector<TH1F*> hBkgs;
    int nBkgTotal = 0;
    for ( size_t bkgIdx=0; bkgIdx<bkgFiles_.size(); ++bkgIdx )
    {
      hBkgs.push_back((TH1F*)(bkgFiles_[bkgIdx]->Get(name)));
      nBkgTotal += hBkgs[bkgIdx]->GetEntries()*bkgScale_[bkgIdx];
    }

    // Build up graph
    TGraph* grpSignif = new TGraph(hSig->GetNbinsX());
    TGraph* grpEffSig = new TGraph(hSig->GetNbinsX());
    TGraph* grpEffBkg = new TGraph(hSig->GetNbinsX());

    // Count # of candidates
    int nSig = 0, nBkg = 0;
    int bin = greaterThan ? hSig->GetNbinsX() : 1;
    for ( int i=0; i<hSig->GetNbinsX(); ++i)
    {
      greaterThan ? --bin : ++bin;

      nSig += hSig->GetBinContent(bin)*sigScale_;
      for ( size_t bkgIdx=0; bkgIdx<bkgFiles_.size(); ++bkgIdx )
      {
        nBkg += hBkgs[bkgIdx]->GetBinContent(bin)*bkgScale_[bkgIdx];
      }

      const double x = hSig->GetBinLowEdge(bin);
      const double signif = nSig+nBkg == 0 ? 0 : nSig/sqrt(nSig+nBkg);
      grpSignif->SetPoint(i, x, signif);
      grpEffSig->SetPoint(i, x, double(nSig)/nSigTotal);
      grpEffBkg->SetPoint(i, x, double(nBkg)/nBkgTotal);

      cout << x << ' ' << signif << ' ' << double(nSig)/nSigTotal << ' ' << double(nBkg)/nBkgTotal << endl;
    }

    TCanvas* c = new TCanvas("can_"+name, "can_"+name);
    c->Divide(1,2);
    grpSignif->SetLineWidth(2);
    grpEffSig->SetLineWidth(2);
    grpEffBkg->SetLineWidth(2);

    grpEffSig->SetLineColor(kRed);
    grpEffBkg->SetLineColor(kGreen);

    c->cd(1);
    grpSignif->Draw("ALP");
    grpSignif->SetMinimum(0);

    c->cd(2);
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

void scanCuts(const double mass = 140, const TString sigNameForm = "Hpp%.0f_EMu_10TeV_GEN_HLT")
{
  CutProcessor cutProcessor;

  // Set signal sample
  const TString sigFileName = TString("res/")+Form(sigNameForm, mass)+".root";
  cutProcessor.setSignal(TFile::Open(sigFileName), getHppXSec(mass), 10000);
 
  // Set background sample
  cutProcessor.addBackground(TFile::Open("res/LLBB_4l_10TeV_GEN.root"), 56200*1.66*0.007, 1063204);
  cutProcessor.addBackground(TFile::Open("res/TT_4l_10TeV_GEN.root"), 280900*1.46*0.01091, 1007062);
  cutProcessor.addBackground(TFile::Open("res/ZZ_4l_10TeV_GEN.root"), 189*(1+0.35+0.2)*0.3165, 898940);

  // Start cut scan
  cutProcessor.scanCut("ElectronPt");
  cutProcessor.scanCut("MuonPt");
  cutProcessor.scanCut("HppPt");
  cutProcessor.scanCut("HmmPt");

  cutProcessor.scanCut("ElectronIso", false);
  cutProcessor.scanCut("MuonIso", false);
  //TString varList[nCuts] = {"MassDiff", "ZToMuMuMass", "ZToEEMass"};
  
}

