#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <iostream>
#include <vector>

#include "HppXSec.cpp"

using namespace std;

void candCounting()
{
  TF1 resolFtn("resolFtn", "pol2", 140, 800);
  resolFtn.SetParameter("p0", 1.42607e+00);
  resolFtn.SetParameter("p1", 6.51756e-03);
  resolFtn.SetParameter("p2", 1.53268e-05);

  std::vector<double> massPoints;
  massPoints.push_back(140);
  massPoints.push_back(160);
  massPoints.push_back(180);
  massPoints.push_back(200);
  massPoints.push_back(220);
  massPoints.push_back(240);
  massPoints.push_back(260);
  massPoints.push_back(280);
  massPoints.push_back(300);
  massPoints.push_back(350);
  massPoints.push_back(400);
  massPoints.push_back(500);
  massPoints.push_back(600);
  massPoints.push_back(800);

  const double lumi = 1;

  const double bkgXSec_TT_4l = 280900*1.46*0.01091;
  const double bkgXSec_ZZ_4l = 189*(1+0.35+0.2)*0.3165;
  const double bkgXSec_LLBB_4l = 56200*1.66*0.007;

  const double bkgNEvents_TT_4l = 1007062;
  const double bkgNEvents_ZZ_4l = 898940;
  const double bkgNEvents_LLBB_4l = 1063204;

  const double sigNEvents = 10000;

  TGraph* grpSig = new TGraph(massPoints.size());
  TGraph* grpBkg = new TGraph(massPoints.size());
  TGraph* grpSignif = new TGraph(massPoints.size());

  for ( unsigned int i=0; i<massPoints.size(); ++i )
  {
    const double mass = massPoints[i];
    const double sigma = resolFtn.Eval(mass);
    const double sigXSec = getHppXSec(mass);

    const double massMin0 = mass-2*sigma;
    const double massMax0 = mass+2*sigma;

    TString sigFileName = Form("res/Hpp%.0f_EMu_10TeV_GEN_HLT.root", mass);
    TFile* sigFile = TFile::Open(sigFileName);
    TH1F* hSigMass = (TH1F*)sigFile->Get("HppMassNearest");
    hSigMass->Scale(lumi*sigXSec/sigNEvents);

    TFile* bkgFile_TT_4l = TFile::Open("res/TT_4l_10TeV_GEN.root");
    TFile* bkgFile_ZZ_4l = TFile::Open("res/ZZ_4l_10TeV_GEN.root");
    TFile* bkgFile_LLBB_4l = TFile::Open("res/LLBB_4l_10TeV_GEN.root");

    TH1F* hBkgMass_TT_4l = (TH1F*)bkgFile_TT_4l->Get("HppMassNearest");
    TH1F* hBkgMass_ZZ_4l = (TH1F*)bkgFile_ZZ_4l->Get("HppMassNearest");
    TH1F* hBkgMass_LLBB_4l = (TH1F*)bkgFile_LLBB_4l->Get("HppMassNearest");

    hBkgMass_TT_4l->Scale(lumi*bkgXSec_TT_4l/bkgNEvents_TT_4l);
    hBkgMass_ZZ_4l->Scale(lumi*bkgXSec_ZZ_4l/bkgNEvents_ZZ_4l);
    hBkgMass_LLBB_4l->Scale(lumi*bkgXSec_LLBB_4l/bkgNEvents_LLBB_4l);

    const int binMin = hSigMass->FindBin(massMin0);
    const int binMax = hSigMass->FindBin(massMax0);
 
    const double nSigCand = hSigMass->Integral(binMin, binMax);
    const double nBkgCand_TT_4l = hBkgMass_TT_4l->Integral(binMin, binMax);
    const double nBkgCand_ZZ_4l = hBkgMass_ZZ_4l->Integral(binMin, binMax);
    const double nBkgCand_LLBB_4l = hBkgMass_LLBB_4l->Integral(binMin, binMax);

    const double nBkgCand = nBkgCand_TT_4l+nBkgCand_ZZ_4l+nBkgCand_LLBB_4l;

    grpSig->SetPoint(i, mass, nSigCand);
    grpBkg->SetPoint(i, mass, nBkgCand);
    grpSignif->SetPoint(i, mass, nSigCand/sqrt(nSigCand+nBkgCand));

    sigFile->Close();
    bkgFile_TT_4l->Close();
    bkgFile_ZZ_4l->Close();
    bkgFile_LLBB_4l->Close();
  }

  grpSig->SetTitle("Number of Signal");
  grpBkg->SetTitle("Number of background");
  grpSignif->SetTitle("Signal significance S/#sqrt{S+B}");

  grpSig->SetLineColor(kBlue);
  grpBkg->SetLineColor(kRed);

  grpSig->SetMarkerStyle(2);
  grpBkg->SetMarkerStyle(2);
  grpSignif->SetMarkerStyle(2);

  TCanvas* cNCand = new TCanvas("cNCand", "cNCand");
  grpSig->Draw("ALP");
  grpBkg->Draw("LP");

  TCanvas* cSignif = new TCanvas("cSignif", "cSignif");
  grpSignif->Draw("ALP");
}
