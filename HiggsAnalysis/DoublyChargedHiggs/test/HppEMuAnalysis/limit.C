#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph2D.h"

#include "TMath.h"
#include "TLimit.h"
#include "TConfidenceLevel.h"
#include "TLimitDataSource.h"

#include <vector>
#include <iostream>

#include "HppXSec.cpp"

using namespace std;

void limit()
{
  // Define constants
  // MC samples's cross-sections and number of events
  const double bkgXSec_TT_4l = 280900*1.46*0.01091;
  const double bkgXSec_ZZ_4l = 189*(1+0.35+0.2)*0.3165;
  const double bkgXSec_LLBB_4l = 56200*1.66*0.007;

  const double bkgNEvents_TT_4l = 1007062;
  const double bkgNEvents_ZZ_4l = 898940;
  const double bkgNEvents_LLBB_4l = 1063204;

  const double sigNEvents = 10000;

  // Number of pseudo-experiments during the limit calculation
  const int nMC = 10000;

  // Luminosity steps
  const int nLumiPoints = 14;
  const double lumiPoints[nLumiPoints] = {0.01, 0.02, 0.05, 0.07, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 5.0, 10.0, 50.0, 100.};

  // Set mass points
  const int nMassPoints = 14;
  double massPoints[nMassPoints] = {140, 160, 180, 200, 220, 240, 260, 280, 300, 350, 400, 500, 600, 800};

  TFile* file_output = TFile::Open("limit.root", "RECREATE");
  TDirectory* dir_CLs = file_output->mkdir("CLs", "CLs");

  TGraph2D* grpCL_Bkg_MassVsLumi = new TGraph2D(nMassPoints*nLumiPoints);
  grpCL_Bkg_MassVsLumi->SetName("grpCL_Bkg_MassVsLumi");
  grpCL_Bkg_MassVsLumi->SetTitle("Exclusion limit CL = 1-CLs (B only hypothesis)");

  TGraph2D* grpCL_SigBkg_MassVsLumi = new TGraph2D(nMassPoints*nLumiPoints);
  grpCL_SigBkg_MassVsLumi->SetName("grpCL_SigBkg_MassVsLumi");
  grpCL_SigBkg_MassVsLumi->SetTitle("Exclusion limit CL = 1-CLs (S+B hypothesis)");

  // Read background samples first
  TFile* file_TT_4l = TFile::Open("res/TT_4l_10TeV_GEN.root");
  TFile* file_ZZ_4l = TFile::Open("res/ZZ_4l_10TeV_GEN.root");
  TFile* file_LLBB_4l = TFile::Open("res/LLBB_4l_10TeV_GEN.root");

  TH1F* hMass_TT_4l = (TH1F*)file_TT_4l->Get("HppMassNearest");
  TH1F* hMass_ZZ_4l = (TH1F*)file_ZZ_4l->Get("HppMassNearest");
  TH1F* hMass_LLBB_4l = (TH1F*)file_LLBB_4l->Get("HppMassNearest");

  hMass_TT_4l->Sumw2();
  hMass_ZZ_4l->Sumw2();
  hMass_LLBB_4l->Sumw2();

  for ( int massIdx=0; massIdx<nMassPoints; ++massIdx )
  {
    const double mass = massPoints[massIdx];
    const double sigXSec = getHppXSec(mass);

    // Read signal sample
    TString fileName_Sig = Form("res/Hpp%.0f_EMu_10TeV_GEN_HLT.root", mass);
    TFile* file_Sig = TFile::Open(fileName_Sig);

    TH1F* hMass_Sig = (TH1F*)file_Sig->Get("HppMassNearest");
    hMass_Sig->Sumw2();

    for ( int lumiIdx=0; lumiIdx<nLumiPoints; ++lumiIdx )
    {
      const double lumi = lumiPoints[lumiIdx];

      // Add up all of backgrounds with consideration of normailzation factor
      // Normalized signal sample
      TH1F* hMass_NormSig = (TH1F*)hMass_Sig->Clone("HppMassNearest_NormSig");
      hMass_NormSig->Scale(lumi*sigXSec/sigNEvents);

      // Sum up backgrounds with normalization
      TH1F* hMass_NormBkg = (TH1F*)hMass_Sig->Clone("HppMassNearest_AllBkg");
      hMass_NormBkg->Reset();
      hMass_NormBkg->Sumw2();
      hMass_NormBkg->Add(hMass_TT_4l, lumi*bkgXSec_TT_4l/bkgNEvents_TT_4l);
      hMass_NormBkg->Add(hMass_ZZ_4l, lumi*bkgXSec_ZZ_4l/bkgNEvents_ZZ_4l);
      hMass_NormBkg->Add(hMass_LLBB_4l, lumi*bkgXSec_LLBB_4l/bkgNEvents_LLBB_4l);

      // Make Pseudo data
      TH1F* hMass_DataSigBkg = (TH1F*)hMass_Sig->Clone("HppMassNearest_DataSigBkg");
      hMass_DataSigBkg->Add(hMass_NormBkg);

      // Now compute exclusion limits
      TLimitDataSource dataSource_SigBkg(hMass_NormSig, hMass_NormBkg, hMass_DataSigBkg);
      TLimitDataSource dataSource_Bkg(hMass_NormSig, hMass_NormBkg, hMass_NormBkg);

      dir_CLs->cd();
      TConfidenceLevel* cl_SigBkg = TLimit::ComputeLimit(&dataSource_SigBkg, nMC);
      TConfidenceLevel* cl_Bkg = TLimit::ComputeLimit(&dataSource_Bkg, nMC);

      cl_SigBkg->Write(Form("cl_SigBkg_M%.1f_L%.1f", mass, lumi));
      cl_Bkg->Write(Form("cl_Bkg_M%.1f_L%.1f", mass, lumi));

      const double cls_SigBkg = cl_SigBkg->CLs();
      const double cls_Bkg = cl_Bkg->CLs();
      //const double clb = cl_Bkg->CLb();
      //const double expCLs = cl->GetExpectedCLs_b();
      //const double expCLb = cl->GetExpectedCLb_b();

      const int idx = lumiIdx + massIdx*nLumiPoints;
      grpCL_SigBkg_MassVsLumi->SetPoint(idx, mass, lumi, 1-cls_SigBkg);
      grpCL_Bkg_MassVsLumi->SetPoint(idx, mass, lumi, 1-cls_Bkg);

      hMass_NormSig->Delete();
      hMass_NormBkg->Delete();
      hMass_DataSigBkg->Delete();
      //cl_SigBkg->Delete();
      //cl_Bkg->Delete();
    }
  }

  grpCL_Bkg_MassVsLumi->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  grpCL_Bkg_MassVsLumi->GetYaxis()->SetTitle("Luminosity [fb^{-1}]");

  grpCL_SigBkg_MassVsLumi->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  grpCL_SigBkg_MassVsLumi->GetYaxis()->SetTitle("Luminosity [fb^{-1}]");
  if ( !gROOT->IsBatch() )
  {
    gStyle->SetPalette(1);

    TCanvas* cCLs_Bkg = new TCanvas("cCLs_Bkg", "cCLs_Bkg");
    cCLs_Bkg->SetLogy();
    grpCL_Bkg_MassVsLumi->Draw("CONT4Z");

    TCanvas* cCLs_SigBkg = new TCanvas("cLCs_SigBkg", "cCLs_SigBkg");
    cCLs_SigBkg->SetLogy();
    grpCL_SigBkg_MassVsLumi->Draw("CONT4Z");

    file_output->Write();
  }
  else
  {
    file_output->Write();
    file_output->Close();
  }

}
