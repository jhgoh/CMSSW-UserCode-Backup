#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"

void massResol(bool drawFit=false)
{
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

  TCanvas* cGrp = new TCanvas("cGrp", "cGrp");
  cGrp->SetGridx();
  cGrp->SetGridy();

  TGraphErrors* grp = new TGraphErrors(massPoints.size());

  for ( unsigned int i=0; i<massPoints.size(); ++i )
  {
    const double mass = massPoints[i];

    TString fileName = Form("res/Hpp%.0f_EMu_10TeV_GEN_HLT.root", mass);
    TFile* f = TFile::Open(fileName);
    TH1F* hMass = (TH1F*)f->Get("HppMassNearest");

    TF1 gausFtn("gausFtn", "gaus", 
                mass-(hMass->GetRMS())/3, hMass->GetXaxis()->GetXmax());
    gausFtn.SetParameters(1, mass, 1);
    gausFtn.FixParameter(1, mass);

    hMass->Fit(&gausFtn, "RBQ");
    if ( drawFit )
    {
      TCanvas* c = new TCanvas;
      c->SetLogy();
      hMass->Draw();
    }
    else
    {
      f->Close();
    }
    const double sigma = gausFtn.GetParameter(2);
    const double sigmaErr = gausFtn.GetParError(2);

    grp->SetPoint(i, mass, fabs(sigma));
    grp->SetPointError(i, 0, sigmaErr);

  }

  TF1* fResol = new TF1("fResol", "pol2", 140, 800);
  fResol->SetLineColor(kRed);

  grp->Fit(fResol, "RB+");
  grp->SetMarkerStyle(1);

  cGrp->cd();
  grp->Draw("ALP");

  cGrp->Print("Resol.gif");
}
