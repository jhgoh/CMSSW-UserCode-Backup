#include "HLTrigger/TPGAnalysis/interface/EfficiencyView.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>

#include "PhysicsTools/RooStatsCms/interface/ClopperPearsonBinomialInterval.h"

#include <vector>
#include <iostream>

using namespace std;

void calculateEfficiency(TH1F* hAccept, TH1F* hTrials, TGraphAsymmErrors* grp, TString fitFtn)
{
  if ( hAccept == 0 || hTrials == 0 || grp == 0 ) return;

  const int nBinsAccept = hAccept->GetNbinsX();
  const int nBinsTrials = hTrials->GetNbinsX();
  if ( nBinsAccept != nBinsTrials )
  {
    cout << "Bins are not consistent" << endl;
    return;
  }

  ClopperPearsonBinomialInterval effCalculator;
  effCalculator.init(1.-0.682);
  for ( int i=1, nPoint = 0; i<nBinsTrials; ++i )
  {
    const double nTrials = hTrials->GetBinContent(i);
    const double nAccept = hAccept->GetBinContent(i);

    if ( nTrials == 0 || nAccept > nTrials ) continue;

    effCalculator.calculate(nAccept, nTrials);

    const double x = hTrials->GetBinCenter(i);
    const double dxLo = x-hTrials->GetXaxis()->GetBinLowEdge(i);
    const double dxUp = hTrials->GetXaxis()->GetBinUpEdge(i)-x;

    const double y = nAccept/nTrials;
    const double dyLo = y - effCalculator.lower();
    const double dyUp = effCalculator.upper() - y;

    grp->SetPoint(nPoint, x, y);
    grp->SetPointError(nPoint, dxLo, dxUp, dyLo, dyUp);

    ++nPoint;
  }

  if ( fitFtn != "" ) grp->Fit(fitFtn, "REX0");
}

void drawOverlayPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT, bool doLogY)
{
  if ( !hOff && !hL1T && !hHLT )
  {
    cout << "Null" << endl;
    return;
  }

  double minY = doLogY == 0 ? 1 : 0;
  double maxY = minY;

  if ( !pad ) pad = new TCanvas;

  pad->SetBorderMode(0);
  pad->SetGridx();
  pad->SetGridy();

  TString drawOpt = "";
  if ( hOff )
  {
    if ( maxY < hOff->GetMaximum() ) maxY = hOff->GetMaximum();
    hOff->SetLineColor(kBlack);
    hOff->Draw(drawOpt);
    if ( drawOpt == "" ) drawOpt = "sames";
  }
  if ( hL1T )
  {
    if ( maxY < hL1T->GetMaximum() ) maxY = hL1T->GetMaximum();
    hL1T->SetLineColor(kRed);
    hL1T->Draw(drawOpt);
    if ( drawOpt == "" ) drawOpt = "sames";
  }
  if ( hHLT )
  {
    if ( maxY < hHLT->GetMaximum() ) maxY = hHLT->GetMaximum();
    hHLT->SetLineColor(kBlue);
    hHLT->Draw(drawOpt);
    if ( drawOpt == "" ) drawOpt = "sames";
  }

  if ( doLogY ) pad->SetLogy();
  pad->Update();
  TPaveStats* stats;

  TLegend* legend = new TLegend(0.50, 0.75, 0.75, 0.95);

  double yNDC = 0.95;

  if ( hOff )
  {
    if ( (stats = (TPaveStats*)hOff->FindObject("stats")) != 0)
    {
      stats->SetLineColor(hOff->GetLineColor());
      stats->SetY1NDC(yNDC-0.2); stats->SetY2NDC(yNDC);

      yNDC -= 0.25;
    }

    legend->AddEntry(hOff, "Offline object", "lp");
  }

  if ( hL1T )
  {
    if ( (stats = (TPaveStats*)hL1T->FindObject("stats")) != 0 )
    {
      stats->SetLineColor(hL1T->GetLineColor());
      stats->SetY1NDC(yNDC-0.2); stats->SetY2NDC(yNDC);

      yNDC -= 0.25;
    }

    legend->AddEntry(hL1T, "L1 object", "lp");
  }

  if ( hHLT )
  {
    if ( (stats = (TPaveStats*)hHLT->FindObject("stats")) != 0 )
    {
      stats->SetLineColor(hHLT->GetLineColor());
      stats->SetY1NDC(yNDC-0.2); stats->SetY2NDC(yNDC);

      yNDC -= 0.25;
    }

    legend->AddEntry(hHLT, "HLT object", "lp");
  }

  legend->Draw();
}

void drawEfficiencyPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT, TString fitFunction)
{
  if ( !hOff && !hL1T && !hHLT )
  {
    cout << "Null" << endl;
    return;
  }

  if ( !pad ) pad = new TCanvas;
  TLegend* effLegend = new TLegend(0.75, 0.15, 0.95, 0.4);

  pad->SetBorderMode(0);
  pad->SetGridx();
  pad->SetGridy();

  TGraphAsymmErrors* grpL1TEff = 0;
  TGraphAsymmErrors* grpHLTEff = 0;
  TGraphAsymmErrors* grpGlbEff = 0;

  if (hL1T && hOff) grpL1TEff = new TGraphAsymmErrors;
  if (hHLT && hL1T) grpHLTEff = new TGraphAsymmErrors;
  if (hHLT && hOff) grpGlbEff = new TGraphAsymmErrors;

  TString drawOpt = "AP";

  TPaveStats* stats;
  double yNDCEff = 0.95;

  if ( grpL1TEff )
  {
    calculateEfficiency(hL1T, hOff, grpL1TEff, fitFunction);
    grpL1TEff->SetName("grpL1TEff");
    grpL1TEff->SetTitle("Efficiency");
    grpL1TEff->GetYaxis()->SetTitle("Efficiency");
    grpL1TEff->GetXaxis()->SetTitle(hOff->GetXaxis()->GetTitle());

    const double xmin = hOff->GetXaxis()->GetXmin();
    const double xmax = hOff->GetXaxis()->GetXmax();
    grpL1TEff->GetXaxis()->SetLimits(xmin, xmax);

    grpL1TEff->SetLineColor(kRed);
    grpL1TEff->SetMinimum(0); grpL1TEff->SetMaximum(1.1);
    grpL1TEff->Draw(drawOpt);
    if ( fitFunction != "" )
    {
      TF1* ftn = grpL1TEff->GetFunction(fitFunction);
      //ftn->SetLineColor(grpL1TEff->GetLineColor());
      ftn->Draw("same");

      if ( (stats = (TPaveStats*)grpL1TEff->FindObject("stats")) != 0)
      {
        stats->SetLineColor(grpL1TEff->GetLineColor());
        stats->SetY1NDC(yNDCEff-0.2); stats->SetY2NDC(yNDCEff);

        yNDCEff -= 0.25;
      }

    }

    if ( drawOpt == "AP" ) drawOpt = "P";
    effLegend->AddEntry(grpL1TEff, "L1 efficiency", "lp");
  }

  if ( grpHLTEff )
  {
    calculateEfficiency(hHLT, hL1T, grpHLTEff, fitFunction);
    grpHLTEff->SetName("grpHLTEff");
    grpHLTEff->SetTitle("Efficiency");
    grpHLTEff->GetYaxis()->SetTitle("Efficiency");
    grpHLTEff->GetXaxis()->SetTitle(hOff->GetXaxis()->GetTitle());

    const double xmin = hHLT->GetXaxis()->GetXmin();
    const double xmax = hHLT->GetXaxis()->GetXmax();
    grpHLTEff->GetXaxis()->SetLimits(xmin, xmax);

    grpHLTEff->SetLineColor(kBlue);
    grpHLTEff->SetMinimum(0); grpHLTEff->SetMaximum(1.1);
    grpHLTEff->Draw(drawOpt);
    if ( fitFunction != "" )
    {
      TF1* ftn = grpHLTEff->GetFunction(fitFunction);
      //ftn->SetLineColor(grpHLTEff->GetLineColor());
      ftn->Draw("same");

      if ( (stats = (TPaveStats*)grpHLTEff->FindObject("stats")) != 0)
      {
        stats->SetLineColor(grpHLTEff->GetLineColor());
        stats->SetY1NDC(yNDCEff-0.2); stats->SetY2NDC(yNDCEff);

        yNDCEff -= 0.25;
      }
    }

    if ( drawOpt == "AP" ) drawOpt = "P";
    effLegend->AddEntry(grpHLTEff, "HLT efficiency", "lp");
  }

  if ( grpGlbEff )
  {
    calculateEfficiency(hHLT, hOff, grpGlbEff, fitFunction);
    grpGlbEff->SetName("grpGlbEff");
    grpGlbEff->SetTitle("Efficiency");
    grpGlbEff->GetYaxis()->SetTitle("Efficiency");
    grpGlbEff->GetXaxis()->SetTitle(hOff->GetXaxis()->GetTitle());

    const double xmin = hHLT->GetXaxis()->GetXmin();
    const double xmax = hHLT->GetXaxis()->GetXmax();
    grpGlbEff->GetXaxis()->SetLimits(xmin, xmax);

    grpGlbEff->SetLineColor(kGreen);
    grpGlbEff->SetMinimum(0); grpGlbEff->SetMaximum(1.1);
    grpGlbEff->Draw(drawOpt);
    if ( fitFunction != "" )
    {
      TF1* ftn = grpGlbEff->GetFunction(fitFunction);
      //ftn->SetLineColor(grpGlbEff->GetLineColor());
      ftn->Draw("same");

      if ( (stats = (TPaveStats*)grpGlbEff->FindObject("stats")) != 0)
      {
        stats->SetLineColor(grpGlbEff->GetLineColor());
        stats->SetY1NDC(yNDCEff-0.2); stats->SetY2NDC(yNDCEff);

        yNDCEff -= 0.25;
      }
    }

    if ( drawOpt == "AP" ) drawOpt = "P";
    effLegend->AddEntry(grpGlbEff, "Global efficiency", "lp");
  }

  effLegend->Draw();
}

