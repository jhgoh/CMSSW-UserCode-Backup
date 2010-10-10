#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>

#include "PhysicsTools/RooStatsCms/interface/ClopperPearsonBinomialInterval.h"

#include <vector>
#include <iostream>

using namespace std;

void view();
void viewJet(TString triggerName, TString selectionName);
void viewMuon(TString triggerName, TString selectionName);
void drawOverlayPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT);
void drawEfficiencyPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT);

void view()
{
/*
  viewJet("HLT_Jet50U_MinimumBias", "jetHLTAnalyzer/All");
  viewJet("HLT_Jet50U_MinimumBias", "jetHLTAnalyzer/Central");
  viewJet("HLT_Jet50U_MinimumBias", "jetHLTAnalyzer/Overlap");
  viewJet("HLT_Jet50U_MinimumBias", "jetHLTAnalyzer/Forward");

  viewJet("HLT_Jet50U_MinimumBias", "jetHLTAnalyzer/AllLeading");
  viewJet("HLT_Jet50U_MinimumBias", "jetHLTAnalyzer/CentralLeading");
  viewJet("HLT_Jet50U_MinimumBias", "jetHLTAnalyzer/OverlapLeading");
  viewJet("HLT_Jet50U_MinimumBias", "jetHLTAnalyzer/ForwardLeading");
*/

  TString sampleName = "JetAnalysis_Mu_Run2010B-PromptReco-v2";
  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/AllJetNoL1");
  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/CentralJetNoL1");
  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/OverlapJetNoL1");
  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/ForwardJetNoL1");

/*
  viewMuon("HLT_Mu9_MinimumBias", "muonHLTAnalyzer/All");
  viewMuon("HLT_Mu9_MinimumBias", "muonHLTAnalyzer/Barrel");
  viewMuon("HLT_Mu9_MinimumBias", "muonHLTAnalyzer/Overlap");
  viewMuon("HLT_Mu9_MinimumBias", "muonHLTAnalyzer/Endcap");

  viewMuon("HLT_Mu9_MinimumBias", "muonHLTAnalyzer/AllLeading");
  viewMuon("HLT_Mu9_MinimumBias", "muonHLTAnalyzer/BarrelLeading");
  viewMuon("HLT_Mu9_MinimumBias", "muonHLTAnalyzer/OverlapLeading");
  viewMuon("HLT_Mu9_MinimumBias", "muonHLTAnalyzer/EndcapLeading");
*/
}

void calculateEfficiency(TH1F* hAccept, TH1F* hTrials, TGraphAsymmErrors* grp)
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
}

void viewJet(TString triggerName, TString baseDir)
{
  TString selectionName = baseDir;
  selectionName.ReplaceAll("/", "_");

  TCanvas* c;

  TFile* f = TFile::Open("result/"+triggerName+".root");

  // Plots for Transverse energy
  TH1F* hEtL1 = 0;//(TH1F*)f->Get(baseDir+"/hEtL1T");
  TH1F* hEtHLT = (TH1F*)f->Get(baseDir+"/hEtHLT");
  TH1F* hRecoEt = (TH1F*)f->Get(baseDir+"/hEtReco");

  c = new TCanvas(triggerName+selectionName+"Et", triggerName+" "+selectionName+" Et", 600, 700); c->Divide(1,2);
  drawEfficiencyPlots(c->cd(1), hRecoEt, hEtL1, hEtHLT);
  drawOverlayPlots(c->cd(2), hRecoEt, hEtL1, hEtHLT);

  c->Print("result/"+triggerName+"/"+selectionName+"_Et.png");

  // Plots for Pseudorapidity
  TH1F* hEtaL1 = 0;//(TH1F*)f->Get(baseDir+"/hEtaL1T");
  TH1F* hEtaHLT = (TH1F*)f->Get(baseDir+"/hEtaHLT");
  TH1F* hEtaReco = (TH1F*)f->Get(baseDir+"/hEtaReco");

  c = new TCanvas(triggerName+selectionName+"_Eta", triggerName+" "+selectionName+" Eta", 600, 700); c->Divide(1,2);
  drawEfficiencyPlots(c->cd(1), hEtaReco, hEtaL1, hEtaHLT);
  drawOverlayPlots(c->cd(2), hEtaReco, hEtaL1, hEtaHLT);

  c->Print("result/"+triggerName+"/"+selectionName+"_Eta.png");
}

void viewMuon(TString triggerName, TString baseDir)
{
  TString selectionName = baseDir;
  selectionName.ReplaceAll("/", "_");

  TCanvas* c;

  TFile* f = TFile::Open("result/"+triggerName+".root");

  TH1F* hEtL1 = (TH1F*)f->Get(baseDir+"/hEtL1T");
  TH1F* hEtHLT = (TH1F*)f->Get(baseDir+"/hEtHLT");
  TH1F* hEtReco = (TH1F*)f->Get(baseDir+"/hEtReco");

  TH1F* hEtaL1 = (TH1F*)f->Get(baseDir+"/hEtaL1T");
  TH1F* hEtaHLT = (TH1F*)f->Get(baseDir+"/hEtaHLT");
  TH1F* hEtaReco = (TH1F*)f->Get(baseDir+"/hEtaReco");

  c = new TCanvas("c"+triggerName+selectionName+"Et", triggerName+" "+selectionName+" Et", 600, 700); c->Divide(1,2);
  drawEfficiencyPlots(c->cd(1), hEtReco, hEtL1, hEtHLT);
  drawOverlayPlots(c->cd(2), hEtReco, hEtL1, hEtHLT);
  c->Print("result/"+triggerName+"/"+selectionName+"_Et.png");

  c = new TCanvas("c"+triggerName+selectionName+"Eta", triggerName+" "+selectionName+" Eta", 600, 700);  c->Divide(1,2);
  drawEfficiencyPlots(c->cd(1), hEtaReco, hEtaL1, hEtaHLT);
  drawOverlayPlots(c->cd(2), hEtaReco, hEtaL1, hEtaHLT);
  c->Print("result/"+triggerName+"/"+selectionName+"_Eta.png");
}

void drawOverlayPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT)
{
  if ( !hOff && !hL1T && !hHLT )
  {
    cout << "Null" << endl;
    return;
  }

  if ( !pad ) pad = new TCanvas;

  pad->SetBorderMode(0);
  pad->SetGridx();
  pad->SetGridy();

  TString drawOpt = "";
  if ( hOff ) 
  {
    hOff->SetLineColor(kBlack);
    hOff->Draw(drawOpt);
    if ( drawOpt == "" ) drawOpt = "sames";
  }
  if ( hL1T )
  {
    hL1T->SetLineColor(kRed);
    hL1T->Draw(drawOpt);
    if ( drawOpt == "" ) drawOpt = "sames";
  }
  if ( hHLT )
  {
    hHLT->SetLineColor(kBlue);
    hHLT->Draw(drawOpt);
    if ( drawOpt == "" ) drawOpt = "sames";
  }

  gPad->Update();
  TPaveStats* stats;

  TLegend* legend = new TLegend(0.50, 0.75, 0.75, 0.95);

  double yNDC = 0.95;

  if ( hOff )
  {
    stats = (TPaveStats*)hOff->FindObject("stats");
    stats->SetLineColor(hOff->GetLineColor());
    stats->SetY1NDC(yNDC-0.2); stats->SetY2NDC(yNDC);

    yNDC -= 0.25;

    legend->AddEntry(hOff, "Offline object", "lp");
  }

  if ( hL1T )
  {
    stats = (TPaveStats*)hL1T->FindObject("stats");
    stats->SetLineColor(hL1T->GetLineColor());
    stats->SetY1NDC(yNDC-0.2); stats->SetY2NDC(yNDC);

    yNDC -= 0.25;

    legend->AddEntry(hL1T, "L1 object", "lp");
  }

  if ( hHLT )
  {
    stats = (TPaveStats*)hHLT->FindObject("stats");
    stats->SetLineColor(hHLT->GetLineColor());
    stats->SetY1NDC(yNDC-0.2); stats->SetY2NDC(yNDC);

    yNDC -= 0.25;

    legend->AddEntry(hHLT, "HLT object", "lp");
  }

  legend->Draw();
}

void drawEfficiencyPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT)
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

  if ( grpL1TEff )
  {
    calculateEfficiency(hL1T, hOff, grpL1TEff);
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

    if ( drawOpt == "AP" ) drawOpt = "P";
    effLegend->AddEntry(grpL1TEff, "L1 efficiency", "lp");
  }

  if ( grpHLTEff )
  {
    calculateEfficiency(hHLT, hL1T, grpHLTEff);
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

    if ( drawOpt == "AP" ) drawOpt = "P";
    effLegend->AddEntry(grpHLTEff, "HLT efficiency", "lp");
  }

  if ( grpGlbEff )
  {
    calculateEfficiency(hHLT, hOff, grpGlbEff);
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

    if ( drawOpt == "AP" ) drawOpt = "P";
    effLegend->AddEntry(grpGlbEff, "Global efficiency", "lp");
  }

  effLegend->Draw();
}
