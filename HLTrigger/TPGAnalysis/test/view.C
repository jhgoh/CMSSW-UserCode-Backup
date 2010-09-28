#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>

#include <vector>

using namespace std;

void view();
void viewJet(TString triggerName, TString selectionName);
void viewMuon(TString triggerName, TString selectionName);
void drawOverlayPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT);
void drawEfficiencyPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT);

void view()
{
  viewJet("HLT_Jet50U_MinimumBias", "All");
  viewJet("HLT_Jet50U_MinimumBias", "Central");
  viewJet("HLT_Jet50U_MinimumBias", "Forward");

  viewJet("HLT_Jet50U_MinimumBias", "AllLeading");
  viewJet("HLT_Jet50U_MinimumBias", "CentralLeading");
  viewJet("HLT_Jet50U_MinimumBias", "ForwardLeading");

  viewMuon("HLT_Mu9_MinimumBias", "All");
  viewMuon("HLT_Mu9_MinimumBias", "Barrel");
  viewMuon("HLT_Mu9_MinimumBias", "Overlap");
  viewMuon("HLT_Mu9_MinimumBias", "Endcap");

  viewMuon("HLT_Mu9_MinimumBias", "AllLeading");
  viewMuon("HLT_Mu9_MinimumBias", "BarrelLeading");
  viewMuon("HLT_Mu9_MinimumBias", "OverlapLeading");
  viewMuon("HLT_Mu9_MinimumBias", "EndcapLeading");
}

void viewJet(TString triggerName, TString selectionName)
{
  TCanvas* c;

  TFile* f = TFile::Open("result/"+triggerName+".root");

  // Plots for Transverse energy
  TH1F* hEtL1 = (TH1F*)f->Get(selectionName+"/hEtL1T");
  TH1F* hEtHLT = (TH1F*)f->Get(selectionName+"/hEtHLT");
  TH1F* hRecoEt = (TH1F*)f->Get(selectionName+"/hEtReco");

  c = new TCanvas(triggerName+selectionName+"Et", triggerName+" "+selectionName+" Et", 600, 700); c->Divide(1,2);
  drawEfficiencyPlots(c->cd(1), hRecoEt, hEtL1, hEtHLT);
  drawOverlayPlots(c->cd(2), hRecoEt, hEtL1, hEtHLT);

  c->Print("result/"+triggerName+"/"+selectionName+"_Et.png");

  // Plots for Pseudorapidity
  TH1F* hEtaL1 = (TH1F*)f->Get(selectionName+"/hEtaL1T");
  TH1F* hEtaHLT = (TH1F*)f->Get(selectionName+"/hEtaHLT");
  TH1F* hEtaReco = (TH1F*)f->Get(selectionName+"/hEtaReco");

  c = new TCanvas(triggerName+selectionName+"_Eta", triggerName+" "+selectionName+" Eta", 600, 700); c->Divide(1,2);
  drawEfficiencyPlots(c->cd(1), hEtaReco, hEtaL1, hEtaHLT);
  drawOverlayPlots(c->cd(2), hEtaReco, hEtaL1, hEtaHLT);

  c->Print("result/"+triggerName+"/"+selectionName+"_Eta.png");
}

void viewMuon(TString triggerName, TString selectionName)
{
  TCanvas* c;

  TFile* f = TFile::Open("result/"+triggerName+".root");

  TH1F* hEtL1 = (TH1F*)f->Get(selectionName+"/hEtL1T");
  TH1F* hEtHLT = (TH1F*)f->Get(selectionName+"/hEtHLT");
  TH1F* hEtReco = (TH1F*)f->Get(selectionName+"/hEtReco");

  TH1F* hEtaL1 = (TH1F*)f->Get(selectionName+"/hEtaL1T");
  TH1F* hEtaHLT = (TH1F*)f->Get(selectionName+"/hEtaHLT");
  TH1F* hEtaReco = (TH1F*)f->Get(selectionName+"/hEtaReco");

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
    grpL1TEff->BayesDivide(hL1T, hOff);
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
    grpHLTEff->BayesDivide(hHLT, hL1T);
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
    grpGlbEff->BayesDivide(hHLT, hOff);
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
