#include <TCanvas.h>
#include <TGraphAsymmErrors.h>

#include <vector>

using namespace std;

void drawOverlayPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT);
void drawEfficiencyPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT);

void view()
{
  viewJet("HLT_Jet50U_MinimumBias", "Central");
//  viewJet("HLT_Jet50U_TPGSkim", "Central");
//  viewMuon("HLT_Mu9_VBTFSkim", "All");
//  viewMuon("HLT_Mu9_VBTFSkim", "Barrel");
//    viewMuon("HLT_Mu9_MinimumBias", "All");
//    viewMuon("HLT_Mu9_MinimumBias", "Barrel");
//    viewMuon("HLT_Mu9_MinimumBias", "Overlap");
//    viewMuon("HLT_Mu9_MinimumBias", "Endcap");
}

void viewJet(TString triggerName, TString selectionName)
{
  TCanvas* c;

  TFile* f = TFile::Open("result/"+triggerName+".root");

  // Plots for Transverse energy
  TH1F* hL1Et = (TH1F*)f->Get(selectionName+"/hEtL1T");
  TH1F* hHLTEt = (TH1F*)f->Get(selectionName+"/hEtHLT");
  TH1F* hRecoEt = (TH1F*)f->Get(selectionName+"/hEtReco");

  c = new TCanvas(triggerName+selectionName+"Et", triggerName+" "+selectionName+" Et", 600, 700); c->Divide(1,2);
  drawEfficiencyPlots(c->cd(1), hRecoEt, hL1Et, hHLTEt);
  drawOverlayPlots(c->cd(2), hRecoEt, hL1Et, hHLTEt);

  c->Print("result/"+triggerName+"/"+selectionName+"_Et.png");

  // Plots for Pseudorapidity
  TH1F* hEtaL1 = (TH1F*)f->Get(selectionName+"/hEtaL1T");
  TH1F* hEtaHLT = (TH1F*)f->Get(selectionName+"/hEtaHLT");
  TH1F* hRecoEta = (TH1F*)f->Get(selectionName+"/hEtaReco");

  c = new TCanvas(triggerName+"_Eta", triggerName+" "+selectionName+" Eta", 600, 700); c->Divide(1,2);
  drawEfficiencyPlots(c->cd(1), hRecoEta, hEtaL1, hEtaHLT);
  drawOverlayPlots(c->cd(2), hRecoEta, hEtaL1, hEtaHLT);

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
  TH1F* hRecoEta = (TH1F*)f->Get(selectionName+"/hEtaReco");

  c = new TCanvas("c"+triggerName+selectionName+"Pt", triggerName+" "+selectionName+" Pt", 600, 700); c->Divide(1,2);
  drawEfficiencyPlots(c->cd(1), hEtReco, hEtL1, hEtHLT);
  drawOverlayPlots(c->cd(2), hEtReco, hEtL1, hEtHLT);
  c->Print("result/"+triggerName+"/"+selectionName+"_Pt.png");

  c = new TCanvas("c"+triggerName+selectionName+"Eta", triggerName+" "+selectionName+" Eta", 600, 700);  c->Divide(1,2);
  drawEfficiencyPlots(c->cd(1), hRecoEta, hEtaL1, hEtaHLT);
  drawOverlayPlots(c->cd(2), hRecoEta, hEtaL1, hEtaHLT);
  c->Print("result/"+triggerName+"/"+selectionName+"_Eta.png");
}

void drawOverlayPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT)
{
  if ( !hOff || !hL1T || !hHLT ) return;
  if ( !pad ) pad = new TCanvas;

  pad->SetBorderMode(0);
  pad->SetGridx();
  pad->SetGridy();

  hOff->SetLineColor(kBlack);
  hL1T->SetLineColor(kRed);
  hHLT->SetLineColor(kBlue);

  hOff->Draw();
  hL1T->Draw("sames");
  hHLT->Draw("sames");

  gPad->Update();
  TPaveStats* stats;

  stats = (TPaveStats*)hOff->FindObject("stats");
  stats->SetLineColor(hOff->GetLineColor());
  stats->SetY1NDC(0.75); stats->SetY2NDC(0.95);

  stats = (TPaveStats*)hL1T->FindObject("stats");
  stats->SetLineColor(hL1T->GetLineColor());
  stats->SetY1NDC(0.50); stats->SetY2NDC(0.70);

  stats = (TPaveStats*)hHLT->FindObject("stats");
  stats->SetLineColor(hHLT->GetLineColor());
  stats->SetY1NDC(0.25); stats->SetY2NDC(0.45);

  TLegend* legend = new TLegend(0.50, 0.75, 0.75, 0.95);
  legend->AddEntry(hOff, "Offline object", "lp");
  legend->AddEntry(hL1T, "L1 object", "lp");
  legend->AddEntry(hHLT, "HLT object", "lp");
  legend->Draw();
}

void drawEfficiencyPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT)
{
  if ( !hOff || !hL1T || !hHLT ) return;
  if ( !pad ) pad = new TCanvas;

  pad->SetBorderMode(0);
  pad->SetGridx();
  pad->SetGridy();

  TGraphAsymmErrors* grpL1TEff =  new TGraphAsymmErrors;
  TGraphAsymmErrors* grpHLTEff =  new TGraphAsymmErrors;
  TGraphAsymmErrors* grpGlbEff =  new TGraphAsymmErrors;

  grpL1TEff->BayesDivide(hL1T, hOff);
  grpHLTEff->BayesDivide(hHLT, hL1T);
  grpGlbEff->BayesDivide(hHLT, hOff);

  grpL1TEff->SetName("grpL1TEff");
  grpHLTEff->SetName("grpHLTEff");
  grpGlbEff->SetName("grpGlbEff");

  grpL1TEff->SetTitle("Efficiency");
  grpHLTEff->SetTitle("Efficiency");
  grpGlbEff->SetTitle("Efficiency");

  grpL1TEff->GetYaxis()->SetTitle("Efficiency");
  grpHLTEff->GetYaxis()->SetTitle("Efficiency");
  grpGlbEff->GetYaxis()->SetTitle("Efficiency");

  grpL1TEff->GetXaxis()->SetTitle(hOff->GetXaxis()->GetTitle());
  grpHLTEff->GetXaxis()->SetTitle(hOff->GetXaxis()->GetTitle());
  grpGlbEff->GetXaxis()->SetTitle(hOff->GetXaxis()->GetTitle());

  const double xmin = hOff->GetXaxis()->GetXmin();
  const double xmax = hOff->GetXaxis()->GetXmax();
  grpL1TEff->GetXaxis()->SetLimits(xmin, xmax);
  grpHLTEff->GetXaxis()->SetLimits(xmin, xmax);
  grpGlbEff->GetXaxis()->SetLimits(xmin, xmax);

  grpL1TEff->SetLineColor(kRed);
  grpHLTEff->SetLineColor(kBlue);
  grpGlbEff->SetLineColor(kGreen);

  grpL1TEff->SetMinimum(0); grpL1TEff->SetMaximum(1.1);
  grpHLTEff->SetMinimum(0); grpHLTEff->SetMaximum(1.1);
  grpGlbEff->SetMinimum(0); grpGlbEff->SetMaximum(1.1);

  grpL1TEff->Draw("AL");
  grpHLTEff->Draw("L");
  grpGlbEff->Draw("L");

  TLegend* effLegend = new TLegend(0.75, 0.15, 0.95, 0.4);
  effLegend->AddEntry(grpL1TEff, "L1 efficiency", "lp");
  effLegend->AddEntry(grpHLTEff, "HLT efficiency", "lp");
  effLegend->AddEntry(grpGlbEff, "Global efficiency", "lp");
  effLegend->Draw();
}
