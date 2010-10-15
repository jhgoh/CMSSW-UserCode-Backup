#include "HLTrigger/TPGAnalysis/interface/EfficiencyView.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>

#include <iostream>

using namespace std;

void view();
void viewJet(TString triggerName, TString selectionName);
void viewMuon(TString triggerName, TString selectionName);

void view()
{
  gStyle->SetOptFit(11111);

  TF1* jetTurnOnCurve = new TF1("jetTurnOnCurve", "[0]/2*(1+TMath::Erf((x-[1])/sqrt(2.)/[2]))", 0, 200);
  jetTurnOnCurve->SetParameters(0.5, 100, 1);
  jetTurnOnCurve->SetParLimits(0, 0, 1);
  jetTurnOnCurve->SetParLimits(1, 0, 200);
  jetTurnOnCurve->SetParLimits(2, 0, 2);

  TString sampleName;
  
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

  sampleName = "JetAnalysis_Mu_Run2010B-PromptReco-v2";
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/AllJetNoL1");
  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/CentralJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/OverlapJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/ForwardJetNoL1");

  sampleName = "JetAnalysis_MinimumBias_Run2010A-Sep17ReReco_v2";
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/AllJetNoL1");
  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/CentralJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/OverlapJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/ForwardJetNoL1");
  
  sampleName = "JetAnalysis_MinimumBias_Run2010B-PromptReco-v2";
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/AllJetNoL1");
  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/CentralJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/OverlapJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/ForwardJetNoL1");
  

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
  drawEfficiencyPlots(c->cd(1), hRecoEt, hEtL1, hEtHLT, "jetTurnOnCurve");
  drawOverlayPlots(c->cd(2), hRecoEt, hEtL1, hEtHLT, true);

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
  drawOverlayPlots(c->cd(2), hEtReco, hEtL1, hEtHLT, true);
  c->Print("result/"+triggerName+"/"+selectionName+"_Et.png");

  c = new TCanvas("c"+triggerName+selectionName+"Eta", triggerName+" "+selectionName+" Eta", 600, 700);  c->Divide(1,2);
  drawEfficiencyPlots(c->cd(1), hEtaReco, hEtaL1, hEtaHLT);
  drawOverlayPlots(c->cd(2), hEtaReco, hEtaL1, hEtaHLT);
  c->Print("result/"+triggerName+"/"+selectionName+"_Eta.png");
}

