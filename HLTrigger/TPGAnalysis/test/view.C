#include "HLTrigger/TPGAnalysis/interface/EfficiencyView.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TRegexp.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>

#include <iostream>

using namespace std;

void view(int runNumber);
void viewJet(TString triggerName, TString selectionName);
void viewMuon(TString triggerName, TString selectionName);
void viewFV(TString fileName, const int runNumber, TString numHLTPath, TString denHLTPath);
int runNumberFromFileName(TString fileName);

void view(int runNumber)
{
  gStyle->SetOptFit(11111);

  TF1* turnOnCurve = new TF1("turnOnCurve", "[0]/2*(1+TMath::Erf((x-[1])/sqrt(2.)/[2]))", 0, 200);
  turnOnCurve->SetLineWidth(1);
  turnOnCurve->SetParameters(0.5, 100, 1);
  turnOnCurve->SetParLimits(0, 0, 1);
  turnOnCurve->SetParLimits(1, 0, 200);
  turnOnCurve->SetParLimits(2, 0, 200);

  TString dataDir = "/home/jhgoh/data/TPGAnalysis/FourVector/";
  TString dqmFile = "";

  dqmFile = dataDir+Form("DQM_V0001_R%09d__Offline__MinimumBias-Run2010B__DQM.root", runNumber);
  //viewFV(dqmFile, runNumber, "HLT_Jet15U", "MinBias");
  viewFV(dqmFile, runNumber, "HLT_Jet30U", "MinBias");
  viewFV(dqmFile, runNumber, "HLT_Jet70U_v2", "MinBias");
  viewFV(dqmFile, runNumber, "HLT_Jet100U_v2", "MinBias");
  viewFV(dqmFile, runNumber, "HLT_Jet140U_v1", "MinBias");

  dqmFile = dataDir+Form("DQM_V0001_R%09d__Offline__Mu-Run2010B__DQM.root", runNumber);
  //viewFV(dqmFile, runNumber, "HLT_Jet15U", "HLT_Mu");
  viewFV(dqmFile, runNumber, "HLT_Jet30U", "HLT_Mu");
  viewFV(dqmFile, runNumber, "HLT_Jet70U_v2", "HLT_Mu");
  viewFV(dqmFile, runNumber, "HLT_Jet100U_v2", "HLT_Mu");
  viewFV(dqmFile, runNumber, "HLT_Jet140U_v1", "HLT_Mu");

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

/*
  TString sampleName;

  sampleName = "JetAnalysis_Mu_Run2010B-PromptReco-v2";
  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/AllJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/CentralJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/OverlapJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/ForwardJetNoL1");

//  sampleName = "JetAnalysis_MinimumBias_Run2010A-Sep17ReReco_v2";
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/AllJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/CentralJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/OverlapJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/ForwardJetNoL1");
  
//  sampleName = "JetAnalysis_MinimumBias_Run2010B-PromptReco-v2";
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/AllJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/CentralJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/OverlapJetNoL1");
//  viewJet(sampleName, "jetHLTAnalyzerHLTJet70U/ForwardJetNoL1");
*/
  

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
  drawEfficiencyPlots(c->cd(1), hRecoEt, hEtL1, hEtHLT, "turnOnCurve");
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

void viewFV(TString fileName, const int runNumber, TString numHLTPath, TString denHLTPath)
{
  TFile* f = TFile::Open(fileName);
  if ( !f || f->IsZombie() )
  {
    cout << "Cannot open file " << fileName << endl;
    return;
  }

  TString histDir = Form("DQMData/Run %d/HLT/Run summary/FourVector/paths/", runNumber);
  histDir += numHLTPath;
  TH1F* hL1TEt = (TH1F*)f->Get(Form(histDir+"/"+numHLTPath+"_wrt_"+denHLTPath+"_offEtL1Off"));
  TH1F* hHLTEt = (TH1F*)f->Get(Form(histDir+"/"+numHLTPath+"_wrt_"+denHLTPath+"_offEtOnOff"));
  TH1F* hOffEt = (TH1F*)f->Get(Form(histDir+"/"+numHLTPath+"_wrt_"+denHLTPath+"_offEtOff"));

  if ( !hL1TEt || !hHLTEt || !hOffEt )
  {
    cout << "Cannot find histogram\n";
    return;
  }

  TCanvas* c = new TCanvas(TString("c")+Form("_Run_%d_", runNumber)+numHLTPath+"_"+denHLTPath+"Et", 
                           Form("Run %d", runNumber)+numHLTPath+" "+denHLTPath+" Et", 600, 700);
  c->Divide(1,2);

  drawEfficiencyPlots(c->cd(1), hOffEt, hL1TEt, hHLTEt);
  drawOverlayPlots(c->cd(2), hOffEt, hL1TEt, hHLTEt, true);

  c->Print(TString(c->GetName())+".png");
}

int runNumberFromFileName(TString dqmFileName)
{
  TString dqmSubStr = dqmFileName(TRegexp("DQM_V[0-9]+_R[0-9]+__"));

  if ( dqmSubStr == "" )
  {
    cout << "Filename is not compatible with DQM format\n";
    return -1;
  }

  TString runStr = dqmSubStr(TRegexp("R[0-9]+"));
  runStr.Remove(0,1);

  int runNumber = atoi((const char*)runStr);
  return runNumber;
}
