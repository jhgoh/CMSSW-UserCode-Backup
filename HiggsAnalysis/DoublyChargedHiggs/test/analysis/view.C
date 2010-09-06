#include <THStack.h>
#include <TH1F.h>
#include <TFile.h>
#include <TList.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TString.h>

#include <iostream>

using namespace std;

TList* makePlots(TDirectory* baseDir, int colorTable[], bool doFillColor = false);
void addPlot(TObject* obj, TString title, int color, THStack* hStack, bool doFillColor);

void view()
{
  TFile* f = TFile::Open("result.root");

  int signalColorTable[20], backgroundColorTable[20];

  for ( int i=0; i<20; ++i )
  {
    signalColorTable[i] = kAzure+10-i;
    backgroundColorTable[i] = kOrange+10-i;
  }

  TList* signalPlots = makePlots((TDirectory*)f->Get("MC_Signal_EMEM"), signalColorTable);
  TList* backgroundPlots = makePlots((TDirectory*)f->Get("MC_Background_EMEM"), backgroundColorTable, true);

  if ( signalPlots == 0 || backgroundPlots == 0 ) return;

  const int nPlots = signalPlots->GetSize();
  for ( int i=0; i<nPlots; ++i )
  {
    THStack* hSignal = (THStack*)signalPlots->At(i);
    THStack* hBackground = (THStack*)backgroundPlots->At(i);

    TString histName = hSignal->GetName();
    bool doLog = histName.Contains("Pt");// || histName.Contains("RelIso");

    TCanvas* c = new TCanvas(TString("c")+hSignal->GetName(), hSignal->GetTitle(), 1200, 600);
    TPad* pad;

    c->Divide(2,1);

    TString xTitle, yTitle;

    pad = (TPad*)c->cd(1);
    if ( doLog ) pad->SetLogy();
    pad->SetBorderSize(0);
    pad->SetBorderMode(0);
    hBackground->Draw();

    xTitle = ((TH1*)hBackground->GetHists()->At(0))->GetXaxis()->GetTitle();
    yTitle = ((TH1*)hBackground->GetHists()->At(0))->GetYaxis()->GetTitle();
    hBackground->GetXaxis()->SetTitle(xTitle);
    hBackground->GetYaxis()->SetTitle(yTitle);

    pad->BuildLegend(0.6, 0.6, 0.98, 0.98);

    pad = (TPad*)c->cd(2);
    if ( doLog ) pad->SetLogy();
    pad->SetBorderSize(0);
    pad->SetBorderMode(0);
    hSignal->Draw("nostack");

    xTitle = ((TH1*)hSignal->GetHists()->At(0))->GetXaxis()->GetTitle();
    yTitle = ((TH1*)hSignal->GetHists()->At(0))->GetYaxis()->GetTitle();
    hSignal->GetXaxis()->SetTitle(xTitle);
    hSignal->GetYaxis()->SetTitle(yTitle);

    pad->BuildLegend(0.6, 0.7, 0.98, 0.98);

    c->Print(TString(c->GetName())+".png");
  }
}

TList* makePlots(TDirectory* baseDir, int colorTable[], bool doFillColor)
{
  if ( !baseDir ) return 0;
  TList* dirList = baseDir->GetListOfKeys();

  if ( !dirList ) return 0;
  const int nKeys = dirList->GetSize();

  THStack* hStackMuonPt = new THStack("hMuonPt", "Muon p_{T};Transverse momentum [GeV/c]");
  THStack* hStackMuonEta = new THStack("hMuonEta", "Muon #eta;Pseudorapidity");
  THStack* hStackMuonPhi = new THStack("hMuonPhi", "Muon #phi;Azimuthal angle");
  THStack* hStackMuonRelIso = new THStack("hMuonRelIso", "Muon relative isolation;Relative isolation");

  THStack* hStackMuon1Pt = new THStack("hMuon1Pt", "Leading Muon p_{T};Transverse momentum [GeV/c]");
  THStack* hStackMuon1Eta = new THStack("hMuon1Eta", "Leading Muon #eta;Pseudorapidity");
  THStack* hStackMuon1Phi = new THStack("hMuon1Phi", "Leading Muon #phi;Azimuthal angle");
  THStack* hStackMuon1RelIso = new THStack("hMuon1RelIso", "Leading Muon relative isolation;Relative isolation");

  THStack* hStackElectronPt = new THStack("hElectronPt", "Electron p_{T};Transverse momentum [GeV/c]");
  THStack* hStackElectronEta = new THStack("hElectronEta", "Electron #eta;Pseudorapidity");
  THStack* hStackElectronPhi = new THStack("hElectronPhi", "Electron #phi;Azimuthal angle");
  THStack* hStackElectronRelIso = new THStack("hElectronRelIso", "Electron relative isolation;Relative isolation");

  THStack* hStackElectron1Pt = new THStack("hElectron1Pt", "Leading Electron p_{T};Transverse momentum [GeV/c]");
  THStack* hStackElectron1Eta = new THStack("hElectron1Eta", "Leading Electron #eta;Pseudorapidity");
  THStack* hStackElectron1Phi = new THStack("hElectron1Phi", "Leading Electron #phi;Azimuthal angle");
  THStack* hStackElectron1RelIso = new THStack("hElectron1RelIso", "Leading Electron relative isolation;Relative isolation");

  THStack* hStackPosEMuMass = new THStack("hPosEMuMass", "H^{++} #rightarrow e#mu mass;Mass [GeV/c^{2}]");
  THStack* hStackPosEMuPt = new THStack("hPosEMuPt", "H^{++} #rightarrow e#mu p_{T};Transverse momentum [GeV/c]");
  THStack* hStackPosEMuEta = new THStack("hPosEMuEta", "H^{++} #rightarrow e#mu #eta;Pseudorapidity");
  THStack* hStackPosEMuPhi = new THStack("hPosEMuPhi", "H^{++} #rightarrow e#mu #phi;Azimuthal angle");

  THStack* hStackNegEMuMass = new THStack("hNegEMuMass", "H^{--} #rightarrow e#mu mass;Mass [GeV/c^{2}]");
  THStack* hStackNegEMuPt = new THStack("hNegEMuPt", "H^{--} #rightarrow e#mu p_{T};Transverse momentum [GeV/c]");
  THStack* hStackNegEMuEta = new THStack("hNegEMuEta", "H^{--} #rightarrow e#mu #eta;Pseudorapidity");
  THStack* hStackNegEMuPhi = new THStack("hNegEMuPhi", "H^{--} #rightarrow e#mu #phi;Azimuthal angle");

  THStack* hStackPosEMuZVetoMass = new THStack("hPosEMuZVetoMass", "H^{++} #rightarrow e#mu mass after Z-veto cut;Mass [GeV/c^{2}]");
  THStack* hStackNegEMuZVetoMass = new THStack("hNegEMuZVetoMass", "H^{--} #rightarrow e#mu mass after Z-veto cut;Mass [GeV/c^{2}]");

  THStack* hStackPosEMuBestMass = new THStack("hPosEMuBestMass", "Best matching H^{++} #rightarrow e#mu mass;Mass [GeV/c^{2}]");
  THStack* hStackNegEMuBestMass = new THStack("hNegEMuBestMass", "Best matching H^{--} #rightarrow e#mu mass;Mass [GeV/c^{2}]");

  THStack* hStackZMuMuMass = new THStack("hZMuMuMass", "Opposite signed dimuon mass;Mass [GeV/c^{2}]");
  THStack* hStackZEEMass = new THStack("hZEEMass", "Opposite signed dielectron mass;Mass [GeV/c^{2}]");
  THStack* hStackFourLeptonSumPt = new THStack("hFourLeptonSumPt", "Scalar sum of four lepton p_{T};#Sigma p_{T} [GeV/c]");
  THStack* hStackDHDeltaPhi = new THStack("hDHDeltaPhi", "Azimuthal angle difference between two H^{#pm#pm} candidates;#Delta#phi [Radian]");

  for ( int i=0; i<nKeys; ++i )
  {
    const int color = colorTable[i];

    TString keyName = dirList->At(i)->GetName();
    TObject* obj = baseDir->Get(keyName);

    if ( !obj || !obj->IsA()->InheritsFrom("TDirectory") ) continue;
    TDirectory* dir = (TDirectory*)obj;

    // Muon histograms
    addPlot(dir->Get("m/hPt"), keyName, color, hStackMuonPt, doFillColor);
    addPlot(dir->Get("m/hEta"), keyName, color, hStackMuonEta, doFillColor);
    addPlot(dir->Get("m/hPhi"), keyName, color, hStackMuonPhi, doFillColor);
    addPlot(dir->Get("m/hRelIso"), keyName, color, hStackMuonRelIso, doFillColor);

    addPlot(dir->Get("m/m1/hPt"), keyName, color, hStackMuon1Pt, doFillColor);
    addPlot(dir->Get("m/m1/hEta"), keyName, color, hStackMuon1Eta, doFillColor);
    addPlot(dir->Get("m/m1/hPhi"), keyName, color, hStackMuon1Phi, doFillColor);
    addPlot(dir->Get("m/m1/hRelIso"), keyName, color, hStackMuon1RelIso, doFillColor);

    // Electron histograms
    addPlot(dir->Get("e/hPt"), keyName, color, hStackElectronPt, doFillColor);
    addPlot(dir->Get("e/hEta"), keyName, color, hStackElectronEta, doFillColor);
    addPlot(dir->Get("e/hPhi"), keyName, color, hStackElectronPhi, doFillColor);
    addPlot(dir->Get("e/hRelIso"), keyName, color, hStackElectronRelIso, doFillColor);

    addPlot(dir->Get("e/e1/hPt"), keyName, color, hStackElectron1Pt, doFillColor);
    addPlot(dir->Get("e/e1/hEta"), keyName, color, hStackElectron1Eta, doFillColor);
    addPlot(dir->Get("e/e1/hPhi"), keyName, color, hStackElectron1Phi, doFillColor);
    addPlot(dir->Get("e/e1/hRelIso"), keyName, color, hStackElectron1RelIso, doFillColor);

    // H++ candidates
    addPlot(dir->Get("posEMuCand/hMass"), keyName, color, hStackPosEMuMass, doFillColor);
    addPlot(dir->Get("posEMuCand/hPt"), keyName, color, hStackPosEMuPt, doFillColor);
    addPlot(dir->Get("posEMuCand/hEta"), keyName, color, hStackPosEMuEta, doFillColor);
    addPlot(dir->Get("posEMuCand/hPhi"), keyName, color, hStackPosEMuPhi, doFillColor);

    // H-- candidates
    addPlot(dir->Get("negEMuCand/hMass"), keyName, color, hStackNegEMuMass, doFillColor);
    addPlot(dir->Get("negEMuCand/hPt"), keyName, color, hStackNegEMuPt, doFillColor);
    addPlot(dir->Get("negEMuCand/hEta"), keyName, color, hStackNegEMuEta, doFillColor);
    addPlot(dir->Get("negEMuCand/hPhi"), keyName, color, hStackNegEMuPhi, doFillColor);

    // Four-lepton cut variables
    addPlot(dir->Get("fourLeptons/hZMuMuMass"), keyName, color, hStackZMuMuMass, doFillColor);
    addPlot(dir->Get("fourLeptons/hZEEMass"), keyName, color, hStackZEEMass, doFillColor);
    addPlot(dir->Get("fourLeptons/hFourLeptonSumPt"), keyName, color, hStackFourLeptonSumPt, doFillColor);
    addPlot(dir->Get("fourLeptons/hDHDeltaPhi"), keyName, color, hStackDHDeltaPhi, doFillColor);

    // DH z-veto and best matching
    addPlot(dir->Get("posEMuCand/posEMuZVeto/hMass"), keyName, color, hStackPosEMuZVetoMass, doFillColor);
    addPlot(dir->Get("negEMuCand/negEMuZVeto/hMass"), keyName, color, hStackNegEMuZVetoMass, doFillColor);
    addPlot(dir->Get("posEMuCand/posEMuBest/hMass"), keyName, color, hStackPosEMuBestMass, doFillColor);
    addPlot(dir->Get("negEMuCand/negEMuBest/hMass"), keyName, color, hStackNegEMuBestMass, doFillColor);
  }

  TList* plotList = new TList;

  plotList->Add(hStackMuonPt);
  plotList->Add(hStackMuonEta);
  plotList->Add(hStackMuonPhi);
  plotList->Add(hStackMuonRelIso);

  plotList->Add(hStackMuon1Pt);
  plotList->Add(hStackMuon1Eta);
  plotList->Add(hStackMuon1Phi);
  plotList->Add(hStackMuon1RelIso);

  plotList->Add(hStackElectronPt);
  plotList->Add(hStackElectronEta);
  plotList->Add(hStackElectronPhi);
  plotList->Add(hStackElectronRelIso);

  plotList->Add(hStackElectron1Pt);
  plotList->Add(hStackElectron1Eta);
  plotList->Add(hStackElectron1Phi);
  plotList->Add(hStackElectron1RelIso);

  plotList->Add(hStackPosEMuMass);
  plotList->Add(hStackPosEMuPt);
  plotList->Add(hStackPosEMuEta);
  plotList->Add(hStackPosEMuPhi);

  plotList->Add(hStackNegEMuMass);
  plotList->Add(hStackNegEMuPt);
  plotList->Add(hStackNegEMuEta);
  plotList->Add(hStackNegEMuPhi);

  plotList->Add(hStackPosEMuZVetoMass);
  plotList->Add(hStackNegEMuZVetoMass);
                       
  plotList->Add(hStackPosEMuBestMass);
  plotList->Add(hStackNegEMuBestMass);
                       
  plotList->Add(hStackZMuMuMass);
  plotList->Add(hStackZEEMass);
  plotList->Add(hStackFourLeptonSumPt);
  plotList->Add(hStackDHDeltaPhi);

  return plotList;
}

void addPlot(TObject* obj, TString title, int color, THStack* hStack, bool doFillColor)
{
  TH1F* h = (TH1F*)obj;
  h->SetTitle(title);
  if ( doFillColor )
  {
    h->SetFillColor(color);
    //h->SetLineWidth(2);
  }
  else
  {
    h->SetLineColor(color);
    h->SetLineWidth(2);
  }
  hStack->Add(h);
}
