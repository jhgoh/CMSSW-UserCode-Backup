void viewResult(TString fileName="hist_beamSplash_120026.root")
{
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1);

  TCanvas* c = new TCanvas("cOverview", "Overview", 800, 1200);

  // Set up canvas
  c->Divide(1,4);
  c->cd(2)->Divide(2,1);
  c->cd(3)->Divide(2,1);
  c->cd(4)->Divide(2,1);

  TPad* cTop = c->cd(1);
  TPad* cREP1 = c->cd(2)->cd(1);
  TPad* cREP2 = c->cd(3)->cd(1);
  TPad* cREP3 = c->cd(4)->cd(1);
  TPad* cREN1 = c->cd(2)->cd(2);
  TPad* cREN2 = c->cd(3)->cd(2);
  TPad* cREN3 = c->cd(4)->cd(2);

  // Retrieve from the root file
  TFile* file = new TFile(fileName);
  TDirectory* dir = (TDirectory*)file->Get("muonTimingAnalyzer");

  TH1F* hREP1 = (TH1F*)dir->Get("DigiBx_RE+1");
  TH1F* hREP2 = (TH1F*)dir->Get("DigiBx_RE+2");
  TH1F* hREP3 = (TH1F*)dir->Get("DigiBx_RE+3");
  TH1F* hREN1 = (TH1F*)dir->Get("DigiBx_RE-1");
  TH1F* hREN2 = (TH1F*)dir->Get("DigiBx_RE-2");
  TH1F* hREN3 = (TH1F*)dir->Get("DigiBx_RE-3");
  TH1F* hBxNumber = (TH1F*)dir->Get("hBxNumber");
/*
  hREP1->SetAxisRange(330, 350);
  hREP2->SetAxisRange(330, 350);
  hREP3->SetAxisRange(330, 350);
  hREN1->SetAxisRange(330, 350);
  hREN2->SetAxisRange(330, 350);
  hREN3->SetAxisRange(330, 350);
*/

  // This is correct ranges
  hREP1->SetAxisRange(2590, 2620);
  hREP2->SetAxisRange(2590, 2620);
  hREP3->SetAxisRange(2590, 2620);
  hREN1->SetAxisRange(2590, 2620);
  hREN2->SetAxisRange(2590, 2620);
  hREN3->SetAxisRange(2590, 2620);


  hREP1->Fit("gaus");
  hREP2->Fit("gaus");
  hREP3->Fit("gaus");
  hREN1->Fit("gaus");
  hREN2->Fit("gaus");
  hREN3->Fit("gaus");
  
  cREP1->cd(); hREP1->Draw();
  cREP2->cd(); hREP2->Draw();
  cREP3->cd(); hREP3->Draw();
  cREN1->cd(); hREN1->Draw();
  cREN2->cd(); hREN2->Draw();
  cREN3->cd(); hREN3->Draw();
  cTop->cd(); hBxNumber->Draw();
}
