void viewResult()
{
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
  TFile* file = new TFile("hist_beamSplash_120026.root");
  TDirectory* dir = (TDirectory*)file->Get("muonTimingAnalyzer");

  TH1F* hREP1 = (TH1F*)dir->Get("DigiBx_REP1");
  TH1F* hREP2 = (TH1F*)dir->Get("DigiBx_REP2");
  TH1F* hREP3 = (TH1F*)dir->Get("DigiBx_REP3");
  TH1F* hREN1 = (TH1F*)dir->Get("DigiBx_REN1");
  TH1F* hREN2 = (TH1F*)dir->Get("DigiBx_REN2");
  TH1F* hREN3 = (TH1F*)dir->Get("DigiBx_REN3");
  TH2F* h2BxVsNDigi = (TH2F*)dir->Get("h2BxVsNDigi");
  
  cREP1->cd(); hREP1->Draw();
  cREP2->cd(); hREP2->Draw();
  cREP3->cd(); hREP3->Draw();
  cREN1->cd(); hREN1->Draw();
  cREN2->cd(); hREN2->Draw();
  cREN3->cd(); hREN3->Draw();
  cTop->cd(); h2BxVsNDigi->Draw();
}
