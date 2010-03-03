#include "TH1F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"

void viewLimit()
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile* limitResult = TFile::Open("limit.root");

  TGraph2D* grpCL_Exclusion = dynamic_cast<TGraph2D*>(limitResult->Get("grpCL_Exclusion_MassVsLumi"));
  TGraph2D* grpCL_Discovery = dynamic_cast<TGraph2D*>(limitResult->Get("grpCL_Discovery_MassVsLumi"));
  if ( !grpCL_Exclusion || !grpCL_Discovery ) return;

  TCanvas* c1 = new TCanvas;
  c1->SetLogy();
  TCanvas* c2 = new TCanvas;
  c2->SetLogy();

  TH2F* hGrpCL_Discovery = new TH2F("hGrpCL_Discovery", "Confidence Level for Exclusion CL=1-CLs;"
                                    "Mass [GeV/c^{2}];Integrated Luminosity [fb^{-1}]", 100, 140, 800, 5000, 1e-2, 1);
  grpCL_Discovery->SetHistogram(hGrpCL_Discovery);

  TH2F* hGrpCL_Exclusion = new TH2F("hGrpCL_Exclusion", "Confidence Level for Discovery CL=1-CLb;"
                                  "Mass [GeV/c^{2}];Integrated Luminosity [fb^{-1}]", 100, 140, 800, 5000, 1e-2, 1);
  grpCL_Exclusion->SetHistogram(hGrpCL_Exclusion);

  c1->cd();
  grpCL_Exclusion->Draw("CONT4Z");

  c2->cd();
  grpCL_Discovery->Draw("CONT4Z");
  //grpCL_Exclusion->Draw("COLZ");
}
