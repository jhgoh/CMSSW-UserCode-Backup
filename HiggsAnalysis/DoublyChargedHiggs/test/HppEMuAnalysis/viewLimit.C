#include "TH1F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"

void viewLimit()
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile* fLimitResult = TFile::Open("limit.root");

  // Luminosity steps
  const int nLumiPoints = 12;
  const double lumiPoints[nLumiPoints] = {0.01, 0.02, 0.05, 0.07, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 5.0, 10.0};
  const double lumiProjection = 1;

  // Set mass points
  const int nMassPoints = 14;
  const double massPoints[nMassPoints] = {140, 160, 180, 200, 220, 240, 260, 280, 300, 350, 400, 500, 600, 800};

  gROOT->cd();

  const int nBinMass = 100;
  const int nBinLumi = 5000;

  const double minMass = massPoints[0], maxMass = massPoints[nMassPoints-1];
  const double minLumi = lumiPoints[0], maxLumi = lumiPoints[nLumiPoints-1];

  TGraph2D* grpCL_Exclusion = new TGraph2D(nMassPoints*nLumiPoints);
  grpCL_Exclusion->SetName("grpCL_Exclusion");
  grpCL_Exclusion->SetTitle("Exclusion limit CL = 1-CLs (B only hypothesis)");

  TGraph2D* grpCL_Discovery = new TGraph2D(nMassPoints*nLumiPoints);
  grpCL_Discovery->SetName("grpCL_Discovery");
  grpCL_Discovery->SetTitle("Discovery limit CL = 1-CLb (S+B hypothesis)");

  // Exclusion, lumi projection
  TGraph* grpCL_LumiProjection_Exclusion_CL = new TGraph(nMassPoints);
  TGraph* grpCL_LumiProjection_Exclusion = new TGraph(nMassPoints);
  TGraphAsymmErrors* grpCL_LumiProjection_Exclusion_1Sigma = new TGraphAsymmErrors(nMassPoints);
  TGraphAsymmErrors* grpCL_LumiProjection_Exclusion_2Sigma = new TGraphAsymmErrors(nMassPoints);

  grpCL_LumiProjection_Exclusion_CL->SetTitle("Exclusion limit");
  grpCL_LumiProjection_Exclusion_1Sigma->SetFillColor(kGreen);
  grpCL_LumiProjection_Exclusion_2Sigma->SetFillColor(kYellow);
  grpCL_LumiProjection_Exclusion->SetLineColor(kRed);

  // Discovery, lumi projection
  TGraph* grpCL_LumiProjection_Discovery_CL = new TGraph(nMassPoints);
  TGraph* grpCL_LumiProjection_Discovery = new TGraph(nMassPoints);
  TGraphAsymmErrors* grpCL_LumiProjection_Discovery_1Sigma = new TGraphAsymmErrors(nMassPoints);
  TGraphAsymmErrors* grpCL_LumiProjection_Discovery_2Sigma = new TGraphAsymmErrors(nMassPoints);

  grpCL_LumiProjection_Discovery->SetTitle("Discovery limit");
  grpCL_LumiProjection_Discovery_1Sigma->SetFillColor(kGreen);
  grpCL_LumiProjection_Discovery_2Sigma->SetFillColor(kYellow);
  grpCL_LumiProjection_Discovery->SetLineColor(kBlue);

  for ( int massIdx=0; massIdx<nMassPoints; ++massIdx )
  {
    const double mass = massPoints[massIdx];
    
    for ( int lumiIdx=0; lumiIdx<nLumiPoints; ++lumiIdx )
    {
      const double lumi = lumiPoints[lumiIdx];

      TConfidenceLevel* cl_Exclusion = fLimitResult->Get(Form("CLs/cl_Exclusion_M%.1f_L%.1f", mass, lumi));
      TConfidenceLevel* cl_Discovery = fLimitResult->Get(Form("CLb/cl_Discovery_M%.1f_L%.1f", mass, lumi));
      
      const int idx = lumiIdx + massIdx*nLumiPoints;
      grpCL_Exclusion->SetPoint(idx, mass, lumi, 1-cl_Exclusion->CLs());
      grpCL_Discovery->SetPoint(idx, mass, lumi, 1-cl_Discovery->CLb());

      if ( lumi == lumiProjection )
      {
        // Block for exclusion
        {
          const double excl_expCL = 1-cl_Exclusion->GetExpectedCLsb_b();
          const double excl_expCL_errLo1 = fabs(cl_Exclusion->GetExpectedCLsb_b(-1)-excl_expCL);
          const double excl_expCL_errLo2 = fabs(cl_Exclusion->GetExpectedCLsb_b(-2)-excl_expCL);
          const double excl_expCL_errUp1 = fabs(cl_Exclusion->GetExpectedCLsb_b( 1)-excl_expCL);
          const double excl_expCL_errUp2 = fabs(cl_Exclusion->GetExpectedCLsb_b( 2)-excl_expCL);

          grpCL_LumiProjection_Exclusion_CL->SetPoint(massIdx, mass, 1-cl_Discovery->CLs());

          grpCL_LumiProjection_Exclusion->SetPoint(massIdx, mass, excl_expCL);
          grpCL_LumiProjection_Exclusion_1Sigma->SetPoint(massIdx, mass, excl_expCL);
          grpCL_LumiProjection_Exclusion_2Sigma->SetPoint(massIdx, mass, excl_expCL);

          grpCL_LumiProjection_Exclusion_1Sigma->SetPointError(massIdx, 0, 0, excl_expCL_errLo1, excl_expCL_errUp1);
          grpCL_LumiProjection_Exclusion_2Sigma->SetPointError(massIdx, 0, 0, excl_expCL_errLo2, excl_expCL_errUp2);
        }

        // Block for discovery
        {
          const double discover_expCL = 1-cl_Discovery->GetExpectedCLb_sb();
          const double discover_expCL_errLo1 = fabs(cl_Discovery->GetExpectedCLb_b(-1)-discover_expCL);
          const double discover_expCL_errLo2 = fabs(cl_Discovery->GetExpectedCLb_b(-2)-discover_expCL);
          const double discover_expCL_errUp1 = fabs(cl_Discovery->GetExpectedCLb_b( 1)-discover_expCL);
          const double discover_expCL_errUp2 = fabs(cl_Discovery->GetExpectedCLb_b( 2)-discover_expCL);

          grpCL_LumiProjection_Discovery_CL->SetPoint(massIdx, mass, 1-cl_Discovery->CLb());

          grpCL_LumiProjection_Discovery->SetPoint(massIdx, mass, discover_expCL);
          grpCL_LumiProjection_Discovery_1Sigma->SetPoint(massIdx, mass, discover_expCL);
          grpCL_LumiProjection_Discovery_2Sigma->SetPoint(massIdx, mass, discover_expCL);

          grpCL_LumiProjection_Discovery_1Sigma->SetPointError(massIdx, 0, 0, discover_expCL_errLo1, discover_expCL_errUp1);
          grpCL_LumiProjection_Discovery_2Sigma->SetPointError(massIdx, 0, 0, discover_expCL_errLo2, discover_expCL_errUp2);

          cout << discover_expCL_errLo1 << ' ' << discover_expCL_errUp1 << endl;
        }
   
      }
    }
  }

  TH2D* hGrpCL_Exclusion = new TH2D("hGrpCL_Exclusion", "Confidence Level for Discovery CL=1-CLs;"
                                    "H^{#pm#pm} mass [GeV/c^{2}];Integrated Luminosity [fb^{-1}]", 
                                    nBinMass, minMass, maxMass, nBinLumi, minLumi, maxLumi);
  double contLevels[] = {0.90, 0.95};
  hGrpCL_Exclusion->SetContour(2, contLevels);
  grpCL_Exclusion->SetHistogram(hGrpCL_Exclusion);

  TCanvas* c1 = new TCanvas;
  c1->SetLogy();
  grpCL_Exclusion->Draw("CONT4");

/*
  TH2F* hGrpCL_Discovery = new TH2F("hGrpCL_Discovery", "Confidence Level for Exclusion CL=1-CLs;"
                                    "H^{#pm#pm} mass [GeV/c^{2}];Integrated Luminosity [fb^{-1}]", 
                                    nBinMass, minMass, maxMass, nBinLumi, minLumi, maxLumi);
  grpCL_Discovery->SetHistogram(hGrpCL_Discovery);
  
  TCanvas* c2 = new TCanvas;
  c1->SetLogy();
  grpCL_Discovery->Draw("CONT4Z");
*/

/*
  TCanvas* c3 = new TCanvas;
  c3->SetGridx();
  c3->SetGridy();
  grpCL_LumiProjection_Exclusion_CL->SetTitle("Exclusion limit;H^{#pm#pm} mass [GeV/c^{2};1-CL_{s}");
  grpCL_LumiProjection_Exclusion_CL->SetMinimum(0);
  grpCL_LumiProjection_Exclusion_CL->SetMaximum(1.1);

  grpCL_LumiProjection_Exclusion_CL->Draw("AL");
  grpCL_LumiProjection_Exclusion_2Sigma->Draw("3");
  grpCL_LumiProjection_Exclusion_1Sigma->Draw("3");
  grpCL_LumiProjection_Exclusion->Draw("LP");
  grpCL_LumiProjection_Exclusion_CL->Draw("LP");

  TCanvas* c4 = new TCanvas;
  c4->SetGridx();
  c4->SetGridy();
  grpCL_LumiProjection_Discovery_CL->SetTitle("Discovery potential;H^{#pm#pm} mass [GeV/c^{2};1-CL_{b}");
  grpCL_LumiProjection_Discovery_CL->SetMinimum(0);
  grpCL_LumiProjection_Discovery_CL->SetMaximum(1.1);

  grpCL_LumiProjection_Discovery_CL->Draw("AL");
  grpCL_LumiProjection_Discovery_2Sigma->Draw("3");
  grpCL_LumiProjection_Discovery_1Sigma->Draw("3");
  grpCL_LumiProjection_Discovery->Draw("LP");
  grpCL_LumiProjection_Discovery_CL->Draw("LP");
*/
}
