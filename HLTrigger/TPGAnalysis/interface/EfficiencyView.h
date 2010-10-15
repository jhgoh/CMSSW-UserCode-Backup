#ifndef HLTrigger_TPGAnalysis_EfficiencyView_H
#define HLTrigger_TPGAnalysis_EfficiencyView_H

#include <TString.h>
#include <TH1F.h>
#include <TPad.h>
#include <TGraphAsymmErrors.h>

void calculateEfficiency(TH1F* hAccept, TH1F* hTrials, TGraphAsymmErrors* grp, TString fitFtn = "");
void drawOverlayPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT, const bool doLogY = false);
void drawEfficiencyPlots(TVirtualPad* pad, TH1F* hOff, TH1F* hL1T, TH1F* hHLT, TString fitFunction = "");

#endif

