{
  using namespace std;

  gSystem->CompileMacro("HppXSec.cpp", "k");

  // Read signal sample files
  const double hppMass = 140;

  TString sigName = Form("Hpp%.0f_EMu_10TeV_GEN_HLT", hppMass);

  TChain* sigChain = new TChain("ntp");
  sigChain->Add("res/"+sigName+".root");

  // Read background sample files
  TChain* bkgChain_TT4l = new TChain("ntp");
  TChain* bkgChain_ZZ4l = new TChain("ntp");
  TChain* bkgChain_LLBB = new TChain("ntp");

  bkgChain_TT4l->Add("res/TT_4l_10TeV_GEN.root");
  bkgChain_ZZ4l->Add("res/ZZ_4l_10TeV_GEN.root");
  bkgChain_LLBB->Add("res/LLBB_4l_10TeV_GEN.root");

  // Build MVA trainer
  gSystem->CompileMacro("LearnMVA.cpp", "k");
  //gInterpreter->LoadMacro("LearnMVA.cpp");

  TString outFileName = "mva/"+sigName+".root";
  LearnMVA learnMVA(outFileName);

  learnMVA.SetSignal(sigChain, getHppXSec(hppMass));
  learnMVA.SetBackground(bkgChain_TT4l, (280900*1.46*0.01091)/(1007062));  
  learnMVA.SetBackground(bkgChain_ZZ4l, (189*(1+0.35+0.2)*0.3165)/(898940));
  learnMVA.SetBackground(bkgChain_LLBB, (56200*1.66*0.007)/(1063204));

  learnMVA.Run();

/*
  if ( !gROOT->IsBatch() ) 
  {
    #include "TMVA/macros/TMVAGui.C"
    gROOT->SetMacroPath(TString(gROOT->GetMacroPath())+":TMVA/macros");
    TMVAGui(outFileName);
  }
*/
}

