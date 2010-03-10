void run(TString sampleName = "")
{
  gSystem->AddLinkedLibs("-L$ROOTSYS/lib -lMathCore");
  gSystem->CompileMacro("AnalyzeHppToEMu.cpp", "k");

  if ( sampleName.Length() == 0 )
  {
    AnalyzeHppToEMu("Hpp140_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp160_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp180_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp200_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp220_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp240_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp260_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp280_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp300_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp350_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp400_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp500_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp600_EMu_10TeV_GEN_HLT");
    AnalyzeHppToEMu("Hpp800_EMu_10TeV_GEN_HLT");

    AnalyzeHppToEMu("LLBB_4l_10TeV_GEN");
    AnalyzeHppToEMu("TT_4l_10TeV_GEN");
    AnalyzeHppToEMu("ZZ_4l_10TeV_GEN");
  }
  else
  {
    AnalyzeHppToEMu(sampleName);
  }
}
