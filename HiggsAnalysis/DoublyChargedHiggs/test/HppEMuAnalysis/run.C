{
  gSystem->AddLinkedLibs("-L$ROOTSYS/lib -lMathCore");
  gSystem->CompileMacro("AnalyzeHppToEMu.cpp", "k");

/*
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M1000_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M140_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M160_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M180_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M200_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M220_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M240_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M260_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M280_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M300_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M350_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M400_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M500_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M600_10TeV_FastSim");
  AnalyzeHppToEMu("PYTHIA6_HppEMu_M800_10TeV_FastSim");
*/

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
