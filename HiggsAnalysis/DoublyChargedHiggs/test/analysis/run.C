void run()
{
  gSystem->CompileMacro("DHCrossSections.h", "k");
  gSystem->CompileMacro("FWLiteAnalyzerBase.cc", "k");
  gSystem->CompileMacro("FWLiteAnalyzer.cc", "k");

  DHCrossSections crossSectionAt7TeV(DHCrossSections::_7TeV);
  FWLiteAnalyzerEMu emuAna("result.root", true);

  // Add signal samples
  //string signalBase = "rfio:/castor/cern.ch/user/j/jhgoh/HiggsAnalysis/DoublyChargedHiggs/CMSSW_3_5_6_patch1/MC_Signal_EMEM/";
  string signalBase = "~/data/HiggsAnalysis/DoublyChargedHiggs/CMSSW_3_5_6_patch1/MC_Signal_EMEM/";
  const int nMassPoint = 10;
  const double massPoints[nMassPoint] = {
    70, 80, 90, 100, 110, 120, 130, 140, 150, 160
  };

  for ( int i=0; i<nMassPoint; ++i )
  {
    string datasetName = Form("MC_Signal_EMEM/DBLH-m%.0f", massPoints[i]);
    string datasetPath = Form("%sDBLH-m%.0f", signalBase, massPoints[i]);
    const double xsec = crossSectionAt7TeV.Eval(massPoints[i]);

    emuAna.AddMCSample(datasetName, datasetPath, xsec/36, 10000);
  }

  // Add background samples
  //string backgroundBase = "rfio:/castor/cern.ch/user/j/jhgoh/HiggsAnalysis/DoublyChargedHiggs/CMSSW_3_5_6_patch1/MC_Background_EMEM/";
  string backgroundBase = "~/data/HiggsAnalysis/DoublyChargedHiggs/CMSSW_3_5_6_patch1/MC_Background_EMEM/";

  emuAna.AddMCSample("MC_Background_EMEM/TTbar", backgroundBase+"TTbar", 94.3*1000, 500000);
  //emuAna.AddMCSample("TTbar_2l", backgroundBase+"TTbar_2l", 17.3*1000, 1000000);
  emuAna.AddMCSample("MC_Background_EMEM/ZZ", backgroundBase+"ZZ", 4.3*1000, 100000);

  emuAna.ListDataFiles();
  emuAna.ProcessEvent();
}
