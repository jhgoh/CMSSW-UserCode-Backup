void run()
{
  gSystem->CompileMacro("FWLiteAnalyzerBase.cc", "k");
  gSystem->CompileMacro("FWLiteAnalyzer.cc", "k");
  FWLiteAnalyzerEMu emuAna("result.root", false);

//  gSystem.CompileMacro("FWLiteAnalyzer.cc");
//  FWLiteAnalyzer* emuAna = new FWLiteAnalyzer();

//  TString signalBase = "rfio:/castor/cern.ch/user/j/jhgoh/HiggsAnalysis/DoublyChargedHiggs/CMSSW_3_5_6_patch1/MC_Signal_EMEM/";
  string signalBase = "~/data/HiggsAnalysis/DoublyChargedHiggs/CMSSW_3_5_6_patch1/candSelection/";
  const int nMassPoint = 10;
  const double massPoints[nMassPoint] = {
    70, 80, 90, 100, 110, 120, 130, 140, 150, 160
  };

  for ( int i=0; i<nMassPoint; ++i )
  {
    string datasetName = Form("DBLH-m%.0f", massPoints[i]);
    string datasetPath = Form("%sDBLH-m%.0f", signalBase, massPoints[i]);
    const double xsec = 1;

    emuAna.AddSignal(datasetName, datasetPath, xsec);
  }

  //emuAna.AddBackground("ZZ", false);

  emuAna.ListDataFiles();
  emuAna.ProcessEvent();
}
