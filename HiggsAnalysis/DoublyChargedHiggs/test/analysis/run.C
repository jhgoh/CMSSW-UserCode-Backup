void run()
{
  gSystem->CompileMacro("PatAnalyzerBase.cc", "k");
  gSystem->CompileMacro("PatAnalyzer.cc", "k");
  PatAnalyzerEMu emuAna("result.root", false);

//  gSystem.CompileMacro("FWLiteAnalyzer.cc");
//  FWLiteAnalyzer* emuAna = new FWLiteAnalyzer();

//  TString signalBase = "rfio:/castor/cern.ch/user/j/jhgoh/HiggsAnalysis/DoublyChargedHiggs/CMSSW_3_5_6_patch1/MC_Signal_EMEM/";
  TString signalBase = "~/data/HiggsAnalysis/DoublyChargedHiggs/CMSSW_3_5_6_patch1/candSelection/";
  emuAna.AddSignal("DBLH-m80", string(signalBase+"DBLH-m80"));
  emuAna.AddSignal("DBLH-m70", string(signalBase+"DBLH-m70"));
  emuAna.AddSignal("DBLH-m90", string(signalBase+"DBLH-m90"));
  emuAna.AddSignal("DBLH-m100", string(signalBase+"DBLH-m100"));
  emuAna.AddSignal("DBLH-m110", string(signalBase+"DBLH-m110"));
  emuAna.AddSignal("DBLH-m120", string(signalBase+"DBLH-m120"));
  emuAna.AddSignal("DBLH-m130", string(signalBase+"DBLH-m130"));
  emuAna.AddSignal("DBLH-m140", string(signalBase+"DBLH-m140"));
  emuAna.AddSignal("DBLH-m150", string(signalBase+"DBLH-m150"));
  emuAna.AddSignal("DBLH-m160", string(signalBase+"DBLH-m160"));

  //emuAna.AddBackground("ZZ", false);

  emuAna.ListDataFiles();
  emuAna.ProcessEvent();
}
