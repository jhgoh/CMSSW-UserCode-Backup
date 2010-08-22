void run()
{
  gSystem->CompileMacro("DHCrossSections.h", "k");
  gSystem->CompileMacro("FWLiteAnalyzerBase.cc", "k");
  gSystem->CompileMacro("FWLiteAnalyzer.cc", "k");

  DHCrossSections crossSectionAt7TeV(DHCrossSections::_7TeV);
  FWLiteAnalyzerEMu emuAna("result.root", false);

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

    //emuAna.AddMCSample(datasetName, datasetPath, xsec/36, 10000/36);
    emuAna.AddMCSample(datasetName, datasetPath, xsec, 10000);
  }

  // Add background samples
  //string backgroundBase = "rfio:/castor/cern.ch/user/j/jhgoh/HiggsAnalysis/DoublyChargedHiggs/CMSSW_3_5_6_patch1/MC_Background_EMEM/";
  string backgroundBase = "~/data/HiggsAnalysis/DoublyChargedHiggs/CMSSW_3_5_6_patch1/MC_Background_EMEM/";

  emuAna.AddMCSample("MC_Background_EMEM/TTbar", backgroundBase+"TTbar", 94.3*1000, 500000);
  //emuAna.AddMCSample("TTbar_2l", backgroundBase+"TTbar_2l", 17.3*1000, 1000000);
  emuAna.AddMCSample("MC_Background_EMEM/ZZ", backgroundBase+"ZZ", 4.3*1000, 100000);

  //emuAna.AddMCSample("MC_Background_EMEM/ZeeJet_Pt0to15"   , backgroundBase+"ZeeJet_Pt0to15"   , 4.434e+03*1000, 200000);
  //emuAna.AddMCSample("MC_Background_EMEM/ZeeJet_Pt15to20"  , backgroundBase+"ZeeJet_Pt15to20"  , 1.454e+02*1000, 200000);
  //emuAna.AddMCSample("MC_Background_EMEM/ZeeJet_Pt20to30"  , backgroundBase+"ZeeJet_Pt20to30"  , 1.318e+02*1000, 150000);
  emuAna.AddMCSample("MC_Background_EMEM/ZeeJet_Pt30to50"  , backgroundBase+"ZeeJet_Pt30to50"  , 8.438e+01*1000, 150000);
  emuAna.AddMCSample("MC_Background_EMEM/ZeeJet_Pt50to80"  , backgroundBase+"ZeeJet_Pt50to80"  , 3.235e+01*1000, 100000);
  emuAna.AddMCSample("MC_Background_EMEM/ZeeJet_Pt80to120" , backgroundBase+"ZeeJet_Pt80to120" , 9.981e+00*1000, 100000);
  emuAna.AddMCSample("MC_Background_EMEM/ZeeJet_Pt120to170", backgroundBase+"ZeeJet_Pt120to170", 2.760e+00*1000, 100000);
  emuAna.AddMCSample("MC_Background_EMEM/ZeeJet_Pt170to230", backgroundBase+"ZeeJet_Pt170to230", 7.241e-01*1000, 100000);
  emuAna.AddMCSample("MC_Background_EMEM/ZeeJet_Pt230to300", backgroundBase+"ZeeJet_Pt230to300", 1.946e-01*1000, 100000);
  emuAna.AddMCSample("MC_Background_EMEM/ZeeJet_Pt300toInf", backgroundBase+"ZeeJet_Pt300toInf", 7.627e-02*1000, 100000);

  //emuAna.AddMCSample("MC_Background_EMEM/ZmumuJet_Pt0to15"   , backgroundBase+"ZmumuJet_Pt0to15"   , 4.434e+03*1000, 200000);
  //emuAna.AddMCSample("MC_Background_EMEM/ZmumuJet_Pt15to20"  , backgroundBase+"ZmumuJet_Pt15to20"  , 1.454e+02*1000, 200000);
  //emuAna.AddMCSample("MC_Background_EMEM/ZmumuJet_Pt20to30"  , backgroundBase+"ZmumuJet_Pt20to30"  , 1.318e+02*1000, 200000);
  emuAna.AddMCSample("MC_Background_EMEM/ZmumuJet_Pt30to50"  , backgroundBase+"ZmumuJet_Pt30to50"  , 8.438e+01*1000, 150000);
  emuAna.AddMCSample("MC_Background_EMEM/ZmumuJet_Pt50to80"  , backgroundBase+"ZmumuJet_Pt50to80"  , 3.235e+01*1000, 150000);
  emuAna.AddMCSample("MC_Background_EMEM/ZmumuJet_Pt80to120" , backgroundBase+"ZmumuJet_Pt80to120" , 9.981e+00*1000, 100000);
  emuAna.AddMCSample("MC_Background_EMEM/ZmumuJet_Pt120to170", backgroundBase+"ZmumuJet_Pt120to170", 2.760e+00*1000, 100000);
  //emuAna.AddMCSample("MC_Background_EMEM/ZmumuJet_Pt170to230", backgroundBase+"ZmumuJet_Pt170to230", 7.241e-01*1000, 100000);
  //emuAna.AddMCSample("MC_Background_EMEM/ZmumuJet_Pt230to300", backgroundBase+"ZmumuJet_Pt230to300", 1.946e-01*1000, 100000);
  //emuAna.AddMCSample("MC_Background_EMEM/ZmumuJet_Pt300toInf", backgroundBase+"ZmumuJet_Pt300toInf", 7.627e-02*1000, 100000);

//  emuAna.ListDataFiles();
  emuAna.ProcessEvent();
}
