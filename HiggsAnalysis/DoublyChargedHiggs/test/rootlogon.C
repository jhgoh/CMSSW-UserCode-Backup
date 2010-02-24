{
  if ( TString(gSystem->Getenv("CMSSW_VERSION")).Length() > 0 )
  {
    cout << "============================" << endl;
    cout << " Detected " << gSystem->Getenv("CMSSW_VERSION") << endl;
    cout << " Setting up ROOT for FWLite " << endl;
    cout << "============================" << endl;

    gSystem->Load("libFWCoreFWLite.so");
    gSystem->Load("libDataFormatsFWLite.so");
    gSystem->Load("libDataFormatsPatCandidates.so");

    AutoLibraryLoader::enable();
  }
}
