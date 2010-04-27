{
  if ( gSystem->Getenv("CMSSW_VERSION") != "" )
  {
    cout << "============================" << endl;
    cout << " Detected " << gSystem->Getenv("CMSSW_VERSION") << endl;
    cout << " Setting up ROOT for FWLite " << endl;
    cout << "============================" << endl;

    gSystem->Load("libFWCoreFWLite.so");
    AutoLibraryLoader::enable();

    gSystem->Load("libDataFormatsFWLite.so");
    gSystem->Load("libDataFormatsPatCandidates.so");
  }
}
