{
  #include <vector>
  #include <string>
  #include <set>

  using namespace std;

  gSystem->CompileMacro("TCutScanner.cpp", "k");

  vector<string> signal;
  signal.push_back("res/HppMuMu.root");

  vector<string> bkg_TTBar;
  bkg_TTBar.push_back("res/TTBar.root");

  vector<string> bkg_QCD;
  bkg_QCD.push_back("res/QCD.root");

  TCutScanner cs;

  cs.SetSignal(signal, 1);
//  cs.SetBkg(bkg_TTBar, 20);
//  cs.SetBkg(bkg_QCD, 30);

  cs.Run();
}
