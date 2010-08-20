#ifndef HiggsAnalysis_DoublyChargedHiggs_FWLiteAnalyzerBase_H
#define HiggsAnalysis_DoublyChargedHiggs_FWLiteAnalyzerBase_H

#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
//#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//#include "DataFormats/Candidate/interface/OverlapChecker.h"

#endif

#include "TFile.h"
#include "TH1F.h"
#include "TString.h"

#include <string>
#include <vector>
#include <map>

namespace pat
{
  class Muon;
  class Electron;
  class CompositeCandidate;
}

struct HMuon
{
  HMuon(TDirectory* baseDir, const double scale = 1, TString prefix = "")
  {
    baseDir->cd();
    if ( prefix.Length() != 0 ) prefix += " ";
    
    const double ptMin = 0;
    const double ptMax = 500;
    const double ptBinWidth = 5;

    hPt = new TH1F("hPt", prefix+Form("Muon tansverse momentum;Transverse momentum p_{T} [GeV/c];Entries per %.0f GeV", ptBinWidth), TMath::Nint((ptMax-ptMin)/ptBinWidth), ptMin, ptMax);
    hEta = new TH1F("hEta", prefix+"Muon pseudorapidity;Pseudorapidity #eta", 100, -2.5, 2.5);
    hPhi = new TH1F("hPhi", prefix+"Muon azimuthal angle;Azimuthal angle #phi [Radian]", 100, -3.15, 3.15); 
    hQ = new TH1F("hQ", prefix+"Muon charge;Electric charge", 5, -2.5, 2.5);

    hRelIso = new TH1F("hRelIso", prefix+"Muon relative isolation;Relative isolation", 100, 0, 1);

    hPt->SetMinimum(0);
    hEta->SetMinimum(0);
    hPhi->SetMinimum(0);
    hQ->SetMinimum(0);
    hRelIso->SetMinimum(0);
  };

  void Fill(const pat::Muon& muon)
  {
    hPt->Fill(muon.pt());
    hEta->Fill(muon.eta());
    hPhi->Fill(muon.phi());
    hQ->Fill(muon.charge());

    hRelIso->Fill((muon.trackIso()+muon.caloIso())/muon.pt());
  };

  void Scale(const double scale)
  {
    hPt->Scale(scale);
    hEta->Scale(scale);
    hPhi->Scale(scale);
    hQ->Scale(scale);
  };

  TH1F* hPt, * hEta, * hPhi, * hQ;
  TH1F* hRelIso;
};

struct HElectron
{
  HElectron(TDirectory* baseDir, const double scale = 1, TString prefix = "")
  {
    baseDir->cd();
    if ( prefix.Length() != 0 ) prefix += " ";
    
    const double ptMin = 0;
    const double ptMax = 500;
    const double ptBinWidth = 5;
    hPt = new TH1F("hPt", prefix+Form("Electron tansverse momentum;Transverse momentum p_{T} [GeV/c];Entries per %.0f GeV/c", ptBinWidth), TMath::Nint((ptMax-ptMin)/ptBinWidth), ptMin, ptMax);
    hEta = new TH1F("hEta", prefix+"Electron pseudorapidity;Pseudorapidity #eta", 100, -2.5, 2.5);
    hPhi = new TH1F("hPhi", prefix+"Electron azimuthal angle;Azimuthal angle #phi [Radian]", 100, -3.15, 3.15); 
    hQ = new TH1F("hQ", prefix+"Electron charge;Electric charge", 5, -2.5, 2.5);

    hRelIso = new TH1F("hRelIso", prefix+"Electron relative isolation;Relative isolation", 100, 0, 1);

    hPt->SetMinimum(0);
    hEta->SetMinimum(0);
    hPhi->SetMinimum(0);
    hQ->SetMinimum(0);
    hRelIso->SetMinimum(0);
  };

  void Fill(const pat::Electron& electron)
  {
    hPt->Fill(electron.pt());
    hEta->Fill(electron.eta());
    hPhi->Fill(electron.phi());
    hQ->Fill(electron.charge());
    hRelIso->Fill((electron.trackIso()+electron.caloIso())/electron.pt());
  };

  void Scale(const double scale)
  {
    hPt->Scale(scale);
    hEta->Scale(scale);
    hPhi->Scale(scale);
    hQ->Scale(scale);
    hRelIso->Scale(scale);
  };

  TH1F* hPt, * hEta, * hPhi, * hQ;
  TH1F* hRelIso;
};

struct HComposite
{
  HComposite(TDirectory* baseDir, const double scale = 1, TString prefix = "")
  {
    baseDir->cd();
    if ( prefix.Length() != 0 ) prefix += " ";

    const double massMin = 0;
    const double massMax = 200;
    const double massBinWidth = 2;

    const double ptMin = 0;
    const double ptMax = 500;
    const double ptBinWidth = 5;

    hMass = new TH1F("hMass", prefix+Form("Candidate invariant mass;Invariant mass [GeV/c^{2}];Entries per %.0f GeV/c^{2}", massBinWidth), TMath::Nint((massMax-massMin)/massBinWidth), massMin, massMax);
    hPt = new TH1F("hPt", prefix+Form("Candidate transverse momentum;Transverse momentum p_{T} [GeV/c];Entries per %.0f GeV/c", ptBinWidth), TMath::Nint((ptMax-ptMin)/ptBinWidth), ptMin, ptMax);
    hEta = new TH1F("hEta", prefix+"Candidate pseudorapidity;Pseudorapidity #eta", 100, -2.5, 2.5);
    hPhi = new TH1F("hPhi", prefix+"Candidate azimuthal angle;Azimuthal angle #phi [Radian]", 100, -3.15, 3.15);
    hQ = new TH1F("hQ", prefix+"Candidate charge;Electric charge", 5, -2.5, 2.5);

    hPt->SetMinimum(0);
    hEta->SetMinimum(0);
    hPhi->SetMinimum(0);
    hQ->SetMinimum(0);

    hPt->Scale(scale);
    hEta->Scale(scale);
    hPhi->Scale(scale);
    hQ->Scale(scale);
  };

  void Fill(const pat::CompositeCandidate& cand)
  {
    hMass->Fill(cand.mass());
    hPt->Fill(cand.pt());
    hEta->Fill(cand.eta());
    hPhi->Fill(cand.phi());
    hQ->Fill(cand.charge());
  };

  void Scale(const double scale)
  {
    hMass->Scale(scale);
    hPt->Scale(scale);
    hEta->Scale(scale);
    hPhi->Scale(scale);
    hQ->Scale(scale);
  };

  TH1F* hMass, * hPt, * hEta, * hPhi, * hQ;
};

class FWLiteAnalyzerBase
{
public:
  FWLiteAnalyzerBase(const std::string outFileName = "result.root", const bool verbose = true);
  ~FWLiteAnalyzerBase()
  {
    outFile_->Write();
  };

  void SetLumi(const double lumi);
  void AddMCSample(const std::string name, const std::string inputFile, const double xsec, const long nGenEvent);
  void ProcessEvent();
  void ListDataFiles();

  typedef std::map<const std::string, std::vector<std::string> > FileMap;

protected:
  virtual void Analyze(const std::string& channelName, const std::vector<std::string>& files) = 0;

  TDirectory* MakeDirectory(const std::string& path);

  TFile* outFile_;
  FileMap mcSamples_;

  double lumi_;
  std::map<const std::string, double> mcScaleFactors_;

  bool verbose_;
  bool isEventLoaded_;
};

#endif

