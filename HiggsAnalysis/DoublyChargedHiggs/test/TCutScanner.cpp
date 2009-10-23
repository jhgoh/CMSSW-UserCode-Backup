#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"

#include <iostream>
#include <vector>
#include <string>
#include <utility>

#endif

using namespace std;

class TCutScanner
{
public:
  TCutScanner(const string options = "DRAW")
  {
    cout << "Loading cut scanner" << endl;

    //options = options.lower();
    optDraw_ = options.find("draw") > 0;
  };

  ~TCutScanner()
  {
  };

  void SetSignal(vector<string> files, double scale)
  {
    sigSample_ = files;
    sigScale_ = scale;
  };

  void SetBkg(vector<string> files, double scale)
  {
    bkgSamples_.push_back(files);
    bkgScales_.push_back(scale);
  };

  void Run()
  {
    hDimuonPt_ = new TH1F("hDimuonPt", "Dimuon pt", 10, 0, 1000);

    cout << "Starting Run" << endl;

    // Set files to be linked with chains
    fwlite::ChainEvent sigEvent(sigSample_);
    std::vector<fwlite::ChainEvent> bkgEvents;
    for(unsigned int i=0; i<bkgSamples_.size(); ++i)
    {
      fwlite::ChainEvent bkgEvent(bkgSamples_[i]);
      bkgEvents.push_back(bkgEvent);
    }

    cout << "Read up data files" << endl;

    // Scan over muon quality

    cout << sigEvent.size() << endl;
    for(sigEvent.toBegin(); !sigEvent.atEnd(); ++sigEvent)
    {
      fwlite::Handle<pat::CompositeCandidateCollection> dimuonHandle;
      dimuonHandle.getByLabel(sigEvent, "deltaToMuMu");

      const pat::CompositeCandidateCollection* dimuonCands = dimuonHandle.ptr();

      for(unsigned int i=0; i<dimuonCands->size(); ++i)
      {
        const pat::CompositeCandidate& dimuon = dimuonCands->at(i);
        hDimuonPt_->Fill(dimuon.pt());
      }
    }

    hDimuonPt_->Draw();
  };

private:
  std::vector<std::string> sigSample_;
  double sigScale_;

  std::vector<std::vector<std::string> > bkgSamples_;
  std::vector<double> bkgScales_;

  // Drawing histograms
  bool optDraw_;

  TH1F* hDimuonPt_;
};
