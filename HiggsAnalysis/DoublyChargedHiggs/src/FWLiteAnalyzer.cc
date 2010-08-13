#include "HiggsAnalysis/DoublyChargedHiggs/interface/FWLiteAnalyzerBase.h"

#include "CommonTools/Utils/interface/PtComparator.h"

#include <list>

using namespace std;

class FWLiteAnalyzerEMu : public FWLiteAnalyzerBase
{
public:
  FWLiteAnalyzerEMu(const string outFileName, const bool verbose):
    FWLiteAnalyzerBase(outFileName, verbose)
  {
  };

  void Analyze(const string& channelName, const vector<string>& files)
  {
    // Process signals
    TDirectory* dir = outFile_->mkdir(channelName.c_str(), channelName.c_str());

    TDirectory* muonDir = dir->mkdir("m", "Muons");
    muonDir->cd();
    hNMuon_ = new TH1F("hNMuon", "Number of muons;Number of muons", 6, -0.5, 5.5);
    hNMuon_->SetMinimum(0);

    hMuon_ = new HMuon(muonDir);
    hMuon1_ = new HMuon(muonDir->mkdir("m1", "Leading muons"), "Leading");
    hMuon2_ = new HMuon(muonDir->mkdir("m2", "2nd leading muons"), "2nd leading");
    hMuon3_ = new HMuon(muonDir->mkdir("m3", "3rd leading muons"), "3rd leading");

    TDirectory* electronDir = dir->mkdir("e", "Electrons");
    electronDir->cd();
    hNElectron_ = new TH1F("hNElectron", "Number of electrons;Number of electrons", 6, -0.5, 5.5);
    hNElectron_->SetMinimum(0);

    hElectron_ = new HElectron(electronDir);
    hElectron1_ = new HElectron(electronDir->mkdir("e1", "Leading electrons"), "Leading");
    hElectron2_ = new HElectron(electronDir->mkdir("e2", "2nd leading electrons"), "2nd leading");
    hElectron3_ = new HElectron(electronDir->mkdir("e3", "3nd leading electrons"), "3rd leading");

    fwlite::ChainEvent event(files);

    for ( event.toBegin(); !event.atEnd(); ++event )
    {
      fwlite::Handle<std::vector<pat::Muon> > muonHandle;
      muonHandle.getByLabel(event, "goodPatMuons");

      fwlite::Handle<std::vector<pat::Electron> > electronHandle;
      electronHandle.getByLabel(event, "goodPatElectrons");

      const int nMuon = muonHandle->size();
      const int nElectron = electronHandle->size();

      // Sort by pT
      vector<pat::Muon> muons(nMuon);
      vector<pat::Electron> electrons(nElectron);

      copy(muonHandle->begin(), muonHandle->end(), muons.begin());
      copy(electronHandle->begin(), electronHandle->end(), electrons.begin());

      sort(muons.begin(), muons.end(), GreaterByPt<pat::Muon>());
      sort(electrons.begin(), electrons.end(), GreaterByPt<pat::Electron>());

      hNMuon_->Fill(nMuon);
      hNElectron_->Fill(nElectron);

      if ( nMuon > 2 ) hMuon3_->Fill(muons[2]);
      if ( nMuon > 1 ) hMuon2_->Fill(muons[1]);
      if ( nMuon > 0 ) hMuon1_->Fill(muons[0]);

      if ( nElectron > 2 ) hElectron3_->Fill(electrons[2]);
      if ( nElectron > 1 ) hElectron2_->Fill(electrons[1]);
      if ( nElectron > 0 ) hElectron1_->Fill(electrons[0]);

      for ( vector<pat::Muon>::const_iterator muon = muons.begin();
            muon != muons.end(); ++muon )
      {
        hMuon_->Fill(*muon);
      }

      for ( vector<pat::Electron>::const_iterator electron = electrons.begin();
            electron != electrons.end(); ++electron )
      {
        hElectron_->Fill(*electron);
      }
    }
  };

private:
  TH1F* hNMuon_, * hNElectron_;

  HMuon* hMuon_, * hMuon1_, * hMuon2_, * hMuon3_;
  HElectron* hElectron_, * hElectron1_, * hElectron2_, * hElectron3_;

};

