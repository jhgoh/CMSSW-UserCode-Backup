#include "HiggsAnalysis/DoublyChargedHiggs/interface/FWLiteAnalyzerBase.h"

#include "CommonTools/Utils/interface/PtComparator.h"

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
    const double scale = mcScaleFactors_[channelName];
    if ( verbose_ ) cout << "@@@@ Scale factor is " << scale << endl;

    TDirectory* dir = MakeDirectory(channelName);

    TDirectory* muonDir = dir->mkdir("m", "Muons");
    muonDir->cd();
    hNMuon_ = new TH1F("hNMuon", "Number of muons;Number of muons", 6, -0.5, 5.5);
    hNMuon_->SetMinimum(0);

    hMuon_ = new HMuon(muonDir, scale);
    hMuon1_ = new HMuon(muonDir->mkdir("m1", "Leading muons"), scale, "Leading");
    hMuon2_ = new HMuon(muonDir->mkdir("m2", "2nd leading muons"), scale, "2nd leading");
    hMuon3_ = new HMuon(muonDir->mkdir("m3", "3rd leading muons"), scale, "3rd leading");

    TDirectory* electronDir = dir->mkdir("e", "Electrons");
    electronDir->cd();
    hNElectron_ = new TH1F("hNElectron", "Number of electrons;Number of electrons", 6, -0.5, 5.5);
    hNElectron_->SetMinimum(0);

    hElectron_ = new HElectron(electronDir, scale);
    hElectron1_ = new HElectron(electronDir->mkdir("e1", "Leading electrons"), scale, "Leading");
    hElectron2_ = new HElectron(electronDir->mkdir("e2", "2nd leading electrons"), scale, "2nd leading");
    hElectron3_ = new HElectron(electronDir->mkdir("e3", "3nd leading electrons"), scale, "3rd leading");

    TDirectory* posEMuCandDir = dir->mkdir("posEMuCand", "H^{++} #rightarrow e#mu candidates");
    posEMuCandDir->cd();
    hNPosEMuCand_ = new TH1F("hNCand", "Number of H^{++} #rightarrow e#mu candidates;Number of candidates", 4, -0.5, 3.5);
    hNPosEMuCand_->SetMinimum(0);

    hPosEMuCand_ = new HComposite(posEMuCandDir, scale);
    hPosEMuCand1_ = new HComposite(posEMuCandDir->mkdir("posEMu1", "Leading H^{++} candidate"), scale, "Leading");
    hPosEMuCand2_ = new HComposite(posEMuCandDir->mkdir("posEMu2", "2nd leading H^{++} candidate"), scale, "2nd leading");

    TDirectory* negEMuCandDir = dir->mkdir("negEMuCand", "H^{--} #rightarrow e#mu candidates");
    negEMuCandDir->cd();
    hNNegEMuCand_ = new TH1F("hNCand", "Number of H^{--} #rightarrow e#mu candidates;Number of candidates", 4, -0.5, 3.5);
    hNNegEMuCand_->SetMinimum(0);

    hNegEMuCand_ = new HComposite(negEMuCandDir, scale);
    hNegEMuCand1_ = new HComposite(negEMuCandDir->mkdir("negEMu1", "Leading H^{--} candidate"), scale, "Leading");
    hNegEMuCand2_ = new HComposite(negEMuCandDir->mkdir("negEMu2", "2nd leading H^{--} candidate"), scale, "2nd leading");

    fwlite::ChainEvent event(files);
    const int nTotalEvent = event.size();
    int nProcessedEvent = 0;

    for ( event.toBegin(); !event.atEnd(); ++event )
    {
      const int eventFraction = 1000*(++nProcessedEvent)/nTotalEvent;
      if ( eventFraction % 10 == 0 ) 
      {
        cout << "@@@@ Processing " << channelName << " : " 
             << eventFraction/10. << "% done\r";
      }

      // First, scan for the leptons
      fwlite::Handle<std::vector<pat::Muon> > muonHandle;
      muonHandle.getByLabel(event, "goodPatMuons");

      fwlite::Handle<std::vector<pat::Electron> > electronHandle;
      electronHandle.getByLabel(event, "goodPatElectrons");

      fwlite::Handle<std::vector<pat::CompositeCandidate> > emuHandle;
      emuHandle.getByLabel(event, "dhCandProducerToEM");

      const int nMuon = muonHandle->size();
      const int nElectron = electronHandle->size();
      const int nEmuCand = emuHandle->size();

      // Basic 4-lepton cut
      if ( nMuon < 2 || nElectron < 2 || nEmuCand < 2 ) continue;

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

      // Scan for the DH candidates
      // Copy and sort by pT
      vector<pat::CompositeCandidate> posEMuCands, negEMuCands;
      for ( vector<pat::CompositeCandidate>::const_iterator emuCand = emuHandle->begin();
            emuCand != emuHandle->end(); ++emuCand )
      {
        if ( emuCand->pt() < 20 ) continue;
        const reco::Candidate* dau1 = emuCand->daughter(0);
        const reco::Candidate* dau2 = emuCand->daughter(1);

        if ( dau1->pt() < 10 || dau2->pt() < 10 ) continue;

        const pat::Muon* muon = dynamic_cast<const pat::Muon*>(dau1->isMuon() ? dau1 : dau2);
        const pat::Electron* electron = dynamic_cast<const pat::Electron*>(dau1->isElectron() ? dau1 : dau2);
        if ( !muon || !electron ) continue;

        const double muonRelIso = (muon->trackIso()+muon->caloIso())/muon->pt();
        const double electronRelIso = (electron->trackIso()+electron->caloIso())/electron->pt();
        if ( muonRelIso > 0.3 || electronRelIso > 0.2 ) continue;

        if ( emuCand->charge() == +2 )
        {
          hPosEMuCand_->Fill(*emuCand);
          posEMuCands.push_back(*emuCand);
        }
        else if ( emuCand->charge() == -2 )
        {
          hNegEMuCand_->Fill(*emuCand);
          negEMuCands.push_back(*emuCand);
        }
      }
      sort(posEMuCands.begin(), posEMuCands.end(), GreaterByPt<pat::CompositeCandidate>());
      sort(negEMuCands.begin(), negEMuCands.end(), GreaterByPt<pat::CompositeCandidate>());

      const int nPosEMuCand = posEMuCands.size();
      const int nNegEMuCand = negEMuCands.size();

      hNPosEMuCand_->Fill(nPosEMuCand);
      hNNegEMuCand_->Fill(nNegEMuCand);

      if ( nPosEMuCand > 1 ) hPosEMuCand2_->Fill(posEMuCands[1]);
      if ( nPosEMuCand > 0 ) hPosEMuCand1_->Fill(posEMuCands[0]);

      if ( nNegEMuCand > 1 ) hNegEMuCand2_->Fill(negEMuCands[1]);
      if ( nNegEMuCand > 0 ) hNegEMuCand1_->Fill(negEMuCands[0]);
    }
    cout << "@@@@ Finished to process " << channelName << endl;

    // Apply scales
    hNMuon_->Scale(scale);
    hNElectron_->Scale(scale);
    hNPosEMuCand_->Scale(scale);
    hNNegEMuCand_->Scale(scale);

    hMuon_->Scale(scale);
    hMuon1_->Scale(scale);
    hMuon2_->Scale(scale);
    hMuon3_->Scale(scale);

    hElectron_->Scale(scale);
    hElectron1_->Scale(scale);
    hElectron2_->Scale(scale);
    hElectron3_->Scale(scale);

    hPosEMuCand_->Scale(scale);
    hPosEMuCand1_->Scale(scale);
    hPosEMuCand2_->Scale(scale);

    hNegEMuCand_->Scale(scale);
    hNegEMuCand1_->Scale(scale);
    hNegEMuCand2_->Scale(scale);
  };

private:
  TH1F* hNMuon_, * hNElectron_;
  TH1F* hNPosEMuCand_, * hNNegEMuCand_;

  HMuon* hMuon_, * hMuon1_, * hMuon2_, * hMuon3_;
  HElectron* hElectron_, * hElectron1_, * hElectron2_, * hElectron3_;

  HComposite* hPosEMuCand_, * hPosEMuCand1_, * hPosEMuCand2_;
  HComposite* hNegEMuCand_, * hNegEMuCand1_, * hNegEMuCand2_;

};

