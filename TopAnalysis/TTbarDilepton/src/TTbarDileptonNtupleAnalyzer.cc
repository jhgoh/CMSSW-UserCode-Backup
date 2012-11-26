#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TMath.h"

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/VectorUtil.h"
//#include "Math/VectorUtil_Cint.h"
using namespace ROOT::Math::VectorUtil;

using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

void printEntryFraction(int i, int n)
{
  static int digit = 1;
  if ( i%max(1, digit/10) == 0 )
  {
    cout << Form(TString("% ")+Form("%d", int(log10(n))+2)+"d/%d processed\r", i, n);
    cout.flush();
  }
  if ( i%digit == 0 ) digit *= 10;

  if ( i==(n-1) )
  {
    cout << Form(TString("% ")+Form("%d", int(log10(n))+2)+"d/%d processed", i, n) << endl;
    digit = 1;
  }
}

class TTbarDileptonNtupleAnalyzer
{
public:
  struct Event;
  struct Selector;
  struct HistSet;

  TTbarDileptonNtupleAnalyzer(const std::string inputFileName, const std::string mode, const std::string outputFileName);
  ~TTbarDileptonNtupleAnalyzer() {};
  void analyze(int verboseLevel=1);
  void endJob(int verboseLevel=1);

private:

public:
  struct Event
  { 
    Event(TDirectory* dir);
    ~Event(); 

    bool isValid_;

    double nEvent_;
    TTree* eventTree_;
    TTree* genTree_;

    int run_, lumi_, event_;
    double weight_;
    int nVertex_;
    std::vector<LorentzVector>* electrons_, * muons_;
    std::vector<int>* electrons_Q_, * muons_Q_;
    std::vector<double>* electrons_Iso_, * muons_Iso_;
    LorentzVector* met_;
    std::vector<LorentzVector>* jets_;
    std::vector<double>* jets_bTag_;

    std::vector<LorentzVector>* genMuons_, * genElectrons_;
    std::vector<LorentzVector>* genMuonNus_, * genElectronNus_;
    std::vector<int>* jets_motherId_;
  };

  struct HistSet
  {
    HistSet(TDirectory* dir);
    void write();

    typedef TH1F* H1;
    TDirectory* dir_;

    H1 hNEvent_;
    H1 hNVertex_;
  
    H1 hNMuon_;
    H1 hMuon1_pt_, hMuon1_eta_, hMuon1_phi_, hMuon1_iso_;
    H1 hMuon2_pt_, hMuon2_eta_, hMuon2_phi_, hMuon2_iso_;

    H1 hNElectron_;
    H1 hElectron1_pt_, hElectron1_eta_, hElectron1_phi_, hElectron1_iso_;
    H1 hElectron2_pt_, hElectron2_eta_, hElectron2_phi_, hElectron2_iso_;

    H1 hNJets_;
    H1 hNTbjets_, hNMbjets_, hNLbjets_;
    H1 hJet1_pt_, hJet1_eta_, hJet1_phi_, hJet1_btag_;
    H1 hJet2_pt_, hJet2_eta_, hJet2_phi_, hJet2_btag_;
    H1 hJet3_pt_, hJet3_eta_, hJet3_phi_, hJet3_btag_;
    H1 hJet4_pt_, hJet4_eta_, hJet4_phi_, hJet4_btag_;

    H1 hMet_pt_, hMet_phi_;

    H1 hLL_q_, hLL_m_, hLL_pt_, hLL_eta_, hLL_phi_;
    H1 hMM_q_, hMM_m_, hMM_pt_, hMM_eta_, hMM_phi_;
    H1 hEE_q_, hEE_m_, hEE_pt_, hEE_eta_, hEE_phi_;
    H1 hME_q_, hME_m_, hME_pt_, hME_eta_, hME_phi_;

    H1 hTTbarVsum_m_, hTTbarVsum_pt_, hTTbarVsum_eta_, hTTbarVsum_phi_;

    H1 hJJ_m_, hJJ_pt_, hJJ_eta_, hJJ_phi_;
    H1 hBB_m_, hBB_pt_, hBB_eta_, hBB_phi_;

    void fill(Selector& selector, int cutStep);
  };

  struct Selector
  {
    Selector(Event* event, int mode);

    enum MODE {MM = 0, EE, ME, END};
    int mode_; 

    bool isGoodEvent(int cutStep) { return isGoodEvent_.at(cutStep); };

    Event* event_;
    std::vector<bool> isGoodEvent_;

    const static int nCutSteps;

    const static double cut_btagTight_;
    const static double cut_btagMedium_;
    const static double cut_btagLoose_;
  };

public:
  Selector* selector_;

private:
  TFile* inputFile_, * outputFile_;
  Event* event_;
  std::vector<HistSet*> hists_;

  int mode_;
};

const int TTbarDileptonNtupleAnalyzer::Selector::nCutSteps = 7;
const double TTbarDileptonNtupleAnalyzer::Selector::cut_btagTight_ = 0.898;
const double TTbarDileptonNtupleAnalyzer::Selector::cut_btagMedium_ = 0.679;
const double TTbarDileptonNtupleAnalyzer::Selector::cut_btagLoose_ = 0.244;

TTbarDileptonNtupleAnalyzer::TTbarDileptonNtupleAnalyzer(const std::string inputFileName, const std::string mode, const std::string outputFileName)
{
  event_ = 0;

  inputFile_ = TFile::Open(inputFileName.c_str());
  if ( !inputFile_ or !inputFile_->IsOpen() or inputFile_->IsZombie() )
  {
    cout << "File open error, " << inputFileName << endl;
    return;
  }
  TDirectory* eventDir = inputFile_->GetDirectory(mode.c_str());
  if ( !eventDir )
  {
    cout << "Incorrect directory structure, " << inputFileName << ":/" << mode << endl;
    return;
  }

  outputFile_ = new TFile(outputFileName.c_str(), "RECREATE");
  if ( !outputFile_ or !outputFile_->IsOpen() or outputFile_->IsZombie() )
  {
    cout << "Output file open error, " << outputFileName << endl;
    return;
  }

  event_ = new Event(eventDir);

  for ( int i=0; i<Selector::nCutSteps; ++i )
  {
    hists_.push_back(new HistSet(outputFile_->mkdir(Form("Step%2d", i))));
  }

  if ( mode == "mm" ) mode_ = Selector::MM;
  else if ( mode == "ee" ) mode_ = Selector::EE;
  else if ( mode == "me" ) mode_ = Selector::ME;
  else mode_ = Selector::END;

  if ( mode_ == Selector::END ) event_ = 0;
}

void TTbarDileptonNtupleAnalyzer::endJob(int verboseLevel)
{
  if ( !event_ ) return;

  if ( verboseLevel >= 1 )
  {
    cout << "Finishing job for " << inputFile_->GetName() << endl;
    cout << "  writing to " << outputFile_->GetName() << endl;
  }
  delete event_;
  for ( int i=0; i<Selector::nCutSteps; ++i )
  {
    hists_.at(i)->write();
  }
  outputFile_->Close();
  for ( int i=0; i<Selector::nCutSteps; ++i )
  {
    delete hists_.at(i);
  }
}

TTbarDileptonNtupleAnalyzer::Event::Event(TDirectory* dir)
{
  nEvent_ = ((TH1F*)dir->Get("hEventCounter"))->GetBinContent(1);

  electrons_     = new std::vector<LorentzVector>();
  electrons_Q_   = new std::vector<int>();
  electrons_Iso_ = new std::vector<double>();

  muons_     = new std::vector<LorentzVector>();
  muons_Q_   = new std::vector<int>();
  muons_Iso_ = new std::vector<double>();

  met_       = new LorentzVector();
  jets_      = new std::vector<LorentzVector>();
  jets_bTag_ = new std::vector<double>();

  genMuons_       = new std::vector<LorentzVector>();
  genMuonNus_     = new std::vector<LorentzVector>();
  genElectrons_   = new std::vector<LorentzVector>();
  genElectronNus_ = new std::vector<LorentzVector>();
  jets_motherId_  = new std::vector<int>();

  eventTree_ = (TTree*)dir->Get("event");
  genTree_ = (TTree*)dir->Get("gen");

  eventTree_->SetBranchAddress("run", &run_);
  eventTree_->SetBranchAddress("lumi", &lumi_);
  eventTree_->SetBranchAddress("event", &event_);
  eventTree_->SetBranchAddress("weight", &weight_);
  eventTree_->SetBranchAddress("nVertex", &nVertex_);

  eventTree_->SetBranchAddress("electrons", &electrons_);
  eventTree_->SetBranchAddress("electrons_Q", &electrons_Q_);
  eventTree_->SetBranchAddress("electrons_Iso", &electrons_Iso_);
  eventTree_->SetBranchAddress("muons", &muons_);
  eventTree_->SetBranchAddress("muons_Q", &muons_Q_);
  eventTree_->SetBranchAddress("muons_Iso", &muons_Iso_);

  eventTree_->SetBranchAddress("met", &met_);

  eventTree_->SetBranchAddress("jets", &jets_);
  eventTree_->SetBranchAddress("jets_bTag", &jets_bTag_);

  genTree_->SetBranchAddress("electrons", &genElectrons_);
  genTree_->SetBranchAddress("electronNus", &genElectronNus_);
  genTree_->SetBranchAddress("muons", &genMuons_);
  genTree_->SetBranchAddress("muonNus", &genMuonNus_);
  genTree_->SetBranchAddress("jets_motherId", &jets_motherId_);

  isValid_ = false;
  do
  {
    if ( eventTree_->GetEntries() != genTree_->GetEntries() )
    {
      cout << "Reco/Gen event tree nEntries does not match\n";
      break;
    }

    isValid_ = true;
  }
  while ( false );

  if ( !isValid_ ) cout << "Event tree is not valid\n";

}

TTbarDileptonNtupleAnalyzer::Event::~Event()
{
  if ( electrons_     ) delete electrons_    ;
  if ( electrons_Q_   ) delete electrons_Q_  ;
  if ( electrons_Iso_ ) delete electrons_Iso_;

  if ( muons_     ) delete muons_    ;
  if ( muons_Q_   ) delete muons_Q_  ;
  if ( muons_Iso_ ) delete muons_Iso_;

  if ( met_       ) delete met_      ;
  if ( jets_      ) delete jets_     ;
  if ( jets_bTag_ ) delete jets_bTag_;

  if ( genMuons_       ) delete genMuons_      ;
  if ( genMuonNus_     ) delete genMuonNus_    ;
  if ( genElectrons_   ) delete genElectrons_  ;
  if ( genElectronNus_ ) delete genElectronNus_;
  if ( jets_motherId_  ) delete jets_motherId_ ;
}

TTbarDileptonNtupleAnalyzer::HistSet::HistSet(TDirectory* dir)
{
  dir_ = dir;
  dir->cd();

  hNEvent_ = new TH1F("hNEvent", "Number of event", 3, 1, 4);
  hNEvent_->GetXaxis()->SetBinLabel(1, "nTotal");
  hNEvent_->GetXaxis()->SetBinLabel(2, "nPass");
  hNEvent_->GetXaxis()->SetBinLabel(3, "nPass_Weighted");
  hNVertex_ = new TH1F("hNVertex", "Number of vertex", 50, 0, 50);

  hNMuon_     = new TH1F("hNMuon"    , "Number of muons;;Events", 5, 0, 5);
  hMuon1_pt_  = new TH1F("hMuon1_pt" , "Muon1 p_{T};p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hMuon1_eta_ = new TH1F("hMuon1_eta", "Muon1 #eta;#eta;Events", 100, -2.5, 2.5);
  hMuon1_phi_ = new TH1F("hMuon1_phi", "Muon1 #phi;#phi;Events", 100, -TMath::Pi(), TMath::Pi());
  hMuon1_iso_ = new TH1F("hMuon1_iso", "Muon1 isolation;Relative isolation;Events", 100, 0, 0.5);
  hMuon2_pt_  = new TH1F("hMuon2_pt" , "Muon2 p_{T};p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hMuon2_eta_ = new TH1F("hMuon2_eta", "Muon2 #eta;#eta;Events", 100, -2.5, 2.5);
  hMuon2_phi_ = new TH1F("hMuon2_phi", "Muon2 #phi;#phi;Events", 100, -TMath::Pi(), TMath::Pi());
  hMuon2_iso_ = new TH1F("hMuon2_iso", "Muon2 isolation;Relative isolation;Events", 100, 0, 0.5);

  hNElectron_     = new TH1F("hNElectron"    , "Number of muons;;Events", 5, 0, 5);
  hElectron1_pt_  = new TH1F("hElectron1_pt" , "Electron1 p_{T};p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hElectron1_eta_ = new TH1F("hElectron1_eta", "Electron1 #eta;#eta;Events", 100, -2.5, 2.5);
  hElectron1_phi_ = new TH1F("hElectron1_phi", "Electron1 #phi;#phi;Events", 100, -TMath::Pi(), TMath::Pi());
  hElectron1_iso_ = new TH1F("hElectron1_iso", "Electron1 isolation;Relative isolation;Events", 100, 0, 0.5);
  hElectron2_pt_  = new TH1F("hElectron2_pt" , "Electron2 p_{T};p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hElectron2_eta_ = new TH1F("hElectron2_eta", "Electron2 #eta;#eta;Events", 100, -2.5, 2.5);
  hElectron2_phi_ = new TH1F("hElectron2_phi", "Electron2 #phi;#phi;Events", 100, -TMath::Pi(), TMath::Pi());
  hElectron2_iso_ = new TH1F("hElectron2_iso", "Electron2 isolation;Relative isolation;Events", 100, 0, 0.5);

  hNJets_    = new TH1F("hNJets"   , "Number of jets;;Events", 10, 0, 10);
  hNTbjets_  = new TH1F("hNTbjets" , "Number of Tight b jets;;Events", 10, 0, 10);
  hNMbjets_  = new TH1F("hNMbjets" , "NUmber of Medium b jets;;Events", 10, 0, 10);
  hNLbjets_  = new TH1F("hNLbjets" , "NUmber of Loose b jets;;Events", 10, 0, 10);
  hJet1_pt_   = new TH1F("hJet1_pt"  , "Jet1 p_{T};p_{T} (GeV/c;Events per 5GeV/c", 100, 0, 500);
  hJet1_eta_  = new TH1F("hJet1_eta" , "Jet1 #eta;#eta;Events", 100, -2.5, 2.5);
  hJet1_phi_  = new TH1F("hJet1_phi" , "Jet1 #phi;#phi;Events", 100, -TMath::Pi(), TMath::Pi());
  hJet1_btag_ = new TH1F("hJet1_btag", "Jet1 b tag;b discriminator;Events", 100, 0, 1);
  hJet2_pt_   = new TH1F("hJet2_pt"  , "Jet2 p_{T};p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hJet2_eta_  = new TH1F("hJet2_eta" , "Jet2 #eta;#eta;Events", 100, -2.5, 2.5);
  hJet2_phi_  = new TH1F("hJet2_phi" , "Jet2 #phi;#phi;Events", 100, -TMath::Pi(), TMath::Pi());
  hJet2_btag_ = new TH1F("hJet2_btag", "Jet1 b tag;b discriminator;Events", 100, 0, 1);
  hJet3_pt_   = new TH1F("hJet3_pt"  , "Jet3 p_{T};p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hJet3_eta_  = new TH1F("hJet3_eta" , "Jet3 #eta;#eta;Events", 100, -2.5, 2.5);
  hJet3_phi_  = new TH1F("hJet3_phi" , "Jet3 #phi;#phi;Events", 100, -TMath::Pi(), TMath::Pi());
  hJet3_btag_ = new TH1F("hJet3_btag", "Jet1 b tag;b discriminator;Events", 100, 0, 1);
  hJet4_pt_   = new TH1F("hJet4_pt"  , "Jet4 p_{T};p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hJet4_eta_  = new TH1F("hJet4_eta" , "Jet4 #eta;#eta;Events", 100, -2.5, 2.5);
  hJet4_phi_  = new TH1F("hJet4_phi" , "Jet4 #phi;#phi;Events", 100, -TMath::Pi(), TMath::Pi());
  hJet4_btag_ = new TH1F("hJet4_btag", "Jet1 b tag;b discriminator;Events", 100, 0, 1);

  hMet_pt_  = new TH1F("hMet_pt" , "Met p_{T};p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hMet_phi_ = new TH1F("hMet_phi", "Met #phi;#phi;Events", 100, -TMath::Pi(), TMath::Pi());

  hLL_q_ = new TH1F("hLL_q", "Lepton pair charge;Lepton pair charge;Events", 3, -1.5, 1.5);
  hLL_m_ = new TH1F("hLL_m", "Lepton pair mass;Lepton pair mass (GeV/c^{2});Events per 5GeV/c^{2}", 100, 0, 500);
  hLL_pt_  = new TH1F("hLL_pt" , "Lepton pair p_{T};Lepton pair p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hLL_eta_ = new TH1F("hLL_eta", "Lepton pair #eta;Lepton pair #eta (GeV/c);Events", 100, -2.5, 2.5);
  hLL_phi_ = new TH1F("hLL_phi", "Lepton pair #phi;Lepton pair #phi (GeV/c);Events", 100, -TMath::Pi(), TMath::Pi());

  hMM_q_ = new TH1F("hMM_q", "Muon pair charge;Muon pair charge;Events", 3, -1.5, 1.5);
  hMM_m_ = new TH1F("hMM_m", "Muon pair mass;Muon pair mass (GeV/c^{2});Events per 5GeV/c^{2}", 100, 0, 500);
  hMM_pt_  = new TH1F("hMM_pt" , "Muon pair p_{T};Muon pair p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hMM_eta_ = new TH1F("hMM_eta", "Muon pair #eta;Muon pair #eta (GeV/c);Events", 100, -2.5, 2.5);
  hMM_phi_ = new TH1F("hMM_phi", "Muon pair #phi;Muon pair #phi (GeV/c);Events", 100, -TMath::Pi(), TMath::Pi());

  hEE_q_ = new TH1F("hEE_q", "Electron pair charge;Electron pair charge;Events", 3, -1.5, 1.5);
  hEE_m_ = new TH1F("hEE_m", "Electron pair mass;Electron pair mass (GeV/c^{2});Events per 5GeV/c^{2}", 100, 0, 500);
  hEE_pt_  = new TH1F("hEE_pt" , "Electron pair p_{T};Electron pair p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hEE_eta_ = new TH1F("hEE_eta", "Electron pair #eta;Electron pair #eta (GeV/c);Events", 100, -2.5, 2.5);
  hEE_phi_ = new TH1F("hEE_phi", "Electron pair #phi;Electron pair #phi (GeV/c);Events", 100, -TMath::Pi(), TMath::Pi());
 
  hME_q_ = new TH1F("hME_q", "Muon-Electron pair charge;Muon-Electron pair charge;Events", 3, -1.5, 1.5);
  hME_m_ = new TH1F("hME_m", "Muon-Electron pair mass;Muon-Electron pair mass (GeV/c^{2});Events per 5GeV/c^{2}", 100, 0, 500);
  hME_pt_  = new TH1F("hME_pt" , "Muon-Electron pair p_{T};Muon-Electron pair p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hME_eta_ = new TH1F("hME_eta", "Muon-Electron pair #eta;Muon-Electron pair #eta (GeV/c);Events", 100, -2.5, 2.5);
  hME_phi_ = new TH1F("hME_phi", "Muon-Electron pair #phi;Muon-Electron pair #phi (GeV/c);Events", 100, -TMath::Pi(), TMath::Pi());

  hTTbarVsum_m_   = new TH1F("hTTbarVsum_m"  , "t#bar{t} candidate vector sum mass;t#bar{t} mass (GeV/c^{2});Events per 20GeV/c^{2}", 100, 0, 2000);
  hTTbarVsum_pt_  = new TH1F("hTTbarVsum_pt" , "t#bar{t} candidate vector sum p_{T};t#bar{t} candidate vector sum p_{T} (GeV/c);Events per 5GeV/c^{2}", 100, 0, 500);
  hTTbarVsum_eta_ = new TH1F("hTTbarVsum_eta", "t#bar{t} candidate vector sum #eta;t#bar{t} candidate vector sum #eta;Events", 100, -2.5, 2.5);
  hTTbarVsum_phi_ = new TH1F("hTTbarVsum_phi", "t#bar{t} candidate vector sum #phi;t#bar{t} candidate vector sum #phi;Events", 100, -TMath::Pi(), TMath::Pi());

  hJJ_m_   = new TH1F("hJJ_m"  , "Dijet mass;Dijet mass (GeV/c^{2});Events per 10GeV/c^{2}", 100, 0, 1000);
  hJJ_pt_  = new TH1F("hJJ_pt" , "Dijet p_{T};Dijet p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hJJ_eta_ = new TH1F("hJJ_eta", "Dijet #eta;Dijet #eta;Events per 5GeV/c", 100, -2.5, 2.5);
  hJJ_phi_ = new TH1F("hJJ_phi", "Dijet #phi;Dijet #phi;Events per 5GeV/c", 100, -TMath::Pi(), TMath::Pi());

  hBB_m_   = new TH1F("hBB_m"  , "B tagged dijet mass;Dijet mass(GeV/c^{2});Events per 10GeV/c^{2}", 100, 0, 1000);
  hBB_pt_  = new TH1F("hBB_pt" , "B tagged dijet p_{T};Dijet p_{T} (GeV/c);Events per 5GeV/c", 100, 0, 500);
  hBB_eta_ = new TH1F("hBB_eta", "B tagged dijet #eta;Dijet #eta;Events per 5GeV/c", 100, -2.5, 2.5);
  hBB_phi_ = new TH1F("hBB_phi", "B tagged dijet #phi;Dijet #phi;Events per 5GeV/c", 100, -TMath::Pi(), TMath::Pi());

}

void TTbarDileptonNtupleAnalyzer::HistSet::write()
{
  dir_->cd();

  hNEvent_->Write();
  hNVertex_->Write();

  hNMuon_->Write();
  hMuon1_pt_->Write();
  hMuon1_eta_->Write();
  hMuon1_phi_->Write();
  hMuon1_iso_->Write();
  hMuon2_pt_->Write();
  hMuon2_eta_->Write();
  hMuon2_phi_->Write();
  hMuon2_iso_->Write();

  hNElectron_->Write();
  hElectron1_pt_->Write();
  hElectron1_eta_->Write();
  hElectron1_phi_->Write();
  hElectron1_iso_->Write();
  hElectron2_pt_->Write();
  hElectron2_eta_->Write();
  hElectron2_phi_->Write();
  hElectron2_iso_->Write();

  hNJets_->Write();
  hNTbjets_->Write();
  hNMbjets_->Write();
  hNLbjets_->Write();
  hJet1_pt_->Write();
  hJet1_eta_->Write();
  hJet1_phi_->Write();
  hJet1_btag_->Write();
  hJet2_pt_->Write();
  hJet2_eta_->Write();
  hJet2_phi_->Write();
  hJet2_btag_->Write();
  hJet3_pt_->Write();
  hJet3_eta_->Write();
  hJet3_phi_->Write();
  hJet3_btag_->Write();
  hJet4_pt_->Write();
  hJet4_eta_->Write();
  hJet4_phi_->Write();
  hJet4_btag_->Write();

  hMet_pt_->Write();
  hMet_phi_->Write();

  hLL_q_->Write();
  hLL_m_->Write();
  hLL_pt_->Write();
  hLL_eta_->Write();
  hLL_phi_->Write();

  hMM_q_->Write();
  hMM_m_->Write();
  hMM_pt_->Write();
  hMM_eta_->Write();
  hMM_phi_->Write();

  hEE_q_->Write();
  hEE_m_->Write();
  hEE_pt_->Write();
  hEE_eta_->Write();
  hEE_phi_->Write();

  hME_q_->Write();
  hME_m_->Write();
  hME_pt_->Write();
  hME_eta_->Write();
  hME_phi_->Write();

  hTTbarVsum_m_  ->Write();
  hTTbarVsum_pt_ ->Write();
  hTTbarVsum_eta_->Write();
  hTTbarVsum_phi_->Write();

  hJJ_m_  ->Write();
  hJJ_pt_ ->Write();
  hJJ_eta_->Write();
  hJJ_phi_->Write();

  hBB_m_  ->Write();
  hBB_pt_ ->Write();
  hBB_eta_->Write();
  hBB_phi_->Write();

}

void TTbarDileptonNtupleAnalyzer::HistSet::fill(TTbarDileptonNtupleAnalyzer::Selector& selector, int cutStep)
{
  if ( !selector.isGoodEvent(cutStep) ) return;

  Event& event = *selector.event_;
  const int mode = selector.mode_;

  hNEvent_->Fill(2);
  hNEvent_->Fill(3, event.weight_);

  hNVertex_->Fill(event.nVertex_, event.weight_);
  hNMuon_->Fill(event.muons_->size(), event.weight_);
  hNElectron_->Fill(event.electrons_->size(), event.weight_);
  hNJets_->Fill(event.jets_->size(), event.weight_);

  int zQ = -999;
  LorentzVector zLVec;

  if ( event.muons_->size() > 0 )
  {
    const LorentzVector& lv1 = event.muons_->at(0);
    const int q1 = event.muons_Q_->at(0);
    const double& iso1 = event.muons_Iso_->at(0);

    hMuon1_pt_ ->Fill(lv1.pt() , event.weight_);
    hMuon1_eta_->Fill(lv1.eta(), event.weight_);
    hMuon1_phi_->Fill(lv1.phi(), event.weight_);
    hMuon1_iso_->Fill(iso1     , event.weight_);

    if ( event.muons_->size() > 1 )
    {
      const LorentzVector& lv2 = event.muons_->at(1);
      const int q2 = event.muons_Q_->at(1);
      const double& iso2 = event.muons_Iso_->at(1);

      hMuon2_pt_ ->Fill(lv2.pt() , event.weight_);
      hMuon2_eta_->Fill(lv2.eta(), event.weight_);
      hMuon2_phi_->Fill(lv2.phi(), event.weight_);
      hMuon2_iso_->Fill(iso2     , event.weight_);

      const LorentzVector llLVec = lv1+lv2;
      hMM_q_->Fill(q1+q2, event.weight_);
      hMM_m_->Fill(llLVec.mass(), event.weight_);
      hMM_pt_ ->Fill(llLVec.pt() , event.weight_);
      hMM_eta_->Fill(llLVec.eta(), event.weight_);
      hMM_phi_->Fill(llLVec.phi(), event.weight_);

      if ( mode == Selector::MM )
      {
        zQ = q1+q2;
        zLVec = llLVec;
      }
    }
  }

  if ( event.electrons_->size() > 0 )
  {
    const LorentzVector& lv1 = event.electrons_->at(0);
    const int q1 = event.electrons_Q_->at(0);
    const double& iso1 = event.electrons_Iso_->at(0);

    hElectron1_pt_ ->Fill(lv1.pt() , event.weight_);
    hElectron1_eta_->Fill(lv1.eta(), event.weight_);
    hElectron1_phi_->Fill(lv1.phi(), event.weight_);
    hElectron1_iso_->Fill(iso1     , event.weight_);

    if ( event.electrons_->size() > 1 )
    {
      const LorentzVector& lv2 = event.electrons_->at(1);
      const int q2 = event.electrons_Q_->at(1);
      const double& iso2 = event.electrons_Iso_->at(1);

      hElectron2_pt_ ->Fill(lv2.pt() , event.weight_);
      hElectron2_eta_->Fill(lv2.eta(), event.weight_);
      hElectron2_phi_->Fill(lv2.phi(), event.weight_);
      hElectron2_iso_->Fill(iso2     , event.weight_);

      const LorentzVector llLVec = lv1+lv2;
      hEE_q_->Fill(q1+q2, event.weight_);
      hEE_m_->Fill(llLVec.mass(), event.weight_);
      hEE_pt_ ->Fill(llLVec.pt() , event.weight_);
      hEE_eta_->Fill(llLVec.eta(), event.weight_);
      hEE_phi_->Fill(llLVec.phi(), event.weight_);

      if ( mode == Selector::EE )
      {
        zQ = q1+q2;
        zLVec = llLVec;
      }
    }
  }

  if ( event.muons_->size() > 0 and event.electrons_->size() > 0 )
  {
    const LorentzVector& lv1 = event.muons_->at(0);
    const LorentzVector& lv2 = event.electrons_->at(0);

    const int q1 = event.muons_Q_->at(0);
    const int q2 = event.electrons_Q_->at(0);

    const LorentzVector llLVec = lv1+lv2;
    hME_q_->Fill(q1+q2, event.weight_);
    hME_m_->Fill(llLVec.mass(), event.weight_);
    hME_pt_ ->Fill(llLVec.pt() , event.weight_);
    hME_eta_->Fill(llLVec.eta(), event.weight_);
    hME_phi_->Fill(llLVec.phi(), event.weight_);

    if ( mode == Selector::ME )
    {
      zQ = q1+q2;
      zLVec = llLVec;
    }
  }

  if ( zQ != -999 )
  {
    hLL_q_->Fill(zQ, event.weight_);
    hLL_m_->Fill(zLVec.mass(), event.weight_);
    hLL_pt_ ->Fill(zLVec.pt() , event.weight_);
    hLL_eta_->Fill(zLVec.eta(), event.weight_);
    hLL_phi_->Fill(zLVec.phi(), event.weight_);
    
    if ( event.jets_->size() >= 2 )
    {
      const LorentzVector ttbarVsumLVec = zLVec + event.jets_->at(0) + event.jets_->at(1) + *event.met_;
      hTTbarVsum_m_  ->Fill(ttbarVsumLVec.mass(), event.weight_);
      hTTbarVsum_pt_ ->Fill(ttbarVsumLVec.pt()  , event.weight_);
      hTTbarVsum_eta_->Fill(ttbarVsumLVec.eta() , event.weight_);
      hTTbarVsum_phi_->Fill(ttbarVsumLVec.phi() , event.weight_);
    }
  }

  hMet_pt_->Fill((event.met_->pt()), event.weight_);
  hMet_phi_->Fill((event.met_->phi()), event.weight_);

  if ( event.jets_->size() > 0 )
  {
    const LorentzVector& lv = event.jets_->at(0);
    const double bTag = event.jets_bTag_->at(0);

    hJet1_pt_  ->Fill((lv.pt() ), event.weight_);
    hJet1_eta_ ->Fill((lv.eta()), event.weight_);
    hJet1_phi_ ->Fill((lv.phi()), event.weight_);
    hJet1_btag_->Fill((bTag    ), event.weight_);
  }
  if ( event.jets_->size() > 1 )
  {
    const LorentzVector& lv = event.jets_->at(1);
    const double bTag = event.jets_bTag_->at(1);

    hJet2_pt_  ->Fill((lv.pt() ), event.weight_);
    hJet2_eta_ ->Fill((lv.eta()), event.weight_);
    hJet2_phi_ ->Fill((lv.phi()), event.weight_);
    hJet2_btag_->Fill((bTag    ), event.weight_);
  }
  if ( event.jets_->size() > 2 )
  {
    const LorentzVector& lv = event.jets_->at(2);
    const double bTag = event.jets_bTag_->at(2);

    hJet3_pt_  ->Fill((lv.pt() ), event.weight_);
    hJet3_eta_ ->Fill((lv.eta()), event.weight_);
    hJet3_phi_ ->Fill((lv.phi()), event.weight_);
    hJet3_btag_->Fill((bTag    ), event.weight_);
  }
  if ( event.jets_->size() > 3 )
  {
    const LorentzVector& lv = event.jets_->at(3);
    const double bTag = event.jets_bTag_->at(3);

    hJet4_pt_  ->Fill((lv.pt() ), event.weight_);
    hJet4_eta_ ->Fill((lv.eta()), event.weight_);
    hJet4_phi_ ->Fill((lv.phi()), event.weight_);
    hJet4_btag_->Fill((bTag    ), event.weight_);
  }

  int nTbjet = 0, nMbjet = 0, nLbjet = 0;
  for ( int i=0, n=event.jets_->size(); i<n; ++i )
  {
    //const LorentzVector& lv = event.jets_->at(i);
    const double bTag = event.jets_bTag_->at(i);

    if ( bTag > selector.cut_btagTight_  ) ++nTbjet;
    if ( bTag > selector.cut_btagMedium_ ) ++nMbjet;
    if ( bTag > selector.cut_btagLoose_  ) ++nLbjet;
  }
  hNTbjets_->Fill(nTbjet, event.weight_);
  hNMbjets_->Fill(nMbjet, event.weight_);
  hNLbjets_->Fill(nLbjet, event.weight_);

  if ( event.jets_->size() >= 2 )
  {
    const LorentzVector jjLVec = event.jets_->at(0)+event.jets_->at(1);
    hJJ_m_  ->Fill(jjLVec.mass(), event.weight_);
    hJJ_pt_ ->Fill(jjLVec.pt()  , event.weight_);
    hJJ_eta_->Fill(jjLVec.eta() , event.weight_);
    hJJ_phi_->Fill(jjLVec.phi() , event.weight_);

    int bidx1 = -1, bidx2 = -1;
    for ( int i=0, n=event.jets_->size(); i<n; ++i )
    {
      if ( event.jets_bTag_->at(i) < selector.cut_btagMedium_ ) continue;

      if ( bidx1 < 0 ) bidx1 = i;
      else if ( bidx2 < 0 ) bidx2 = i;
      else break;
    }
    if ( bidx1 >= 0 and bidx2 >= 0 )
    {
      const LorentzVector bbLVec = event.jets_->at(bidx1)+event.jets_->at(bidx2);
      hBB_m_  ->Fill(bbLVec.mass(), event.weight_);
      hBB_pt_ ->Fill(bbLVec.pt()  , event.weight_);
      hBB_eta_->Fill(bbLVec.eta() , event.weight_);
      hBB_phi_->Fill(bbLVec.phi() , event.weight_);
    }
  }

}

TTbarDileptonNtupleAnalyzer::Selector::Selector(Event* event, int mode)
{
  event_ = event;
  mode_ = mode;

  isGoodEvent_.resize(nCutSteps);
  for ( int i=0; i<nCutSteps; ++i ) isGoodEvent_[i] = false;

  int q1, q2;
  LorentzVector lep1, lep2;

  if ( mode == EE )
  {
    if ( event->electrons_->size() < 2 ) return;
    lep1 = event->electrons_->at(0);
    lep2 = event->electrons_->at(1);
    q1 = event->electrons_Q_->at(0);
    q2 = event->electrons_Q_->at(1);
  }
  else if ( mode == MM )
  {
    if ( event->muons_->size() < 2 ) return;
    lep1 = event->muons_->at(0);
    lep2 = event->muons_->at(1);
    q1 = event->muons_Q_->at(0);
    q2 = event->muons_Q_->at(1);
  }
  else if ( mode == ME )
  {
    if ( event->electrons_->size() < 1 or event->muons_->size() < 1 ) return;
    lep1 = event->muons_->at(0);
    lep2 = event->electrons_->at(0);
    q1 = event->muons_Q_->at(0);
    q2 = event->electrons_Q_->at(0);
  }
  else
  {
    cout << "Wrong channel mode, " << mode << endl;
    event_ = 0;
    return;
  }

  do
  {
    const double mLL = (lep1+lep2).mass();

    // Cut step 0 : low mass cut, opposite signed pair
    if ( mLL < 12 or q1+q2 != 0 ) break;
    isGoodEvent_[0] = true;

    // Z veto
    if ( mode != ME and abs(mLL-91.2) < 15 ) break;
    isGoodEvent_[1] = true;

    // Jet multiplicity
    if ( event->jets_->size() < 2 ) break;
    isGoodEvent_[2] = true;

    // Missing ET cut
    if ( mode != ME and event->met_->pt() < 30 ) break;
    else if ( event->met_->pt() < 20 ) break;
    isGoodEvent_[3] = true;

    int nBjets = 0;
    for ( int i=0, n=event->jets_bTag_->size(); i<n; ++i )
    {
      if ( event->jets_bTag_->at(i) > cut_btagMedium_ ) ++nBjets;
    }
    // B jet multiplicity >= 2
    if ( nBjets < 2 ) break;
    isGoodEvent_[4] = true;

    // Jet multiplicity >= 4 
    if ( event->jets_->size() < 4 ) break;
    isGoodEvent_[5] = true;

    if ( nBjets < 3 ) break;
    isGoodEvent_[6] = true;
  } while ( false );

}

void TTbarDileptonNtupleAnalyzer::analyze(int verboseLevel)
{
  if ( !event_ or !event_->isValid_ )
  {
    cout << "Invalid event tree\n";
    return;
  }

  for ( int iCutStep = 0, nCutStep = hists_.size(); iCutStep < nCutStep; ++iCutStep )
  {
    hists_[iCutStep]->hNEvent_->Fill(1, event_->nEvent_);
  }

  if ( verboseLevel >= 1 ) cout << "Beginning processing sample " << inputFile_->GetName() << endl;
  for ( int iEntry = 0, nEntry = event_->eventTree_->GetEntries(); iEntry<nEntry; ++iEntry )
  {
    if ( verboseLevel >= 2 ) printEntryFraction(iEntry, nEntry);

    event_->eventTree_->GetEntry(iEntry);
    event_->genTree_->GetEntry(iEntry);

    Selector selector(event_, mode_);
    for ( int iCutStep = 0, nCutStep = hists_.size(); iCutStep < nCutStep; ++iCutStep )
    {
      hists_[iCutStep]->fill(selector, iCutStep);
    }
  }
}

