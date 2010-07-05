#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"

#include "TH1F.h"

class MuonHLTAnalyzer : public edm::EDAnalyzer
{
public:
  explicit MuonHLTAnalyzer(const edm::ParameterSet& pset);
  ~MuonHLTAnalyzer();

  virtual void beginJob();
  virtual void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup);
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void endJob();

private:
  edm::ESHandle<MagneticField> theBField;
  edm::ESHandle<GlobalTrackingGeometry> theGeometry;
  edm::ESHandle<Propagator> thePropagator;
  PropagateToMuon theL1Matcher;

  edm::InputTag l1Label_, l2SeedLabel_, l3SeedLabel_, l2Label_, l3Label_;
  edm::InputTag muonLabel_;
  edm::InputTag trigLabel_;

  // Histograms
  typedef TH1F* TH1FP;
  typedef TH2F* TH2FP;

  TH1FP hGlbMu_N_, hGlbMu_Pt_, hGlbMu_Eta_, hGlbMu_Phi_, hGlbMu_Charge_;
  TH1FP hL1Mu_N_, hL1Mu_Pt_, hL1Mu_Eta_, hL1Mu_Phi_;
};
