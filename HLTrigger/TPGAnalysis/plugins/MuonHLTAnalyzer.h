#ifndef HLTrigger_TPGAnalysis_MuonHLTAnalyzer_H
#define HLTrigger_TPGAnalysis_MuonHLTAnalyzer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MuonAnalysis/MuonAssociators/interface/L1MuonMatcherAlgo.h"

#include <TH1F.h>
#include <TH2F.h>

#include <vector>

namespace reco
{
  class Muon;
}

class TrajectoryStateOnSurface;
class MagneticField;
class Propagator;
class GlobalTrackingGeometry;
class Histograms;

class MuonHLTAnalyzer : public edm::EDAnalyzer
{
public:
  MuonHLTAnalyzer(const edm::ParameterSet& pset);
  ~MuonHLTAnalyzer();

  void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup);
  void endRun();
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  typedef std::vector<std::string> VString;
  typedef TH1F* TH1FP;

  edm::InputTag l1MuonTag_, triggerEventTag_, recoMuonTag_;
  VString muonL1TNames_;
  
  L1MuonMatcherAlgo l1Matcher_;
  edm::ESHandle<MagneticField> bField_;
  edm::ESHandle<GlobalTrackingGeometry> geometry_;
  edm::ESHandle<Propagator> propagator_;

  // Cut variables
  double minPt_, maxRelIso_;

  // List of histograms
  TH1F* hNEvent_;

  std::vector<Histograms*> histograms_;

  // Run by run histograms
  std::map<int, TH1F*> hNEvent_ByRun_;
  std::map<int, std::vector<Histograms*> > histograms_ByRun_;
};

#endif

