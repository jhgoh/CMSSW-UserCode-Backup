#ifndef CosmicAnalysis_CosmicAnalyzer_MuonTrackAnalyzer_H
#define CosmicAnalysis_CosmicAnalyzer_MuonTrackAnalyzer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"


namespace edm {
  class InputTag;
}

class MuonServiceProxy;

class HistogramGroup;

class MuonTrackAnalyzer : public edm::EDAnalyzer
{
  public:
    MuonTrackAnalyzer(const edm::ParameterSet& pset);
    ~MuonTrackAnalyzer();

    virtual void beginJob(const edm::EventSetup& eventSetup);
    virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);
    virtual void endJob();

  private:
    edm::InputTag trkLabel_;
    edm::InputTag rpcHitLabel_;

    HistogramGroup* hTrk_;

    MuonServiceProxy* theMuonService;
    edm::ESHandle<Propagator> thePropagator;

    TrackDetectorAssociator theTrkDetAssociator;
    TrackAssociatorParameters trkDetAssocParams_;
};

#endif
