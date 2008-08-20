#include "CosmicAnalysis/CosmicAnalyzer/src/MuonTrackAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/InputTag.h"

//#include "DataFormats/Common/interface/Ref.h"
//#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"

//#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
//#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
//#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

//
/*
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
*/
//

#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>

#include <memory>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

string analyzerName = "MuonTrackAnalyzer";

class HistogramGroup
{
  public:
    HistogramGroup(const string& name, Service<TFileService>& fs):
      dir_(fs->mkdir(name))
    {
      const int nBins = 50;
      const int nBinsEta = 100;
      const double minEta = -2.5, maxEta = 2.5;
      const double minPhi = -3.15, maxPhi = 3.15;

      hN_ = dir_.make<TH1F>("nTracks", "Number of tracks", 20, 0, 20);
      hAlgo_ = dir_.make<TH1F>("algo", "Tracking algorithm", 10, -1, 8);

      hQ_ = dir_.make<TH1F>("charge", "Track charge", 5, -2, 2);
      hEta_Phi_ = dir_.make<TH2F>("Eta_Phi", "#eta vs #phi", nBinsEta, minEta, maxEta, nBins, minPhi, maxPhi);

      hNValidHits_ = dir_.make<TH1F>("nHits", "Number of hits", 100, 0, 100); //FIXME : set # of bins correctly
      hChi2NDof_ = dir_.make<TH1F>("Chi2NDof", "Number of DoF", nBins, 0, 20);
      hChi2Norm_ = dir_.make<TH1F>("Chi2Norm", "Normalized #Chi^{2}", nBins, 0, 10);
      hChi2Prob_ = dir_.make<TH1F>("Chi2Prob", "#Chi^{2} probability", nBins, 0, 1);

      const int nBinsPos = 100;
      const double minX = -1000, maxX = 1000;
      const double minY = -1000, maxY = 1000;
      const double minZ = -750, maxZ = 750;

      hInnerX_Y_ = dir_.make<TH2F>("innerX_Y", "y vs x of inner position", nBinsPos, minX, maxX, nBinsPos, minY, maxY);
      hOuterX_Y_ = dir_.make<TH2F>("outerX_Y", "y vs x of outer position", nBinsPos, minX, maxX, nBinsPos, minY, maxY);

      hInnerZ_ = dir_.make<TH1F>("innerZ", "z of inner position", nBinsPos, minZ, maxZ);
      hOuterZ_ = dir_.make<TH1F>("outerZ", "z of outer position", nBinsPos, minZ, maxZ);

      hRecHitsX_Y_ = dir_.make<TH2F>("RecHitsX_Y", "y vs x of recHits", nBinsPos, minX, maxX, nBinsPos, minY, maxY);
      hRecHitsZ_ = dir_.make<TH1F>("RecHitsZ", "z of recHits", nBinsPos, minZ, maxZ);

      hRecHitsDTX_Y_ = dir_.make<TH2F>("RecHitsDTX_Y", "y vs x of DT recHits", nBinsPos, minX, maxX, nBinsPos, minY, maxY);
      hRecHitsDTZ_ = dir_.make<TH1F>("RecHitsDTZ", "z of DT recHits", nBinsPos, minZ, maxZ);

      hRecHitsCSCX_Y_ = dir_.make<TH2F>("RecHitsCSCX_Y", "y vs x of CSC recHits", nBinsPos, minX, maxX, nBinsPos, minY, maxY);
      hRecHitsCSCZ_ = dir_.make<TH1F>("RecHitsCSCZ", "z of CSC recHits", nBinsPos, minZ, maxZ);

      hRecHitsRPCX_Y_ = dir_.make<TH2F>("RecHitsRPCX_Y", "y vs x of RPC recHits", nBinsPos, minX, maxX, nBinsPos, minY, maxY);
      hRecHitsRPCZ_ = dir_.make<TH1F>("RecHitsRPCZ", "z of RPC recHits", nBinsPos, minZ, maxZ);

      const int maxClusterSize = 5;
      hRPCClusterSize_ = dir_.make<TH1F>("RPCClusterSize", "cluster size of recHits", maxClusterSize, 0, maxClusterSize);
      hRPCClusterSize_Eta_ = dir_.make<TH2F>("RPCClusterSize_Eta", "#eta vs cluster size of recHits", maxClusterSize, 0, maxClusterSize, nBinsEta, minEta, maxEta);
      hRPCClusterSize_Phi_ = dir_.make<TH2F>("RPCClusterSize_Phi", "#phi vs cluster size of recHits", maxClusterSize, 0, maxClusterSize, nBins, minPhi, maxPhi);
    };

    TFileDirectory dir_;

    typedef TH1F* TH1FP;
    typedef TH2F* TH2FP;

    TH1FP hN_;
    TH1FP hQ_;
    TH1FP hAlgo_;

    TH2FP hEta_Phi_;
    TH1FP hChi2NDof_, hChi2Norm_, hChi2Prob_;

    TH2FP hInnerX_Y_, hOuterX_Y_;
    TH1FP hInnerZ_, hOuterZ_;

    TH1FP hNValidHits_; 

    TH2FP hRecHitsX_Y_;
    TH1FP hRecHitsZ_;

    TH2FP hRecHitsDTX_Y_, hRecHitsRPCX_Y_, hRecHitsCSCX_Y_;
    TH1FP hRecHitsDTZ_, hRecHitsCSCZ_, hRecHitsRPCZ_;

    TH1FP hRPCClusterSize_;
    TH2FP hRPCClusterSize_Eta_, hRPCClusterSize_Phi_;
};

MuonTrackAnalyzer::MuonTrackAnalyzer(const ParameterSet& pset)
{
  trkLabel_ = pset.getUntrackedParameter<InputTag>("track");
  //rpcHitLabel_ = pset.getUntrackedParameter<InputTag>("rpcHitLabel");

  edm::Service<TFileService> fs;
  hTrk_ = new HistogramGroup("CosmicBarrel", fs);

  ParameterSet muonServParam = pset.getParameter<ParameterSet>("ServiceParameters");
  theMuonService = new MuonServiceProxy(muonServParam);

  ParameterSet trkDetAssocParams = pset.getParameter<ParameterSet>("TrackAssociatorParameters");
  trkDetAssocParams_.loadParameters(trkDetAssocParams);
  theTrkDetAssociator.useDefaultPropagator();

}

MuonTrackAnalyzer::~MuonTrackAnalyzer()
{
//  LogDebug("MuonTrackAnalyzer") << "Destroy MuonTrackAnalyzer" << endl;
  if ( hTrk_ ) delete hTrk_;
  if ( theMuonService ) delete theMuonService;
}

void MuonTrackAnalyzer::beginJob(const edm::EventSetup& eventSetup)
{
}

void MuonTrackAnalyzer::endJob()
{
}

void MuonTrackAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
//  LogDebug("MuonTrackAnalyzer") << "MuonTrackAnalyzer::analyze()" << endl;

  Handle<TrackCollection> trkColl;
  event.getByLabel(trkLabel_, trkColl);

  hTrk_->hN_->Fill(trkColl->size());

  Handle<RPCRecHitCollection> rpcHitColl;
  event.getByLabel(rpcHitLabel_, rpcHitColl);

  cout << rpcHitColl->size() << endl;

//  ESHandle<MagneticField> bField;
//  eventSetup.get<IdealMagneticFieldRecord>().get(bField);

  ESHandle<GlobalTrackingGeometry> trkGeometry;
  eventSetup.get<GlobalTrackingGeometryRecord>().get(trkGeometry);

  ESHandle<TransientTrackBuilder> transTrkBuilder;
  eventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrkBuilder);

  for(TrackCollection::const_iterator iTrk = trkColl->begin();
      iTrk != trkColl->end(); ++iTrk) 
  {
    hTrk_->hAlgo_->Fill(int(iTrk->algo()));

    hTrk_->hQ_->Fill(iTrk->charge());
    hTrk_->hEta_Phi_->Fill(iTrk->eta(), iTrk->phi());

    hTrk_->hNValidHits_->Fill(iTrk->numberOfValidHits());
    hTrk_->hChi2NDof_->Fill(iTrk->ndof());
    hTrk_->hChi2Norm_->Fill(iTrk->normalizedChi2());
    hTrk_->hChi2Prob_->Fill(TMath::Prob(iTrk->chi2(), int(iTrk->ndof())));

    const math::XYZPoint innerPosition = iTrk->innerPosition();
    const math::XYZPoint outerPosition = iTrk->outerPosition();

    hTrk_->hInnerX_Y_->Fill(innerPosition.X(), innerPosition.Y());
    hTrk_->hOuterX_Y_->Fill(outerPosition.X(), outerPosition.Y());

    hTrk_->hInnerZ_->Fill(outerPosition.Z());
    hTrk_->hOuterZ_->Fill(outerPosition.Z());

    //for(TrackingRecHitRefVector::const_iterator iHit = iTrk->recHitsBegin();
    for(trackingRecHit_iterator iHit = iTrk->recHitsBegin();
        iHit != iTrk->recHitsEnd(); ++iHit) 
    {
      if ( !((*iHit)->isValid()) ) continue;

      const TrackingRecHit* hit = iHit->get();

      const DetId& detId = hit->geographicalId();
      const GlobalPoint& point = trkGeometry->idToDet(detId)->surface().toGlobal(hit->localPosition());

      hTrk_->hRecHitsX_Y_->Fill(point.x(), point.y());
      hTrk_->hRecHitsZ_->Fill(point.z());

      if ( detId.det() == DetId::Muon ) 
      {
        if ( detId.subdetId() == 1 ) // DT hits
        {
          hTrk_->hRecHitsDTX_Y_->Fill(point.x(), point.y());
          hTrk_->hRecHitsDTZ_->Fill(point.z());
        }
        else if ( detId.subdetId() == 2 ) // CSC hits
        {
          hTrk_->hRecHitsCSCX_Y_->Fill(point.x(), point.y());
          hTrk_->hRecHitsCSCZ_->Fill(point.z());
        }
        else if ( detId.subdetId() == 3 ) // RPC hits
        {
          hTrk_->hRecHitsRPCX_Y_->Fill(point.x(), point.y());
          hTrk_->hRecHitsRPCZ_->Fill(point.z());

/*
          const RPCRecHit* rpcRecHit = dynamic_cast<const RPCRecHit*>(hit);
          if ( rpcRecHit ) {
            const int clusterSize = rpcRecHit->clusterSize();
            hTrk_->hRPCClusterSize_->Fill(clusterSize);
            hTrk_->hRPCClusterSize_Eta_->Fill(clusterSize, iTrk->eta());
            hTrk_->hRPCClusterSize_Phi_->Fill(clusterSize, iTrk->phi());
          }
          else {
            LogDebug(analyzerName) << "Hit came from RPC, but cannot convert TrackingRecHit to RPCRecHit.\n";
          }
*/
        }
      }
    }

/*
    const reco::TransientTrack transTrk = transTrkBuilder->build(*iTrk);
    const FreeTrajectoryState fts = transTrk.initialFreeState();
    const GlobalTrajectoryParameters trajParam = fts.parameters();
*/
/*
    TrackDetMatchInfo trkDetMatchInfo = theTrkDetAssociator.associate(event, eventSetup, fts, trkDetAssocParams_);

    for(vector<TAMuonChamberMatch>::const_iterator iChamber = trkDetMatchInfo.chambers.begin();
        iChamber != trkDetMatchInfo.chambers.end(); ++iChamber) {
      const GeomDet* detector = trkGeometry->idToDet(iChamber->id);

      const DTLayer* dtLayer = dynamic_cast<const DTLayer*>(detector);
      const RPCRoll* rpcRoll = dynamic_cast<const RPCRoll*>(detector);
    }
*/
  }
}
