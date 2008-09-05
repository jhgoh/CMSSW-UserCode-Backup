#include "CosmicAnalysis/CosmicAnalyzer/src/MuonTrackAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>

#include <memory>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

string analyzerName = "MuonTrackAnalyzer";

typedef TrajectoryStateOnSurface TSOS;

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
      hAlgo_ = dir_.make<TH1F>("algo", "Tracking algorithm", 8, 0, 8);

      hQ_ = dir_.make<TH1F>("charge", "Track charge", 5, -2, 3);
      hEta_Phi_ = dir_.make<TH2F>("Eta_Phi", "#eta vs #phi", nBinsEta, minEta, maxEta, nBins, minPhi, maxPhi);

      hNValidHits_ = dir_.make<TH1F>("nHits", "Number of hits", 60, 0, 60);
      hChi2NDof_ = dir_.make<TH1F>("Chi2NDof", "Number of DoF", 60, 0, 60);
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

      hAllHits_X_Y_ = dir_.make<TH2F>("AllHitsX_Y", "y vs x of recHits", nBinsPos, minX, maxX, nBinsPos, minY, maxY);
      hAllHits_Z_ = dir_.make<TH1F>("AllHitsZ", "z of recHits", nBinsPos, minZ, maxZ);

      hDTHits_X_Y_ = dir_.make<TH2F>("DTHits_X_Y", "y vs x of DT recHits", nBinsPos, minX, maxX, nBinsPos, minY, maxY);
      hDTHits_Z_ = dir_.make<TH1F>("DTHits_Z", "z of DT recHits", nBinsPos, minZ, maxZ);

      hCSCHits_X_Y_ = dir_.make<TH2F>("CSCHits_X_Y", "y vs x of CSC recHits", nBinsPos, minX, maxX, nBinsPos, minY, maxY);
      hCSCHits_Z_ = dir_.make<TH1F>("CSCHits_Z", "z of CSC recHits", nBinsPos, minZ, maxZ);

      hRPCHits_X_Y_ = dir_.make<TH2F>("RPCHits_X_Y", "y vs x of RPC recHits", nBinsPos, minX, maxX, nBinsPos, minY, maxY);
      hRPCHits_Z_ = dir_.make<TH1F>("RPCHits_Z", "z of RPC recHits", nBinsPos, minZ, maxZ);

      hMatchedRPCHits_X_Y_ = dir_.make<TH2F>("MatchedRPCHits_X_Y", "y vs x of matched RPC recHits", nBinsPos, minX, maxX, nBinsPos, minY, maxY);
      hMatchedRPCHits_Z_ = dir_.make<TH1F>("MatchedRPCHits_Z", "z of matched RPC recHits", nBinsPos, minZ, maxZ);

      const int maxClusterSize = 10;
      hRPCClusterSize_ = dir_.make<TH1F>("RPCClusterSize", "cluster size of recHits", maxClusterSize, 0, maxClusterSize);
      hRPC_Eta_ClusterSize_ = dir_.make<TH2F>("RPC_Eta_ClusterSize", "#eta vs cluster size of recHits", nBinsEta, minEta, maxEta, maxClusterSize, 0, maxClusterSize);
      hRPC_Phi_ClusterSize_ = dir_.make<TH2F>("RPC_Phi_ClusterSize", "#phi vs cluster size of recHits",  nBins, minPhi, maxPhi, maxClusterSize, 0, maxClusterSize);

      hMatchedRPC_Eta_ClusterSize_ = dir_.make<TH2F>("MatchedRPC_Eta_ClusterSize", "matched cluster size of recHits vs eta", nBinsEta, minEta, maxEta, maxClusterSize, 0, maxClusterSize);
      hMatchedRPC_Phi_ClusterSize_ = dir_.make<TH2F>("MatchedRPC_Phi_ClusterSize", "matched cluster size of recHits vs phi", nBins, minPhi, maxPhi, maxClusterSize, 0, maxClusterSize);

      hMatchedRPC_Eta_LocalAngle_ = dir_.make<TH2F>("MatchedRPC_Eta_LocalAngle", "matched RPC recHits' eta vs local angle", nBinsEta, minEta, maxEta, nBins, -1, 1);
      hMatchedRPC_Phi_LocalAngle_ = dir_.make<TH2F>("MatchedRPC_Phi_LocalAngle", "matched RPC recHits' phi vs local angle", nBins, minPhi, maxPhi, nBins, -1, 1);

      hMatchedRPC_CosLocalAngle_ClusterSize_ = dir_.make<TH2F>("MatchedRPC_CosLocalAngle_ClusterSize", "matched cluster size of rechits vs cos(local angle)", nBins, 0, 1, maxClusterSize, 0, maxClusterSize);
      hMatchedRPC_CosLocalAngleX_ClusterSize_ = dir_.make<TH2F>("MatchedRPC_CosLocalAngleX_ClusterSize", "matched cluster size of rechits vsx/(x^2+z^2)", nBins, -1, 1, maxClusterSize, 0, maxClusterSize);
      hMatchedRPC_CosLocalAngleY_ClusterSize_ = dir_.make<TH2F>("MatchedRPC_CosLocalAngleY_ClusterSize", "matched cluster size of rechits vsy/(y^2+z^2)", nBins, -1, 1, maxClusterSize, 0, maxClusterSize);

      hMatchedRPC_NStrips_ = dir_.make<TH1F>("MatchedRPC_NStrip", "Reminder of # of strip", nBins, 0, 1);
      hMatchedRPC_NStrips_ClusterSize_ = dir_.make<TH2F>("MatchedRPC_NStrips_ClusterSize", "matched cluster size of rechits vs Reminder of # of strip", nBins, 0, 1, maxClusterSize, 0, maxClusterSize);

      hMatchedRPC_PredictedStrip_ = dir_.make<TH1F>("MatchedRPC_PredictedStrip", "Predicted strip", nBins, 0, 1);
      hMatchedRPC_PredictedStrip_ClusterSize_ = dir_.make<TH2F>("MatchedRPC_PredictedStrip_ClusterSize", "matched cluster size of rechits vs Reminder of predicted strip", nBins, 0, 1, maxClusterSize, 0, maxClusterSize);

      hN_->GetXaxis()->SetTitle("Number of tracks");
      hAlgo_->GetXaxis()->SetTitle("Tracking algorithm #");

      hQ_->GetXaxis()->SetTitle("Track charge");
      hEta_Phi_->GetXaxis()->SetTitle("#eta");
      hEta_Phi_->GetYaxis()->SetTitle("#phi");

      hNValidHits_->GetXaxis()->SetTitle("Number of valid hits");
      hChi2NDof_->GetXaxis()->SetTitle("Number of D.o.F");
      hChi2Norm_->GetXaxis()->SetTitle("Normalized #Chi^{2}");
      hChi2Prob_->GetXaxis()->SetTitle("#Chi^{2} probability");

      hInnerX_Y_->GetXaxis()->SetTitle("X position (cm)"); 
      hInnerX_Y_->GetYaxis()->SetTitle("Y position (cm)"); 
      hInnerZ_->GetXaxis()->SetTitle("Z position (cm)");

      hOuterX_Y_->GetXaxis()->SetTitle("X position (cm)"); 
      hOuterX_Y_->GetYaxis()->SetTitle("Y position (cm)"); 
      hOuterZ_->GetXaxis()->SetTitle("Z position (cm)");

      hAllHits_X_Y_->GetXaxis()->SetTitle("X position (cm)");
      hAllHits_X_Y_->GetYaxis()->SetTitle("Y position (cm)");
      hAllHits_Z_->GetXaxis()->SetTitle("Z position (cm)");

      hDTHits_X_Y_->GetXaxis()->SetTitle("X position (cm)");
      hDTHits_X_Y_->GetYaxis()->SetTitle("Y position (cm)");
      hDTHits_Z_->GetXaxis()->SetTitle("Z position (cm)");

      hCSCHits_X_Y_->GetXaxis()->SetTitle("X position (cm)");
      hCSCHits_X_Y_->GetYaxis()->SetTitle("Y position (cm)");
      hCSCHits_Z_->GetXaxis()->SetTitle("Z position (cm)");

      hRPCHits_X_Y_->GetXaxis()->SetTitle("X position (cm)");
      hRPCHits_X_Y_->GetYaxis()->SetTitle("Y position (cm)");
      hRPCHits_Z_->GetXaxis()->SetTitle("Z position (cm)");

      hRPC_Phi_ClusterSize_->GetXaxis()->SetTitle("#phi");
      hRPC_Phi_ClusterSize_->GetYaxis()->SetTitle("Cluster size");

      hMatchedRPCHits_X_Y_->GetXaxis()->SetTitle("X position (cm)");
      hMatchedRPCHits_X_Y_->GetYaxis()->SetTitle("Y position (cm)");
      hMatchedRPCHits_Z_->GetXaxis()->SetTitle("Z position (cm)");

      hMatchedRPC_Eta_ClusterSize_->GetXaxis()->SetTitle("#eta");
      hMatchedRPC_Eta_ClusterSize_->GetYaxis()->SetTitle("Cluster size");

      hMatchedRPC_Phi_ClusterSize_->GetXaxis()->SetTitle("#phi");
      hMatchedRPC_Phi_ClusterSize_->GetYaxis()->SetTitle("Cluster size");

      hMatchedRPC_CosLocalAngle_ClusterSize_->GetXaxis()->SetTitle("cosine of local direction");
      hMatchedRPC_CosLocalAngle_ClusterSize_->GetYaxis()->SetTitle("Cluster size");  

      hMatchedRPC_CosLocalAngleX_ClusterSize_->GetXaxis()->SetTitle("cosine of local direction");
      hMatchedRPC_CosLocalAngleX_ClusterSize_->GetYaxis()->SetTitle("Cluster size");  

      hMatchedRPC_CosLocalAngleY_ClusterSize_->GetXaxis()->SetTitle("cosine of local direction");
      hMatchedRPC_CosLocalAngleY_ClusterSize_->GetYaxis()->SetTitle("Cluster size");  

      hMatchedRPC_Eta_LocalAngle_->GetXaxis()->SetTitle("#eta");
      hMatchedRPC_Eta_LocalAngle_->GetYaxis()->SetTitle("cosine of local angle");

      hMatchedRPC_Phi_LocalAngle_->GetXaxis()->SetTitle("#phi");
      hMatchedRPC_Phi_LocalAngle_->GetYaxis()->SetTitle("cosine of local angle");

      hMatchedRPC_NStrips_->GetXaxis()->SetTitle("nStrip");
      hMatchedRPC_NStrips_ClusterSize_->GetXaxis()->SetTitle("nStrip");
      hMatchedRPC_NStrips_ClusterSize_->GetYaxis()->SetTitle("cluster size");

      hMatchedRPC_PredictedStrip_->GetXaxis()->SetTitle("Reminder of predicted strip number");
      hMatchedRPC_PredictedStrip_ClusterSize_->GetXaxis()->SetTitle("Reminder of predicted strip number");
      hMatchedRPC_PredictedStrip_ClusterSize_->GetYaxis()->SetTitle("cluster size");
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

    TH2FP hAllHits_X_Y_;
    TH1FP hAllHits_Z_;

    TH2FP hDTHits_X_Y_, hRPCHits_X_Y_, hCSCHits_X_Y_;
    TH1FP hDTHits_Z_, hCSCHits_Z_, hRPCHits_Z_;

    TH2FP hMatchedRPCHits_X_Y_;
    TH1FP hMatchedRPCHits_Z_;

    TH1FP hRPCClusterSize_;
    TH2FP hRPC_Eta_ClusterSize_, hRPC_Phi_ClusterSize_;

    TH2FP hMatchedRPC_Eta_ClusterSize_, hMatchedRPC_Phi_ClusterSize_;

    TH2FP hMatchedRPC_CosLocalAngle_ClusterSize_;
    TH2FP hMatchedRPC_CosLocalAngleX_ClusterSize_;
    TH2FP hMatchedRPC_CosLocalAngleY_ClusterSize_;
    TH2FP hMatchedRPC_Eta_LocalAngle_, hMatchedRPC_Phi_LocalAngle_;

    TH1FP hMatchedRPC_NStrips_;
    TH2FP hMatchedRPC_NStrips_ClusterSize_;

    TH1FP hMatchedRPC_PredictedStrip_;
    TH2FP hMatchedRPC_PredictedStrip_ClusterSize_;
};

MuonTrackAnalyzer::MuonTrackAnalyzer(const ParameterSet& pset)
{
  trkLabel_ = pset.getUntrackedParameter<InputTag>("track");
//  cscHitLabel_ = pset.getUntrackedParameter<InputTag>("CSCHits");
//  rpcHitLabel_ = pset.getUntrackedParameter<InputTag>("RPCHits");

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

//  Handle<CSCRecHit2DCollection> cscHitColl;
//  event.getByLabel(cscHitLabel_, cscHitColl);

//  Handle<RPCRecHitCollection> rpcHitColl;
//  event.getByLabel(rpcHitLabel_, rpcHitColl);

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

    for(TrackingRecHitRefVector::const_iterator iHit = iTrk->recHitsBegin();
        iHit != iTrk->recHitsEnd(); ++iHit) 
    {
      if ( !((*iHit)->isValid()) ) continue;

      const TrackingRecHit* hit = iHit->get();

      const DetId& detId = hit->geographicalId();
      const GeomDet* detector = trkGeometry->idToDet(detId);
      const GlobalPoint& point = detector->surface().toGlobal(hit->localPosition());

      const double trkX = point.x(), trkY = point.y();
      const double trkZ = point.z();

      hTrk_->hAllHits_X_Y_->Fill(trkX, trkY);
      hTrk_->hAllHits_Z_->Fill(trkZ);

      const reco::TransientTrack transTrk = transTrkBuilder->build(*iTrk);
      const FreeTrajectoryState fts = transTrk.initialFreeState();

      if ( detId.det() == DetId::Muon ) 
      {
        if ( detId.subdetId() == 1 ) // DT hits
        {
          hTrk_->hDTHits_X_Y_->Fill(trkX, trkY);
          hTrk_->hDTHits_Z_->Fill(trkZ);
        }
        else if ( detId.subdetId() == 2 ) // CSC hits
        {
          hTrk_->hCSCHits_X_Y_->Fill(trkX, trkY);
          hTrk_->hCSCHits_Z_->Fill(trkZ);
        }
        else if ( detId.subdetId() == 3 ) // RPC hits
        {
          // Find RPC hits from RPCRecHitCollection
          const RPCRecHit* rpcHit = dynamic_cast<const RPCRecHit*>(hit);
          if ( !rpcHit ) continue;

          hTrk_->hRPCHits_X_Y_->Fill(trkX, trkY);
          hTrk_->hRPCHits_Z_->Fill(trkZ);

          const int clusterSize = rpcHit->clusterSize();
          hTrk_->hRPCClusterSize_->Fill(clusterSize);
          hTrk_->hRPC_Eta_ClusterSize_->Fill(iTrk->eta(), clusterSize);
          hTrk_->hRPC_Phi_ClusterSize_->Fill(iTrk->phi(), clusterSize);

          hTrk_->hMatchedRPCHits_X_Y_->Fill(trkX, trkY);
          hTrk_->hMatchedRPCHits_Z_->Fill(trkZ);

          hTrk_->hMatchedRPC_Eta_ClusterSize_->Fill(iTrk->eta(), clusterSize);
          hTrk_->hMatchedRPC_Phi_ClusterSize_->Fill(iTrk->phi(), clusterSize);

          TSOS tsosAtDet(fts, detector->surface());
          LocalVector trjDir = tsosAtDet.localDirection();

          const double lx = trjDir.x();
          const double ly = trjDir.y();
          const double lz = trjDir.z();
          const double lr = hypot(lx, ly);

          hTrk_->hMatchedRPC_CosLocalAngle_ClusterSize_->Fill(lr/hypot(lr, lz), clusterSize);
          hTrk_->hMatchedRPC_CosLocalAngleX_ClusterSize_->Fill(lx/hypot(lx, lz), clusterSize);
          hTrk_->hMatchedRPC_CosLocalAngleY_ClusterSize_->Fill(ly/hypot(ly, lz), clusterSize);

          const RPCRoll* rpcRoll = dynamic_cast<const RPCRoll*>(detector);
          if ( rpcRoll ) 
          {
            const float predictedStrip = rpcRoll->strip(tsosAtDet.localPosition());
            const float nStrips = rpcRoll->nstrips();

            const float remPredictedStrip = predictedStrip - floor(predictedStrip);
            const float remNStrips = nStrips - floor(nStrips);

            hTrk_->hMatchedRPC_NStrips_->Fill(remNStrips);
            hTrk_->hMatchedRPC_NStrips_ClusterSize_->Fill(remNStrips, clusterSize);

            hTrk_->hMatchedRPC_PredictedStrip_->Fill(remPredictedStrip);
            hTrk_->hMatchedRPC_PredictedStrip_ClusterSize_->Fill(remPredictedStrip, clusterSize);
          }
        }
      }
    }
  }
}

/* vim:set ts=2 sts=2 sw=2 expandtab: */
