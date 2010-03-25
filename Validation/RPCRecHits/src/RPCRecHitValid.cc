#include "Validation/RPCRecHits/interface/RPCRecHitValid.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

using namespace std;

struct HName
{
  enum
  {
    ClusterSize, Res, Pull,
    NSimRoll_Wheel, NSimRoll_Disk,
    NRecRoll_Wheel, NRecRoll_Disk,
    NLostRoll_Wheel, NLostRoll_Disk,
    NNoisyRoll_Wheel, NNoisyRoll_Disk,
    NSimHitXY_Barrel, NSimHitXY_Endcap, NSimHitRZ,
    NRecHitXY_Barrel, NRecHitXY_Endcap, NRecHitRZ,
    Res_WM2, Res_WM1, Res_W00, Res_WP1, Res_WP2,
    Res_DM3, Res_DM2, Res_DM1, Res_DP1, Res_DP2, Res_DP3,
    Pull_WM2, Pull_WM1, Pull_W00, Pull_WP1, Pull_WP2,
    Pull_DM3, Pull_DM2, Pull_DM1, Pull_DP1, Pull_DP2, Pull_DP3,
    END
  };
};

RPCRecHitValid::RPCRecHitValid(const edm::ParameterSet& pset)
{
  rootFileName_ = pset.getUntrackedParameter<string>("rootFileName", "");
  simHitLabel_ = pset.getParameter<edm::InputTag>("simHit");
  recHitLabel_ = pset.getParameter<edm::InputTag>("recHit");

  isStandAloneMode_ = pset.getUntrackedParameter<bool>("standAloneMode", false);

  dbe_ = edm::Service<DQMStore>().operator->();
  if ( !dbe_ ) 
  {
    edm::LogError("RPCRecHitValid") << "No DQMStore instance\n";
    return;
  }

  // Book MonitorElements
  string subDir("RPCRecHitsV");
  dbe_->setCurrentFolder(subDir);

  // Global plots
  h_[HName::ClusterSize] = dbe_->book1D("ClusterSize", "Cluster size;Cluster size", 11, -10.5, 10.5);

  h_[HName::Res] = dbe_->book1D("Res", "Global residuals;Residual [cm]", 100, -8, 8);
  h_[HName::Pull] = dbe_->book1D("Pull", "Global pulls;Pull", 100, -5, 5);

  h_[HName::NSimRoll_Wheel] = dbe_->book1D("NSimRoll_Wheel", "Number of rolls containing SimHits;Wheel", 5, -2.5, 2.5);
  h_[HName::NSimRoll_Disk] = dbe_->book1D("NSimRoll_Disk", "Number of rolls containing SimHits;Disk", 7, -3.5, 3.5);

  h_[HName::NRecRoll_Wheel] = dbe_->book1D("NRecRoll_Wheel", "Number of rolls containing RecHits;Wheel", 5, -2.5, 2.5);
  h_[HName::NRecRoll_Disk] = dbe_->book1D("NRecRoll_Disk", "Number of rolls containing RecHits;Disk", 7, -3.5, 3.5);

  h_[HName::NLostRoll_Wheel] = dbe_->book1D("NLostRoll_Wheel", "Number of rolls with hits lost;Wheel", 5, -2.5, 2.5);
  h_[HName::NLostRoll_Disk] = dbe_->book1D("NLostRoll_Disk", "Number of rolls with hits lost;Disk", 7, -3.5, 3.5);

  h_[HName::NNoisyRoll_Wheel] = dbe_->book1D("NNoisyRoll_Wheel", "Number of rolls with hits noisy;Wheel", 5, -2.5, 2.5);
  h_[HName::NNoisyRoll_Disk] = dbe_->book1D("NNoisyRoll_Disk", "Number of rolls with hits noisy;Disk", 7, -3.5, 3.5);

  // XY overview
  if ( isStandAloneMode_ )
  {
    h_[HName::NSimHitXY_Barrel] = dbe_->book2D("NSimHitXY_Barrel", "Number of SimHits in Barrel, global XY view;X;Y", 1000, -700, 700, 1000, -700, 700);
    h_[HName::NSimHitXY_Endcap] = dbe_->book2D("NSimHitXY_Endcap", "Number of SimHits in Endcap, global XY view;X;Y", 1000, -700, 700, 1000, -700, 700);
    h_[HName::NSimHitRZ] = dbe_->book2D("NSimHitRZ", "Number of SimHits in global RZ view;R;Z", 1000, -1100, 1100, 1000, 0, 700);

    h_[HName::NRecHitXY_Barrel] = dbe_->book2D("NRecHitXY_Barrel", "Number of RecHits in global XY view;X;Y", 1000, -700, 700, 1000, -700, 700);
    h_[HName::NRecHitXY_Endcap] = dbe_->book2D("NRecHitXY_Endcap", "Number of RecHits in global XY view;X;Y", 1000, -700, 700, 1000, -700, 700);
    h_[HName::NRecHitRZ] = dbe_->book2D("NRecHitRZ", "Number of RecHits in global RZ view;R;Z", 1000, -1100, 1100, 1000, 0, 700);
  }

  // Residuals and pulls
  h_[HName::Res_WM2] = dbe_->book1D("Res_WM2", "Residuals for Wheel -2;Residual [cm]", 100, -8, 8);
  h_[HName::Res_WM1] = dbe_->book1D("Res_WM1", "Residuals for Wheel -1;Residual [cm]", 100, -8, 8);
  h_[HName::Res_W00] = dbe_->book1D("Res_W00", "Residuals for Wheel  0;Residual [cm]", 100, -8, 8);
  h_[HName::Res_WP1] = dbe_->book1D("Res_WP1", "Residuals for Wheel +1;Residual [cm]", 100, -8, 8);
  h_[HName::Res_WP2] = dbe_->book1D("Res_WP2", "Residuals for Wheel +2;Residual [cm]", 100, -8, 8);

  h_[HName::Res_DM3] = dbe_->book1D("Res_DM3", "Residuals for Disk -3;Residual [cm]", 100, -8, 8);
  h_[HName::Res_DM2] = dbe_->book1D("Res_DM2", "Residuals for Disk -2;Residual [cm]", 100, -8, 8);
  h_[HName::Res_DM1] = dbe_->book1D("Res_DM1", "Residuals for Disk -1;Residual [cm]", 100, -8, 8);
  h_[HName::Res_DP1] = dbe_->book1D("Res_DP1", "Residuals for Disk +1;Residual [cm]", 100, -8, 8);
  h_[HName::Res_DP2] = dbe_->book1D("Res_DP2", "Residuals for Disk +2;Residual [cm]", 100, -8, 8);
  h_[HName::Res_DP3] = dbe_->book1D("Res_DP3", "Residuals for Disk +3;Residual [cm]", 100, -8, 8);

  h_[HName::Pull_WM2] = dbe_->book1D("Pull_WM2", "Pull for Wheel -2;Pull", 100, -5, 5);
  h_[HName::Pull_WM1] = dbe_->book1D("Pull_WM1", "Pull for Wheel -1;Pull", 100, -5, 5);
  h_[HName::Pull_W00] = dbe_->book1D("Pull_W00", "Pull for Wheel  0;Pull", 100, -5, 5);
  h_[HName::Pull_WP1] = dbe_->book1D("Pull_WP1", "Pull for Wheel +1;Pull", 100, -5, 5);
  h_[HName::Pull_WP2] = dbe_->book1D("Pull_WP2", "Pull for Wheel +2;Pull", 100, -5, 5);

  h_[HName::Pull_DM3] = dbe_->book1D("Pull_DM3", "Pull for Disk -3;Pull", 100, -5, 5);
  h_[HName::Pull_DM2] = dbe_->book1D("Pull_DM2", "Pull for Disk -2;Pull", 100, -5, 5);
  h_[HName::Pull_DM1] = dbe_->book1D("Pull_DM1", "Pull for Disk -1;Pull", 100, -5, 5);
  h_[HName::Pull_DP1] = dbe_->book1D("Pull_DP1", "Pull for Disk +1;Pull", 100, -5, 5);
  h_[HName::Pull_DP2] = dbe_->book1D("Pull_DP2", "Pull for Disk +2;Pull", 100, -5, 5);
  h_[HName::Pull_DP3] = dbe_->book1D("Pull_DP3", "Pull for Disk +3;Pull", 100, -5, 5);
}

RPCRecHitValid::~RPCRecHitValid()
{
}

void RPCRecHitValid::beginJob()
{
}

void RPCRecHitValid::endJob()
{
  if ( dbe_ )
  {
    if ( !rootFileName_.empty() ) dbe_->save(rootFileName_);
  }
}

void RPCRecHitValid::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  if ( !dbe_ ) 
  {
    edm::LogError("RPCRecHitValid") << "No DQMStore instance\n";
    return;
  }

  // Get the RPC Geometry
  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  // Retrieve SimHits from the event
  edm::Handle<edm::PSimHitContainer> simHitHandle;
  event.getByLabel(simHitLabel_, simHitHandle);

  // Retrieve RecHits from the event
  edm::Handle<RPCRecHitCollection> recHitHandle;
  event.getByLabel(recHitLabel_, recHitHandle);

  typedef edm::PSimHitContainer::const_iterator SimHitIter;
  typedef RPCRecHitCollection::const_iterator RecHitIter;

  typedef std::map<RPCDetId, vector<SimHitIter> > RPCDetToSimHitMap;
  typedef std::map<RPCDetId, vector<RecHitIter> > RPCDetToRecHitMap;
  RPCDetToSimHitMap detToSimHits;
  RPCDetToRecHitMap detToRecHits;

  for ( SimHitIter simHitIter = simHitHandle->begin();
        simHitIter != simHitHandle->end(); ++simHitIter )
  {
    //if ( abs(simHit.particleType()) != 13 ) continue;

    const RPCDetId rpcDetId = static_cast<const RPCDetId>(simHitIter->detUnitId());
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(rpcDetId));
    if ( !roll ) continue;

    detToSimHits[rpcDetId].push_back(simHitIter);

    if ( isStandAloneMode_ )
    {
      const GlobalPoint pos = roll->toGlobal(simHitIter->localPosition());
      h_[HName::NSimHitXY_Endcap]->Fill(pos.x(), pos.y());
      h_[HName::NSimHitRZ]->Fill(pos.z(), pos.perp());
    }
  }

  for ( RecHitIter recHitIter = recHitHandle->begin();
        recHitIter != recHitHandle->end(); ++recHitIter )
  {
    const RPCDetId rpcDetId = static_cast<const RPCDetId>(recHitIter->rpcId());
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(rpcDetId));
    if ( !roll ) continue;

    detToRecHits[rpcDetId].push_back(recHitIter);

    if ( isStandAloneMode_ )
    {
      const GlobalPoint pos = roll->toGlobal(recHitIter->localPosition());
      h_[HName::NRecHitXY_Endcap]->Fill(pos.x(), pos.y());
      h_[HName::NRecHitRZ]->Fill(pos.z(), pos.perp());
    }

    const int clusterSize = recHitIter->clusterSize();
    h_[HName::ClusterSize]->Fill(clusterSize);

    //const int firstClusterStrip = recHitIter->firstClusterStrip();
  }

  // Sim to Reco association
  for ( RPCDetToSimHitMap::const_iterator detToSimHitsIter = detToSimHits.begin();
        detToSimHitsIter != detToSimHits.end(); ++detToSimHitsIter )
  {
    const RPCDetId simRPCDetId = detToSimHitsIter->first;
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(simRPCDetId));

    const int region = roll->id().region();
    const int ring = roll->id().ring();
    //const int sector = roll->id().sector();
    const int station = abs(roll->id().station());
    //const int layer = roll->id().layer();
    //const int subsector = roll->id().subsector();

    if ( region == 0 ) h_[HName::NSimRoll_Wheel]->Fill(ring);
    else h_[HName::NSimRoll_Disk]->Fill(region*station);

    RPCDetToRecHitMap::const_iterator detToRecHitsIter = detToRecHits.find(simRPCDetId);
    if ( detToRecHitsIter == detToRecHits.end() )
    {
      // No Sim to Rec matching : could be dead chambers, "LostRolls"
      if ( region == 0 ) h_[HName::NLostRoll_Wheel]->Fill(ring);
      else h_[HName::NLostRoll_Disk]->Fill(region*station);
    }
    else
    {
      // Sim to Rec matching exists, try to associate simHit and recHit
      std::vector<SimHitIter> simHitIters = detToSimHitsIter->second;
      std::vector<RecHitIter> recHitIters = detToRecHitsIter->second;

      const size_t nSimHits = simHitIters.size();
      const size_t nRecHits = recHitIters.size();

      // Start sim-rec matching
      typedef map<SimHitIter, RecHitIter> SimToRecHitMap;
      SimToRecHitMap simToRecHitMap;

      for ( size_t i=0; i<nSimHits; ++i )
      {
        SimHitIter simHitIter = simHitIters[i];
        const double simHitX = simHitIter->localPosition().x();

        for ( size_t j=0; j<nRecHits; ++j )
        {
          RecHitIter recHitIter = recHitIters[j];
          const double recHitX = recHitIter->localPosition().x();
          const double newDx = fabs(recHitX - simHitX);

          SimToRecHitMap::iterator prevMatch = simToRecHitMap.find(simHitIter);
          if ( prevMatch == simToRecHitMap.end() )
          {
            simToRecHitMap.insert(std::make_pair(simHitIter, recHitIter));
          }
          else
          {
            const double oldRecHitX = prevMatch->second->localPosition().x();
            const double oldDx = fabs(oldRecHitX - simHitX);

            if ( newDx < oldDx ) prevMatch->second = recHitIter;
          }
        }
      }

      for ( SimToRecHitMap::const_iterator simToRecHitIter = simToRecHitMap.begin();
            simToRecHitIter != simToRecHitMap.end(); ++simToRecHitIter )
      {
        const SimHitIter simHit = simToRecHitIter->first;
        const RecHitIter recHit = simToRecHitIter->second;

        const double simHitX = simHit->localPosition().x();
        const double recHitX = recHit->localPosition().x();
        const double errX = recHit->localPositionError().xx();
        const double dX = recHitX - simHitX;
        const double pull = dX/errX;

        h_[HName::Res]->Fill(dX);
        h_[HName::Pull]->Fill(pull);

        if ( region == 0 )
        {
          h_[HName::Res_W00+ring]->Fill(dX);
          h_[HName::Pull_W00+ring]->Fill(pull);
        }
        else if ( region == -1 and station < 4 )
        {
          h_[HName::Res_DM1-(station-1)]->Fill(dX);
          h_[HName::Pull_DM1-(station-1)]->Fill(pull);
        }
        else if ( region == 1 and station < 4 )
        {
          h_[HName::Res_DP1+(station-1)]->Fill(dX);
          h_[HName::Pull_DP1+(station-1)]->Fill(pull);
        }
      }
    }
  }

  // Reco to Sim association
  for ( RPCDetToRecHitMap::const_iterator detToRecHitsIter = detToRecHits.begin();
        detToRecHitsIter != detToRecHits.end(); ++detToRecHitsIter )
  {
    const RPCDetId recRPCDetId = detToRecHitsIter->first;
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(recRPCDetId));

    const int region = roll->id().region();
    const int ring = roll->id().ring();
    //const int sector = roll->id().sector();
    const int station = abs(roll->id().station());
    //const int layer = roll->id().layer();
    //const int subSector = roll->id().subsector();

    if ( region == 0 ) h_[HName::NRecRoll_Wheel]->Fill(ring);
    else h_[HName::NRecRoll_Disk]->Fill(region*station);

    RPCDetToSimHitMap::const_iterator detToSimHitsIter = detToSimHits.find(recRPCDetId);
    if ( detToSimHitsIter == detToSimHits.end() )
    {
      // Noisy rolls
      if ( region == 0 ) h_[HName::NNoisyRoll_Wheel]->Fill(ring);
      else h_[HName::NNoisyRoll_Disk]->Fill(region*station);
    }
    // Matched case already analyzed in Sim to Reco association loop
  }

}

