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

    NSimHit_Wheel, NSimHit_Disk,
    NRecHit_Wheel, NRecHit_Disk,

    NLostHit_Wheel, NLostHit_Disk,
    NNoisyHit_Wheel, NNoisyHit_Disk,

    NSimHitMatch_Wheel, NSimHitMatch_Disk,
    NRecHitMatch_Wheel, NRecHitMatch_Disk,

    NSimHitRZ, NRecHitRZ,
    NSimHitMatchRZ, NRecHitMatchRZ,

    NSimHitXY_WM2, NSimHitXY_WM1, NSimHitXY_W00, NSimHitXY_WP1, NSimHitXY_WP2,
    NRecHitXY_WM2, NRecHitXY_WM1, NRecHitXY_W00, NRecHitXY_WP1, NRecHitXY_WP2,
    NSimHitXY_DM3, NSimHitXY_DM2, NSimHitXY_DM1, NSimHitXY_DP1, NSimHitXY_DP2, NSimHitXY_DP3,
    NRecHitXY_DM3, NRecHitXY_DM2, NRecHitXY_DM1, NRecHitXY_DP1, NRecHitXY_DP2, NRecHitXY_DP3,

    NSimHitMatchXY_WM2, NSimHitMatchXY_WM1, NSimHitMatchXY_W00, NSimHitMatchXY_WP1, NSimHitMatchXY_WP2,
    NRecHitMatchXY_WM2, NRecHitMatchXY_WM1, NRecHitMatchXY_W00, NRecHitMatchXY_WP1, NRecHitMatchXY_WP2,
    NSimHitMatchXY_DM3, NSimHitMatchXY_DM2, NSimHitMatchXY_DM1, NSimHitMatchXY_DP1, NSimHitMatchXY_DP2, NSimHitMatchXY_DP3,
    NRecHitMatchXY_DM3, NRecHitMatchXY_DM2, NRecHitMatchXY_DM1, NRecHitMatchXY_DP1, NRecHitMatchXY_DP2, NRecHitMatchXY_DP3,

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
  h_[HName::ClusterSize] = dbe_->book1D("ClusterSize", "Cluster size;Cluster size", 21, -10.5, 10.5);

  h_[HName::Res] = dbe_->book1D("Res", "Global residuals;Residual [cm]", 100, -8, 8);
  h_[HName::Pull] = dbe_->book1D("Pull", "Global pulls;Pull", 100, -5, 5);

  h_[HName::NSimHit_Wheel] = dbe_->book1D("NSimHit_Wheel", "Number of SimHits;Wheel", 5, -2.5, 2.5);
  h_[HName::NSimHit_Disk] = dbe_->book1D("NSimHit_Disk", "Number of SimHits;Disk", 7, -3.5, 3.5);

  h_[HName::NRecHit_Wheel] = dbe_->book1D("NRecHit_Wheel", "Number of RecHits;Wheel", 5, -2.5, 2.5);
  h_[HName::NRecHit_Disk] = dbe_->book1D("NRecHit_Disk", "Number of RecHits;Disk", 7, -3.5, 3.5);

  h_[HName::NLostHit_Wheel] = dbe_->book1D("NLostHit_Wheel", "Number of lost hits;Wheel", 5, -2.5, 2.5);
  h_[HName::NLostHit_Disk] = dbe_->book1D("NLostHit_Disk", "Number of lost hits;Disk", 7, -3.5, 3.5);

  h_[HName::NNoisyHit_Wheel] = dbe_->book1D("NNoisyHit_Wheel", "Number of noisy hits;Wheel", 5, -2.5, 2.5);
  h_[HName::NNoisyHit_Disk] = dbe_->book1D("NNoisyHit_Disk", "Number of noisy hits;Disk", 7, -3.5, 3.5);

  h_[HName::NSimHitMatch_Wheel] = dbe_->book1D("NSimHitMatch_Wheel", "Number of SimHitMatchs;Wheel", 5, -2.5, 2.5);
  h_[HName::NSimHitMatch_Disk] = dbe_->book1D("NSimHitMatch_Disk", "Number of SimHitMatchs;Disk", 7, -3.5, 3.5);

  h_[HName::NRecHitMatch_Wheel] = dbe_->book1D("NRecHitMatch_Wheel", "Number of RecHitMatchs;Wheel", 5, -2.5, 2.5);
  h_[HName::NRecHitMatch_Disk] = dbe_->book1D("NRecHitMatch_Disk", "Number of RecHitMatchs;Disk", 7, -3.5, 3.5);

  // XY overview
  if ( isStandAloneMode_ )
  {
    const int nBin = 1000;
    const double xmin = -1000, xmax = 1000;
    const double ymin = -1000, ymax = 1000;

    h_[HName::NSimHitRZ] = dbe_->book2D("NSimHitRZ", "Number of SimHits;Z;R", nBin, -1100, 1100, nBin, 0, xmax);
    h_[HName::NRecHitRZ] = dbe_->book2D("NRecHitRZ", "Number of RecHits;Z;R", nBin, -1100, 1100, nBin, 0, xmax);

    h_[HName::NSimHitMatchRZ] = dbe_->book2D("NSimHitMatchRZ", "Number of Matched SimHits;Z;R", nBin, -1100, 1100, nBin, 0, xmax);
    h_[HName::NRecHitMatchRZ] = dbe_->book2D("NRecHitMatchRZ", "Number of Matched RecHits;Z;R", nBin, -1100, 1100, nBin, 0, xmax);

    h_[HName::NSimHitXY_WM2] = dbe_->book2D("NSimHitXY_WM2", "Number of SimHits WM2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitXY_WM1] = dbe_->book2D("NSimHitXY_WM1", "Number of SimHits WM1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitXY_W00] = dbe_->book2D("NSimHitXY_W00", "Number of SimHits W00;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitXY_WP1] = dbe_->book2D("NSimHitXY_WP1", "Number of SimHits WP1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitXY_WP2] = dbe_->book2D("NSimHitXY_WP2", "Number of SimHits WP2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
                                                                                   
    h_[HName::NSimHitXY_DM3] = dbe_->book2D("NSimHitXY_DM3", "Number of SimHits DM3;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitXY_DM2] = dbe_->book2D("NSimHitXY_DM2", "Number of SimHits DM2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitXY_DM1] = dbe_->book2D("NSimHitXY_DM1", "Number of SimHits DM1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitXY_DP1] = dbe_->book2D("NSimHitXY_DP1", "Number of SimHits DP1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitXY_DP2] = dbe_->book2D("NSimHitXY_DP2", "Number of SimHits DP2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitXY_DP3] = dbe_->book2D("NSimHitXY_DP3", "Number of SimHits DP3;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
                                                                                   
    h_[HName::NRecHitXY_WM2] = dbe_->book2D("NRecHitXY_WM2", "Number of RecHits WM2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitXY_WM1] = dbe_->book2D("NRecHitXY_WM1", "Number of RecHits WM1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitXY_W00] = dbe_->book2D("NRecHitXY_W00", "Number of RecHits W00;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitXY_WP1] = dbe_->book2D("NRecHitXY_WP1", "Number of RecHits WP1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitXY_WP2] = dbe_->book2D("NRecHitXY_WP2", "Number of RecHits WP2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
                                                                                   
    h_[HName::NRecHitXY_DM3] = dbe_->book2D("NRecHitXY_DM3", "Number of RecHits DM3;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitXY_DM2] = dbe_->book2D("NRecHitXY_DM2", "Number of RecHits DM2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitXY_DM1] = dbe_->book2D("NRecHitXY_DM1", "Number of RecHits DM1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitXY_DP1] = dbe_->book2D("NRecHitXY_DP1", "Number of RecHits DP1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitXY_DP2] = dbe_->book2D("NRecHitXY_DP2", "Number of RecHits DP2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitXY_DP3] = dbe_->book2D("NRecHitXY_DP3", "Number of RecHits DP3;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);

    h_[HName::NSimHitMatchXY_WM2] = dbe_->book2D("NSimHitMatchXY_WM2", "Number of Matched SimHits WM2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitMatchXY_WM1] = dbe_->book2D("NSimHitMatchXY_WM1", "Number of Matched SimHits WM1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitMatchXY_W00] = dbe_->book2D("NSimHitMatchXY_W00", "Number of Matched SimHits W00;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitMatchXY_WP1] = dbe_->book2D("NSimHitMatchXY_WP1", "Number of Matched SimHits WP1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitMatchXY_WP2] = dbe_->book2D("NSimHitMatchXY_WP2", "Number of Matched SimHits WP2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
                                                                                   
    h_[HName::NSimHitMatchXY_DM3] = dbe_->book2D("NSimHitMatchXY_DM3", "Number of Matched SimHits DM3;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitMatchXY_DM2] = dbe_->book2D("NSimHitMatchXY_DM2", "Number of Matched SimHits DM2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitMatchXY_DM1] = dbe_->book2D("NSimHitMatchXY_DM1", "Number of Matched SimHits DM1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitMatchXY_DP1] = dbe_->book2D("NSimHitMatchXY_DP1", "Number of Matched SimHits DP1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitMatchXY_DP2] = dbe_->book2D("NSimHitMatchXY_DP2", "Number of Matched SimHits DP2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NSimHitMatchXY_DP3] = dbe_->book2D("NSimHitMatchXY_DP3", "Number of Matched SimHits DP3;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
                                                                                   
    h_[HName::NRecHitMatchXY_WM2] = dbe_->book2D("NRecHitMatchXY_WM2", "Number of Matched RecHits WM2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitMatchXY_WM1] = dbe_->book2D("NRecHitMatchXY_WM1", "Number of Matched RecHits WM1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitMatchXY_W00] = dbe_->book2D("NRecHitMatchXY_W00", "Number of Matched RecHits W00;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitMatchXY_WP1] = dbe_->book2D("NRecHitMatchXY_WP1", "Number of Matched RecHits WP1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitMatchXY_WP2] = dbe_->book2D("NRecHitMatchXY_WP2", "Number of Matched RecHits WP2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
                                                                                   
    h_[HName::NRecHitMatchXY_DM3] = dbe_->book2D("NRecHitMatchXY_DM3", "Number of Matched RecHits DM3;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitMatchXY_DM2] = dbe_->book2D("NRecHitMatchXY_DM2", "Number of Matched RecHits DM2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitMatchXY_DM1] = dbe_->book2D("NRecHitMatchXY_DM1", "Number of Matched RecHits DM1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitMatchXY_DP1] = dbe_->book2D("NRecHitMatchXY_DP1", "Number of Matched RecHits DP1;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitMatchXY_DP2] = dbe_->book2D("NRecHitMatchXY_DP2", "Number of Matched RecHits DP2;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);
    h_[HName::NRecHitMatchXY_DP3] = dbe_->book2D("NRecHitMatchXY_DP3", "Number of Matched RecHits DP3;X;Y", nBin, xmin, xmax, nBin, ymin, ymax);

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
  if ( dbe_ )
  {
    if ( !rootFileName_.empty() ) dbe_->save(rootFileName_);
  }
}

void RPCRecHitValid::beginJob()
{
}

void RPCRecHitValid::endJob()
{
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

  // Loop over simHits, fill histograms which does not need associations
  for ( SimHitIter simHitIter = simHitHandle->begin();
        simHitIter != simHitHandle->end(); ++simHitIter )
  {
    const RPCDetId detId = static_cast<const RPCDetId>(simHitIter->detUnitId());
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(detId()));
    if ( !roll ) continue;

    const int region = roll->id().region();
    const int ring = roll->id().ring();
    //const int sector = roll->id().sector();
    const int station = abs(roll->id().station());
    //const int layer = roll->id().layer();
    //const int subSector = roll->id().subsector();

    if ( region == 0 ) h_[HName::NSimHit_Wheel]->Fill(ring);
    else h_[HName::NSimHit_Disk]->Fill(region*station);

    if ( isStandAloneMode_ )
    {
      const GlobalPoint pos = roll->toGlobal(simHitIter->localPosition());

      h_[HName::NSimHitRZ]->Fill(pos.z(), pos.perp());
      if ( region == 0 )
      {
        h_[HName::NSimHitXY_W00+ring]->Fill(pos.x(), pos.y());
      }
      else if ( region == -1 and station < 4 )
      {
        h_[HName::NSimHitXY_DM1-(station-1)]->Fill(pos.x(), pos.y());
      }
      else if ( region == 1 and station < 4 )
      {
        h_[HName::NSimHitXY_DP1+(station-1)]->Fill(pos.x(), pos.y());
      }
    }
  }

  // Loop over recHits, fill histograms which does not need associations
  for ( RecHitIter recHitIter = recHitHandle->begin();
        recHitIter != recHitHandle->end(); ++recHitIter )
  {
    const RPCDetId detId = static_cast<const RPCDetId>(recHitIter->rpcId());
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(detId()));
    if ( !roll ) continue;

    const int region = roll->id().region();
    const int ring = roll->id().ring();
    //const int sector = roll->id().sector();
    const int station = abs(roll->id().station());
    //const int layer = roll->id().layer();
    //const int subSector = roll->id().subsector();

    h_[HName::ClusterSize]->Fill(recHitIter->clusterSize());

    if ( region == 0 ) h_[HName::NRecHit_Wheel]->Fill(ring);
    else h_[HName::NRecHit_Disk]->Fill(region*station);

    if ( isStandAloneMode_ )
    {
      const GlobalPoint pos = roll->toGlobal(recHitIter->localPosition());

      h_[HName::NRecHitRZ]->Fill(pos.z(), pos.perp());
      if ( region == 0 )
      {
        h_[HName::NRecHitXY_W00+ring]->Fill(pos.x(), pos.y());
      }
      else if ( region == -1 and station < 4 )
      {
        h_[HName::NRecHitXY_DM1-(station-1)]->Fill(pos.x(), pos.y());
      }
      else if ( region == 1 and station < 4 )
      {
        h_[HName::NRecHitXY_DP1+(station-1)]->Fill(pos.x(), pos.y());
      }
    }
  }

  // Start matching SimHits to RecHits
  // Trick used : To remove duplicated matching to the same recHit, 
  //               we use RecHit->SimHit map 

  //typedef std::map<SimHitIter, RecHitIter> SimToRecHitMap;
  typedef std::map<RecHitIter, SimHitIter> RecToSimHitMap;
  //SimToRecHitMap simToRecHitMap;
  RecToSimHitMap recToSimHitMap;

  for ( SimHitIter simHitIter = simHitHandle->begin();
        simHitIter != simHitHandle->end(); ++simHitIter )
  {
    const RPCDetId simDetId = static_cast<const RPCDetId>(simHitIter->detUnitId());
    const RPCRoll* simRoll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(simDetId));
    if ( !simRoll ) continue;

    const double simX = simHitIter->localPosition().x();

    for ( RecHitIter recHitIter = recHitHandle->begin();
          recHitIter != recHitHandle->end(); ++recHitIter )
    {
      const RPCDetId recDetId = static_cast<const RPCDetId>(recHitIter->rpcId());
      const RPCRoll* recRoll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(recDetId));
      if ( !recRoll ) continue;

      if ( simDetId != recDetId ) continue;

      const double recX = recHitIter->localPosition().x();
      const double newDx = fabs(recX - simX);

      // Associate SimHit to RecHit
      RecToSimHitMap::const_iterator prevMatch = recToSimHitMap.find(recHitIter);
      if ( prevMatch == recToSimHitMap.end() )
      {
        recToSimHitMap.insert(std::make_pair(recHitIter, simHitIter));
      }
      else
      {
        const double prevX = prevMatch->first->localPosition().x();
        const double prevDx = fabs(prevX - simX);
        
        if ( newDx < prevDx )
        {
          recToSimHitMap[recHitIter] = simHitIter;
        }
      }
    }
  }

  // Now we have simHit-recHit mapping
  // So we can fill up relavant histograms
  for ( RecToSimHitMap::const_iterator match = recToSimHitMap.begin();
        match != recToSimHitMap.end(); ++match )
  {
    RecHitIter recHitIter = match->first;
    SimHitIter simHitIter = match->second;

    const RPCDetId detId = static_cast<const RPCDetId>(simHitIter->detUnitId());
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(detId));

    const int region = roll->id().region();
    const int ring = roll->id().ring();
    //const int sector = roll->id().sector();
    const int station = abs(roll->id().station());
    //const int layer = roll->id().layer();
    //const int subsector = roll->id().subsector();

    const double simX = simHitIter->localPosition().x();
    const double recX = recHitIter->localPosition().x();
    const double errX = recHitIter->localPositionError().xx();
    const double dX = recX - simX;
    const double pull = errX == 0 ? -999 : dX/errX;
    
    h_[HName::Res]->Fill(dX);
    h_[HName::Pull]->Fill(pull);

    if ( isStandAloneMode_ )
    {
      const GlobalPoint simPos = roll->toGlobal(simHitIter->localPosition());
      const GlobalPoint recPos = roll->toGlobal(recHitIter->localPosition());

      h_[HName::NSimHitMatchRZ]->Fill(simPos.z(), simPos.perp());
      h_[HName::NRecHitMatchRZ]->Fill(recPos.z(), recPos.perp());

      if ( region == 0 )
      {
        h_[HName::NSimHitMatchXY_W00+ring]->Fill(simPos.x(), simPos.y());
        h_[HName::NRecHitMatchXY_W00+ring]->Fill(recPos.x(), recPos.y());
      }
      else if ( region == -1 and station < 4 )
      {
        h_[HName::NSimHitMatchXY_DM1-(station-1)]->Fill(simPos.x(), simPos.y());
        h_[HName::NRecHitMatchXY_DM1-(station-1)]->Fill(recPos.x(), recPos.y());
      }
      else if ( region == 1 and station < 4 )
      {
        h_[HName::NSimHitMatchXY_DP1+(station-1)]->Fill(simPos.x(), simPos.y());
        h_[HName::NRecHitMatchXY_DP1+(station-1)]->Fill(recPos.x(), recPos.y());
      }
    }

    if ( region == 0 )
    {
      h_[HName::NRecHitMatch_Wheel]->Fill(ring);
      h_[HName::Res_W00+ring]->Fill(dX);
      h_[HName::Pull_W00+ring]->Fill(pull);
    }
    else if ( region == -1 and station < 4 )
    {
      h_[HName::NRecHitMatch_Disk]->Fill(region*station);
      h_[HName::Res_DM1-(station-1)]->Fill(dX);
      h_[HName::Pull_DM1-(station-1)]->Fill(pull);
    }
    else if ( region == 1 and station < 4 )
    {
      h_[HName::NRecHitMatch_Disk]->Fill(region*station);
      h_[HName::Res_DP1+(station-1)]->Fill(dX);
      h_[HName::Pull_DP1+(station-1)]->Fill(pull);
    }
  }

  // Find Lost hits
  for ( SimHitIter simHitIter = simHitHandle->begin();
        simHitIter != simHitHandle->end(); ++simHitIter )
  {
    const RPCDetId detId = static_cast<const RPCDetId>(simHitIter->detUnitId());
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(detId));

    const int region = roll->id().region();
    const int ring = roll->id().ring();
    //const int sector = roll->id().sector();
    const int station = abs(roll->id().station());
    //const int layer = roll->id().layer();
    //const int subsector = roll->id().subsector();

    bool matched = false;
    for ( RecToSimHitMap::const_iterator match = recToSimHitMap.begin();
          match != recToSimHitMap.end(); ++match )
    {
      if ( simHitIter == match->second )
      {
        matched = true;
        break;
      }
    }

    if ( !matched )
    {
      if ( region == 0 ) h_[HName::NLostHit_Wheel]->Fill(ring);
      else h_[HName::NLostHit_Disk]->Fill(region*station);
    }
  }

  // Find Noisy hits
  for ( RecHitIter recHitIter = recHitHandle->begin();
        recHitIter != recHitHandle->end(); ++recHitIter )
  {
    const RPCDetId detId = static_cast<const RPCDetId>(recHitIter->rpcId());
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(detId));

    const int region = roll->id().region();
    const int ring = roll->id().ring();
    //const int sector = roll->id().sector();
    const int station = abs(roll->id().station());
    //const int layer = roll->id().layer();
    //const int subsector = roll->id().subsector();

    bool matched = false;
    for ( RecToSimHitMap::const_iterator match = recToSimHitMap.begin();
          match != recToSimHitMap.end(); ++match )
    {
      if ( recHitIter == match->first )
      {
        matched = true;
        break;
      }
    }

    if ( !matched )
    {
      if ( region == 0 ) h_[HName::NNoisyHit_Wheel]->Fill(ring);
      else h_[HName::NNoisyHit_Disk]->Fill(region*station);
    }
  }
}

