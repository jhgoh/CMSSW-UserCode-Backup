#include "CosmicAnalysis/CosmicAnalyzer/interface/MuonRPCDigisAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"

#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile2D.h>

#include <memory>
#include <string>

using namespace std;

namespace H1
{
  enum Histograms
  {
    Strip, Bx, nDigi
  };
}

namespace HPrf2
{
  enum Histograms
  {
    Strip_XY_Rm2, Strip_XY_Rm1, Strip_XY_R00, Strip_XY_Rp1, Strip_XY_Rp2,
    Bx_XY_Rm2, Bx_XY_Rm1, Bx_XY_R00, Bx_XY_Rp1, Bx_XY_Rp2
  };
}

MuonRPCDigisAnalyzer::MuonRPCDigisAnalyzer(const edm::ParameterSet& pset)
{
  digiLabel_ = pset.getParameter<edm::InputTag>("digiLabel");

  edm::Service<TFileService> fs;

  h1_[H1::Strip] = fs->make<TH1F>("hStrip", "Strip profile", 100, 0, 100);
  h1_[H1::Bx] = fs->make<TH1F>("hBx", "Bunch crossing", 11, -5.5, 5.5);
  h1_[H1::nDigi] = fs->make<TH1F>("hNDigi", "Number of digi", 100, 0, 100);

  const double xMin = -1000, xMax = 1000, yMin = -1000, yMax = 1000;

  hPrf2_[HPrf2::Strip_XY_Rm2] = fs->make<TProfile2D>("hStrip_XY_Rm2", "Strip profile on XY plane, Ring -2", 100, xMin, xMax, 100, yMin, yMax);
  hPrf2_[HPrf2::Strip_XY_Rm1] = fs->make<TProfile2D>("hStrip_XY_Rm1", "Strip profile on XY plane, Ring -1", 100, xMin, xMax, 100, yMin, yMax);
  hPrf2_[HPrf2::Strip_XY_R00] = fs->make<TProfile2D>("hStrip_XY_R00", "Strip profile on XY plane, Ring  0", 100, xMin, xMax, 100, yMin, yMax);
  hPrf2_[HPrf2::Strip_XY_Rp1] = fs->make<TProfile2D>("hStrip_XY_Rp1", "Strip profile on XY plane, Ring -1", 100, xMin, xMax, 100, yMin, yMax);
  hPrf2_[HPrf2::Strip_XY_Rp2] = fs->make<TProfile2D>("hStrip_XY_Rp2", "Strip profile on XY plane, Ring -2", 100, xMin, xMax, 100, yMin, yMax);

  hPrf2_[HPrf2::Bx_XY_Rm2] = fs->make<TProfile2D>("hBx_XY_Rm2", "Bunch crossing profile on XY plane, Ring -2", 100, xMin, xMax, 100, yMin, yMax);
  hPrf2_[HPrf2::Bx_XY_Rm1] = fs->make<TProfile2D>("hBx_XY_Rm1", "Bunch crossing profile on XY plane, Ring -2", 100, xMin, xMax, 100, yMin, yMax);
  hPrf2_[HPrf2::Bx_XY_R00] = fs->make<TProfile2D>("hBx_XY_R00", "Bunch crossing profile on XY plane, Ring -2", 100, xMin, xMax, 100, yMin, yMax);
  hPrf2_[HPrf2::Bx_XY_Rp1] = fs->make<TProfile2D>("hBx_XY_Rp1", "Bunch crossing profile on XY plane, Ring -2", 100, xMin, xMax, 100, yMin, yMax);
  hPrf2_[HPrf2::Bx_XY_Rp2] = fs->make<TProfile2D>("hBx_XY_Rp2", "Bunch crossing profile on XY plane, Ring -2", 100, xMin, xMax, 100, yMin, yMax);
}

MuonRPCDigisAnalyzer::~MuonRPCDigisAnalyzer()
{
  
}

void MuonRPCDigisAnalyzer::beginJob(const edm::EventSetup& eventSetup)
{
}

void MuonRPCDigisAnalyzer::endJob()
{
}

void MuonRPCDigisAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  edm::Handle<RPCDigiCollection> rpcDigis;
  event.getByLabel(digiLabel_, rpcDigis);

  for ( RPCDigiCollection::DigiRangeIterator detUnitIt = rpcDigis->begin();
        detUnitIt != rpcDigis->end(); ++detUnitIt )
  {
    const RPCDetId Rsid = (*detUnitIt).first;
    const RPCDigiCollection::Range& range = (*detUnitIt).second;

    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(Rsid));

    if ( !roll ) continue;

    int nDigi = 0;
    for ( RPCDigiCollection::const_iterator digiIt = range.first;
          digiIt != range.second; ++digiIt )
    {
      nDigi++;

      const int strip = digiIt->strip();
      const int bx = digiIt->bx();

      h1_[H1::Strip]->Fill(strip);
      h1_[H1::Bx]->Fill(bx);

      const GlobalPoint gPoint = roll->toGlobal(roll->centreOfStrip(strip));
      const double globalX = gPoint.x();
      const double globalY = gPoint.y();
      
      if ( abs(Rsid.ring()) > 2 ) continue;

      hPrf2_[HPrf2::Strip_XY_R00+Rsid.ring()]->Fill(globalX, globalY, strip);
      hPrf2_[HPrf2::Bx_XY_R00+Rsid.ring()]->Fill(globalX, globalY, bx);
    }
    h1_[H1::nDigi]->Fill(nDigi);
  }
}

