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

MuonRPCDigisAnalyzer::MuonRPCDigisAnalyzer(const edm::ParameterSet& pset)
{
  digiLabel_ = pset.getParameter<edm::InputTag>("digiLabel");

  edm::Service<TFileService> fs;

  h1_[H1::Strip] = fs->make<TH1F>("hStrip", "Strip profile", 100, 0, 100);
  h1_[H1::Bx] = fs->make<TH1F>("hBx", "Bunch crossing", 10, -5.5, -5.5+10);
  h1_[H1::nDigi] = fs->make<TH1F>("hNDigi", "Number of digi", 100, 0, 100);
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
      h1_[H1::Strip]->Fill(digiIt->strip());
      h1_[H1::Bx]->Fill(digiIt->bx());
      nDigi++;
    }
    h1_[H1::nDigi]->Fill(nDigi);
  }
}

