#include "CosmicAnalysis/CosmicAnalyzer/interface/MuonTimingAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"

#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "CondTools/RPC/interface/RPCRunIOV.h"

#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TString.h>
#include <TDatime.h>

#include <memory>

using namespace std;

MuonTimingAnalyzer::MuonTimingAnalyzer(const edm::ParameterSet& pset)
{
  digisLabel_ = pset.getParameter<edm::InputTag>("digisLabel");
  recHitsLabel_ = pset.getParameter<edm::InputTag>("recHitsLabel");

  // Book profile histograms for Time vs Conditions
  prf_["BxVsNDigi"] = fs_->make<TProfile>("prfBxVsNDigi", "Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0,3600);
  prf_["BxVsNRecHits"] = fs_->make<TProfile>("prfBxVsNRecHits", "Bx number vs Number of RecHits;Bx number;Number of RecHits", 3600, 0, 3600);

  h2_["BxVsNDigi"] = fs_->make<TH2F>("h2BxVsNDigi", "Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 250, 0, 250);
  h2_["BxVsNRecHits"] = fs_->make<TH2F>("h2BxVsNRecHits", "Bx number vs Number of RecHits;Bx number;Number of RecHits", 3600, 0, 3600, 250, 0, 250);

  minBxNumber_ = 0;
  maxBxNumber_ = 0;
}

MuonTimingAnalyzer::~MuonTimingAnalyzer()
{
  
}

void MuonTimingAnalyzer::beginJob(const edm::EventSetup& eventSetup)
{
}

void MuonTimingAnalyzer::endJob()
{
  cout << "MuonTimingAnalyzer::endJob() called" << endl
       << " Minimum bx number = " << minBxNumber_ << endl
       << " Maximum bx number = " << maxBxNumber_ << endl;
}

void MuonTimingAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  const int bxNumber = event.bunchCrossing();
  const int orbitNumber = event.orbitNumber();
  //const int storeNumber = event.storeNumber();

  if ( bxNumber > maxBxNumber_ ) maxBxNumber_ = bxNumber;
  if ( bxNumber < minBxNumber_ ) minBxNumber_ = bxNumber;

  edm::Handle<RPCDigiCollection> digisHandle;
  event.getByLabel(digisLabel_, digisHandle);

  edm::Handle<RPCRecHitCollection> recHitsHandle;
  event.getByLabel(recHitsLabel_, recHitsHandle);

  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  unsigned int nDigis = 0;
  for ( RPCDigiCollection::DigiRangeIterator detUnitIt = digisHandle->begin();
        detUnitIt != digisHandle->end(); ++detUnitIt )
  {
    const RPCDetId Rsid = (*detUnitIt).first;
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(Rsid));
    if ( !roll ) continue;
    const RPCDigiCollection::Range& range = (*detUnitIt).second;

    for ( RPCDigiCollection::const_iterator digiIt = range.first; 
          digiIt != range.second; ++digiIt )
    {
      nDigis++;
    }
  }

  unsigned int nRecHits = 0;
  for ( RPCRecHitCollection::const_iterator recHitsIter = recHitsHandle->begin();
        recHitsIter != recHitsHandle->end(); ++recHitsIter )
  {
    nRecHits++;
  }
  
  prf_["BxVsNDigi"]->Fill(bxNumber, nDigis);
  prf_["BxVsNRecHits"]->Fill(bxNumber, nRecHits);

  h2_["BxVsNDigi"]->Fill(bxNumber, nDigis);
  h2_["BxVsNRecHits"]->Fill(bxNumber, nRecHits);
}

