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
  digiLabel_ = pset.getParameter<edm::InputTag>("digiLabel");
  recHitLabel_ = pset.getParameter<edm::InputTag>("recHitLabel");

  // Book profile histograms for Time vs Conditions
  prf_["BxVsNDigi"] = fs_->make<TProfile>("prfBxVsNDigi", "Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0,3600);
  prf_["BxVsNRecHit"] = fs_->make<TProfile>("prfBxVsNRecHit", "Bx number vs Number of RecHit;Bx number;Number of RecHit", 3600, 0, 3600);

  h1_["NDigi"] = fs_->make<TH1F>("hNDigi", "Number of Digis;Number of Digis", 250, 0, 250);
  h1_["NRecHit"] = fs_->make<TH1F>("hNRecHit", "Number of RecHit;Number of RecHit", 250, 0, 250);

  h1_["BxNumber"] = fs_->make<TH1F>("hBxNumber", "Bx number;Bx number;Number of Events", 3600, 0, 3600);

  h1_["DigiBx_RE+1"] = fs_->make<TH1F>("DigiBx_RE+1", "Bx number from all Digis in the Endcap+1;Bx number", 3600, 0, 3600);
  h1_["DigiBx_RE+2"] = fs_->make<TH1F>("DigiBx_RE+2", "Bx number from all Digis in the Endcap+2;Bx number", 3600, 0, 3600);
  h1_["DigiBx_RE+3"] = fs_->make<TH1F>("DigiBx_RE+3", "Bx number from all Digis in the Endcap+3;Bx number", 3600, 0, 3600);
  h1_["DigiBx_RE-1"] = fs_->make<TH1F>("DigiBx_RE-1", "Bx number from all Digis in the Endcap-1;Bx number", 3600, 0, 3600);
  h1_["DigiBx_RE-2"] = fs_->make<TH1F>("DigiBx_RE-2", "Bx number from all Digis in the Endcap-2;Bx number", 3600, 0, 3600);
  h1_["DigiBx_RE-3"] = fs_->make<TH1F>("DigiBx_RE-3", "Bx number from all Digis in the Endcap-3;Bx number", 3600, 0, 3600);

  h1_["DigiRelBx_RE+1"] = fs_->make<TH1F>("DigiRelBx_RE+1", "RelBx number from all Digis in the Endcap+1;RelBx number", 10, -5, -5);
  h1_["DigiRelBx_RE+2"] = fs_->make<TH1F>("DigiRelBx_RE+2", "RelBx number from all Digis in the Endcap+2;RelBx number", 10, -5, -5);
  h1_["DigiRelBx_RE+3"] = fs_->make<TH1F>("DigiRelBx_RE+3", "RelBx number from all Digis in the Endcap+3;RelBx number", 10, -5, -5);
  h1_["DigiRelBx_RE-1"] = fs_->make<TH1F>("DigiRelBx_RE-1", "RelBx number from all Digis in the Endcap-1;RelBx number", 10, -5, -5);
  h1_["DigiRelBx_RE-2"] = fs_->make<TH1F>("DigiRelBx_RE-2", "RelBx number from all Digis in the Endcap-2;RelBx number", 10, -5, -5);
  h1_["DigiRelBx_RE-3"] = fs_->make<TH1F>("DigiRelBx_RE-3", "RelBx number from all Digis in the Endcap-3;RelBx number", 10, -5, -5);

  h1_["RecHitBx_RE+1"] = fs_->make<TH1F>("RecHitBx_RE+1", "Bx number from all RecHit in the Endcap+1;Bx number", 3600, 0, 3600);
  h1_["RecHitBx_RE+2"] = fs_->make<TH1F>("RecHitBx_RE+2", "Bx number from all RecHit in the Endcap+2;Bx number", 3600, 0, 3600);
  h1_["RecHitBx_RE+3"] = fs_->make<TH1F>("RecHitBx_RE+3", "Bx number from all RecHit in the Endcap+3;Bx number", 3600, 0, 3600);
  h1_["RecHitBx_RE-1"] = fs_->make<TH1F>("RecHitBx_RE-1", "Bx number from all RecHit in the Endcap-1;Bx number", 3600, 0, 3600);
  h1_["RecHitBx_RE-2"] = fs_->make<TH1F>("RecHitBx_RE-2", "Bx number from all RecHit in the Endcap-2;Bx number", 3600, 0, 3600);
  h1_["RecHitBx_RE-3"] = fs_->make<TH1F>("RecHitBx_RE-3", "Bx number from all RecHit in the Endcap-3;Bx number", 3600, 0, 3600);

  h1_["RecHitRelBx_RE+1"] = fs_->make<TH1F>("RecHitRelBx_RE+1", "RelBx number from all RecHit in the Endcap+1;RelBx number", 10, -5, -5);
  h1_["RecHitRelBx_RE+2"] = fs_->make<TH1F>("RecHitRelBx_RE+2", "RelBx number from all RecHit in the Endcap+2;RelBx number", 10, -5, -5);
  h1_["RecHitRelBx_RE+3"] = fs_->make<TH1F>("RecHitRelBx_RE+3", "RelBx number from all RecHit in the Endcap+3;RelBx number", 10, -5, -5);
  h1_["RecHitRelBx_RE-1"] = fs_->make<TH1F>("RecHitRelBx_RE-1", "RelBx number from all RecHit in the Endcap-1;RelBx number", 10, -5, -5);
  h1_["RecHitRelBx_RE-2"] = fs_->make<TH1F>("RecHitRelBx_RE-2", "RelBx number from all RecHit in the Endcap-2;RelBx number", 10, -5, -5);
  h1_["RecHitRelBx_RE-3"] = fs_->make<TH1F>("RecHitRelBx_RE-3", "RelBx number from all RecHit in the Endcap-3;RelBx number", 10, -5, -5);

  h2_["BxVsNDigi"] = fs_->make<TH2F>("h2BxVsNDigi", "Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 250, 0, 250);
  h2_["BxVsNRecHit"] = fs_->make<TH2F>("h2BxVsNRecHit", "Bx number vs Number of RecHit;Bx number;Number of RecHit", 3600, 0, 3600, 250, 0, 250);

  h2_["BxVsNDigi_RE+1"] = fs_->make<TH2F>("h2BxVsNDigi_RE+1", "Endcap+1 Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 250, 0, 250);
  h2_["BxVsNDigi_RE+2"] = fs_->make<TH2F>("h2BxVsNDigi_RE+2", "Endcap+2 Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 250, 0, 250);
  h2_["BxVsNDigi_RE+3"] = fs_->make<TH2F>("h2BxVsNDigi_RE+3", "Endcap+3 Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 250, 0, 250);
  h2_["BxVsNDigi_RE-1"] = fs_->make<TH2F>("h2BxVsNDigi_RE-1", "Endcap-1 Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 250, 0, 250);
  h2_["BxVsNDigi_RE-2"] = fs_->make<TH2F>("h2BxVsNDigi_RE-2", "Endcap-2 Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 250, 0, 250);
  h2_["BxVsNDigi_RE-3"] = fs_->make<TH2F>("h2BxVsNDigi_RE-3", "Endcap-3 Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 250, 0, 250);

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
  //const int orbitNumber = event.orbitNumber();
  //const int storeNumber = event.storeNumber();

  if ( bxNumber < 0 )
  {
    cout << "DEBUG : Negative Bx number!!!" << bxNumber << endl;
  }

  if ( bxNumber > maxBxNumber_ ) maxBxNumber_ = bxNumber;
  if ( bxNumber < minBxNumber_ ) minBxNumber_ = bxNumber;

  edm::Handle<RPCDigiCollection> digiHandle;
  event.getByLabel(digiLabel_, digiHandle);

  edm::Handle<RPCRecHitCollection> recHitHandle;
  event.getByLabel(recHitLabel_, recHitHandle);

  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  h1_["BxNumber"]->Fill(bxNumber);

  unsigned int nDigis = 0;
  for ( RPCDigiCollection::DigiRangeIterator detUnitIt = digiHandle->begin();
        detUnitIt != digiHandle->end(); ++detUnitIt )
  {
    const RPCDetId Rsid = (*detUnitIt).first;
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(Rsid));
    if ( !roll ) continue;
    const RPCDigiCollection::Range& range = (*detUnitIt).second;

    // Determine the histogram name
    TString detName;
    if ( Rsid.region() == 1 ) detName = "RE+";
    else if ( Rsid.region() == -1 ) detName = "RE-";
    else continue;
    detName += Form("%d", Rsid.station());

    double sumBx = 0;
    for ( RPCDigiCollection::const_iterator digiIt = range.first; 
          digiIt != range.second; ++digiIt )
    {
      sumBx += digiIt->bx();
      nDigis++;

      h1_[string("DigiBx_"+detName)]->Fill(bxNumber+digiIt->bx());
      h1_[string("DigiRelBx_"+detName)]->Fill(digiIt->bx());
    }
    h1_["NDigi"]->Fill(nDigis);
    h2_[string("BxVsNDigi_"+detName)]->Fill(bxNumber+sumBx/nDigis, nDigis);

  }

  unsigned int nRecHit = 0;
  for ( RPCRecHitCollection::const_iterator recHitIter = recHitHandle->begin();
        recHitIter != recHitHandle->end(); ++recHitIter )
  {
    const RPCDetId detId = recHitIter->rpcId();

    TString detName;
    if ( detId.region() == 1 ) detName = "RE+";
    else if ( detId.region() == -1 ) detName = "RE-";
    else continue;
    detName += Form("%d", detId.station());

    nRecHit++;

    h1_[string("RecHitBx_"+detName)]->Fill(bxNumber+recHitIter->BunchX());
    h1_[string("RecHitRelBx_"+detName)]->Fill(recHitIter->BunchX());
  }
  h1_["NRecHit"]->Fill(nRecHit);

  prf_["BxVsNDigi"]->Fill(bxNumber, nDigis);
  prf_["BxVsNRecHit"]->Fill(bxNumber, nRecHit);

  h2_["BxVsNDigi"]->Fill(bxNumber, nDigis);
  h2_["BxVsNRecHit"]->Fill(bxNumber, nRecHit);
}

