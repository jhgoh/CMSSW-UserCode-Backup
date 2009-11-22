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

#include "DataFormats/RPCDigi/interface/QualityTrigBx.h"
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
  minNDigiCut_ = pset.getUntrackedParameter<int>("minNDigiCut", 100);

  digiLabel_ = pset.getParameter<edm::InputTag>("digiLabel");
  recHitLabel_ = pset.getParameter<edm::InputTag>("recHitLabel");

  // Book profile histograms for Time vs Conditions
  TString histTitle;

  prf_["BxVsNDigi"] = fs_->make<TProfile>("prfBxVsNDigi", "Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0,3600);
  prf_["BxVsNRecHit"] = fs_->make<TProfile>("prfBxVsNRecHit", "Bx number vs Number of RecHit;Bx number;Number of RecHit", 3600, 0, 3600);

  h1_["NDigi"] = fs_->make<TH1F>("hNDigi", "Number of Digis;Number of Digis", 250, 0, 250);
  h1_["NRecHit"] = fs_->make<TH1F>("hNRecHit", "Number of RecHit;Number of RecHit", 250, 0, 250);

  h1_["BxNumber"] = fs_->make<TH1F>("hBxNumber", "Bx number;Bx number;Number of Events", 3600, 0, 3600);

  TFileDirectory dirREP1 = fs_->mkdir("RE+1");
  TFileDirectory dirREP2 = fs_->mkdir("RE+2");
  TFileDirectory dirREP3 = fs_->mkdir("RE+3");
  TFileDirectory dirREN1 = fs_->mkdir("RE-1");
  TFileDirectory dirREN2 = fs_->mkdir("RE-2");
  TFileDirectory dirREN3 = fs_->mkdir("RE-3");

  h1_["RE+1_DigiBx"] = dirREP1.make<TH1F>("DigiBx", "Bx number from all Digis in the Endcap+1;Bx number", 3600, 0, 3600);
  h1_["RE+2_DigiBx"] = dirREP2.make<TH1F>("DigiBx", "Bx number from all Digis in the Endcap+2;Bx number", 3600, 0, 3600);
  h1_["RE+3_DigiBx"] = dirREP3.make<TH1F>("DigiBx", "Bx number from all Digis in the Endcap+3;Bx number", 3600, 0, 3600);
  h1_["RE-1_DigiBx"] = dirREN1.make<TH1F>("DigiBx", "Bx number from all Digis in the Endcap-1;Bx number", 3600, 0, 3600);
  h1_["RE-2_DigiBx"] = dirREN2.make<TH1F>("DigiBx", "Bx number from all Digis in the Endcap-2;Bx number", 3600, 0, 3600);
  h1_["RE-3_DigiBx"] = dirREN3.make<TH1F>("DigiBx", "Bx number from all Digis in the Endcap-3;Bx number", 3600, 0, 3600);

  h1_["RE+1_DigiRelBx"] = dirREP1.make<TH1F>("DigiRelBx", "RelBx number from all Digis in the Endcap+1;RelBx number", 10, -5, 5);
  h1_["RE+2_DigiRelBx"] = dirREP2.make<TH1F>("DigiRelBx", "RelBx number from all Digis in the Endcap+2;RelBx number", 10, -5, 5);
  h1_["RE+3_DigiRelBx"] = dirREP3.make<TH1F>("DigiRelBx", "RelBx number from all Digis in the Endcap+3;RelBx number", 10, -5, 5);
  h1_["RE-1_DigiRelBx"] = dirREN1.make<TH1F>("DigiRelBx", "RelBx number from all Digis in the Endcap-1;RelBx number", 10, -5, 5);
  h1_["RE-2_DigiRelBx"] = dirREN2.make<TH1F>("DigiRelBx", "RelBx number from all Digis in the Endcap-2;RelBx number", 10, -5, 5);
  h1_["RE-3_DigiRelBx"] = dirREN3.make<TH1F>("DigiRelBx", "RelBx number from all Digis in the Endcap-3;RelBx number", 10, -5, 5);

  h1_["RE+1_DigiBxWithNDigiCut"] = dirREP1.make<TH1F>("DigiBxWithNDigiCut", 
                                                      Form("Bx number with nDigi > %d in the Endcap+1;Bx number", minNDigiCut_), 3600, 0, 3600);
  h1_["RE+2_DigiBxWithNDigiCut"] = dirREP2.make<TH1F>("DigiBxWithNDigiCut",
                                                      Form("Bx number with nDigi > %d in the Endcap+2;Bx number", minNDigiCut_), 3600, 0, 3600);
  h1_["RE+3_DigiBxWithNDigiCut"] = dirREP3.make<TH1F>("DigiBxWithNDigiCut",
                                                      Form("Bx number with nDigi > %d in the Endcap+3;Bx number", minNDigiCut_), 3600, 0, 3600);
  h1_["RE-1_DigiBxWithNDigiCut"] = dirREN1.make<TH1F>("DigiBxWithNDigiCut",
                                                      Form("Bx number with nDigi > %d in the Endcap-1;Bx number", minNDigiCut_), 3600, 0, 3600);
  h1_["RE-2_DigiBxWithNDigiCut"] = dirREN2.make<TH1F>("DigiBxWithNDigiCut",
                                                      Form("Bx number with nDigi > %d in the Endcap-2;Bx number", minNDigiCut_), 3600, 0, 3600);
  h1_["RE-3_DigiBxWithNDigiCut"] = dirREN3.make<TH1F>("DigiBxWithNDigiCut",
                                                      Form("Bx number with nDigi > %d in the Endcap-3;Bx number", minNDigiCut_), 3600, 0, 3600);

  h1_["RE+1_DigiRelBxWithNDigiCut"] = dirREP1.make<TH1F>("DigiRelBxWithNDigiCut", 
                                                         Form("RelBx number with nDigi > %d in the Endcap+1;RelBx number", minNDigiCut_), 10, -5, 5);
  h1_["RE+2_DigiRelBxWithNDigiCut"] = dirREP2.make<TH1F>("DigiRelBxWithNDigiCut",                                                  
                                                         Form("RelBx number with nDigi > %d in the Endcap+2;RelBx number", minNDigiCut_), 10, -5, 5);
  h1_["RE+3_DigiRelBxWithNDigiCut"] = dirREP3.make<TH1F>("DigiRelBxWithNDigiCut",                                                  
                                                         Form("RelBx number with nDigi > %d in the Endcap+3;RelBx number", minNDigiCut_), 10, -5, 5);
  h1_["RE-1_DigiRelBxWithNDigiCut"] = dirREN1.make<TH1F>("DigiRelBxWithNDigiCut",                                                  
                                                         Form("RelBx number with nDigi > %d in the Endcap-1;RelBx number", minNDigiCut_), 10, -5, 5);
  h1_["RE-2_DigiRelBxWithNDigiCut"] = dirREN2.make<TH1F>("DigiRelBxWithNDigiCut",                                                  
                                                         Form("RelBx number with nDigi > %d in the Endcap-2;RelBx number", minNDigiCut_), 10, -5, 5);
  h1_["RE-3_DigiRelBxWithNDigiCut"] = dirREN3.make<TH1F>("DigiRelBxWithNDigiCut",                                                  
                                                         Form("RelBx number with nDigi > %d in the Endcap-3;RelBx number", minNDigiCut_), 10, -5, 5);

  h1_["RE+1_RecHitBx"] = dirREP1.make<TH1F>("RecHitBx", "Bx number from all RecHit in the Endcap+1;Bx number", 3600, 0, 3600);
  h1_["RE+2_RecHitBx"] = dirREP2.make<TH1F>("RecHitBx", "Bx number from all RecHit in the Endcap+2;Bx number", 3600, 0, 3600);
  h1_["RE+3_RecHitBx"] = dirREP3.make<TH1F>("RecHitBx", "Bx number from all RecHit in the Endcap+3;Bx number", 3600, 0, 3600);
  h1_["RE-1_RecHitBx"] = dirREN1.make<TH1F>("RecHitBx", "Bx number from all RecHit in the Endcap-1;Bx number", 3600, 0, 3600);
  h1_["RE-2_RecHitBx"] = dirREN2.make<TH1F>("RecHitBx", "Bx number from all RecHit in the Endcap-2;Bx number", 3600, 0, 3600);
  h1_["RE-3_RecHitBx"] = dirREN3.make<TH1F>("RecHitBx", "Bx number from all RecHit in the Endcap-3;Bx number", 3600, 0, 3600);

  h1_["RE+1_RecHitRelBx"] = dirREP1.make<TH1F>("RecHitRelBx", "RelBx number from all RecHit in the Endcap+1;RelBx number", 10, -5, 5);
  h1_["RE+2_RecHitRelBx"] = dirREP2.make<TH1F>("RecHitRelBx", "RelBx number from all RecHit in the Endcap+2;RelBx number", 10, -5, 5);
  h1_["RE+3_RecHitRelBx"] = dirREP3.make<TH1F>("RecHitRelBx", "RelBx number from all RecHit in the Endcap+3;RelBx number", 10, -5, 5);
  h1_["RE-1_RecHitRelBx"] = dirREN1.make<TH1F>("RecHitRelBx", "RelBx number from all RecHit in the Endcap-1;RelBx number", 10, -5, 5);
  h1_["RE-2_RecHitRelBx"] = dirREN2.make<TH1F>("RecHitRelBx", "RelBx number from all RecHit in the Endcap-2;RelBx number", 10, -5, 5);
  h1_["RE-3_RecHitRelBx"] = dirREN3.make<TH1F>("RecHitRelBx", "RelBx number from all RecHit in the Endcap-3;RelBx number", 10, -5, 5);

  h2_["BxVsNDigi"] = fs_->make<TH2F>("h2BxVsNDigi", "Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 500, 0, 500);
  h2_["BxVsNRecHit"] = fs_->make<TH2F>("h2BxVsNRecHit", "Bx number vs Number of RecHit;Bx number;Number of RecHit", 3600, 0, 3600, 500, 0, 500);

  h2_["RE+1_BxVsNDigi"] = dirREP1.make<TH2F>("h2BxVsNDigi", "Endcap+1 Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 500, 0, 500);
  h2_["RE+2_BxVsNDigi"] = dirREP2.make<TH2F>("h2BxVsNDigi", "Endcap+2 Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 500, 0, 500);
  h2_["RE+3_BxVsNDigi"] = dirREP3.make<TH2F>("h2BxVsNDigi", "Endcap+3 Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 500, 0, 500);
  h2_["RE-1_BxVsNDigi"] = dirREN1.make<TH2F>("h2BxVsNDigi", "Endcap-1 Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 500, 0, 500);
  h2_["RE-2_BxVsNDigi"] = dirREN2.make<TH2F>("h2BxVsNDigi", "Endcap-2 Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 500, 0, 500);
  h2_["RE-3_BxVsNDigi"] = dirREN3.make<TH2F>("h2BxVsNDigi", "Endcap-3 Bx number vs Number of Digis;Bx number;Number of Digis", 3600, 0, 3600, 500, 0, 500);

  h2_["NDigisVsEvtNum_SplashOnly"] = fs_->make<TH2F>("h2BxNumberVsEvtNum_SplashOnly", "Bx number vs Event number;Event number;Bx number", 1, 0, 1, 3600, 0, 3600);

  h2_["RE+1_BxOfStationVsEvtNum_SplashOnly"] = dirREN1.make<TH2F>("hBxOfStationVsEvtNum_SplashOnly", "Bx of stations vs Event number in RE+1;Event number;Bx of stations", 1, 0, 1, 1, 0, 1);
  h2_["RE+2_BxOfStationVsEvtNum_SplashOnly"] = dirREN2.make<TH2F>("hBxOfStationVsEvtNum_SplashOnly", "Bx of stations vs Event number in RE+2;Event number;Bx of stations", 1, 0, 1, 1, 0, 1);
  h2_["RE+3_BxOfStationVsEvtNum_SplashOnly"] = dirREN3.make<TH2F>("hBxOfStationVsEvtNum_SplashOnly", "Bx of stations vs Event number in RE+3;Event number;Bx of stations", 1, 0, 1, 1, 0, 1);
  h2_["RE-1_BxOfStationVsEvtNum_SplashOnly"] = dirREP1.make<TH2F>("hBxOfStationVsEvtNum_SplashOnly", "Bx of stations vs Event number in RE-1;Event number;Bx of stations", 1, 0, 1, 1, 0, 1);
  h2_["RE-2_BxOfStationVsEvtNum_SplashOnly"] = dirREP2.make<TH2F>("hBxOfStationVsEvtNum_SplashOnly", "Bx of stations vs Event number in RE-2;Event number;Bx of stations", 1, 0, 1, 1, 0, 1);
  h2_["RE-3_BxOfStationVsEvtNum_SplashOnly"] = dirREP3.make<TH2F>("hBxOfStationVsEvtNum_SplashOnly", "Bx of stations vs Event number in RE-3;Event number;Bx of stations", 1, 0, 1, 1, 0, 1);

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
  const int eventNumber = event.id().event();
  const int bxNumber = event.bunchCrossing();
  //const int orbitNumber = event.orbitNumber();
  //const int storeNumber = event.storeNumber();

  if ( bxNumber > maxBxNumber_ ) maxBxNumber_ = bxNumber;
  if ( bxNumber < minBxNumber_ ) minBxNumber_ = bxNumber;

  edm::Handle<RPCDigiCollection> digiHandle;
  event.getByLabel(digiLabel_, digiHandle);

  // Additional information to check wrong bx number
  edm::Handle<QualityTrigBx> qtbHandle;
  event.getByLabel(digiLabel_, qtbHandle);
  int headerTriggerBx = qtbHandle->bx();
  //cout << "DEBUG : eventNumber = " << eventNumber << " event.bx=" << bxNumber << " header.bx= " << headerTriggerBx << endl;

  edm::Handle<RPCRecHitCollection> recHitHandle;
  event.getByLabel(recHitLabel_, recHitHandle);

  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  h1_["BxNumber"]->Fill(headerTriggerBx);

  // First, collect list of valid RPC digis, without duplication
  int nDigi = 0;
  std::map<RPCDetId, std::set<RPCDigi> > digiToDetIdMap;
  for ( RPCDigiCollection::DigiRangeIterator detUnitIt = digiHandle->begin();
        detUnitIt != digiHandle->end(); ++detUnitIt )
  {
    const RPCDetId rpcDetId = (*detUnitIt).first;
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(rpcDetId));
    if ( !roll ) continue;
    const RPCDigiCollection::Range& range = (*detUnitIt).second;

    std::set<RPCDigi> digiSet;
    for ( RPCDigiCollection::const_iterator digiIt = range.first;
          digiIt != range.second; ++digiIt )
    {
      const RPCDigi& digi = *digiIt;

      // Keep non-duplicated list of digis
      if ( digiSet.find(digi) == digiSet.end() )
      {
        digiSet.insert(digi);
        ++nDigi;
      }
    }

    digiToDetIdMap.insert(std::make_pair(rpcDetId, digiSet));
  }

  h1_["NDigi"]->Fill(nDigi);
  prf_["BxVsNDigi"]->Fill(headerTriggerBx, nDigi);
  h2_["BxVsNDigi"]->Fill(headerTriggerBx, nDigi);

  for ( std::map<RPCDetId, std::set<RPCDigi> >::const_iterator digiToDetIdMapIter = digiToDetIdMap.begin();
        digiToDetIdMapIter != digiToDetIdMap.end(); ++digiToDetIdMapIter )
  {
    const RPCDetId& detId = digiToDetIdMapIter->first;
    const std::set<RPCDigi>& digiSet = digiToDetIdMapIter->second;

    for ( std::set<RPCDigi>::const_iterator rpcDigiIter = digiSet.begin();
          rpcDigiIter != digiSet.end(); ++rpcDigiIter )
    {
      const RPCDigi& digi = *rpcDigiIter;
      
      // Determine the histogram name
      // And skip for Barrels now (Here we are insterested in the splash events)
      TString detName;
      if ( detId.region() == 1 ) detName = "RE+";
      else if ( detId.region() == -1 ) detName = "RE-";
      else continue;
      detName += Form("%d", detId.station());

      h1_[string(detName+"_DigiBx")]->Fill(headerTriggerBx+digi.bx());
      h1_[string(detName+"_DigiRelBx")]->Fill(digi.bx());

      if ( nDigi > minNDigiCut_ ) // This should be the actual beam splash
      {
        h1_[string(detName+"_DigiBxWithNDigiCut")]->Fill(headerTriggerBx+digi.bx());
        h1_[string(detName+"_DigiRelBxWithNDigiCut")]->Fill(digi.bx());

        h2_["NDigisVsEvtNum_SplashOnly"]->Fill(Form("Event=%d", eventNumber), headerTriggerBx+digi.bx(), 1);
      }
    }
  }

  // Analyze for the RecHits
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

    h1_[string(detName+"_RecHitBx")]->Fill(headerTriggerBx+recHitIter->BunchX());
    h1_[string(detName+"_RecHitRelBx")]->Fill(recHitIter->BunchX());
  }

  h1_["NRecHit"]->Fill(nRecHit);
  prf_["BxVsNRecHit"]->Fill(headerTriggerBx, nRecHit);
  h2_["BxVsNRecHit"]->Fill(headerTriggerBx, nRecHit);
}

