#include "CosmicAnalysis/CosmicAnalyzer/interface/MuonRPCAnalyzer.h"

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
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "CondTools/RPC/interface/RPCRunIOV.h"

#include <TMath.h>
#include <TH1F.h>
#include <TString.h>

#include <memory>
#include <string>

using namespace std;

std::string getSubDetName(const int region, const int ring, const int station)
{
  return Form("%d_%d_%d", region, ring, station);
}

MuonRPCAnalyzer::MuonRPCAnalyzer(const edm::ParameterSet& pset)
{
  digiLabel_ = pset.getParameter<edm::InputTag>("digiLabel");
  
  h1_["strip"] = fs_->make<TH1F>("hStrip", "Strip profile", 100, 0, 100);
  h1_["bx"] = fs_->make<TH1F>("hBx", "Bunch crossing", 11, -5.5, 5.5);
  h1_["nDigi"] = fs_->make<TH1F>("hNDigi", "Number of digi", 100, 0, 100);

  h1_["T"] = fs_->make<TH1F>("hT", "Temperature;Temperature", 100, 15, 30);
  h1_["I"] = fs_->make<TH1F>("hI", "Current;Current", 100, 0, 1e5);
  h1_["V"] = fs_->make<TH1F>("hV", "Voltage;Voltage", 100, 7000, 10000);

  for ( int region = RPCDetId::minRegionId; region <= RPCDetId::maxRegionId; ++region )
  {
    const int minRingId = (region==0) ? RPCDetId::minRingBarrelId : RPCDetId::minRingForwardId;
    const int maxRingId = (region==0) ? RPCDetId::maxRingBarrelId : RPCDetId::maxRingForwardId;
    string regionName;
    switch(region)
    {
      case  0: regionName = "Barr"; break;
      case  1: regionName = "FwdP"; break;
      case -1: regionName = "FwdN"; break;
    }
    TFileDirectory regionDir = fs_->mkdir(regionName, regionName);

    // Summary plot for each regions
    h1_[Form("%d_strip", region)] = regionDir.make<TH1F>("hStrip", "Strip profile", 100, 0, 100);
    h1_[Form("%d_bx", region)] = regionDir.make<TH1F>("hBx", "Bunch crossing", 11, -5.5, 5.5);
    h1_[Form("%d_nDigi", region)] = regionDir.make<TH1F>("hNDigi", "Number of digi", 100, 0, 100);

    h1_[Form("%d_T", region)] = regionDir.make<TH1F>("hT", "Temperature;Temperature", 100, 15, 30);
    h1_[Form("%d_I", region)] = regionDir.make<TH1F>("hI", "Current;Current", 100, 0, 1e5);
    h1_[Form("%d_V", region)] = regionDir.make<TH1F>("hV", "Voltage;Voltage", 100, 7000, 10000);

    for ( int ring = minRingId; ring <= maxRingId; ++ring )
    {
      for ( int station = RPCDetId::minStationId; station <= RPCDetId::maxStationId; ++station )
      {
        const string subDetName = getSubDetName(region, ring ,station);
        TFileDirectory subDir = regionDir.mkdir(subDetName, subDetName);

        h1_[subDetName+"_T"] = subDir.make<TH1F>("hT", "Temperature;Temperature", 100, 15, 30);
        h1_[subDetName+"_I"] = subDir.make<TH1F>("hI", "Current;Current", 100, 0, 1e5);
        h1_[subDetName+"_V"] = subDir.make<TH1F>("hV", "Voltage;Voltage", 100, 7000, 10000);

        rpcIValues_[subDetName] = std::make_pair(0,0);
        rpcVValues_[subDetName] = std::make_pair(0,0);
        rpcTValues_[subDetName] = std::make_pair(0,0);
      }
    }
  }
}

MuonRPCAnalyzer::~MuonRPCAnalyzer()
{
  
}

void MuonRPCAnalyzer::beginJob(const edm::EventSetup& eventSetup)
{
}

void MuonRPCAnalyzer::endJob()
{
/*
  if ( !eventNumbers_.empty() )
  {
    const float* evtNumArr = &eventNumbers_[0];
    const float* na = 0;

    if ( !rpcAvgI_.empty() )
    {
      const int size = rpcAvgI_.size();
      TGraphErrors* grpAvgI = fs_->make<TGraphErrors>(size, evtNumArr, &rpcAvgI_[0], na, &rpcErrI_[0]);
      grpAvgI->SetName("grpEvtNumVsI");
      grpAvgI->SetTitle("Event number vs Current;Event number;Current");
    }
    else
    {
      cerr << "No rpc current information" << endl;
    }

    if ( !rpcAvgV_.empty() )
    {
      const int size = rpcAvgV_.size();
      TGraphErrors* grpAvgV = fs_->make<TGraphErrors>(size, evtNumArr, &rpcAvgI_[0], na, &rpcErrV_[0]);
      grpAvgV->SetName("grpEvtNumVsV");
      grpAvgV->SetTitle("Event number vs Voltage;Event number;Voltage");
    }
    else
    {
      cerr << "No RPC voltage information" << endl;
    }

    if ( !rpcAvgT_.empty() )
    {
      const int size = rpcAvgT_.size();
      TGraphErrors* grpAvgT = fs_->make<TGraphErrors>(size, evtNumArr, &rpcAvgT_[0], na, &rpcErrT_[0]);
      grpAvgT->SetName("grpEvtNumVsT");
      grpAvgT->SetTitle("Event number vs Temperature;Event number;Temperature");
    }
    else
    {
      cerr << "No RPC temperature information" << endl;
    }
  }
*/
}

void MuonRPCAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  // Initialize the RPC cond values map
  for ( DetIOVMap::iterator iValue = rpcIValues_.begin(); iValue != rpcIValues_.end(); ++iValue )
  {
    iValue->second = std::make_pair(0,0.);
  }
  for ( DetIOVMap::iterator vValue = rpcVValues_.begin(); vValue != rpcVValues_.end(); ++vValue )
  {
    vValue->second = std::make_pair(0,0.);
  }
  for ( DetIOVMap::iterator tValue = rpcTValues_.begin(); tValue != rpcTValues_.end(); ++tValue )
  {
    tValue->second = std::make_pair(0,0.);
  }

  // Analyze IOV information
  RPCRunIOV rpcRunIOV(eventSetup);

  std::map<int, RPCObPVSSmap::Item> pvssMap = rpcRunIOV.getPVSSMap();

  std::vector<RPCObImon::I_Item> rpcImon = rpcRunIOV.getImon();
  std::vector<RPCObVmon::V_Item> rpcVmon = rpcRunIOV.getVmon();
  std::vector<RPCObTemp::T_Item> rpcTmon = rpcRunIOV.getTemp();

  // Retrieve IOV information
  for ( std::vector<RPCObImon::I_Item>::const_iterator imon = rpcImon.begin();
        imon != rpcImon.end(); ++imon )
  {
    const double rpcI = imon->value;
    const RPCObPVSSmap::Item pvss = pvssMap[imon->dpid];

    const string subDetName = getSubDetName(pvss.region, pvss.ring, pvss.station);
    rpcIValues_[subDetName].first++;
    rpcIValues_[subDetName].second += rpcI;
  }

  for ( std::vector<RPCObVmon::V_Item>::const_iterator vmon = rpcVmon.begin();
        vmon != rpcVmon.end(); ++vmon )
  {
    const double rpcV = vmon->value;
    const RPCObPVSSmap::Item pvss = pvssMap[vmon->dpid];

    const string subDetName = getSubDetName(pvss.region, pvss.ring, pvss.station);
    rpcVValues_[subDetName].first++;
    rpcVValues_[subDetName].second += rpcV;
  }

  for ( std::vector<RPCObTemp::T_Item>::const_iterator tmon = rpcTmon.begin();
        tmon != rpcTmon.end(); ++tmon )
  {
    const double rpcT = tmon->value;
    const RPCObPVSSmap::Item pvss = pvssMap[tmon->dpid];

    const string subDetName = getSubDetName(pvss.region, pvss.ring, pvss.station);
    rpcTValues_[subDetName].first++;
    rpcTValues_[subDetName].second += rpcT;
  }

  // Calclulate average values of IOV for each detector cell and fill histograms
  for ( DetIOVMap::const_iterator detIOV = rpcTValues_.begin(); detIOV != rpcTValues_.end(); ++detIOV)
  {
    const string& subDetName = detIOV->first;
    const pair<unsigned int, double>& numberAndSum = detIOV->second;

    const unsigned int n = numberAndSum.first;
    const double sum = numberAndSum.second;

    const double avg = (n == 0) ? 0 : sum/n;

    h1_[subDetName+"_T"]->Fill(avg);
  }

  for ( DetIOVMap::const_iterator detIOV = rpcIValues_.begin(); detIOV != rpcIValues_.end(); ++detIOV)
  {
    const string& subDetName = detIOV->first;
    const pair<unsigned int, double>& numberAndSum = detIOV->second;

    const unsigned int n = numberAndSum.first;
    const double sum = numberAndSum.second;

    const double avg = (n == 0) ? 0 : sum/n;

    h1_[subDetName+"_I"]->Fill(avg);
  }

  for ( DetIOVMap::const_iterator detIOV = rpcVValues_.begin(); detIOV != rpcVValues_.end(); ++detIOV)
  {
    const string& subDetName = detIOV->first;
    const pair<unsigned int, double>& numberAndSum = detIOV->second;

    const unsigned int n = numberAndSum.first;
    const double sum = numberAndSum.second;

    const double avg = (n == 0) ? 0 : sum/n;

    h1_[subDetName+"_V"]->Fill(avg);
  }

///////////// Ignore below //////////////
/*
  unsigned int nI = 0, nV = 0, nT = 0;
  double sumI = 0, sumV = 0, sumT = 0;
  double sumI2 = 0, sumV2 = 0, sumT2 = 0;

  for ( std::vector<RPCObImon::I_Item>::const_iterator imon = rpcImon.begin();
        imon != rpcImon.end(); ++imon )
  {
    const double rpcI = imon->value;
    sumI += rpcI;
    sumI2 += rpcI*rpcI;
    ++nI;

    const RPCObPVSSmap::Item pvss = pvssMap[imon->dpid];
    const TString detName = detNameToForm(pvss.region, pvss.ring, pvss.station);

    h1_[detName+"_I"]->Fill(rpcI);
  }
  const double avgI = sumI/nI;
  const double errI = TMath::Sqrt((sumI2 - sumI*sumI/nI)/nI);
  //rpcAvgI_.push_back(avgI);
  //rpcErrI_.push_back(errI);

  for ( std::vector<RPCObVmon::V_Item>::const_iterator vmon = rpcVmon.begin();
        vmon != rpcVmon.end(); ++vmon )
  {
    const double rpcV = vmon->value;
    sumV += rpcV;
    sumV2 += rpcV*rpcV;
    ++nV;
  }
  const double avgV = sumV/nV;
  const double errV = TMath::Sqrt((sumV2 - sumV*sumV/nV)/nV);
  //rpcAvgV_.push_back(avgV);
  //rpcErrV_.push_back(errV);

  for ( std::vector<RPCObTemp::T_Item>::const_iterator tmon = rpcTmon.begin();
        tmon != rpcTmon.end(); ++tmon )
  {
    const double rpcT = tmon->value;
    sumT += rpcT;
    sumT2 += rpcT*rpcT;
    ++nT;
  }
  const double avgT = sumT/nT;
  const double errT = TMath::Sqrt((sumT2 - sumT*sumT/nT)/nT);
  //rpcAvgT_.push_back(avgT);
  //rpcErrT_.push_back(errT);

  // Analyze Digi information
  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  edm::Handle<RPCDigiCollection> rpcDigis;
  event.getByLabel(digiLabel_, rpcDigis);

  // Loop over Digis
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

      h1_["strip"]->Fill(strip);
      h1_["bx"]->Fill(bx);

      const GlobalPoint gPoint = roll->toGlobal(roll->centreOfStrip(strip));
      const double globalX = gPoint.x();
      const double globalY = gPoint.y();
      
      if ( abs(Rsid.ring()) > 2 ) continue;
    }
    h1_["nDigi"]->Fill(nDigi);
  }
*/
}

