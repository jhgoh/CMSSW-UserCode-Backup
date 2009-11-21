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
#include <TProfile.h>
#include <TString.h>
#include <TDatime.h>

#include <memory>

using namespace std;

std::string getSubDetName(const int region, const int ring, const int station)
{
  return Form("%d_%d_%d", region, ring, station);
}

bool isValidDetId(const RPCObPVSSmap::Item& pvss)
{
  if ( pvss.region < RPCDetId::minRegionId || pvss.region > RPCDetId::maxRegionId ) return false;
  if ( pvss.ring < RPCDetId::minRingBarrelId || pvss.ring > RPCDetId::maxRingBarrelId ) return false;
  if ( pvss.station < RPCDetId::minStationId || pvss.station > RPCDetId::maxStationId ) return false;
  if ( pvss.sector == 0 && pvss.layer == 0 && pvss.subsector == 0 ) return false;

  else return true;
}

MuonRPCAnalyzer::MuonRPCAnalyzer(const edm::ParameterSet& pset)
{
  rpcIMinTime_ = rpcVMinTime_ = rpcTMinTime_ = 0;
  rpcIMaxTime_ = rpcVMaxTime_ = rpcTMaxTime_ = 0;

  digiLabel_ = pset.getParameter<edm::InputTag>("digiLabel");

  edm::ParameterSet histoDimensions = pset.getParameter<edm::ParameterSet>("histoDimensions");
  
  const uint64_t minDateTime = histoDimensions.getUntrackedParameter<uint64_t>("minUTime");
  const uint64_t maxDateTime = histoDimensions.getUntrackedParameter<uint64_t>("maxUTime");
  const uint64_t dDateTime = histoDimensions.getUntrackedParameter<uint64_t>("dTime");
  timeOffset_ = minDateTime;

  h1_["strip"] = fs_->make<TH1F>("hStrip", "Strip profile", 100, 0, 100);
  h1_["bx"] = fs_->make<TH1F>("hBx", "Bunch crossing", 11, -5.5, 5.5);
  h1_["nDigi"] = fs_->make<TH1F>("hNDigi", "Number of digi", 100, 0, 100);

  h1_["T"] = fs_->make<TH1F>("hT", "Temperature;Temperature [#circC]", 100, 10, 25);
  h1_["I"] = fs_->make<TH1F>("hI", "Current;Current [#muA]", 100, 0, 10);
  h1_["V"] = fs_->make<TH1F>("hV", "Voltage;Voltage", 100, 0, 1);

  // Book profile histograms for Time vs Conditions
  prf_["TimeVsT"] = fs_->make<TProfile>("prfTimeVsT", "Time vs Temperature;Time [YY-MM-DD hh:mm:ss];Temperature [#circC]",
                                        TMath::Nint(double(maxDateTime-minDateTime)/dDateTime), 0, maxDateTime-minDateTime);
  prf_["TimeVsI"] = fs_->make<TProfile>("prfTimeVsI", "Time vs Current;Time [YY-MM-DD hh:mm:ss];Current [#muA]",
                                        TMath::Nint(double(maxDateTime-minDateTime)/dDateTime), 0, maxDateTime-minDateTime);
  prf_["TimeVsV"] = fs_->make<TProfile>("prfTimeVsV", "Time vs Voltage;Time [YY-MM-DD hh:mm:ss];Voltage [kV]",
                                        TMath::Nint(double(maxDateTime-minDateTime)/dDateTime), 0, maxDateTime-minDateTime);

  // Set axis to correspond to Timestamp format
  prf_["TimeVsT"]->GetXaxis()->SetTimeDisplay(1);
  prf_["TimeVsI"]->GetXaxis()->SetTimeDisplay(1);
  prf_["TimeVsV"]->GetXaxis()->SetTimeDisplay(1);

  //const TString timeFormat("%d-%m-%y");
  const TString timeFormat("%d %H:%M:%S");
  prf_["TimeVsT"]->GetXaxis()->SetTimeFormat(timeFormat);
  prf_["TimeVsI"]->GetXaxis()->SetTimeFormat(timeFormat);
  prf_["TimeVsV"]->GetXaxis()->SetTimeFormat(timeFormat);

  prf_["TimeVsT"]->GetXaxis()->SetTimeOffset(timeOffset_);
  prf_["TimeVsT"]->GetXaxis()->SetTimeOffset(timeOffset_);
  prf_["TimeVsT"]->GetXaxis()->SetTimeOffset(timeOffset_);

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

    h1_[Form("%d_T", region)] = regionDir.make<TH1F>("hT", "Temperature;Temperature [#circC]", 100, 0, 25);
    h1_[Form("%d_I", region)] = regionDir.make<TH1F>("hI", "Current;Current [#muA]", 100, 0, 10);
    h1_[Form("%d_V", region)] = regionDir.make<TH1F>("hV", "Voltage;Voltage", 100, 0, 1);

    for ( int ring = minRingId; ring <= maxRingId; ++ring )
    {
      for ( int station = RPCDetId::minStationId; station <= RPCDetId::maxStationId; ++station )
      {
        const string subDetName = getSubDetName(region, ring ,station);
        TFileDirectory subDir = regionDir.mkdir(subDetName, subDetName);

        h1_[subDetName+"_T"] = subDir.make<TH1F>("hT", "Temperature;Temperature [#circC]", 100, 0, 25);
        h1_[subDetName+"_I"] = subDir.make<TH1F>("hI", "Current;Current [#muA]", 100, 0, 10);
        h1_[subDetName+"_V"] = subDir.make<TH1F>("hV", "Voltage;Voltage", 100, 0, 1);

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
  cout << "---------------------------\n"
       << " endJob() called           \n"
       << " rpcT TimeRange = " << rpcTMinTime_ << ":" << rpcTMaxTime_ << endl
       << " rpcI TimeRange = " << rpcIMinTime_ << ":" << rpcIMaxTime_ << endl
       << " rpcV TimeRange = " << rpcVMinTime_ << ":" << rpcVMaxTime_ << endl;
}

void MuonRPCAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  unsigned long long eventDAQTime = event.time().value();
  const uint64_t eventTime = RPCRunIOV::DAQtoUNIX(&eventDAQTime);

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

  // Analyze conditional information only when the IOV is changed
  if ( eventTime < rpcIMinTime_ || eventTime > rpcIMaxTime_ )
  {
    std::vector<RPCObImon::I_Item> rpcImon = rpcRunIOV.getImon();

    // Retrieve IOV information
    for ( std::vector<RPCObImon::I_Item>::const_iterator imon = rpcImon.begin();
          imon != rpcImon.end(); ++imon )
    {
      const double rpcI = imon->value;
      const uint64_t rpcITime = RPCRunIOV::toUNIX(imon->day, imon->time);

      if ( rpcIMinTime_ == 0 || rpcIMinTime_ > rpcITime ) rpcIMinTime_ = rpcITime;
      if ( rpcIMaxTime_ < rpcITime ) rpcIMaxTime_ = rpcITime;

      const RPCObPVSSmap::Item pvss = pvssMap[imon->dpid];
      if ( !isValidDetId(pvss) ) continue;

      const string subDetName = getSubDetName(pvss.region, pvss.ring, pvss.station);
      rpcIValues_[subDetName].first++;
      rpcIValues_[subDetName].second += rpcI;

      // Fill condition values averaging over all detector cells
      prf_["TimeVsI"]->Fill(rpcITime-timeOffset_, rpcI);
      cout << "I(" << prf_["TimeVsI"]->GetXaxis()->GetXmin() << " - " << rpcITime << " - " << prf_["TimeVsI"]->GetXaxis()->GetXmax() << endl;
    }

    // Calclulate average condition values in IOV for each detector cell and fill histograms
    for ( DetIOVMap::const_iterator detIOV = rpcIValues_.begin(); detIOV != rpcIValues_.end(); ++detIOV)
    {
      const string& subDetName = detIOV->first;
      const pair<unsigned int, double>& numberAndSum = detIOV->second;

      const unsigned int n = numberAndSum.first;
      const double sum = numberAndSum.second;

      if ( n != 0 ) h1_[subDetName+"_I"]->Fill(sum/n);
      h1_["I"]->Fill(sum/n);
    }
  }

  if ( eventTime < rpcIMinTime_ || eventTime > rpcIMaxTime_ )
  {
    std::vector<RPCObVmon::V_Item> rpcVmon = rpcRunIOV.getVmon();

    for ( std::vector<RPCObVmon::V_Item>::const_iterator vmon = rpcVmon.begin();
          vmon != rpcVmon.end(); ++vmon )
    {
      const double rpcV = vmon->value;
      const uint64_t rpcVTime = RPCRunIOV::toUNIX(vmon->day, vmon->time);

      if ( rpcVMinTime_ ==0 || rpcVMinTime_ > rpcVTime ) rpcVMinTime_ = rpcVTime;
      if ( rpcVMaxTime_ < rpcVTime ) rpcVMaxTime_ = rpcVTime;

      const RPCObPVSSmap::Item pvss = pvssMap[vmon->dpid];
      if ( !isValidDetId(pvss) ) continue;

      const string subDetName = getSubDetName(pvss.region, pvss.ring, pvss.station);
      rpcVValues_[subDetName].first++;
      rpcVValues_[subDetName].second += rpcV;

      // Fill condition values averaging over all detector cells
      prf_["TimeVsV"]->Fill(rpcVTime-timeOffset_, rpcV);
      cout << "V(" << prf_["TimeVsV"]->GetXaxis()->GetXmin() << " - " << rpcVTime << " - " << prf_["TimeVsV"]->GetXaxis()->GetXmax() << endl;
    }

    // Calclulate average condition values in IOV for each detector cell and fill histograms
    for ( DetIOVMap::const_iterator detIOV = rpcVValues_.begin(); detIOV != rpcVValues_.end(); ++detIOV)
    {
      const string& subDetName = detIOV->first;
      const pair<unsigned int, double>& numberAndSum = detIOV->second;

      const unsigned int n = numberAndSum.first;
      const double sum = numberAndSum.second;

      if ( n != 0) h1_[subDetName+"_V"]->Fill(sum/n);
      h1_["V"]->Fill(sum/n);
    }
  }

  if ( eventTime < rpcTMinTime_ || eventTime > rpcTMaxTime_ )
  {
    std::vector<RPCObTemp::T_Item> rpcTmon = rpcRunIOV.getTemp();

    for ( std::vector<RPCObTemp::T_Item>::const_iterator tmon = rpcTmon.begin();
          tmon != rpcTmon.end(); ++tmon )
    {
      const double rpcT = tmon->value;
      const uint64_t rpcTTime = RPCRunIOV::toUNIX(tmon->day, tmon->time);

      if ( rpcTMinTime_ == 0 || rpcTMinTime_ > rpcTTime ) rpcTMinTime_ = rpcTTime;
      if ( rpcTMaxTime_ < rpcTTime ) rpcTMaxTime_ = rpcTTime;

      const RPCObPVSSmap::Item pvss = pvssMap[tmon->dpid];
      if ( !isValidDetId(pvss) ) continue;

      const string subDetName = getSubDetName(pvss.region, pvss.ring, pvss.station);
      rpcTValues_[subDetName].first++;
      rpcTValues_[subDetName].second += rpcT;

      // Fill condition values averaging over all detector cells
      prf_["TimeVsT"]->Fill(rpcTTime-timeOffset_, rpcT);
      cout << "T(" << prf_["TimeVsT"]->GetXaxis()->GetXmin() << " - " << rpcTTime << " - " << prf_["TimeVsT"]->GetXaxis()->GetXmax() << endl;
    }

    // Calclulate average condition values in IOV for each detector cell and fill histograms
    for ( DetIOVMap::const_iterator detIOV = rpcTValues_.begin(); detIOV != rpcTValues_.end(); ++detIOV)
    {
      const string& subDetName = detIOV->first;
      const pair<unsigned int, double>& numberAndSum = detIOV->second;

      const unsigned int n = numberAndSum.first;
      const double sum = numberAndSum.second;

      if ( n != 0 ) h1_[subDetName+"_T"]->Fill(sum/n);
      h1_["T"]->Fill(sum/n);
    }
  }

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
}

