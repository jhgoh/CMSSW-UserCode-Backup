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
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TProfile2D.h>

#include <memory>
#include <string>

namespace H1
{
  enum H
  {
    Bx, Strip, nDigi,
    Temperature
  };
}

using namespace std;

MuonRPCAnalyzer::MuonRPCAnalyzer(const edm::ParameterSet& pset)
{
  digiLabel_ = pset.getParameter<edm::InputTag>("digiLabel");

  h1_[H1::Strip] = fs_->make<TH1F>("hStrip", "Strip profile", 100, 0, 100);
  h1_[H1::Bx] = fs_->make<TH1F>("hBx", "Bunch crossing", 11, -5.5, 5.5);
  h1_[H1::nDigi] = fs_->make<TH1F>("hNDigi", "Number of digi", 100, 0, 100);

  h1_[H1::Temperature] = fs_->make<TH1F>("hTemperature", "Temerature", 100, 0, 30);
}

MuonRPCAnalyzer::~MuonRPCAnalyzer()
{
  
}

void MuonRPCAnalyzer::beginJob(const edm::EventSetup& eventSetup)
{
}

void MuonRPCAnalyzer::endJob()
{
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

    if ( !rpcAvgV_.empty() )
    {
      const int size = rpcAvgV_.size();
      TGraphErrors* grpAvgV = fs_->make<TGraphErrors>(size, evtNumArr, &rpcAvgI_[0], na, &rpcErrV_[0]);
      grpAvgV->SetName("grpEvtNumVsV");
      grpAvgV->SetTitle("Event number vs Voltage;Event number;Voltage");
    }

    if ( !rpcAvgT_.empty() )
    {
      const int size = rpcAvgT_.size();
      TGraphErrors* grpAvgT = fs_->make<TGraphErrors>(size, evtNumArr, &rpcAvgT_[0], na, &rpcErrT_[0]);
      grpAvgT->SetName("grpEvtNumVsT");
      grpAvgT->SetTitle("Event number vs Temperature;Event number;Temperature");
    }
  }
}

void MuonRPCAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  //eventNumbers_.push_back(event.id().event());
  cout << event.id().event() << endl;

  // Analyze IOV information
  RunRunIOV rpcRunIOV(eventSetup);

  std::vector<RPCObImon::I_Item> imon = rpcRunIOV.getImon();
  std::vector<RPCObVmon::V_Item> vmon = rpcRunIOV.getVmon();
  std::vector<RPCObTemp::T_Item> tmon = rpcRunIOV.getTemp();

  unsigned int nI = 0, nV = 0, nT = 0;
  double sumI = 0, sumV = 0, sumT = 0;
  double sumI2 = 0, sumV2 = 0, sumT2 = 0;

  return;
  for ( std::vector<RPCObImon::I_Item>::const_iterator imon = rpcImon.begin();
        imon != rpcImon.end(); ++imon )
  {
    const double rpcI = imon->value;
    sumI += rpcI;
    sumI2 += rpcI*rpcI;
    ++nI;
  }
  const double avgI = sumI/nI;
  const double errI = TMath::Sqrt((sumI2 - sumI*sumI/nI)/nI);
  cout << event.id().event() << ' ' << avgI << ' ' << errI << endl;
  //rpcAvgI_.push_back(avgI);
  //rpcErrI_.push_back(errI);

  return;

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

      h1_[H1::Strip]->Fill(strip);
      h1_[H1::Bx]->Fill(bx);

/*
      const GlobalPoint gPoint = roll->toGlobal(roll->centreOfStrip(strip));
      const double globalX = gPoint.x();
      const double globalY = gPoint.y();
      
      if ( abs(Rsid.ring()) > 2 ) continue;
*/
    }
    h1_[H1::nDigi]->Fill(nDigi);
  }
}

