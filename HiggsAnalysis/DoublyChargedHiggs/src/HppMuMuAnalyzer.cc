// -*- C++ -*-
//
// Package:    DoublyChargedHiggs
// Class:      HppMuMuAnalyzer
// 
/**\class HppMuMuAnalyzer HiggsAnalysis/DoublyChargedHiggs/src/HppMuMuAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jongseok Lee
//         Created:  Mon Jan 26 11:22:51 CET 2009
// $Id$
//
//

#include "HiggsAnalysis/DoublyChargedHiggs/src/HppMuMuAnalyzer.h"

// user include files
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//Jet
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenericJetCollection.h"
#include "DataFormats/JetReco/interface/GenericJet.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
 
//Jet Corrections
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/ChainedJetCorrector.h"

//MET
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CorrMETData.h"

//Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/VertexReconstructor.h"
#include "RecoVertex/VertexPrimitives/interface/VertexFitter.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableAdaptiveReconstructor.h"
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableAdaptiveFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "DataFormats/Common/interface/Handle.h" 
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTS.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h" 

//#include "DataFormats/MuonReco/interface/MuonIsolation.h"
//#include "DataFormats/MuonReco/interface/MuIsoDeposit.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
//
// class decleration
//
using namespace reco;
using namespace edm;
using namespace std;

bool debug=false, vtxon=true;
int events=0, nbackmcmuon=0, nbackmcmuon_ptetacut=0;
int isolPtcut=10;
float dphicut_4mu=2.5;
//float dzcut=0.2;
float Tptr_24etamucut[50]={0,}, Tptr2_24etamucut[50]={0,}, ncorgenjetpt_24etamucut[50]={0,};
float Teta_24etamucut[50]={0,}, Teta2_24etamucut[50]={0,}, ncorgenjeteta_24etamucut[50]={0,};
float Tptr_30pt24etamucut[50]={0,}, Tptr2_30pt24etamucut[50]={0,}, ncorgenjetpt_30pt24etamucut[50]={0,};
float Teta_30pt24etamucut[50]={0,}, Teta2_30pt24etamucut[50]={0,}, ncorgenjeteta_30pt24etamucut[50]={0,};
typedef std::pair<float, int> ptvsindex;

double deltaphi(double phi1, double phi2) 
{
  if (fabs(phi2-phi1)<TMath::Pi()) return phi2-phi1;
  else if (phi2-phi1>TMath::Pi()) return phi2-phi1-2*TMath::Pi();
  else return 2*TMath::Pi()+phi2-phi1;
//  if(phi1<0 && phi2>0) return phi2-phi1-2*TMath::Pi();
//  else if(phi1>0 && phi2<0) return phi2+2*TMath::Pi()-phi1;
//  else return phi2-phi1;
}

double isolPt(TLorentzVector t, edm::Handle<TrackCollection> TrC)
{
  double sumPt=0;
  for (TrackCollection::const_iterator tr=TrC->begin(); tr!=TrC->end(); tr++)
  {
    TLorentzVector ttr(tr->px(),tr->py(),tr->pz(),tr->p());
    if(ttr.DeltaR(t)<0.3 && fabs(tr->eta())<2.4) sumPt+=tr->pt();
  }
  sumPt-=t.Pt();
  return sumPt;
}

double isolPt_gmu(TLorentzVector t, edm::Handle<TrackCollection> TrC)
{
  double sumPt=0, dR=0.3, pt_match=0;
  for (TrackCollection::const_iterator tr=TrC->begin(); tr!=TrC->end(); tr++)
  {
    TLorentzVector ttr(tr->px(),tr->py(),tr->pz(),tr->p());
    if(ttr.DeltaR(t)<0.3 && fabs(tr->eta())<2.4) sumPt+=tr->pt();
    if(ttr.DeltaR(t)<dR && fabs(t.Pt()-ttr.Pt())/t.Pt()<0.5) {dR=ttr.DeltaR(t); pt_match=ttr.Pt();}
  }
  sumPt-=pt_match;
  return sumPt;
}

double isolPt_mc(TLorentzVector t, edm::Handle<GenParticleCollection> GenC)
{
  double sumPt=0;
  for (GenParticleCollection::const_iterator gen=GenC->begin(); gen!=GenC->end(); gen++)
  {
    int ndau = gen->numberOfDaughters();
    if(fabs(gen->eta())<2.4 && gen->charge()!=0 && ndau==0)
    {
      TLorentzVector tgen(gen->px(),gen->py(),gen->pz(),gen->p());
      if(tgen.DeltaR(t)<0.3) sumPt+=gen->pt();
    }
  }
  sumPt-=t.Pt();
  return sumPt;
}

//---------------exclude matched track with global muon---------------//
double isolPt2(TLorentzVector t, edm::Handle<TrackCollection> TrC, edm::Handle<TrackCollection> MuC)
{
  double sumPt=0;
  for (TrackCollection::const_iterator tr=TrC->begin(); tr!=TrC->end(); tr++)
  {
    TLorentzVector ttr(tr->px(),tr->py(),tr->pz(),tr->p());
    if(ttr.DeltaR(t)<0.3 && fabs(tr->eta())<2.4 && ttr.DeltaR(t)!=0 && (ttr.Pt()-t.Pt())!=0)
    {
      sumPt+=tr->pt();
      for (TrackCollection::const_iterator gm=MuC->begin(); gm!=MuC->end(); gm++)
      {
        TLorentzVector tgm(gm->px(),gm->py(),gm->pz(),gm->p());
//        if(ttr.DeltaR(tgm)<0.01 && fabs(gm->eta())<2.4)  sumPt-=ttr.Pt();
        if(ttr.DeltaR(tgm)<0.01 && fabs(gm->eta())<2.4 && gm->pt()>10) sumPt-=ttr.Pt();
      }
    }
  }
  return sumPt;
}

double isolPt2_gmu(TLorentzVector t, edm::Handle<TrackCollection> TrC, edm::Handle<TrackCollection> MuC)
{
  double sumPt=0;
  for (TrackCollection::const_iterator tr=TrC->begin(); tr!=TrC->end(); tr++)
  {
    TLorentzVector ttr(tr->px(),tr->py(),tr->pz(),tr->p());
    if(ttr.DeltaR(t)<0.3 && fabs(tr->eta())<2.4 && ttr.DeltaR(t)>0.01)
    {
      sumPt+=tr->pt();
      for (TrackCollection::const_iterator gm=MuC->begin(); gm!=MuC->end(); gm++)
      {
        TLorentzVector tgm(gm->px(),gm->py(),gm->pz(),gm->p());
//        if(ttr.DeltaR(tgm)<0.01 && fabs(gm->eta())<2.4)  sumPt-=ttr.Pt();
        if(ttr.DeltaR(tgm)<0.01 && fabs(gm->eta())<2.4 && gm->pt()>10) sumPt-=ttr.Pt();
      }
    }
  }
  return sumPt;
}

double isolPt2_mc(TLorentzVector t, edm::Handle<GenParticleCollection> GenC)
{
  double sumPt=0;
  for (GenParticleCollection::const_iterator gen=GenC->begin(); gen!=GenC->end(); gen++)
  {
    int ndau = gen->numberOfDaughters();
    if(fabs(gen->eta())<2.4 && gen->charge()!=0 && ndau==0)
    {
      TLorentzVector tgen(gen->px(),gen->py(),gen->pz(),gen->p());
      if(tgen.DeltaR(t)<0.3 && tgen.DeltaR(t)!=0 && (tgen.Pt()-t.Pt())!=0)
      {
        const Candidate * Mom  = gen->mother();
        const Candidate * Mom2 = 0;
        int Momid=0, Mom2id=0;
        if(Mom!=0) {Momid = Mom->pdgId(); Mom2 = (*Mom).mother();}
        if(Mom!=0 && Mom2!=0) Mom2id = Mom2->pdgId();
        if(Mom2id!=9900041&&Mom2id!=-9900041&&Mom2id!=9900042&&Mom2id!=-9900042) sumPt+=gen->pt();
      }
    }
  }
  return sumPt;
}

//---------------jet isolation---------------//
bool jetisol(TLorentzVector t, edm::Handle<CaloJetCollection> cJets)
{
  bool isol=true;
  for( CaloJetCollection::const_iterator cor = cJets->begin(); cor != cJets->end(); ++cor )
  {
    TLorentzVector tcor(cor->px(),cor->py(),cor->pz(),cor->p());
    if(cor->pt()>70 && t.DeltaR(tcor)<0.5) isol=false;
  }
  return isol;
}

//---------------distance between track line & vertex in xy plane---------------//
double distance(double x0, double y0, double x1, double y1, double phi) 
{
  return fabs((x1-x0)*tan(phi)-y1+y0)/sqrt(tan(phi)*tan(phi)+1);
}

//---------------smallest distance between track line & vertex in xy plane---------------//
double length(const Track *tr, edm::Handle<reco::VertexCollection> VtC)
//double length(const Track *tr, Handle<VertexCollection> VtC)
{
  double d1=100, d2=100;
  for( VertexCollection::const_iterator vt = VtC->begin(); vt != VtC->end();++vt )
  {
    d2=distance(tr->vx(),tr->vy(),vt->x(),vt->y(),tr->phi());
    if(d1>d2) d1=d2;
  }
  return d1;
}

double length_cor(const reco::Track *tr, edm::Handle<VertexCollection> VtC, double beamx, double beamy)
{
  double d1=100, d2=100;
  for( VertexCollection::const_iterator vt = VtC->begin(); vt != VtC->end();++vt )
  {
    d2=distance(tr->vx()+beamx,tr->vy()+beamy,vt->x(),vt->y(),tr->phi());
    if(d1>d2) d1=d2;
  }
  return d1;
}

double lengthz(const Track *tr, edm::Handle<VertexCollection> VtC)
{
  double d1=100, d2=100;
  for( VertexCollection::const_iterator vt = VtC->begin(); vt != VtC->end();++vt )
  {
    d2=tr->vz()-vt->z(); if(d2<0) d2=-d2;
    if(d1>d2) d1=d2;
  }
  return d1;
}

//---------------distance between reference position & vertex in 3D---------------//
TLorentzVector rvd(const Track *tr, edm::Handle<VertexCollection> VtC)
{
  double d1=100;
  TLorentzVector Ttr, ttr(tr->vx(),tr->vy(),tr->vz(),1);
  for( VertexCollection::const_iterator vt = VtC->begin(); vt != VtC->end();++vt )
        {
    TLorentzVector tvt(vt->x(),vt->y(),vt->z(),1);
    if(d1>(ttr-tvt).P()) {d1=(ttr-tvt).P(); Ttr=ttr-tvt;}
  }
  return Ttr;
}

//---------------distance between reference positions in 3D---------------//
TLorentzVector rrd(const Track *tr1, const Track *tr2)
{
  TLorentzVector ttr1(tr1->vx(),tr1->vy(),tr1->vz(),1);
  TLorentzVector ttr2(tr2->vx(),tr2->vy(),tr2->vz(),1);
  return ttr1-ttr2;
}

//---------------matching track with MC muon---------------//
bool matmcmu(TLorentzVector *mcmu, const Track *tr)
//bool matmcmu( &mcmu, const reco::Track *tr)
{
  TLorentzVector ttr(tr->px(),tr->py(),tr->pz(),tr->p());
  bool trmcmu=false;
  for(int i=0;i<4;i++) if(mcmu[i].Pt()>10 && fabs(mcmu[i].Eta())<2.4 && mcmu[i].DeltaR(ttr)<0.01 && ((i<=1&&tr->charge()==1)||(i>=2&&tr->charge()==-1))) trmcmu=true;
  return trmcmu;
}

//---------------vector calulation---------------//
double vdot(TLorentzVector t1, TLorentzVector t2) //vector dot
{
  return t1.Px()*t2.Px()+t1.Py()*t2.Py()+t1.Pz()*t2.Pz();
}
 
TLorentzVector vcross(TLorentzVector t1, TLorentzVector t2) //vector cross
{
        TLorentzVector t3;
        t3.SetPxPyPzE(t1.Py()*t2.Pz()-t1.Pz()*t2.Py(), t1.Pz()*t2.Px()-t1.Px()*t2.Pz(), t1.Px()*t2.Py()-t1.Py()*t2.Px(), 1);
        return t3;
}

//---------------impact point---------------//
TLorentzVector D(const Track *tr, edm::Handle<VertexCollection> VtC) //impact vector
{
        double th=tr->theta(), phi=tr->phi();
        TLorentzVector t(sin(th)*cos(phi),sin(th)*sin(phi),cos(th));
        TLorentzVector p0(tr->vx(),tr->vy(),tr->vz(),1);
        TLorentzVector p2, N, D, Dh;
        double l=100, d=100;
        for( VertexCollection::const_iterator vt = VtC->begin(); vt != VtC->end();++vt )
        {
                TLorentzVector p1(vt->x(),vt->y(),vt->z(),1);
                p2=p1-p0;
                N=vcross(p2,t);
                Dh = vcross(N*(1/N.P()),t);
                d = -vdot(Dh,p2);
                if(l>d) {l=d; D=d*Dh;}
        }
        return D;
}

TLorentzVector ip(const Track *tr, edm::Handle<VertexCollection> VtC) //impact point
{
        double th=tr->theta(), phi=tr->phi();
        TLorentzVector t(sin(th)*cos(phi),sin(th)*sin(phi),cos(th));
        TLorentzVector p0(tr->vx(),tr->vy(),tr->vz(),1);
        TLorentzVector p2, N, D, Dh, tip;
        double l=100, d=100;
        for( VertexCollection::const_iterator vt = VtC->begin(); vt != VtC->end();++vt )
        {
                TLorentzVector p1(vt->x(),vt->y(),vt->z(),1);
                p2=p1-p0;
                N=vcross(p2,t);
                Dh = vcross(N*(1/N.P()),t);
                d = -vdot(Dh,p2);
                D = d*Dh;
                if(l>d) {l=d; tip = p1+D;}
        }
        return tip;
}

//---------------distance between track lines---------------//
TLorentzVector dt(const Track *tr1, const Track *tr2)
{
        double th1=tr1->theta(), phi1=tr1->phi(), th2=tr2->theta(), phi2=tr2->phi();
        TLorentzVector t1(sin(th1)*cos(phi1),sin(th1)*sin(phi1),cos(th1)), t2(sin(th2)*cos(phi2),sin(th2)*sin(phi2),cos(th2));
        TLorentzVector p1(tr1->vx(),tr1->vy(),tr1->vz(),1), p2(tr2->vx(),tr2->vy(),tr2->vz(),1);
  TLorentzVector p3=p2-p1;
  double t3=vdot(t1,t2);
  double a = ( vdot(p3,t1)-vdot(p3,t2)*t3 )/( 1-t3*t3 );
  double b = ( vdot(p3,t1)*t3-vdot(p3,t2) )/( 1-t3*t3 );
  TLorentzVector n = p3+b*t2-a*t1;
  return n;
}

//---------------closest point between track lines---------------//
TLorentzVector cp1(const Track *tr1, const Track *tr2) // track 1
{
        double th1=tr1->theta(), phi1=tr1->phi(), th2=tr2->theta(), phi2=tr2->phi();
        TLorentzVector t1(sin(th1)*cos(phi1),sin(th1)*sin(phi1),cos(th1)), t2(sin(th2)*cos(phi2),sin(th2)*sin(phi2),cos(th2));
        TLorentzVector p1(tr1->vx(),tr1->vy(),tr1->vz(),1), p2(tr2->vx(),tr2->vy(),tr2->vz(),1);
  TLorentzVector p3=p2-p1;
  double t3=vdot(t1,t2);
  double a = ( vdot(p3,t1)-vdot(p3,t2)*t3 )/( 1-t3*t3 );
//  double b = ( vdot(p3,t1)*t3-vdot(p3,t2) )/( 1-t3*t3 );
//  TLorentzVector n = p3+b*t2-a*t1;
  TLorentzVector q1 = p1+a*t1;
  return q1;
}

TLorentzVector cp2(const Track *tr1, const Track *tr2) // track 2
{
        double th1=tr1->theta(), phi1=tr1->phi(), th2=tr2->theta(), phi2=tr2->phi();
        TLorentzVector t1(sin(th1)*cos(phi1),sin(th1)*sin(phi1),cos(th1)), t2(sin(th2)*cos(phi2),sin(th2)*sin(phi2),cos(th2));
        TLorentzVector p1(tr1->vx(),tr1->vy(),tr1->vz(),1), p2(tr2->vx(),tr2->vy(),tr2->vz(),1);
  TLorentzVector p3=p2-p1;
  double t3=vdot(t1,t2);
//  double a = ( vdot(p3,t1)-vdot(p3,t2)*t3 )/( 1-t3*t3 );
  double b = ( vdot(p3,t1)*t3-vdot(p3,t2) )/( 1-t3*t3 );
//  TLorentzVector n = p3+b*t2-a*t1;
  TLorentzVector q2 = p2+b*t2;
  return q2;
}

//---------------distance between recoVertex & MC muon vertex---------------//
TLorentzVector dvt(const Candidate *mc, edm::Handle<VertexCollection> VtC)
{
  TLorentzVector Dvt, temp, mcmuvt(mc->vx(),mc->vy(),mc->vz());
  double d1=100;
  for( VertexCollection::const_iterator vt = VtC->begin(); vt != VtC->end();++vt )
  {
    TLorentzVector tvt(vt->x(),vt->y(),vt->z());
    temp=mcmuvt-tvt;
    if(d1>temp.P()) {d1=temp.P(); Dvt=temp;}
  }
  return Dvt;
}

TLorentzVector dvt_BS(const Candidate *mc, edm::Handle<VertexCollection> VtC, double x0, double y0, double z0) //without samed recoVertex with BeamSpot
{
  TLorentzVector Dvt(0,0,0), temp, mcmuvt(mc->vx(),mc->vy(),mc->vz()), BS(x0,y0,z0);
  double d1=100;
  for( VertexCollection::const_iterator vt = VtC->begin(); vt != VtC->end();++vt )
  {
    TLorentzVector tvt(vt->x(),vt->y(),vt->z());
    temp=mcmuvt-tvt;
    if(d1>temp.P() && !(vt->x()==x0 && vt->y()==y0 && vt->z()==z0)) {d1=temp.P(); Dvt=temp;}
  }
  return Dvt;
}

TLorentzVector dvt_cor(const Candidate *mc, edm::Handle<VertexCollection> VtC, double beamx, double beamy)
{
  TLorentzVector Dvt, temp, tvt, mcmuvt(mc->vx()+beamx,mc->vy()+beamy,mc->vz());
  double d1=100;
  for( VertexCollection::const_iterator vt = VtC->begin(); vt != VtC->end();++vt )
  {
    TLorentzVector tvt(vt->x(),vt->y(),vt->z());
    temp=mcmuvt-tvt;
    if(d1>temp.P()) {d1=temp.P(); Dvt=temp;}
  }
  return Dvt;
}

TLorentzVector tV(edm::Handle<reco::VertexCollection> VtC)
{ 
  TLorentzVector ttV;
  for( reco::VertexCollection::const_iterator vt = VtC->begin(); vt != VtC->end();vt++ )
  {
    TLorentzVector temp(vt->x(),vt->y(),vt->z());
    ttV=temp;
  }
  return ttV;
}
TLorentzVector tT(edm::Handle<TrackCollection> TrC)
{ 
  TLorentzVector ttT;
  for( TrackCollection::const_iterator tr = TrC->begin(); tr != TrC->end();++tr )
  {
    TLorentzVector temp(tr->px(),tr->py(),tr->pz(),tr->p());
    ttT=temp;
  }
  return ttT;
}


HppMuMuAnalyzer::HppMuMuAnalyzer(const edm::ParameterSet& iConfig):
  Tracks_(iConfig.getUntrackedParameter<edm::InputTag>("Tracks")),
  dzcut_(iConfig.getUntrackedParameter<double>("dzcut")),
  min(iConfig.getUntrackedParameter<double>("min")),
  max(iConfig.getUntrackedParameter<double>("max")),
  minvtx(iConfig.getUntrackedParameter<double>("minvtx")),
  maxvtx(iConfig.getUntrackedParameter<double>("maxvtx")),
  minvty(iConfig.getUntrackedParameter<double>("minvty")),
  maxvty(iConfig.getUntrackedParameter<double>("maxvty")),
  minvt2x(iConfig.getUntrackedParameter<double>("minvt2x")),
  maxvt2x(iConfig.getUntrackedParameter<double>("maxvt2x")),
  minvt2y(iConfig.getUntrackedParameter<double>("minvt2y")),
  maxvt2y(iConfig.getUntrackedParameter<double>("maxvt2y")),
  minvt3x(iConfig.getUntrackedParameter<double>("minvt3x")),
  maxvt3x(iConfig.getUntrackedParameter<double>("maxvt3x")),
  minvt3y(iConfig.getUntrackedParameter<double>("minvt3y")),
  maxvt3y(iConfig.getUntrackedParameter<double>("maxvt3y")),
  minmcvtx(iConfig.getUntrackedParameter<double>("minmcvtx")),
  maxmcvtx(iConfig.getUntrackedParameter<double>("maxmcvtx")),
  minmcvty(iConfig.getUntrackedParameter<double>("minmcvty")),
  maxmcvty(iConfig.getUntrackedParameter<double>("maxmcvty")),
  vtx_cor(iConfig.getUntrackedParameter<double>("vtx_cor")),
  vty_cor(iConfig.getUntrackedParameter<double>("vty_cor")),
  //JetCorrectionService_(iConfig.getParameter<std::string>("JetCorrectionService")),
  CaloJetAlgorithm_(iConfig.getParameter<std::string>("CaloJetAlgorithm"))
  //GenJetAlgorithm_(iConfig.getParameter<std::string>("GenJetAlgorithm"))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
//---------------jet---------------//
  hncalojet = fs->make<TH1F>("hncalojet","# of calojet per event",50,0,50);
  hncalojet_ptcut = fs->make<TH1F>("hncalojet_ptcut","# of calojet per event",10,0,10);
  hcalojetpt = fs->make<TH1F>("hcalojetpt","pt of calojet",50,0,100);
  hcalojetpt_ptcut = fs->make<TH1F>("hcalojetpt_ptcut","pt of calojet",50,0,100);
  hcalojetpt_ptcut_ex = fs->make<TH1F>("hcalojetpt_ptcut_ex","pt of calojet",50,0,300);
  hcalojeteta = fs->make<TH1F>("hcalojeteta","#eta of calojet",50,-5,5);
  hcalojeteta_ptcut = fs->make<TH1F>("hcalojeteta_ptcut","#eta of calojet",50,-5,5);

  hncorjet = fs->make<TH1F>("hncorjet","# of corjet per event",50,0,50);
  hncorjet_ptcut = fs->make<TH1F>("hncorjet_ptcut","# of corjet per event",10,0,10);
  hncorjet_ptmucut = fs->make<TH1F>("hncorjet_ptmucut","# of corjet per event",10,0,10);
  hcorjetpt = fs->make<TH1F>("hcorjetpt","pt of corjet",50,0,100);
  hcorjetpt_ptcut = fs->make<TH1F>("hcorjetpt_ptcut","pt of corjet",50,0,100);
  hcorjetpt_ptcut_ex = fs->make<TH1F>("hcorjetpt_ptcut_ex","pt of corjet",50,0,300);
  hcorjetpt_ptmucut = fs->make<TH1F>("hcorjetpt_ptmucut","pt of corjet",50,0,100);
  hcorjetpt_ptmucut_ex = fs->make<TH1F>("hcorjetpt_ptmucut_ex","pt of corjet",50,0,300);
  hcorjeteta = fs->make<TH1F>("hcorjeteta","#eta of corjet",50,-5,5);
  hcorjeteta_ptcut = fs->make<TH1F>("hcorjeteta_ptcut","#eta of corjet",50,-5,5);
  hcorjeteta_ptmucut = fs->make<TH1F>("hcorjeteta_ptmucut","#eta of corjet",50,-5,5);

  hcormcmuptdR_ptetacut = fs->make<TH2F>("hcormcmuptdR_ptetacut","#DeltaR(cor,mu) VS P_{T,cor}/P_{T,mcmu}",200,0,2,200,0,4);
  hcormcmuptdR_30ptetacut = fs->make<TH2F>("hcormcmuptdR_30ptetacut","#DeltaR(cor,mu) VS P_{T,cor}/P_{T,mcmu}",200,0,2,200,0,4);

  hngenjet = fs->make<TH1F>("hngenjet","# of genjet per event",50,0,50);
  hngenjet_ptcut = fs->make<TH1F>("hngenjet_ptcut","# of genjet per event",10,0,10);
  hngenjet_ptmucut = fs->make<TH1F>("hngenjet_ptmucut","# of genjet per event",10,0,10);
  hgenjetpt = fs->make<TH1F>("hgenjetpt","pt of genjet",50,0,100);
  hgenjetpt_ptcut = fs->make<TH1F>("hgenjetpt_ptcut","pt of genjet",50,0,100);
  hgenjetpt_ptcut_ex = fs->make<TH1F>("hgenjetpt_ptcut_ex","pt of genjet",50,0,300);
  hgenjetpt_ptmucut = fs->make<TH1F>("hgenjetpt_ptmucut","pt of genjet",50,0,100);
  hgenjetpt_ptmucut_ex = fs->make<TH1F>("hgenjetpt_ptmucut_ex","pt of genjet",50,0,300);
  hgenjetpt_pt24etamucut = fs->make<TH1F>("hgenjetpt_pt24etamucut","pt of genjet",50,0,100);
  hgenjetpt_20pt24etamucut = fs->make<TH1F>("hgenjetpt_20pt24etamucut","pt of genjet",50,0,100);
  hgenjetpt_30pt24etamucut = fs->make<TH1F>("hgenjetpt_30pt24etamucut","pt of genjet",50,0,100);
  hgenjetpt_ptmucut_calojetmatch = fs->make<TH1F>("hgenjetpt_ptmucut_calojetmatch","pt of genjet",50,0,100);
  hgenjetpt_ptmucut_corjetmatch = fs->make<TH1F>("hgenjetpt_ptmucut_corjetmatch","pt of genjet",50,0,100);
  hgenjetpt_pt24etamucut_corjetmatch = fs->make<TH1F>("hgenjetpt_pt24etamucut_corjetmatch","pt of genjet",50,0,100);
  hgenjetpt_20pt24etamucut_corjetmatch = fs->make<TH1F>("hgenjetpt_20pt24etamucut_corjetmatch","pt of genjet",50,0,100);
  hgenjetpt_30pt24etamucut_corjetmatch = fs->make<TH1F>("hgenjetpt_30pt24etamucut_corjetmatch","pt of genjet",50,0,100);
  hgenjeteta = fs->make<TH1F>("hgenjeteta","#eta of genjet",50,-5,5);
  hgenjeteta_ptcut = fs->make<TH1F>("hgenjeteta_ptcut","#eta of genjet",50,-5,5);
  hgenjeteta_ptmucut = fs->make<TH1F>("hgenjeteta_ptmucut","#eta of genjet",50,-5,5);
  hgenjeteta_pt24etamucut = fs->make<TH1F>("hgenjeteta_pt24etamucut","#eta of genjet",50,-2.4,2.4);
  hgenjeteta_20pt24etamucut = fs->make<TH1F>("hgenjeteta_20pt24etamucut","#eta of genjet",50,-2.4,2.4);
  hgenjeteta_30pt24etamucut = fs->make<TH1F>("hgenjeteta_30pt24etamucut","#eta of genjet",50,-2.4,2.4);
  hgenjeteta_ptmucut_calojetmatch = fs->make<TH1F>("hgenjeteta_ptmucut_calojetmatch","#eta of genjet",50,-5,5);
  hgenjeteta_ptmucut_corjetmatch = fs->make<TH1F>("hgenjeteta_ptmucut_corjetmatch","#eta of genjet",50,-5,5);
  hgenjeteta_pt24etamucut_corjetmatch = fs->make<TH1F>("hgenjeteta_pt24etamucut_corjetmatch","#eta of genjet",50,-2.4,2.4);
  hgenjeteta_20pt24etamucut_corjetmatch = fs->make<TH1F>("hgenjeteta_20pt24etamucut_corjetmatch","#eta of genjet",50,-2.4,2.4);
  hgenjeteta_30pt24etamucut_corjetmatch = fs->make<TH1F>("hgenjeteta_30pt24etamucut_corjetmatch","#eta of genjet",50,-2.4,2.4);

  hgenmcmuptdR_ptetacut = fs->make<TH2F>("hgenmcmuptdR_ptetacut","#DeltaR(gen,mcmu) VS P_{T,gen}/P_{T,mcmu}",200,0,2,200,0,4);
  hgenmcmuptdR_30ptetacut = fs->make<TH2F>("hgenmcmuptdR_30ptetacut","#DeltaR(gen,mcmu) VS P_{T,gen}/P_{T,mcmu}",200,0,2,200,0,4);

  hjetresponseVSpt_24etamucut = fs->make<TH1F>("hjetresponseVSpt_24etamucut","Mean(corPt/genPt) VS genPt",50,0,100);
  hjetresponseVSeta_24etamucut = fs->make<TH1F>("hjetresponseVSeta_24etamucut","Mean(corPt/genPt) VS gen #eta",50,-2.4,2.4);
  hjetresolutionVSpt_24etamucut = fs->make<TH1F>("hjetresolutionVSpt_24etamucut","Sigma(corPt/genPt) VS genPt",50,0,100);
  hjetresolutionVSeta_24etamucut = fs->make<TH1F>("hjetresolutionVSeta_24etamucut","Sigma(corPt/genPt) VS gen #eta",50,-2.4,2.4);
  hjetresponseVSpt_30pt24etamucut = fs->make<TH1F>("hjetresponseVSpt_30pt24etamucut","Mean(corPt/genPt) VS genPt",50,0,100);
  hjetresponseVSeta_30pt24etamucut = fs->make<TH1F>("hjetresponseVSeta_30pt24etamucut","Mean(corPt/genPt) VS gen #eta",50,-2.4,2.4);
  hjetresolutionVSpt_30pt24etamucut = fs->make<TH1F>("hjetresolutionVSpt_30pt24etamucut","Sigma(corPt/genPt) VS genPt",50,0,100);
  hjetresolutionVSeta_30pt24etamucut = fs->make<TH1F>("hjetresolutionVSeta_30pt24etamucut","Sigma(corPt/genPt) VS gen #eta",50,-2.4,2.4);

  htrjetpt = fs->make<TH1F>("htrjetpt","pt of tracker jet",50,0,500);
  htrjetcorjetpt = fs->make<TH2F>("htrjetcorjetpt","tracker jet pt VS corjet pt",200,0,1000,200,0,1000);
  htrjetpt_mucut = fs->make<TH1F>("htrjetpt_mucut","pt of tracker jet except muon",200,0,1000);
  htrjetcorjetpt_mucut = fs->make<TH2F>("htrjetcorjetpt_mucut","tracker jet pt VS corjet pt",200,0,1000,200,0,1000);
//---------------track---------------//
  hntr = fs->make<TH1F>("hntr","number of tracks per event",200,0,200);
  hntr_ptetacut = fs->make<TH1F>("hntr_ptetacut","number of track pers event (P_{T}>10GeV, |#eta|<2.4)",10,0,10);
  hntr_ptetaisolcut = fs->make<TH1F>("hntr_ptetaisolcut",Form("number of tracks per event (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),10,0,10);
  hntr_ptetaisol2cut = fs->make<TH1F>("hntr_ptetaisol2cut",Form("number of tracks per event (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),10,0,10);
  hntr_ptetaisol2jcut = fs->make<TH1F>("hntr_ptetaisol2jcut",Form("number of tracks per event (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),10,0,10);
  hntr_i2j = fs->make<TH1F>("hntr_i2j","number of tracks per event when there are muon candidates",200,0,200);
  hntr_i2j3mu = fs->make<TH1F>("hntr_i2j3mu","number of tracks per event when there are 3 muon candidates",200,0,200);
  hntr_i2j4mu = fs->make<TH1F>("hntr_i2j4mu","number of tracks per event when there are 4 muon candidates",200,0,200);
  hntr_BS = fs->make<TH1F>("hntr_BS","number of tracks per event when recoVertex is equal to the BeamSpot",20,0,20);
  hntr_i2jBS = fs->make<TH1F>("hntr_i2jBS","number of muon candidates per event when recoVertex is equal to the BeamSpot",20,0,20);
  htrackpt = fs->make<TH1F>("htrackpt","P_{T} of track",50,0,1000);
  htrackpt_ptetacut = fs->make<TH1F>("htrackpt_ptetacut","P_{T} of track (P_{T}>10GeV, |#eta|<2.4)",50,0,1000);
  htrackpt_ptetaisolcut = fs->make<TH1F>("htrackpt_ptetaisolcut",Form("P_{T} of track (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,0,1000);
  htrackpt_ptetaisol2cut = fs->make<TH1F>("htrackpt_ptetaisol2cut",Form("P_{T} of track (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,0,1000);
  htrackpt_ptetaisol2jcut = fs->make<TH1F>("htrackpt_ptetaisol2jcut",Form("P_{T} of track (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,0,1000);
  htrackpt_mcmumatch = fs->make<TH1F>("htrackpt_mcmumatch","P_{T} of track (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{tr,MC#mu}<0.01)",50,0,1000);
  htrackpt_mcmumatchisol = fs->make<TH1F>("htrackpt_mcmumatchisol",Form("P_{T} of track (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{tr,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,0,1000);
  htrackpt_mcmumatchisol2 = fs->make<TH1F>("htrackpt_mcmumatchisol2",Form("P_{T} of track (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{tr,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,0,1000);
  htrackpt_mcmumatchisol2j = fs->make<TH1F>("htrackpt_mcmumatchisol2j",Form("P_{T} of track (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{tr,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,0,1000);
  htracketa = fs->make<TH1F>("htracketa","#eta of track",50,-2.4,2.4);
  htracketa_ptetacut = fs->make<TH1F>("htracketa_ptetacut","#eta of track (P_{T}>10GeV, |#eta|<2.4)",50,-2.4,2.4);
  htracketa_ptetaisolcut = fs->make<TH1F>("htracketa_ptetaisolcut",Form("#eta of track (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  htracketa_ptetaisol2cut = fs->make<TH1F>("htracketa_ptetaisol2cut",Form("#eta of track (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  htracketa_ptetaisol2jcut = fs->make<TH1F>("htracketa_ptetaisol2jcut",Form("#eta of track (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  htracketa_mcmumatch = fs->make<TH1F>("htracketa_mcmumatch","#eta of track (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{tr,MC#mu}<0.01)",50,-2.4,2.4);
  htracketa_mcmumatchisol = fs->make<TH1F>("htracketa_mcmumatchisol",Form("#eta of track (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{tr,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  htracketa_mcmumatchisol2 = fs->make<TH1F>("htracketa_mcmumatchisol2",Form("#eta of track (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{tr,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  htracketa_mcmumatchisol2j = fs->make<TH1F>("htracketa_mcmumatchisol2j",Form("#eta of track (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{tr,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  htrackisolPt_ptetacut01 = fs->make<TH1F>("htrackisolPt_ptetacut01","isolP_{T} of track (P_{T}>10GeV, |#eta|<2.4)",50,0,0.1);
  htrackisolPt_ptetacut1 = fs->make<TH1F>("htrackisolPt_ptetacut1","isolP_{T} of track (P_{T}>10GeV, |#eta|<2.4)",50,0,1);
  htrackisolPt_ptetacut10 = fs->make<TH1F>("htrackisolPt_ptetacut10","isolP_{T} of track (P_{T}>10GeV, |#eta|<2.4)",50,0,10);
  htrackisolPt_ptetacut100 = fs->make<TH1F>("htrackisolPt_ptetacut100","isolP_{T} of track (P_{T}>10GeV, |#eta|<2.4)",50,0,100);
  htrackisolPt_ptetacut1000 = fs->make<TH1F>("htrackisolPt_ptetacut1000","isolP_{T} of track (P_{T}>10GeV, |#eta|<2.4)",50,0,1000);
  htrackisolPt2_ptetacut10 = fs->make<TH1F>("htrackisolPt2_ptetacut10","isolP_{T} of track (P_{T}>10GeV, |#eta|<2.4)",50,0,10);
  htrackisolPt2_ptetacut100 = fs->make<TH1F>("htrackisolPt2_ptetacut100","isolP_{T} of track (P_{T}>10GeV, |#eta|<2.4)",50,0,100);
  h2trdeltaphi = fs->make<TH1F>("h2trdeltaphi","#Delta#phi between tr+tr+ & tr-tr-",50,0,TMath::Pi());
  h2trdeltaphi_isol = fs->make<TH1F>("h2trdeltaphi_isol","#Delta#phi between tr+tr+ & tr-tr-",50,0,TMath::Pi());
  h2trdeltaphi_isol2 = fs->make<TH1F>("h2trdeltaphi_isol2","#Delta#phi between tr+tr+ & tr-tr-",50,0,TMath::Pi());
  h2tr2mcmudphi_isol2 = fs->make<TH2F>("h2tr2mcmudphi_isol2","#Delta#phi 2 tr VS 2 MC muon",200,0,TMath::Pi(),200,0,TMath::Pi());

  htrpt_i2junv = fs->make<TH1F>("htrpt_i2junv","P_{T} of track",50,0,1000);

  htrackz = fs->make<TH1F>("htrackz","z position of tracks",200,-15,15);
  htrackz2 = fs->make<TH1F>("htrackz2","z position of tracks",200,-15,15);
  htrackxy = fs->make<TH1F>("htrackxy","xy of tracks",200,-0.1,0.1);
  htracksz = fs->make<TH1F>("htracksz","sz of tracks",200,-15,15);
  htrackz_ptetacut = fs->make<TH1F>("htrackz_ptetacut","z position of tracks",200,-15,15);
  htrackz_ptetaisol2cut = fs->make<TH1F>("htrackz_ptetaisol2cut","z position of tracks",200,-15,15);
  htrackz_ptetaisol2jcut = fs->make<TH1F>("htrackz_ptetaisol2jcut","z position of tracks",200,-15,15);
  htrackz2_ptetaisol2jcut = fs->make<TH1F>("htrackz2_ptetaisol2jcut","z position of tracks",200,-15,15);
  htrackxy_ptetaisol2jcut = fs->make<TH1F>("htrackxy_ptetaisol2jcut","xy of tracks",200,-0.1,0.1);
  htracksz_ptetaisol2jcut = fs->make<TH1F>("htracksz_ptetaisol2jcut","sz of tracks",200,-15,15);
  htrackdz = fs->make<TH1F>("htrackdz","#Deltaz between tracks",200,-0.5,0.5);
  htrackdxy = fs->make<TH1F>("htrackdxy","#Deltaxy between tracks",200,-0.1,0.1);
  htrackdsz = fs->make<TH1F>("htrackdsz","#Deltasz between tracks",200,-10,10);
  htrackdz_ptetaisol2jcut = fs->make<TH1F>("htrackdz_ptetaisol2jcut","#Deltaz between tracks",200,-0.5,0.5);
  htrackdxy_ptetaisol2jcut = fs->make<TH1F>("htrackdxy_ptetaisol2jcut","#Deltaxy between tracks",200,-0.1,0.1);
  htrackdsz_ptetaisol2jcut = fs->make<TH1F>("htrackdsz_ptetaisol2jcut","#Deltasz between tracks",200,-10,10);
  h2trdz_isol2j = fs->make<TH1F>("h2trdz_isol2j","#Deltaz between two highest pt tracks",200,-0.05,0.05);
  h2trdz_isol2j_ex = fs->make<TH1F>("h2trdz_isol2j_ex","#Deltaz between two highest pt tracks",200,-0.25,0.25);
  htrackX = fs->make<TH1F>("htrackX","x position of tracks",200,-0.1,0.1);
  htrackY = fs->make<TH1F>("htrackY","y position of tracks",200,-0.1,0.1);
  htrackZ = fs->make<TH1F>("htrackZ","z position of tracks",200,-15,15);
  htrackXY2D = fs->make<TH2F>("htrackXY2D","x&y position of tracks",500,-0.5,0.5,500,-0.5,0.5);
  htrackdX = fs->make<TH1F>("htrackdX","dx between tracks",200,-0.1,0.1);
  htrackdY = fs->make<TH1F>("htrackdY","dy between tracks",200,-0.1,0.1);
  htrackdZ = fs->make<TH1F>("htrackdZ","dz between tracks",200,-0.5,0.5);
  htrackX_ptetaisol2jcut = fs->make<TH1F>("htrackX_ptetaisol2jcut","x position of tracks",200,-0.1,0.1);
  htrackY_ptetaisol2jcut = fs->make<TH1F>("htrackY_ptetaisol2jcut","y position of tracks",200,-0.1,0.1);
  htrackZ_ptetaisol2jcut = fs->make<TH1F>("htrackZ_ptetaisol2jcut","z position of tracks",200,-15,15);
  htrackdX_ptetaisol2jcut = fs->make<TH1F>("htrackdX_ptetaisol2jcut","dx between tracks",200,-0.1,0.1);
  htrackdY_ptetaisol2jcut = fs->make<TH1F>("htrackdY_ptetaisol2jcut","dy between tracks",200,-0.1,0.1);
  htrackdZ_ptetaisol2jcut = fs->make<TH1F>("htrackdZ_ptetaisol2jcut","dz between tracks",200,-0.05,0.05);
  htrackdZ_ptetaisol2jcut_ex = fs->make<TH1F>("htrackdZ_ptetaisol2jcut_ex","dz between tracks",200,-0.25,0.25);
  h2trdZ_isol2j = fs->make<TH1F>("h2trdZ_isol2j","#DeltaZ between two highest pt tracks",200,-0.05,0.05);
  h2trdZ_isol2j_ex = fs->make<TH1F>("h2trdZ_isol2j_ex","#DeltaZ between two highest pt tracks",200,-0.25,0.25);
  htrackXY2D_ptetaisol2jcut = fs->make<TH2F>("htrackXY2D_ptetaisol2jcut","x&y position of tracks in xy plane",500,-0.1,0.1,500,-0.1,0.1);
  htrackXY2D_ptetaisol2jcut_beamP = fs->make<TH2F>("htrackXY2D_ptetaisol2jcut_beamP","x&y position of tracks in xy plane",500,-0.1,0.1,500,-0.1,0.1);
  htrackXY2D_ptetaisol2jcut_beamM = fs->make<TH2F>("htrackXY2D_ptetaisol2jcut_beamM","x&y position of tracks in xy plane",500,-0.1,0.1,500,-0.1,0.1);
  hbeamXY2D = fs->make<TH2F>("hbeamXY2D","x&y position of beams in xy plane",500,-0.1,0.1,500,-0.1,0.1);
  hbeamz = fs->make<TH1F>("hbeamz","z position of beams",200,-0.5,0.5);
  hbeamsigmaz = fs->make<TH1F>("hbeamsigmaz","sigmaz of beams",200,-5,5);

  hHPPmass_tr = fs->make<TH1F>("hHPPmass_tr","Mass of H++",100,min,max);
  hHMMmass_tr = fs->make<TH1F>("hHMMmass_tr","Mass of H--",100,min,max);
  hHPPHMMmass_tr = fs->make<TH1F>("hHPPHMMmass_tr","Mass of H++&H--",100,0,2000);
  hHPPmass_trisol = fs->make<TH1F>("hHPPmass_trisol","Mass of H++",100,min,max);
  hHMMmass_trisol = fs->make<TH1F>("hHMMmass_trisol","Mass of H--",100,min,max);
  hHPPHMMmass_trisol = fs->make<TH1F>("hHPPHMMmass_trisol","Mass of H++&H--",100,0,2000);
  hHPPmass_trisol2 = fs->make<TH1F>("hHPPmass_trisol2","Mass of H++",100,min,max);
  hHMMmass_trisol2 = fs->make<TH1F>("hHMMmass_trisol2","Mass of H--",100,min,max);
  hHPPHMMmass_trisol2 = fs->make<TH1F>("hHPPHMMmass_trisol2","Mass of H++&H--",100,0,2000);
  hHPPmass_trisol2j = fs->make<TH1F>("hHPPmass_trisol2j","Mass of H++",100,min,max);
  hHMMmass_trisol2j = fs->make<TH1F>("hHMMmass_trisol2j","Mass of H--",100,min,max);
  hHPPHMMmass_trisol2j = fs->make<TH1F>("hHPPHMMmass_trisol2j","Mass of H++&H--",100,0,2000);
  hHPPmass_trisol2dphi = fs->make<TH1F>("hHPPmass_trisol2dphi","Mass of H++",100,min,max);
  hHMMmass_trisol2dphi = fs->make<TH1F>("hHMMmass_trisol2dphi","Mass of H--",100,min,max);
  hHPPHMMmass_trisol2dphi = fs->make<TH1F>("hHPPHMMmass_trisol2dphi","Mass of H++&H--",100,0,2000);
  hHPPmass_trisol2jdz = fs->make<TH1F>("hHPPmass_trisol2jdz","Mass of H++",100,min,max);
  hHMMmass_trisol2jdz = fs->make<TH1F>("hHMMmass_trisol2jdz","Mass of H--",100,min,max);
  hHPPmass_trisol2jdz_ex = fs->make<TH1F>("hHPPmass_trisol2jdz_ex","Mass of H++",200,0,700);
  hHMMmass_trisol2jdz_ex = fs->make<TH1F>("hHMMmass_trisol2jdz_ex","Mass of H--",200,0,700);
  hHPPmass_trisol2jv = fs->make<TH1F>("hHPPmass_trisol2jv","Mass of H++",100,min,max);
  hHMMmass_trisol2jv = fs->make<TH1F>("hHMMmass_trisol2jv","Mass of H--",100,min,max);
  hHPPHMMmass_trisol2jv = fs->make<TH1F>("hHPPHMMmass_trisol2jv","Mass of H++&H--",100,0,2000);
  hHPPmass_trisol2jvc = fs->make<TH1F>("hHPPmass_trisol2jvc","Mass of H++",100,min,max);
  hHMMmass_trisol2jvc = fs->make<TH1F>("hHMMmass_trisol2jvc","Mass of H--",100,min,max);
  hHPPHMMmass_trisol2jvc = fs->make<TH1F>("hHPPHMMmass_trisol2jvc","Mass of H++&H--",100,0,2000);
  hHPPmass_trisol2jvs = fs->make<TH1F>("hHPPmass_trisol2jvs","Mass of H++",100,min,max);
  hHMMmass_trisol2jvs = fs->make<TH1F>("hHMMmass_trisol2jvs","Mass of H--",100,min,max);
  hHPPHMMmass_trisol2jvs = fs->make<TH1F>("hHPPHMMmass_trisol2jvs","Mass of H++&H--",100,0,2000);
  hHPPmass_trisol2jvcs = fs->make<TH1F>("hHPPmass_trisol2jvcs","Mass of H++",100,min,max);
  hHMMmass_trisol2jvcs = fs->make<TH1F>("hHMMmass_trisol2jvcs","Mass of H--",100,min,max);
  hHPPHMMmass_trisol2jvcs = fs->make<TH1F>("hHPPHMMmass_trisol2jvcs","Mass of H++&H--",100,0,2000);

  hHPPmass_trmcmumatch = fs->make<TH1F>("hHPPmass_trmcmumatch","Mass of H++",100,min,max);
  hHMMmass_trmcmumatch = fs->make<TH1F>("hHMMmass_trmcmumatch","Mass of H--",100,min,max);
  hHPPHMMmass_trmcmumatch = fs->make<TH1F>("hHPPHMMmass_trmcmumatch","Mass of H++&H--",100,0,2000);
  hHPPmass_trmcmumatchisol = fs->make<TH1F>("hHPPmass_trmcmumatchisol","Mass of H++",100,min,max);
  hHMMmass_trmcmumatchisol = fs->make<TH1F>("hHMMmass_trmcmumatchisol","Mass of H--",100,min,max);
  hHPPHMMmass_trmcmumatchisol = fs->make<TH1F>("hHPPHMMmass_trmcmumatchisol","Mass of H++&H--",100,0,2000);
  hHPPmass_trmcmumatchisol2 = fs->make<TH1F>("hHPPmass_trmcmumatchisol2","Mass of H++",100,min,max);
  hHMMmass_trmcmumatchisol2 = fs->make<TH1F>("hHMMmass_trmcmumatchisol2","Mass of H--",100,min,max);
  hHPPHMMmass_trmcmumatchisol2 = fs->make<TH1F>("hHPPHMMmass_trmcmumatchisol2","Mass of H++&H--",100,0,2000);
  hHPPmass_trmcmumatchisol2j = fs->make<TH1F>("hHPPmass_trmcmumatchisol2j","Mass of H++",100,min,max);
  hHMMmass_trmcmumatchisol2j = fs->make<TH1F>("hHMMmass_trmcmumatchisol2j","Mass of H--",100,min,max);
  hHPPHMMmass_trmcmumatchisol2j = fs->make<TH1F>("hHPPHMMmass_trmcmumatchisol2j","Mass of H++&H--",100,0,2000);
  hHPPmass_trmcmumatchisol24mu = fs->make<TH1F>("hHPPmass_trmcmumatchisol24mu","Mass of H++",100,min,max);
  hHMMmass_trmcmumatchisol24mu = fs->make<TH1F>("hHMMmass_trmcmumatchisol24mu","Mass of H--",100,min,max);
  hHPPmass_trmcmumatchisol2dphi = fs->make<TH1F>("hHPPmass_trmcmumatchisol2dphi","Mass of H++",100,min,max);
  hHMMmass_trmcmumatchisol2dphi = fs->make<TH1F>("hHMMmass_trmcmumatchisol2dphi","Mass of H--",100,min,max);
  hHPPHMMmass_trmcmumatchisol2dphi = fs->make<TH1F>("hHPPHMMmass_trmcmumatchisol2dphi","Mass of H++&H--",100,0,2000);
  hHPPmass_trmcmumatchisol2jdz = fs->make<TH1F>("hHPPmass_trmcmumatchisol2jdz","Mass of H++",100,min,max);
  hHMMmass_trmcmumatchisol2jdz = fs->make<TH1F>("hHMMmass_trmcmumatchisol2jdz","Mass of H--",100,min,max);
  hHPPmass_trmcmumatchisol2jdz_ex = fs->make<TH1F>("hHPPmass_trmcmumatchisol2jdz_ex","Mass of H++",200,0,700);
  hHMMmass_trmcmumatchisol2jdz_ex = fs->make<TH1F>("hHMMmass_trmcmumatchisol2jdz_ex","Mass of H--",200,0,700);

  htrgmupt_ptetacut = fs->make<TH2F>("htrgmupt_ptetacut","track pt VS global muon pt",200,0,1000,200,0,1000);
  htrgmupt_ptetaisol2cut = fs->make<TH2F>("htrgmupt_ptetaisol2cut","track pt VS global muon pt",200,0,1000,200,0,1000);
  htrgmupt_ptetaisol2jcut = fs->make<TH2F>("htrgmupt_ptetaisol2jcut","track pt VS global muon pt",200,0,1000,200,0,1000);

  htrstapt_ptetacut = fs->make<TH2F>("htrstapt_ptetacut","tr P_{T} VS sta P_{T}",200,0,1000,200,0,1000);
  htrstapt_ptetaisol2notjcut = fs->make<TH2F>("htrstapt_ptetaisol2notjcut","tr P_{T} VS sta P_{T}",200,0,1000,200,0,1000);
//---------------global muon---------------//
  hngmu = fs->make<TH1F>("hngmu","number of global muon per event",10,0,10);
  hngmu_ptetacut = fs->make<TH1F>("hngmu_ptetacut","number of global muon per event (P_{T}>10GeV, |#eta|<2.4)",10,0,10);
  hngmu_ptetaisolcut = fs->make<TH1F>("hngmu_ptetaisolcut",Form("# of global muon per event (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),10,0,10);
  hngmu_ptetaisol2cut = fs->make<TH1F>("hngmu_ptetaisol2cut",Form("# of global muon per event (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),10,0,10);
  hngmu_ptetaisol2jcut = fs->make<TH1F>("hngmu_ptetaisol2jcut",Form("# of global muon per event (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),10,0,10);
  hgmupt = fs->make<TH1F>("hgmupt","P_{T} of global muon",50,0,1000);
  hgmupt_ptetacut = fs->make<TH1F>("hgmupt_ptetacut","P_{T} of global muon (P_{T}>10GeV, |#eta|<2.4)",50,0,1000);
  hgmupt_ptetaisolcut = fs->make<TH1F>("hgmupt_ptetaisolcut",Form("P_{T} of global muon (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,0,1000);
  hgmupt_ptetaisol2cut = fs->make<TH1F>("hgmupt_ptetaisol2cut",Form("P_{T} of global muon (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,0,1000);
  hgmupt_ptetaisol2jcut = fs->make<TH1F>("hgmupt_ptetaisol2jcut",Form("P_{T} of global muon (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,0,1000);
  hgmupt_mcmumatch = fs->make<TH1F>("hgmupt_mcmumatch","P_{T} of global muon (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01)",50,0,1000);
  hgmupt_mcmumatchisol = fs->make<TH1F>("hgmupt_mcmumatchisol",Form("P_{T} of global muon (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,0,1000);
  hgmupt_mcmumatchisol2 = fs->make<TH1F>("hgmupt_mcmumatchisol2",Form("P_{T} of global muon (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,0,1000);
  hgmupt_mcmumatchisol2j = fs->make<TH1F>("hgmupt_mcmumatchisol2j",Form("P_{T} of global muon (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,0,1000);
  hgmueta = fs->make<TH1F>("hgmueta","#eta of global muon",50,-2.4,2.4);
  hgmueta_ptetacut = fs->make<TH1F>("hgmueta_ptetacut","#eta of global muon (P_{T}>10GeV, |#eta|<2.4)",50,-2.4,2.4);
  hgmueta_ptetaisolcut=fs->make<TH1F>("hgmueta_ptetaisolcut",Form("#eta of global muon (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  hgmueta_ptetaisol2cut=fs->make<TH1F>("hgmueta_ptetaisol2cut",Form("#eta of global muon (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  hgmueta_ptetaisol2jcut=fs->make<TH1F>("hgmueta_ptetaisol2jcut",Form("#eta of global muon (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  hgmueta_mcmumatch = fs->make<TH1F>("hgmueta_mcmumatch","#eta of global muon (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01)",50,-2.4,2.4);
  hgmueta_mcmumatchisol = fs->make<TH1F>("hgmueta_mcmumatchisol",Form("#eta of global muon (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  hgmueta_mcmumatchisol2 = fs->make<TH1F>("hgmueta_mcmumatchisol2",Form("#eta of global muon (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  hgmueta_mcmumatchisol2j = fs->make<TH1F>("hgmueta_mcmumatchisol2j",Form("#eta of global muon (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  hgmuisolPt_ptetacut01 = fs->make<TH1F>("hgmuisolPt_ptetacut01","isolP_{T} of global muon (P_{T}>10GeV, |#eta|<2.4)",50,0,0.1);
  hgmuisolPt_ptetacut1 = fs->make<TH1F>("hgmuisolPt_ptetacut1","isolP_{T} of global muon (P_{T}>10GeV, |#eta|<2.4)",50,0,1);
  hgmuisolPt_ptetacut10 = fs->make<TH1F>("hgmuisolPt_ptetacut10","isolP_{T} of global muon (P_{T}>10GeV, |#eta|<2.4)",50,0,10);
  hgmuisolPt_ptetacut100 = fs->make<TH1F>("hgmuisolPt_ptetacut100","isolP_{T} of global muon (P_{T}>10GeV, |#eta|<2.4)",50,0,100);
  hgmuisolPt_ptetacut1000 = fs->make<TH1F>("hgmuisolPt_ptetacut1000","isolP_{T} of global muon (P_{T}>10GeV, |#eta|<2.4)",50,0,1000);
  hgmuisolPt2_ptetacut10 = fs->make<TH1F>("hgmuisolPt2_ptetacut10","isolP_{T} of global muon (P_{T}>10GeV, |#eta|<2.4)",50,0,10);
  hgmuisolPt2_ptetacut100 = fs->make<TH1F>("hgmuisolPt2_ptetacut100","isolP_{T} of global muon (P_{T}>10GeV, |#eta|<2.4)",50,0,100);
  h2gmudeltaphi = fs->make<TH1F>("h2gmudeltaphi","#Delta#phi between gmu+gmu+ & gmu-gmu-",50,0,TMath::Pi());
  h2gmudeltaphi_isol = fs->make<TH1F>("h2gmudeltaphi_isol","#Delta#phi between gmu+gmu+ & gmu-gmu-",50,0,TMath::Pi());
  h2gmudeltaphi_isol2 = fs->make<TH1F>("h2gmudeltaphi_isol2","#Delta#phi between gmu+gmu+ & gmu-gmu-",50,0,TMath::Pi());
  h2gmu2mcmudphi_isol2 = fs->make<TH2F>("h2gmu2mcmudphi_isol2","#Delta#phi 2 gmu VS 2 MC muon",200,0,TMath::Pi(),200,0,TMath::Pi());
  hgmuz = fs->make<TH1F>("hgmuz","z position of global muons",200,-15,15);
  hgmuz_ptetacut = fs->make<TH1F>("hgmuz_ptetacut","z position of global muons",200,-15,15);
  hgmuz_ptetaisol2cut = fs->make<TH1F>("hgmuz_ptetaisol2cut","z position of global muons",200,-15,15);
  hgmuz_ptetaisol2jcut = fs->make<TH1F>("hgmuz_ptetaisol2jcut","z position of global muons",200,-15,15);
  hgmudz = fs->make<TH1F>("hgmudz","#Deltaz between global muons",200,-0.1,0.1);
  hgmudz_ptetaisol2jcut = fs->make<TH1F>("hgmudz_ptetaisol2jcut","#Deltaz between global muons",200,-0.1,0.1);
  h2gmudz_isol2j = fs->make<TH1F>("h2gmudz_isol2j","#Deltaz between two highest pt global muons",200,-0.05,0.05);
  h2gmudz_isol2j_ex = fs->make<TH1F>("h2gmudz_isol2j_ex","#Deltaz between two highest pt global muons",200,-0.25,0.25);

  hHPPmass = fs->make<TH1F>("hHPPmass","Mass of H++",100,min,max);
  hHMMmass = fs->make<TH1F>("hHMMmass","Mass of H--",100,min,max);
  hHPPHMMmass = fs->make<TH1F>("hHPPHMMmass","Mass of H++&H--",100,0,2000);
  hHPPmass_isol = fs->make<TH1F>("hHPPmass_isol","Mass of H++",100,min,max);
  hHMMmass_isol = fs->make<TH1F>("hHMMmass_isol","Mass of H--",100,min,max);
  hHPPHMMmass_isol = fs->make<TH1F>("hHPPHMMmass_isol","Mass of H++&H--",100,0,2000);
  hHPPmass_isol2 = fs->make<TH1F>("hHPPmass_isol2","Mass of H++",100,min,max);
  hHMMmass_isol2 = fs->make<TH1F>("hHMMmass_isol2","Mass of H--",100,min,max);
  hHPPHMMmass_isol2 = fs->make<TH1F>("hHPPHMMmass_isol2","Mass of H++&H--",100,0,2000);
  hHPPmass_isol2j = fs->make<TH1F>("hHPPmass_isol2j","Mass of H++",100,min,max);
  hHMMmass_isol2j = fs->make<TH1F>("hHMMmass_isol2j","Mass of H--",100,min,max);
  hHPPHMMmass_isol2j = fs->make<TH1F>("hHPPHMMmass_isol2j","Mass of H++&H--",100,0,2000);
  hHPPmass_isol2dphi = fs->make<TH1F>("hHPPmass_isol2dphi","Mass of H++",100,min,max);
  hHMMmass_isol2dphi = fs->make<TH1F>("hHMMmass_isol2dphi","Mass of H--",100,min,max);
  hHPPHMMmass_isol2dphi = fs->make<TH1F>("hHPPHMMmass_isol2dphi","Mass of H++&H--",100,0,2000);
  hHPPmass_isol2jdz = fs->make<TH1F>("hHPPmass_isol2jdz","Mass of H++",100,min,max);
  hHMMmass_isol2jdz = fs->make<TH1F>("hHMMmass_isol2jdz","Mass of H--",100,min,max);
  hHPPmass_isol2jdz_ex = fs->make<TH1F>("hHPPmass_isol2jdz_ex","Mass of H++",200,0,700);
  hHMMmass_isol2jdz_ex = fs->make<TH1F>("hHMMmass_isol2jdz_ex","Mass of H--",200,0,700);

  hHPPmass_mcmumatch = fs->make<TH1F>("hHPPmass_mcmumatch","Mass of H++",100,min,max);
  hHMMmass_mcmumatch = fs->make<TH1F>("hHMMmass_mcmumatch","Mass of H--",100,min,max);
  hHPPHMMmass_mcmumatch = fs->make<TH1F>("hHPPHMMmass_mcmumatch","Mass of H++&H--",100,0,2000);
  hHPPmass_mcmumatchisol = fs->make<TH1F>("hHPPmass_mcmumatchisol","Mass of H++",100,min,max);
  hHMMmass_mcmumatchisol = fs->make<TH1F>("hHMMmass_mcmumatchisol","Mass of H--",100,min,max);
  hHPPHMMmass_mcmumatchisol = fs->make<TH1F>("hHPPHMMmass_mcmumatchisol","Mass of H++&H--",100,0,2000);
  hHPPmass_mcmumatchisol2 = fs->make<TH1F>("hHPPmass_mcmumatchisol2","Mass of H++",100,min,max);
  hHMMmass_mcmumatchisol2 = fs->make<TH1F>("hHMMmass_mcmumatchisol2","Mass of H--",100,min,max);
  hHPPHMMmass_mcmumatchisol2 = fs->make<TH1F>("hHPPHMMmass_mcmumatchisol2","Mass of H++&H--",100,0,2000);
  hHPPmass_mcmumatchisol2j = fs->make<TH1F>("hHPPmass_mcmumatchisol2j","Mass of H++",100,min,max);
  hHMMmass_mcmumatchisol2j = fs->make<TH1F>("hHMMmass_mcmumatchisol2j","Mass of H--",100,min,max);
  hHPPHMMmass_mcmumatchisol2j = fs->make<TH1F>("hHPPHMMmass_mcmumatchisol2j","Mass of H++&H--",100,0,2000);
  hHPPmass_mcmumatchisol24mu = fs->make<TH1F>("hHPPmass_mcmumatchisol24mu","Mass of H++",100,min,max);
  hHMMmass_mcmumatchisol24mu = fs->make<TH1F>("hHMMmass_mcmumatchisol24mu","Mass of H--",100,min,max);
  hHPPmass_mcmumatchisol2dphi = fs->make<TH1F>("hHPPmass_mcmumatchisol2dphi","Mass of H++",100,min,max);
  hHMMmass_mcmumatchisol2dphi = fs->make<TH1F>("hHMMmass_mcmumatchisol2dphi","Mass of H--",100,min,max);
  hHPPHMMmass_mcmumatchisol2dphi = fs->make<TH1F>("hHPPHMMmass_mcmumatchisol2dphi","Mass of H++&H--",100,0,2000);
  hHPPmass_mcmumatchisol2jdz = fs->make<TH1F>("hHPPmass_mcmumatchisol2jdz","Mass of H++",100,min,max);
  hHMMmass_mcmumatchisol2jdz = fs->make<TH1F>("hHMMmass_mcmumatchisol2jdz","Mass of H--",100,min,max);
  hHPPmass_mcmumatchisol2jdz_ex = fs->make<TH1F>("hHPPmass_mcmumatchisol2jdz_ex","Mass of H++",200,0,700);
  hHMMmass_mcmumatchisol2jdz_ex = fs->make<TH1F>("hHMMmass_mcmumatchisol2jdz_ex","Mass of H--",200,0,700);

  hgmustapt_ptetacut = fs->make<TH2F>("hgmustapt_ptetacut","gmu P_{T} VS sta P_{T}",200,0,1000,200,0,1000);
  hgmustapt_ptetaisol2notjcut = fs->make<TH2F>("hgmustapt_ptetaisol2notjcut","gmu P_{T} VS sta P_{T}",200,0,1000,200,0,1000);
//---------------MC particle---------------//
  hnmcmu_ptetacut = fs->make<TH1F>("hnmcmu_ptetacut","number of MC muon from H++ per event (P_{T}>10GeV, |#eta|<2.4)",10,0,10);
  hnmcmu_ptetaisolcut = fs->make<TH1F>("hnmcmu_ptetaisolcut",Form("number of MC muon from H++ per event (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),10,0,10);
  hnmcmu_ptetaisol2cut = fs->make<TH1F>("hnmcmu_ptetaisol2cut",Form("number of MC muon from H++ per event (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),10,0,10);
  hnmcmu_ptetaisol2jcut = fs->make<TH1F>("hnmcmu_ptetaisol2jcut",Form("number of MC muon from H++ per event (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d, P_{T,jet}>15, #DeltaR(mu,jet)>0.5)",isolPtcut),10,0,10);
  hmcmupt = fs->make<TH1F>("hmcmupt","P_{T} of MC muon from Higgs",50,0,1000);
  hmcmupt_ptetacut = fs->make<TH1F>("hmcmupt_ptetacut","P_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,1000);
  hmcmupt_ptetaisolcut = fs->make<TH1F>("hmcmupt_ptetaisolcut",Form("P_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,0,1000);
  hmcmupt_ptetaisol2cut = fs->make<TH1F>("hmcmupt_ptetaisol2cut",Form("P_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,0,1000);
  hmcmupt_ptetaisol2jcut = fs->make<TH1F>("hmcmupt_ptetaisol2jcut",Form("P_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,0,1000);
  hmcmupt_gmumatch = fs->make<TH1F>("hmcmupt_gmumatch","P_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01)",50,0,1000);
  hmcmupt_gmumatchisol = fs->make<TH1F>("hmcmupt_gmumatchisol",Form("P_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,0,1000);
  hmcmupt_gmumatchisol2 = fs->make<TH1F>("hmcmupt_gmumatchisol2",Form("P_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,0,1000);
  hmcmupt_trmatch_gmuunmatch = fs->make<TH1F>("hmcmupt_trmatch_gmuunmatch","P_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{track,MC#mu}<0.01, #DeltaR_{glob,MC#mu}>0.01)",50,0,1000);
  hmcmueta = fs->make<TH1F>("hmcmueta","#eta of MC muon from Higgs",50,-5,5);
  hmcmueta_ptetacut = fs->make<TH1F>("hmcmueta_ptetacut","#eta of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,-2.4,2.4);
  hmcmueta_ptetaisolcut = fs->make<TH1F>("hmcmueta_ptetaisolcut",Form("#eta of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  hmcmueta_ptetaisol2cut = fs->make<TH1F>("hmcmueta_ptetaisol2cut",Form("#eta of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  hmcmueta_ptetaisol2jcut = fs->make<TH1F>("hmcmueta_ptetaisol2jcut",Form("#eta of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  hmcmueta_gmumatch = fs->make<TH1F>("hmcmueta_gmumatch","#eta of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01)",50,-2.4,2.4);
  hmcmueta_gmumatchisol = fs->make<TH1F>("hmcmueta_gmumatchisol",Form("#eta of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  hmcmueta_gmumatchisol2 = fs->make<TH1F>("hmcmueta_gmumatchisol2",Form("#eta of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{glob,MC#mu}<0.01, isolP_{T}<%d)",isolPtcut),50,-2.4,2.4);
  hmcmueta_trmatch_gmuunmatch = fs->make<TH1F>("hmcmueta_trmatch_gmuunmatch","#eta of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, #DeltaR_{track,MC#mu}<0.01, #DeltaR_{glob,MC#mu}>0.01)",50,-2.4,2.4);
  hmcisolPt_ptetacut01 = fs->make<TH1F>("hmcisolPt_ptetacut01","isolP_{T} of MC particles (P_{T}>10GeV, |#eta|<2.4)",50,0,0.1);
  hmcisolPt_ptetacut1 = fs->make<TH1F>("hmcisolPt_ptetacut1","isolP_{T} of MC particles (P_{T}>10GeV, |#eta|<2.4)",50,0,1);
  hmcisolPt_ptetacut10 = fs->make<TH1F>("hmcisolPt_ptetacut10","isolP_{T} of MC particles (P_{T}>10GeV, |#eta|<2.4)",50,0,10);
  hmcisolPt_ptetacut100 = fs->make<TH1F>("hmcisolPt_ptetacut100","isolP_{T} of MC particles (P_{T}>10GeV, |#eta|<2.4)",50,0,100);
  hmcisolPt_ptetacut1000 = fs->make<TH1F>("hmcisolPt_ptetacut1000","isolP_{T} of MC particles (P_{T}>10GeV, |#eta|<2.4)",50,0,1000);
  hmcisolPt2_ptetacut10 = fs->make<TH1F>("hmcisolPt2_ptetacut10","isolP_{T} of MC particles (P_{T}>10GeV, |#eta|<2.4)",50,0,10);
  hmcisolPt2_ptetacut100 = fs->make<TH1F>("hmcisolPt2_ptetacut100","isolP_{T} of MC particles (P_{T}>10GeV, |#eta|<2.4)",50,0,100);
  hmcmuisolPt_ptetacut01 = fs->make<TH1F>("hmcmuisolPt_ptetacut01","isolP_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,0.1);
  hmcmuisolPt_ptetacut1 = fs->make<TH1F>("hmcmuisolPt_ptetacut1","isolP_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,1);
  hmcmuisolPt_ptetacut10 = fs->make<TH1F>("hmcmuisolPt_ptetacut10","isolP_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,10);
  hmcmuisolPt_ptetacut100 = fs->make<TH1F>("hmcmuisolPt_ptetacut100","isolP_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,100);
  hmcmuisolPt_ptetacut1000 = fs->make<TH1F>("hmcmuisolPt_ptetacut1000","isolP_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,1000);
  hmcmuisolPt2_ptetacut10 = fs->make<TH1F>("hmcmuisolPt2_ptetacut10","isolP_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,10);
  hmcmuisolPt2_ptetacut100 = fs->make<TH1F>("hmcmuisolPt2_ptetacut100","isolP_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,100);
  hmcmuisolPt_mcmudR03 = fs->make<TH1F>("hmcmuisolPt_mcmudR03","isolP_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, #DeltaR(#mu,#mu)<0.3)",50,0,500);
  hmcmuisolPt2_mcmudR03 = fs->make<TH1F>("hmcmuisolPt2_mcmudR03","isolP_{T} of MC muon from Higgs (P_{T}>10GeV, |#eta|<2.4, #DeltaR(#mu,#mu)<0.3)",50,0,500);

  hmcmudeltaphi = fs->make<TH1F>("hmcmudeltaphi","#Delta#phi between MC muons from Higgs",50,0,TMath::Pi());
  hmcmudeltaphi2 = fs->make<TH1F>("hmcmudeltaphi2","#Delta#phi between MC muons from Z",50,0,TMath::Pi());
  hmcHdeltaphi = fs->make<TH1F>("hmcHdeltaphi","#Delta#phi between MC H++ & H--",50,0,TMath::Pi());
  h2mcmudeltaphi = fs->make<TH1F>("h2mcmudeltaphi","#Delta#phi between MC mu+mu+ & mu-mu-",50,0,TMath::Pi());
  h2mcmudeltaphi_isol = fs->make<TH1F>("h2mcmudeltaphi_isol","#Delta#phi between MC mu+mu+ & mu-mu-",50,0,TMath::Pi());
  h2mcmudeltaphi_isol2 = fs->make<TH1F>("h2mcmudeltaphi_isol2","#Delta#phi between MC mu+mu+ & mu-mu-",50,0,TMath::Pi());
  hmcmudR = fs->make<TH1F>("hmcmudR","#DeltaR between MC muons from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,sqrt(2.4*2.4+TMath::Pi()*TMath::Pi()));
  hmcmudR_ex = fs->make<TH1F>("hmcmudR_ex","#DeltaR between MC muons from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,10);
  hmcmudR_allmu = fs->make<TH1F>("hmcmudR_allmu","#DeltaR between MC muons from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,sqrt(2.4*2.4+TMath::Pi()*TMath::Pi()));
  hmcmudR_upperisolPt_10 = fs->make<TH1F>("hmcmudR_upperisolPt_10","#DeltaR between MC muons from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,sqrt(2.4*2.4+TMath::Pi()*TMath::Pi()));
  hmcmudR_upperisolPt_20 = fs->make<TH1F>("hmcmudR_upperisolPt_20","#DeltaR between MC muons from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,sqrt(2.4*2.4+TMath::Pi()*TMath::Pi()));
  hmcmudR_upperisolPt2_10 = fs->make<TH1F>("hmcmudR_upperisolPt2_10","#DeltaR between MC muons from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,sqrt(2.4*2.4+TMath::Pi()*TMath::Pi()));
  hmcmudR_upperisolPt2_20 = fs->make<TH1F>("hmcmudR_upperisolPt2_20","#DeltaR between MC muons from Higgs (P_{T}>10GeV, |#eta|<2.4)",50,0,sqrt(2.4*2.4+TMath::Pi()*TMath::Pi()));
  hmcmuxy = fs->make<TH2F>("hmcmuxy","vertex of MC muons in xy plane",300,minmcvtx,maxmcvtx,300,minmcvty,maxmcvty);
  hmcmuxy_upper1k = fs->make<TH2F>("hmcmuxy_upper1k","vertex of MC muons in xy plane",300,minmcvtx,maxmcvtx,300,minmcvty,maxmcvty);
  hmcmuxy_10times1k = fs->make<TH2F>("hmcmuxy_10times1k","vertex of MC muons in xy plane",300,minmcvtx,maxmcvtx,300,minmcvty,maxmcvty);

  hHPPmass_mc = fs->make<TH1F>("hHPPmass_mc","Mass of MC H++",100,min,max);
  hHMMmass_mc = fs->make<TH1F>("hHMMmass_mc","Mass of MC H--",100,min,max);
  hHPPHMMmass_mc = fs->make<TH1F>("hHPPHMMmass_mc","Mass of MC H++&H--",100,0,2000);
  hHPPmass_mcisol = fs->make<TH1F>("hHPPmass_mcisol","Mass of MC H++",100,min,max);
  hHMMmass_mcisol = fs->make<TH1F>("hHMMmass_mcisol","Mass of MC H--",100,min,max);
  hHPPHMMmass_mcisol = fs->make<TH1F>("hHPPHMMmass_mcisol","Mass of MC H++&H--",100,0,2000);
  hHPPmass_mcisol2 = fs->make<TH1F>("hHPPmass_mcisol2","Mass of MC H++",100,min,max);
  hHMMmass_mcisol2 = fs->make<TH1F>("hHMMmass_mcisol2","Mass of MC H--",100,min,max);
  hHPPHMMmass_mcisol2 = fs->make<TH1F>("hHPPHMMmass_mcisol2","Mass of MC H++&H--",100,0,2000);
  hHPPmass_mcisol2j = fs->make<TH1F>("hHPPmass_mcisol2j","Mass of MC H++",100,min,max);
  hHMMmass_mcisol2j = fs->make<TH1F>("hHMMmass_mcisol2j","Mass of MC H--",100,min,max);
  hHPPHMMmass_mcisol2j = fs->make<TH1F>("hHPPHMMmass_mcisol2j","Mass of MC H++&H--",100,0,2000);
  hHPPmass_mcisol2dphi = fs->make<TH1F>("hHPPmass_mcisol2dphi","Mass of MC H++",100,min,max);
  hHMMmass_mcisol2dphi = fs->make<TH1F>("hHMMmass_mcisol2dphi","Mass of MC H--",100,min,max);
  hHPPHMMmass_mcisol2dphi = fs->make<TH1F>("hHPPHMMmass_mcisol2dphi","Mass of MC H++&H--",100,0,2000);

  hgmumcmu_ptratio = fs->make<TH2F>("hgmumcmu_ptratio","|(P_{T(glob)}-P_{T(MC#mu)}|/P_{T(MC#mu)}",200,0,1000,200,0,0.2);
  hgmumcmu_ptratio_ex = fs->make<TH2F>("hgmumcmu_ptratio_ex","|P_{T(glob)}-P_{T(MC#mu)}|/P_{T(MC#mu)}",200,0,1000,200,0,1);
  hgmumcmu_ptratio_2nd = fs->make<TH2F>("hgmumcmu_ptratio_2nd","|P_{T(glob)}-P_{T(MC#mu)}|/P_{T(MC#mu)}",200,0,1000,200,0,0.2);
  hgmumcmu_ptratio_2nd_ex = fs->make<TH2F>("hgmumcmu_ptratio_2nd_ex","|P_{T(glob)}-P_{T(MC#mu)}|/P_{T(MC#mu)}",200,0,1000,200,0,1);
  hgmumcmu_ptratio_2ndtrue = fs->make<TH2F>("hgmumcmu_ptratio_2ndtrue","|P_{T(glob)}-P_{T(MC#mu)}|/P_{T(MC#mu)}",200,0,1000,200,0,0.2);
  hgmumcmu_ptratio_2ndtrue_ex = fs->make<TH2F>("hgmumcmu_ptratio_2ndtrue_ex","|P_{T(glob)}-P_{T(MC#mu)}|/P_{T(MC#mu)}",200,0,1000,200,0,1);
  hgmumcmu_ptratioisol = fs->make<TH2F>("hgmumcmu_ptratioisol","|P_{T(glob)}-P_{T(MC#mu)}|/P_{T(MC#mu)}",200,0,1000,200,0,0.2);
  hgmumcmu_ptratioisol_ex = fs->make<TH2F>("hgmumcmu_ptratioisol_ex","|P_{T(glob)}-P_{T(MC#mu)}|/P_{T(MC#mu)}",200,0,1000,200,0,1);
  hgmumcmu_ptratioisol2 = fs->make<TH2F>("hgmumcmu_ptratioisol2","|P_{T(glob)}-P_{T(MC#mu)}|/P_{T(MC#mu)}",200,0,1000,200,0,0.2);
  hgmumcmu_ptratioisol2_ex = fs->make<TH2F>("hgmumcmu_ptratioisol2_ex","|P_{T(glob)}-P_{T(MC#mu)}|/P_{T(MC#mu)}",200,0,1000,200,0,1);
  hgmumcmu_ptratioisol2j = fs->make<TH2F>("hgmumcmu_ptratioisol2j","|P_{T(glob)}-P_{T(MC#mu)}|/P_{T(MC#mu)}",200,0,1000,200,0,0.2);
  hgmumcmu_ptratioisol2j_ex = fs->make<TH2F>("hgmumcmu_ptratioisol2j_ex","|P_{T(glob)}-P_{T(MC#mu)}|/P_{T(MC#mu)}",200,0,1000,200,0,1);

//---------------MET---------------//
  hcalmet = fs->make<TH1F>("hcalmet","CaloMET",200,0,300);
  hcalmet01mu = fs->make<TH1F>("hcalmet01mu","CaloMET",200,0,300);
  hgenmet = fs->make<TH1F>("hgenmet","GenMET",200,0,300);
  hgenmetnonu = fs->make<TH1F>("hgenmetnonu","GenMET",200,0,300);
  hgenmet01mu = fs->make<TH1F>("hgenmet01mu","GenMET",200,0,300);
  hcalgenmet = fs->make<TH2F>("hcalgenmet","CaloMET",300,0,300,300,0,300);
  hcalgenmetnonu = fs->make<TH2F>("hcalgenmetnonu","CaloMET",300,0,300,300,0,300);
  hcalgenmet01mu = fs->make<TH2F>("hcalgenmet01mu","CaloMET",300,0,300,300,0,300);
  hcal01mugenmet = fs->make<TH2F>("hcal01mugenmet","CaloMET",300,0,300,300,0,300);
  hcalmet0gmu = fs->make<TH1F>("hcalmet0gmu","CaloMET",200,0,300);
  hgenmet0gmu = fs->make<TH1F>("hgenmet0gmu","GenMET",200,0,300);
  hcalgenmet0gmu = fs->make<TH2F>("hcalgenmet0gmu","CaloMET",300,0,300,300,0,300);
  hcalmet1gmu = fs->make<TH1F>("hcalmet1gmu","CaloMET",200,0,300);
  hgenmet1gmu = fs->make<TH1F>("hgenmet1gmu","GenMET",200,0,300);
  hcalgenmet1gmu = fs->make<TH2F>("hcalgenmet1gmu","CaloMET",300,0,300,300,0,300);
  hcalmet2gmu = fs->make<TH1F>("hcalmet2gmu","CaloMET",200,0,300);
  hgenmet2gmu = fs->make<TH1F>("hgenmet2gmu","GenMET",200,0,300);
  hcalgenmet2gmu = fs->make<TH2F>("hcalgenmet2gmu","CaloMET",300,0,300,300,0,300);
  hcalmet3gmu = fs->make<TH1F>("hcalmet3gmu","CaloMET",200,0,300);
  hgenmet3gmu = fs->make<TH1F>("hgenmet3gmu","GenMET",200,0,300);
  hcalgenmet3gmu = fs->make<TH2F>("hcalgenmet3gmu","CaloMET",300,0,300,300,0,300);
  hcalmet4gmu = fs->make<TH1F>("hcalmet4gmu","CaloMET",200,0,300);
  hgenmet4gmu = fs->make<TH1F>("hgenmet4gmu","GenMET",200,0,300);
  hcalgenmet4gmu = fs->make<TH2F>("hcalgenmet4gmu","CaloMET",300,0,300,300,0,300);

//---------------ZZ->4mu---------------//
  hPMdphi_trisol2j= fs->make<TH1F>("hPMdphi_trisol2j","#Delta#phi between #mu+ and #mu-",100,0,TMath::Pi());

  hPMmassall_trisol2j = fs->make<TH1F>("hPMmassall_trisol2j","Mass of #mu+#mu-",100,0,1000);
  hPMmassbig_trisol2j = fs->make<TH1F>("hPMmassbig_trisol2j","Mass of #mu+#mu-",100,0,1000);
  hPMmassZ_trisol2j = fs->make<TH1F>("hPMmassZ_trisol2j","Mass of #mu+#mu-",100,0,200);
  hPMmassdphi_trisol2j = fs->make<TH1F>("hPMmassdphi_trisol2j","Mass of #mu+#mu-",100,0,1000);

//---------------Vertex---------------//
  hnvt = fs->make<TH1F>("hnvt","# of vertex per event",6,0,6);
  hvtx = fs->make<TH1F>("hvtx","x position of vertices",100,minvtx,maxvtx);
  hvty = fs->make<TH1F>("hvty","y position of vertices",100,minvty,maxvty);
  hvtz = fs->make<TH1F>("hvtz","z position of vertices",100,-15,15);
  hvtxy = fs->make<TH2F>("hvtxy","x&y position of vertices",300,minvtx,maxvtx,300,minvty,maxvty);
  hvtz_BS = fs->make<TH1F>("hvtz_BS","z position of vertices",100,-15,15);
  hnvt2 = fs->make<TH1F>("hnvt2","# of vertices per event",6,0,6);
  hvt2x = fs->make<TH1F>("hvt2x","x position of vertices",100,minvt2x,maxvt2x);
  hvt2y = fs->make<TH1F>("hvt2y","y position of vertices",100,minvt2y,maxvt2y);
  hvt2z = fs->make<TH1F>("hvt2z","z position of vertices",100,-15,15);
  hvt2xy = fs->make<TH2F>("hvt2xy","x&y position of vertices",300,minvt2x,maxvt2x,300,minvt2y,maxvt2y);
  hvt2z_BS = fs->make<TH1F>("hvt2z_BS","z position of vertices",100,-15,15);
  hnvt3 = fs->make<TH1F>("hnvt3","# of vertices per event",6,0,6);
  hvt3x = fs->make<TH1F>("hvt3x","x position of vertices",100,minvt3x,maxvt3x);
  hvt3y = fs->make<TH1F>("hvt3y","y position of vertices",100,minvt3y,maxvt3y);
  hvt3z = fs->make<TH1F>("hvt3z","z position of vertices",100,-15,15);
  hvt3xy = fs->make<TH2F>("hvt3xy","x&y position of vertices",300,minvt3x,maxvt3x,300,minvt3y,maxvt3y);
  hvt3z_BS = fs->make<TH1F>("hvt3z_BS","z position of vertices",100,-15,15);
  hvtntr = fs->make<TH1F>("hvtntr","number of tracks",120,0,120);
  hvt2ntr = fs->make<TH1F>("hvt2ntr","number of tracks",120,0,120);
  hvt3ntr = fs->make<TH1F>("hvt3ntr","number of tracks",120,0,120);

  hmcmuz = fs->make<TH1F>("hmcmuz","vertex of MC muons",100,-15,15);
  hmcmux_BS = fs->make<TH1F>("hmcmux_BS","Vertex of MC muons (x position)",50,minvtx,maxvtx);
  hmcmuy_BS = fs->make<TH1F>("hmcmuy_BS","Vertex of MC muons (y position)",50,minvty,maxvty);
  hmcmuz_BS = fs->make<TH1F>("hmcmuz_BS","Vertex of MC muons (z position)",50,-20,20);
  hmcmuxy_BS = fs->make<TH2F>("hmcmuxy_BS","Vertex of MC muons (x&y position)",300,minvtx,maxvtx,300,minvty,maxvty);
  hmcmu2x_BS = fs->make<TH1F>("hmcmu2x_BS","Vertex of MC muons (x position)",50,minvtx,maxvtx);
  hmcmu2y_BS = fs->make<TH1F>("hmcmu2y_BS","Vertex of MC muons (y position)",50,minvty,maxvty);
  hmcmu2z_BS = fs->make<TH1F>("hmcmu2z_BS","Vertex of MC muons (z position)",50,-20,20);
  hmcmu2xy_BS = fs->make<TH2F>("hmcmu2xy_BS","Vertex of MC muons (x&y position)",300,minvtx,maxvtx,300,minvty,maxvty);
  hmcmu3x_BS = fs->make<TH1F>("hmcmu3x_BS","Vertex of MC muons (x position)",50,minvtx,maxvtx);
  hmcmu3y_BS = fs->make<TH1F>("hmcmu3y_BS","Vertex of MC muons (y position)",50,minvty,maxvty);
  hmcmu3z_BS = fs->make<TH1F>("hmcmu3z_BS","Vertex of MC muons (z position)",50,-20,20);
  hmcmu3xy_BS = fs->make<TH2F>("hmcmu3xy_BS","Vertex of MC muons (x&y position)",300,minvtx,maxvtx,300,minvty,maxvty);

  htrvtdxy = fs->make<TH1F>("htrvtdxy","distance between track line & vertex in xy plane",100,0,0.15);
  htrvtdxy_isol2j = fs->make<TH1F>("htrvtdxy_isol2j","distance between track line & vertex in xy plane",100,0,0.15);
  htrvtdxy_isol2j_in = fs->make<TH1F>("htrvtdxy_isol2j_in","distance between track line & vertex in xy plane",100,0,0.015);
  htrvtdxy_isol2jmc = fs->make<TH1F>("htrvtdxy_isol2jmc","distance between track line & vertex in xy plane",100,0,0.15);
  htrvtdxy_isol2jmc_in = fs->make<TH1F>("htrvtdxy_isol2jmc_in","distance between track line & vertex in xy plane",100,0,0.015);
  htrvtdxy_isol2j4mu = fs->make<TH1F>("htrvtdxy_isol2j4mu","distance between track line & vertex in xy plane",100,0,0.15);
  htrvtdxy_isol2j4mu_in = fs->make<TH1F>("htrvtdxy_isol2j4mu_in","distance between track line & vertex in xy plane",100,0,0.015);
  htrvtdxy_isol2j4mumc = fs->make<TH1F>("htrvtdxy_isol2j4mumc","distance between track line & vertex in xy plane",100,0,0.15);
  htrvtdxy_isol2j4mumc_in = fs->make<TH1F>("htrvtdxy_isol2j4mumc_in","distance between track line & vertex in xy plane",100,0,0.015);
  htrvtdxy_cor = fs->make<TH1F>("htrvtdxy_cor","distance between track line & vertex in xy plane",100,0,0.15);
  htrvtdxy_isol2j_cor = fs->make<TH1F>("htrvtdxy_isol2j_cor","distance between track line & vertex in xy plane",100,0,0.15);
  htrvtdz = fs->make<TH1F>("htrvtdz","#Deltaz between track & vertex",100,0,0.3);
  htrvtdz2 = fs->make<TH1F>("htrvtdz2","#Deltaz between track & vertex",100,0,0.3);
  htrvtdz3 = fs->make<TH1F>("htrvtdz3","#Deltaz between track & vertex",100,0,0.3);
  htrvtdz_isol2j = fs->make<TH1F>("htrvtdz_isol2j","#Deltaz between track & vertex",100,0,0.1);
  htrvtdz2_isol2j = fs->make<TH1F>("htrvtdz2_isol2j","#Deltaz between track & vertex",100,0,0.1);
  htrvtdz3_isol2j = fs->make<TH1F>("htrvtdz3_isol2j","#Deltaz between track & vertex",100,0,0.1);
  htrvtdz_isol2jmc = fs->make<TH1F>("htrvtdz_isol2jmc","#Deltaz between track & vertex",100,0,0.1);
  htrvtd = fs->make<TH1F>("htrvtd","distance between track line & vertex",100,0,0.3);

  htrvtd_i2j = fs->make<TH1F>("htrvtd_i2j","distance between reference position & vertex",100,0,0.3);
  htrvtd_i2jmc = fs->make<TH1F>("htrvtd_i2jmc","distance between reference position & vertex",100,0,0.3);
  htrvtdz_i2j = fs->make<TH1F>("htrvtdz_i2j","#Deltaz between reference position & vertex",100,0,0.3);
  htrvtdz_i2jmc = fs->make<TH1F>("htrvtdz_i2jmc","#Deltaz between reference position & vertex",100,0,0.3);
  htrvtdxy_i2j = fs->make<TH1F>("htrvtdxy_i2j","#Deltaxy between reference position & vertex",100,0,0.3);
  htrvtdxy_i2jmc = fs->make<TH1F>("htrvtdxy_i2jmc","#Deltaxy between reference position & vertex",100,0,0.3);
  htripvtd = fs->make<TH1F>("htripvtd","distance between track line & vertex",100,0,0.1);
  htripvtd_i2j = fs->make<TH1F>("htripvtd_i2j","distance between track line & vertex",100,0,0.015);
  htripvtd_i2jmc = fs->make<TH1F>("htripvtd_i2jmc","distance between track line & vertex",100,0,0.015);
  htripvtdz_i2j = fs->make<TH1F>("htripvtdz_i2j","#Deltaz between impact point & vertex",100,0,0.015);
  htripvtdz_i2jmc = fs->make<TH1F>("htripvtdz_i2jmc","#Deltaz between impact point & vertex",100,0,0.015);
  htripvtdxy_i2j = fs->make<TH1F>("htripvtdxy_i2j","#Deltaxy between track line & vertex",100,0,0.015);
  htripvtdxy_i2jmc = fs->make<TH1F>("htripvtdxy_i2jmc","#Deltaxy between track line & vertex",100,0,0.015);
  htrtrd_i2j = fs->make<TH1F>("htrtrd_i2j","distance between reference positions",100,0,0.2);
  htrtrd_i2jmc = fs->make<TH1F>("htrtrd_i2jmc","distance between reference positions",100,0,0.2);
  htrtrdz_i2j = fs->make<TH1F>("htrtrdz_i2j","#Deltaz between reference positions",100,0,0.2);
  htrtrdz_i2jmc = fs->make<TH1F>("htrtrdz_i2jmc","#Deltaz between reference positions",100,0,0.2);
  htrtrdxy_i2j = fs->make<TH1F>("htrtrdxy_i2j","#Deltaxy between reference positions",100,0,0.2);
  htrtrdxy_i2jmc = fs->make<TH1F>("htrtrdxy_i2jmc","#Deltaxy between reference positions",100,0,0.2);
  htrtripd_i2j = fs->make<TH1F>("htrtripd_i2j","distance between impact points",100,0,0.015);
  htrtripd_i2jmc = fs->make<TH1F>("htrtripd_i2jmc","distance between impact points",100,0,0.015);
  htrtripdz_i2j = fs->make<TH1F>("htrtripdz_i2j","#Deltaz between impact points",100,0,0.015);
  htrtripdz_i2jmc = fs->make<TH1F>("htrtripdz_i2jmc","#Deltaz between impact points",100,0,0.015);
  htrtripdxy_i2j = fs->make<TH1F>("htrtripdxy_i2j","#Deltaxy between impact points",100,0,0.015);
  htrtripdxy_i2jmc = fs->make<TH1F>("htrtripdxy_i2jmc","#Deltaxy between impact points",100,0,0.015);
  htrtrsd_i2j = fs->make<TH1F>("htrtrsd_i2j","distance between track lines",100,0,0.015);
  htrtrsd_i2jmc = fs->make<TH1F>("htrtrsd_i2jmc","distance between track lines",100,0,0.015);
  htrtrsdz_i2j = fs->make<TH1F>("htrtrsdz_i2j","#Deltaz between track lines",100,0,0.015);
  htrtrsdz_i2jmc = fs->make<TH1F>("htrtrsdz_i2jmc","#Deltaz between track lines",100,0,0.015);
  htrtrsdxy_i2j = fs->make<TH1F>("htrtrsdxy_i2j","#Deltaxy between track lines",100,0,0.015);
  htrtrsdxy_i2jmc = fs->make<TH1F>("htrtrsdxy_i2jmc","#Deltaxy between track lines",100,0,0.015);

  htStSd_i2j4mu = fs->make<TH1F>("htStSd_i2j4mu","distance between same charged tracks",100,0,0.025);
  htStSdz_i2j4mu = fs->make<TH1F>("htStSdz_i2j4mu","#Deltaz between same charged tracks",100,0,0.025);
  htStSdxy_i2j4mu = fs->make<TH1F>("htStSdxy_i2j4mu","#Deltaxy between same charged tracks",100,0,0.025);
  htPtMd_i2j4mu = fs->make<TH1F>("htPtMd_i2j4mu","distance between opposite charged tracks",100,0,0.025);
  htPtMdz_i2j4mu = fs->make<TH1F>("htPtMdz_i2j4mu","#Deltaz between opposite charged tracks",100,0,0.025);
  htPtMdxy_i2j4mu = fs->make<TH1F>("htPtMdxy_i2j4mu","#Deltaxy between opposite charged tracks",100,0,0.025);

  hntr_isol2j_VS_vt = fs->make<TH2F>("hntr_isol2j_VS_vt","# of muon candidate VS # of vertex",10,0,10,10,0,10);

  hnmcvt = fs->make<TH1F>("hnmcvt","# of vertex per event(MC)",50,0,500);
  hmcvtxy10 = fs->make<TH2F>("hmcvtxy10","x&y position of vertex in xy plane(MC)",300,minmcvtx,maxmcvtx,300,minmcvty,maxmcvty);
//  htrvtxy = fs->make<TH2F>("htrvtxy","position of track vertex in xy plane",300,-0.1,0.1,300,-0.1,0.1);
//  htrvtxy_isol2j = fs->make<TH2F>("htrvtxy_isol2j","position of track vertex in xy plane",300,-0.1,0.1,300,-0.1,0.1);
  hmcmurecovtd = fs->make<TH1F>("hmcmurecovtd","distance between recoVertex & MC muon vertex",100,0,0.05);
  hmcmurecovt2d = fs->make<TH1F>("hmcmurecovt2d","distance between recoVertex & MC muon vertex",100,0,0.05);
  hmcmurecovt3d = fs->make<TH1F>("hmcmurecovt3d","distance between recoVertex & MC muon vertex",100,0,0.05);
  hmcmurecovt2d_ex = fs->make<TH1F>("hmcmurecovt2d_ex","distance between recoVertex & MC muon vertex",100,0.0208,0.0408);
  hmcmurecovt3d_ex = fs->make<TH1F>("hmcmurecovt3d_ex","distance between recoVertex & MC muon vertex",100,0.0208,0.0408);
  hmcmurecovtcd = fs->make<TH1F>("hmcmurecovtcd","distance between recoVertex & MC muon vertex",100,0.0208,0.0408);
  hmcmurecovtc2d = fs->make<TH1F>("hmcmurecovtc2d","distance between recoVertex & MC muon vertex",100,0,0.05);
  hmcmurecovtc3d = fs->make<TH1F>("hmcmurecovtc3d","distance between recoVertex & MC muon vertex",100,0,0.05);
  hmcmuxvtx = fs->make<TH2F>("hmcmuxvtx","MC muon vertex VS recoVertex (x position)",300,minmcvtx,maxmcvtx,300,minvtx,maxvtx);
  hmcmuyvty = fs->make<TH2F>("hmcmuyvty","MC muon vertex VS recoVertex (y position)",300,minmcvty,maxmcvty,300,minvty,maxvty);
  hmcmuzvtz = fs->make<TH2F>("hmcmuzvtz","MC muon vertex VS recoVertex (z position)",300,-15,15,300,-15,15);
  hmcmuxvtx_BS = fs->make<TH2F>("hmcmuxvtx_BS","MC muon vertex VS recoVertex (x position)",300,minmcvtx,maxmcvtx,300,minvtx,maxvtx);
  hmcmuyvty_BS = fs->make<TH2F>("hmcmuyvty_BS","MC muon vertex VS recoVertex (y position)",300,minmcvty,maxmcvty,300,minvty,maxvty);
  hmcmuzvtz_BS = fs->make<TH2F>("hmcmuzvtz_BS","MC muon vertex VS recoVertex (z position)",300,-15,15,300,-15,15);
//  hmcmurecovtxy = fs->make<TH2F>("hmcmurecovtxy","Vertex of MC muons & recoVertex (x&y position)",300,minmcvtx,maxmcvtx,300,minmcvty,maxmcvty);
//  hmcmurecovt2xy = fs->make<TH2F>("hmcmurecovt2xy","Vertex of MC muons & recoVertex (x&y position)",300,minmcvtx,maxmcvtx,300,minmcvty,maxmcvty);

  hmcmuxcpx = fs->make<TH2F>("hmcmuxcpx","MC muon vertex VS closest point between two same charged muon candidates (x position)",300,minmcvtx,maxmcvtx,300,minvtx,maxvtx);
  hmcmuycpy = fs->make<TH2F>("hmcmuycpy","MC muon vertex VS closest point between two same charged muon candidates (y position)",300,minmcvty,maxmcvty,300,minvty,maxvty);
  hmcmuzcpz = fs->make<TH2F>("hmcmuzcpz","MC muon vertex VS closest point between two same charged muon candidates (z position)",300,-15,15,300,-15,15);
  hvtxcpx_BS = fs->make<TH2F>("hvtxcpx_BS","recoVertex VS closest point between two same charged muon candidates (x position)",300,minmcvtx,maxmcvtx,300,minvtx,maxvtx);
  hvtycpy_BS = fs->make<TH2F>("hvtycpy_BS","recoVertex VS closest point between two same charged muon candidates (y position)",300,minmcvty,maxmcvty,300,minvty,maxvty);
  hvtzcpz_BS = fs->make<TH2F>("hvtzcpz_BS","recoVertex VS closest point between two same charged muon candidates (z position)",300,-15,15,300,-15,15);

  hvtchi2_BS = fs->make<TH1F>("hvtchi2_BS","normalized chi2 of recoVertex",100,0,2);
  hvtchi2_BS3mu = fs->make<TH1F>("hvtchi2_BS3mu","normalized chi2 of recoVertex",100,0,2);
  hvtchi2_BS4mu = fs->make<TH1F>("hvtchi2_BS4mu","normalized chi2 of recoVertex",100,0,2);
  hvt2chi2_BS = fs->make<TH1F>("hvt2chi2_BS","normalized chi2 of recoVertex",100,0,2);
  hvt3chi2_BS = fs->make<TH1F>("hvt3chi2_BS","normalized chi2 of recoVertex",100,0,2);

  hevents = fs->make<TH1F>("hevents","number of events for crab",2,0,2);

}


HppMuMuAnalyzer::~HppMuMuAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called to for each event  ------------
void HppMuMuAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  events++; hevents->Fill(1); if(int(events)%100==0) cout<<"Event "<<events<<endl;

  Handle<TrackCollection> TrCollection;
  Handle<TrackCollection> MuCollection;
  Handle<TrackCollection> StaCollection;
  Handle<GenParticleCollection> GenCollection;
//  iEvent.getByLabel("ctfWithMaterialTracks", TrCollection);
//  iEvent.getByLabel(Tracks_, TrCollection);
  iEvent.getByLabel("generalTracks", TrCollection);
  iEvent.getByLabel("globalMuons", MuCollection);
  iEvent.getByLabel("standAloneMuons", StaCollection);
  iEvent.getByLabel("genParticles", GenCollection);
  hntr->Fill(TrCollection->size());
  hngmu->Fill(MuCollection->size());

  Handle<CaloJetCollection> caloJets;
  Handle<CaloJetCollection> corJets;
  Handle<GenJetCollection> genJets;
//  Handle<JetCollection> caloJets;
//  Handle<JetCollection> corJets;
//  Handle<JetCollection> genJets;
  iEvent.getByLabel("iterativeCone5CaloJets",caloJets);
  iEvent.getByLabel(CaloJetAlgorithm_,corJets);
  iEvent.getByLabel("iterativeCone5GenJets",genJets);
//  iEvent.getByLabel(GenJetAlgorithm_,genJets);
//  const JetCorrector* corrector = JetCorrector::getJetCorrector (JetCorrectionService_, iSetup);
  hncalojet->Fill(caloJets->size());
  hncorjet->Fill(corJets->size());
  hngenjet->Fill(genJets->size());

  Handle<BeamSpot> Bspot;
  iEvent.getByLabel("offlineBeamSpot", Bspot);

  double x0 = Bspot->x0();
  double y0 = Bspot->y0();
  double z0 = Bspot->z0();
  double sigmaz = Bspot->sigmaZ();
//  double dxdz = Bspot->dxdz();
  double Bwidth = Bspot->BeamWidth();
//  double BwidthX = Bspot->BeamWidthX();
//  double BwidthY = Bspot->BeamWidthY();
  hbeamXY2D->Fill(x0,y0);
  hbeamz->Fill(z0);
  hbeamsigmaz->Fill(sigmaz);
  if(events==1) cout<<"BeamX="<<x0<<", BeamY="<<y0<<", BeamZ="<<z0<<", BeamsigmaZ="<<sigmaz<<", BeamWidth="<<Bwidth<<endl;

//  Handle<VertexCollection> VtCollection;
//  Handle<VertexCollection> VtCollection2;
//  Handle<VertexCollection> VtCollection3;
  edm::Handle<reco::VertexCollection> VtCollection;
  edm::Handle<reco::VertexCollection> VtCollection2;
  edm::Handle<reco::VertexCollection> VtCollection3;
  iEvent.getByLabel("offlinePrimaryVertices",VtCollection);
  iEvent.getByLabel("pixelVertices",VtCollection2);
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",VtCollection3);

//---------------------------MC particles(mc000, line000)-----------------------------//
  TLorentzVector mcmu[4], mcH[2];
  int Hidx=0, Muidx=0, MuidxP=0, MuidxM=0;
  int nmcvt=0, nmcvtmax=500;
//  double  mcvtx[nmcvtmax]={0,}, mcvty[nmcvtmax]={0,};
  double  mcvtx[500]={0,}, mcvty[500]={0,};
  const Candidate * mcmu1=0;
  for (GenParticleCollection::const_iterator gen=GenCollection->begin(); gen!=GenCollection->end(); gen++)
  {
//if(debug) cout<<"1"<<endl;
    const Candidate * mom  = gen->mother();
    const Candidate * mom2 = 0;
    int momid=0, mom2id=0;
    int ndau = gen->numberOfDaughters();
    if(mom!=0) {momid = mom->pdgId(); mom2 = (*mom).mother();}
    if(mom!=0 && mom2!=0) mom2id = mom2->pdgId();
    int id = gen->pdgId();
//    cout<<mom2id<<", "<<momid<<", "<<id<<" : ";
//    for(int i=0;i<ndau;i++) cout<<gen->daughter(i)->pdgId()<<", "; cout<<endl;
    if( (id==9900041||id==-9900041||id==9900042||id==-9900042) && ndau==3) 
    {
      int n=0,m=0;
      if(id<0) n=2;
      for(int i=n;i<(n+2);i++,m++) 
      {
        const Candidate *Hdau = gen->daughter(m);
        mcmu[i].SetPxPyPzE(Hdau->px(),Hdau->py(),Hdau->pz(),Hdau->p());
//        if(Hdau->pt()>10 && fabs(Hdau->eta())<2.4) hmcmuxy->Fill(Hdau->vx(),Hdau->vy());
        if(Muidx==0)
        {
//          hmcmurecovtxy->Fill(Hdau->vx(),Hdau->vy()); hmcmurecovt2xy->Fill(Hdau->vx(),Hdau->vy()); hmcmurecovt3xy->Fill(Hdau->vx(),Hdau->vy());
          mcmu1 = &(*Hdau);

//          if(!((tmcmu-tdvt).Px()==x0 && (tmcmu-tdvt).Py()==y0 && (tmcmu-tdvt).Pz()==z0)) hmcmuzvtz_BS->Fill(tmcmu.Pz(),(tmcmu-tdvt).Pz());
//          if(vtxon) if((tmcmu-tdvt).Px()==x0 && (tmcmu-tdvt).Py()==y0 && (tmcmu-tdvt).Pz()==z0) cout<<"hello"<<endl;
//cout<<dvt(&(*Hdau),VtCollection).P()<<", "<<dvt(&(*Hdau),VtCollection2).P()<<", "<<dvt(&(*Hdau),VtCollection3).P()<<", "<<dvt_cor(&(*Hdau),VtCollection,x0,y0,z0).P()<<", "<<dvt_cor(&(*Hdau),VtCollection2,x0,y0,z0).P()<<", "<<dvt_cor(&(*Hdau),VtCollection3,x0,y0,z0).P()<<endl;
        }
        if(Muidx==0 && events<=1000) hmcmuxy_upper1k->Fill(Hdau->vx(),Hdau->vy());
        if(Muidx==0 && int(events)%1000<=100) hmcmuxy_10times1k->Fill(Hdau->vx(),Hdau->vy());
        Muidx++;
//if(VtCollection->size()==2) cout<<Hdau->vx()<<", "<<Hdau->vy()<<", ";
      }//mcmu[0]=mu+ from H++, mcmu[1]=mu+ from H++, mcmu[2]=mu- from H--, mcmu[3]=mu- from H--
if(debug) cout<<"2"<<endl;
      mcH[Hidx].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->p());
      Hidx++;
if(debug) cout<<"3"<<endl;
//if(VtCollection->size()==2) cout<<endl;
//cout<<gen->vx()<<", "<<gen->vy()<<endl;
    }
    if(abs(id)==13 && momid!=9900041 && momid!=-9900041 && momid!=9900042 && momid!=-9900042 && abs(momid)!=13)
    {
      nbackmcmuon++;
      if(gen->pt()>10 && fabs(gen->eta())<2.4) nbackmcmuon_ptetacut++;
    }
    if(fabs(gen->eta())<2.4 && gen->charge()!=0 && gen->pt()>10 && ndau==0)
    {
      TLorentzVector tgen(gen->px(),gen->py(),gen->pz(),gen->p());
//      hmcisolPt_ptetacut01->Fill(isolPt_mc(tgen,GenCollection));
//      hmcisolPt_ptetacut1->Fill(isolPt_mc(tgen,GenCollection));
      hmcisolPt_ptetacut10->Fill(isolPt_mc(tgen,GenCollection));
      hmcisolPt_ptetacut100->Fill(isolPt_mc(tgen,GenCollection));
//      hmcisolPt_ptetacut1000->Fill(isolPt_mc(tgen,GenCollection));
      hmcisolPt2_ptetacut10->Fill(isolPt2_mc(tgen,GenCollection));
      hmcisolPt2_ptetacut100->Fill(isolPt2_mc(tgen,GenCollection));
    }
    bool vtmatching=true; 
    if(ndau==0 && gen->charge()!=0 && nmcvt<nmcvtmax)
    {
      for(int i=0;i<=nmcvt;i++)
      {
        if(mcvtx[i]!=gen->vx() && mcvty[i]!=gen->vy()) vtmatching=false;
      }
      if(!vtmatching)
      {
        mcvtx[nmcvt]=gen->vx(); mcvty[nmcvt]=gen->vy();
        nmcvt++;
      }
    }
/*    if(id==23 && ndau==3) 
    {
      int n=4,m=0;
      if(Muidx==2) n=3;
//      for(int i=n;i<(n+2);i+=2,m++) 
      for(int i=n;i<4;i-=2,m++) 
      {
        const Candidate *Zdau = gen->daughter(m);
        mcmu[i].SetPxPyPzE(Zdau->px(),Zdau->py(),Zdau->pz(),Zdau->p());
//        if(Hdau->pt()>10 && fabs(Hdau->eta())<2.4) hmcmuxy->Fill(Hdau->vx(),Hdau->vy());
        if(Muidx==0) mcmu1 = &(*gen);
        Muidx++;
      }
    }
*/
/*    if(mom2id==23 && abs(id)==13)
    {
      if(id<0) {mcmu[MuidxP].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->p()); MuidxP++;}
      if(id>0) {mcmu[MuidxM].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->p()); MuidxM++;}
      if(Muidx==0) mcmu1 = &(*gen);
      Muidx++;
    }
*/
    if(id==23 && ndau==3) 
    {
//      int n=4,m=0;
//      if(Muidx==2) n=3;
//      for(int i=n;i<(n+2);i+=2,m++) 
      for(int i=0;i<2;i++) 
      {
        const Candidate *Zdau = gen->daughter(i);
//        mcmu[i].SetPxPyPzE(Zdau->px(),Zdau->py(),Zdau->pz(),Zdau->p());
        if(Zdau->pdgId()<0) {mcmu[MuidxP].SetPxPyPzE(Zdau->px(),Zdau->py(),Zdau->pz(),Zdau->p()); MuidxP++;}
        if(Zdau->pdgId()>0) {mcmu[MuidxM+2].SetPxPyPzE(Zdau->px(),Zdau->py(),Zdau->pz(),Zdau->p()); MuidxM++;}
//        if(Hdau->pt()>10 && fabs(Hdau->eta())<2.4) hmcmuxy->Fill(Hdau->vx(),Hdau->vy());
        if(Muidx==0) mcmu1 = &(*Zdau);
        Muidx++;
//cout<<"id : "<<Zdau->pdgId()<<endl;
      }
      mcH[Hidx].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->p());
      Hidx++;
    }
  }
if(debug) cout<<"4"<<endl;
if(debug) for(int i=0;i<4;i++) cout<<"mcmu pt : "<<mcmu[i].Pt()<<endl;
  hnmcvt->Fill(nmcvt);
  if(events<=10) for(int i=0;i<nmcvt;i++) hmcvtxy10->Fill(mcvtx[i],mcvty[i]);

  hmcmuxy->Fill(mcmu1->vx(),mcmu1->vy());
  hmcmuz->Fill(mcmu1->vz());
  TLorentzVector tmcmu(mcmu1->vx(),mcmu1->vy(),mcmu1->vz(),1);
  TLorentzVector tdvt; if(vtxon) tdvt=dvt(&(*mcmu1),VtCollection);
  TLorentzVector tdvt2; if(vtxon) tdvt2=dvt_BS(&(*mcmu1),VtCollection,x0,y0,z0);
  hmcmurecovtd->Fill(tdvt.P());
  if(vtxon) hmcmurecovt2d->Fill(dvt(&(*mcmu1),VtCollection2).P());
  if(vtxon) hmcmurecovt3d->Fill(dvt(&(*mcmu1),VtCollection3).P());
  if(vtxon) hmcmurecovt2d_ex->Fill(dvt(&(*mcmu1),VtCollection2).P());
  if(vtxon) hmcmurecovt3d_ex->Fill(dvt(&(*mcmu1),VtCollection3).P());
  if(vtxon) hmcmurecovtcd->Fill(dvt_cor(&(*mcmu1),VtCollection,x0,y0).P());
  if(vtxon) hmcmurecovtc2d->Fill(dvt_cor(&(*mcmu1),VtCollection2,x0,y0).P());
  if(vtxon) hmcmurecovtc3d->Fill(dvt_cor(&(*mcmu1),VtCollection3,x0,y0).P());
  hmcmuxvtx->Fill(tmcmu.Px(),(tmcmu-tdvt).Px());
  hmcmuyvty->Fill(tmcmu.Py(),(tmcmu-tdvt).Py());
  hmcmuzvtz->Fill(tmcmu.Pz(),(tmcmu-tdvt).Pz());
  if(tdvt2.P()!=0)
  {
    hmcmuxvtx_BS->Fill(tmcmu.Px(),(tmcmu-tdvt2).Px());
    hmcmuyvty_BS->Fill(tmcmu.Py(),(tmcmu-tdvt2).Py());
    hmcmuzvtz_BS->Fill(tmcmu.Pz(),(tmcmu-tdvt2).Pz());
  }

  int mcmuidx_ptetacut=0, mcmuidx_ptetaisolcut=0, mcmui[4]={0,};
  int mcmuidx_ptetaisol2cut=0, mcmui2[4]={0,}, mcmui2j[4]={0,};
  float mcmuisolPt[4]={0,}, mcmuisolPt2[4]={0,};
  for(int i=0;i<4;i++)
  {
    hmcmupt->Fill(mcmu[i].Pt());
    hmcmueta->Fill(mcmu[i].Eta());
    if(mcmu[i].Pt()>10 && fabs(mcmu[i].Eta())<2.4)
    {
      hmcmupt_ptetacut->Fill(mcmu[i].Pt());
      hmcmueta_ptetacut->Fill(mcmu[i].Eta());
      mcmuidx_ptetacut++;
      mcmuisolPt[i]=isolPt_mc(mcmu[i],GenCollection);
//      hmcmuisolPt_ptetacut01->Fill(mcmuisolPt[i]);
//      hmcmuisolPt_ptetacut1->Fill(mcmuisolPt[i]);
      hmcmuisolPt_ptetacut10->Fill(mcmuisolPt[i]);
      hmcmuisolPt_ptetacut100->Fill(mcmuisolPt[i]);
//      hmcmuisolPt_ptetacut1000->Fill(mcmuisolPt[i]);
      mcmuisolPt2[i]=isolPt2_mc(mcmu[i],GenCollection);
      hmcmuisolPt2_ptetacut10->Fill(mcmuisolPt2[i]);
      hmcmuisolPt2_ptetacut100->Fill(mcmuisolPt2[i]);
      if(mcmuisolPt[i]<isolPtcut)
      {
        hmcmupt_ptetaisolcut->Fill(mcmu[i].Pt());
        hmcmueta_ptetaisolcut->Fill(mcmu[i].Eta());
        mcmuidx_ptetaisolcut++;
        mcmui[i]=1;
      }
      if(mcmuisolPt2[i]<isolPtcut)
      {
        hmcmupt_ptetaisol2cut->Fill(mcmu[i].Pt());
        hmcmueta_ptetaisol2cut->Fill(mcmu[i].Eta());
        mcmuidx_ptetaisol2cut++;
        mcmui2[i]=1;
if(debug) cout<<"10"<<endl;
        if(jetisol(mcmu[i],corJets))
        {
          mcmui2j[i]=1;
          hmcmupt_ptetaisol2jcut->Fill(mcmu[i].Pt());
          hmcmueta_ptetaisol2jcut->Fill(mcmu[i].Eta());
        }
if(debug) cout<<"20"<<endl;
      }
      float dR_mcmuons=10, dR_mcmuons1=10, dR_mcmuons2=10, dR_mcmuons10=10, dR_mcmuons20=10;
      for(int n=0;n<4;n++)
      {
        if(n!=i && mcmu[n].Pt()>10 && fabs(mcmu[n].Eta())<2.4)
        {
          float dR=mcmu[i].DeltaR(mcmu[n]);
          if(mcmuisolPt[i]>10 && dR<dR_mcmuons1) dR_mcmuons1=dR;
          if(dR_mcmuons1<10) hmcmudR_upperisolPt_10->Fill(dR_mcmuons1);
          if(mcmuisolPt[i]>20 && dR<dR_mcmuons2) dR_mcmuons2=dR;
          if(dR_mcmuons2<10) hmcmudR_upperisolPt_20->Fill(dR_mcmuons2);

          if(mcmuisolPt2[i]>10 && dR<dR_mcmuons10) dR_mcmuons10=dR;
          if(dR_mcmuons10<10) hmcmudR_upperisolPt2_10->Fill(dR_mcmuons10);
          if(mcmuisolPt2[i]>20 && dR<dR_mcmuons20) dR_mcmuons20=dR;
          if(dR_mcmuons20<10) hmcmudR_upperisolPt2_20->Fill(dR_mcmuons20);

          if(dR<dR_mcmuons) dR_mcmuons=dR;
        }
      }
      if(dR_mcmuons<0.3) {hmcmuisolPt_mcmudR03->Fill(mcmuisolPt[i]); hmcmuisolPt2_mcmudR03->Fill(mcmuisolPt2[i]);}
      hmcmudR_allmu->Fill(dR_mcmuons);
    }
  }
//for(int i=0;i<4;i++) cout<<mcmu[i].Pt()<<", "<<mcmu[i].Eta()<<", ";
//cout<<endl;
  hnmcmu_ptetacut->Fill(mcmuidx_ptetacut);
  hnmcmu_ptetaisolcut->Fill(mcmuidx_ptetaisolcut);
  hnmcmu_ptetaisol2cut->Fill(mcmuidx_ptetaisol2cut);
  hnmcmu_ptetaisol2jcut->Fill(mcmui2j[0]+mcmui2j[1]+mcmui2j[2]+mcmui2j[3]);
  for(int i=0,h=0;i<3;i+=2,h++)
  {
    if(mcmu[i].Pt()>10 && fabs(mcmu[i].Eta())<2.4 && mcmu[i+1].Pt()>10 && fabs(mcmu[i+1].Eta())<2.4)
    {
      hHPPmass_mc->Fill((mcmu[i]+mcmu[i+1]).M());
      hmcmudeltaphi->Fill(fabs(deltaphi(mcmu[i].Phi(),mcmu[i+1].Phi())));
      hmcmudeltaphi2->Fill(fabs(deltaphi(mcmu[h].Phi(),mcmu[h+2].Phi())));
      hmcmudR->Fill(mcmu[i].DeltaR(mcmu[i+1]));
      hmcmudR_ex->Fill(mcmu[i].DeltaR(mcmu[i+1]));
    }
  }
  if(mcmu[0].Pt()>10 && fabs(mcmu[0].Eta())<2.4 && mcmu[1].Pt()>10 && fabs(mcmu[1].Eta())<2.4)
  if(mcmu[2].Pt()>10 && fabs(mcmu[2].Eta())<2.4 && mcmu[2].Pt()>10 && fabs(mcmu[3].Eta())<2.4)
  {
    hHPPHMMmass_mc->Fill((mcmu[0]+mcmu[1]+mcmu[2]+mcmu[3]).M());
    h2mcmudeltaphi->Fill(fabs(deltaphi((mcmu[0]+mcmu[1]).Phi(),(mcmu[2]+mcmu[3]).Phi())));
    if(mcH[0].P()!=0) hmcHdeltaphi->Fill(fabs(deltaphi(mcH[0].Phi(),mcH[1].Phi())));
  }
  if(mcmui[0]&&mcmui[1]) hHPPmass_mcisol->Fill((mcmu[0]+mcmu[1]).M());
  if(mcmui[2]&&mcmui[3]) hHMMmass_mcisol->Fill((mcmu[2]+mcmu[3]).M());
  if(mcmui[0]&&mcmui[1]&&mcmui[2]&&mcmui[3])
  {
    hHPPHMMmass_mcisol->Fill((mcmu[0]+mcmu[1]+mcmu[2]+mcmu[3]).M());
    h2mcmudeltaphi_isol->Fill(fabs(deltaphi((mcmu[0]+mcmu[1]).Phi(),(mcmu[2]+mcmu[3]).Phi())));
  }
  if(mcmui2[0]&&mcmui2[1]) hHPPmass_mcisol2->Fill((mcmu[0]+mcmu[1]).M());
  if(mcmui2[2]&&mcmui2[3]) hHMMmass_mcisol2->Fill((mcmu[2]+mcmu[3]).M());
  if(mcmui2[0]&&mcmui2[1]&&mcmui2[2]&&mcmui2[3])
  {
    hHPPHMMmass_mcisol2->Fill((mcmu[0]+mcmu[1]+mcmu[2]+mcmu[3]).M());
    h2mcmudeltaphi_isol2->Fill(fabs(deltaphi((mcmu[0]+mcmu[1]).Phi(),(mcmu[2]+mcmu[3]).Phi())));
    if(fabs(deltaphi((mcmu[0]+mcmu[1]).Phi(),(mcmu[2]+mcmu[3]).Phi()))>dphicut_4mu)
    {
      hHPPmass_mcisol2dphi->Fill((mcmu[0]+mcmu[1]).M());
      hHMMmass_mcisol2dphi->Fill((mcmu[2]+mcmu[3]).M());
      hHPPHMMmass_mcisol2dphi->Fill((mcmu[0]+mcmu[1]+mcmu[2]+mcmu[3]).M());
    }
  }
if(debug) cout<<"50"<<endl;
  if(mcmui2j[0]&&mcmui2j[1]) hHPPmass_mcisol2j->Fill((mcmu[0]+mcmu[1]).M());
  if(mcmui2j[2]&&mcmui2j[3]) hHMMmass_mcisol2j->Fill((mcmu[2]+mcmu[3]).M());
  if(mcmui2j[0]&&mcmui2j[1]&&mcmui2j[2]&&mcmui2j[3])
  {
    hHPPHMMmass_mcisol2j->Fill((mcmu[0]+mcmu[1]+mcmu[2]+mcmu[3]).M());
  }
if(debug) cout<<"60"<<endl;

//---------------------------jet000-----------------------------//
  int ncalojet_ptcut=0;
  for( CaloJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); ++cal )
  {
    TLorentzVector tcal(cal->px(),cal->py(),cal->pz(),cal->p());
    hcalojetpt->Fill(cal->pt());
    hcalojeteta->Fill(cal->eta());
    if(cal->pt()>15)
    {
      hcalojetpt_ptcut->Fill(cal->pt()); hcalojetpt_ptcut_ex->Fill(cal->pt());
      hcalojeteta_ptcut->Fill(cal->eta());
      ncalojet_ptcut++;
    }
  }
  hncalojet_ptcut->Fill(ncalojet_ptcut);

if(debug) cout<<"61"<<endl;
  int ncorjet_ptcut=0, ncorjet_ptmucut=0;
  for( CaloJetCollection::const_iterator cor = corJets->begin(); cor != corJets->end(); ++cor )
  {
    TLorentzVector tcor(cor->px(),cor->py(),cor->pz(),cor->p());
    hcorjetpt->Fill(cor->pt());
    hcorjeteta->Fill(cor->eta());
    if(cor->pt()>15)
    {
      hcorjetpt_ptcut->Fill(cor->pt()); hcorjetpt_ptcut_ex->Fill(cor->pt());
      hcorjeteta_ptcut->Fill(cor->eta());
      ncorjet_ptcut++;
      int ncorjet_mcmumatch=0;
      for(int i=0;i<4;i++) if(mcmu[i].Pt()>10 && mcmu[i].DeltaR(tcor)<0.5) ncorjet_mcmumatch++;
      if(ncorjet_mcmumatch==0)
      {
        hcorjetpt_ptmucut->Fill(cor->pt()); hcorjetpt_ptmucut_ex->Fill(cor->pt());
        hcorjeteta_ptmucut->Fill(cor->eta());
        ncorjet_ptmucut++;
      }
    }
    if(cor->pt()>15 && fabs(cor->eta())<2.4)
    {
      float dR_cormcmu=100;
      int mcmuidx=10;
      for(int i=0;i<4;i++)
      {
        if(mcmu[i].Pt()>10 && mcmu[i].DeltaR(tcor)<dR_cormcmu)
        {
          dR_cormcmu=mcmu[i].DeltaR(tcor);
          mcmuidx=i;
        }
      }
      hcormcmuptdR_ptetacut->Fill(cor->pt()/mcmu[mcmuidx].Pt(),tcor.DeltaR(mcmu[mcmuidx]));

    }
    if(cor->pt()>30 && fabs(cor->eta())<2.4)
    {
      float dR_cormcmu=100;
      int mcmuidx=10;
      for(int i=0;i<4;i++)
      {
        if(mcmu[i].Pt()>10 && mcmu[i].DeltaR(tcor)<dR_cormcmu)
        {
          dR_cormcmu=mcmu[i].DeltaR(tcor);
          mcmuidx=i;
        }
      }
      hcormcmuptdR_30ptetacut->Fill(cor->pt()/mcmu[mcmuidx].Pt(),tcor.DeltaR(mcmu[mcmuidx]));
    }
  }
  hncorjet_ptcut->Fill(ncorjet_ptcut);
  hncorjet_ptmucut->Fill(ncorjet_ptmucut);

if(debug) cout<<"62"<<endl;
  int ngenjet_ptcut=0, ngenjet_ptmucut=0;
  for( GenJetCollection::const_iterator gen = genJets->begin(); gen != genJets->end(); ++gen )
  {
    TLorentzVector tgen(gen->px(),gen->py(),gen->pz(),gen->p());
    hgenjetpt->Fill(gen->pt());
    hgenjeteta->Fill(gen->eta());
    if(gen->pt()>15)
    {
      hgenjetpt_ptcut->Fill(gen->pt()); hgenjetpt_ptcut_ex->Fill(gen->pt());
      hgenjeteta_ptcut->Fill(gen->eta());
      ngenjet_ptcut++;
      int ngenjet_mcmumatch=0;
      for(int i=0;i<4;i++) if(mcmu[i].Pt()>10 && mcmu[i].DeltaR(tgen)<0.5) ngenjet_mcmumatch++;
      if(ngenjet_mcmumatch==0)
      {
        hgenjetpt_ptmucut->Fill(gen->pt()); hgenjetpt_ptmucut_ex->Fill(gen->pt());
        hgenjeteta_ptmucut->Fill(gen->eta());
        ngenjet_ptmucut++;
        if(fabs(gen->eta())<2.4)
        {
          hgenjetpt_pt24etamucut->Fill(gen->pt());
          hgenjeteta_pt24etamucut->Fill(gen->eta());
          if(gen->pt()>20)
          {
            hgenjetpt_20pt24etamucut->Fill(gen->pt());
            hgenjeteta_20pt24etamucut->Fill(gen->eta());
          }
          if(gen->pt()>30)
          {
            hgenjetpt_30pt24etamucut->Fill(gen->pt());
            hgenjeteta_30pt24etamucut->Fill(gen->eta());
          }
        }
        bool calojetmatch=false;
        for( CaloJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); ++cal )
        {
          TLorentzVector tcal(cal->px(),cal->py(),cal->pz(),cal->p());
          if(cal->pt()>15 && tgen.DeltaR(tcal)<0.2)
          {
            calojetmatch=true;
          }
        }
        if(calojetmatch)
        {
          hgenjetpt_ptmucut_calojetmatch->Fill(gen->pt());
                                  hgenjeteta_ptmucut_calojetmatch->Fill(gen->eta());
        }

if(debug) cout<<"63"<<endl;
        bool corjetmatch=false, corjetmatch_24etacut=false, corjetmatch_20pt24etacut=false, corjetmatch_30pt24etacut=false;
        float dR_24etamucut=0.2, dR_30pt24etamucut=0.2;
        float ptr_24etamucut=0, ptr_30pt24etamucut=0;
        for( CaloJetCollection::const_iterator cor = corJets->begin(); cor != corJets->end(); ++cor )
        {
          TLorentzVector tcor(cor->px(),cor->py(),cor->pz(),cor->p());
          int ncorjet_mcmumatch=0;
          for(int i=0;i<4;i++) if(mcmu[i].Pt()>10 && mcmu[i].DeltaR(tcor)<0.5 && cor->pt()>15) ncorjet_mcmumatch++;
          if(cor->pt()>15 && tgen.DeltaR(tcor)<0.2 && ncorjet_mcmumatch==0)
          {
            corjetmatch=true;
            if(fabs(gen->eta())<2.4 && fabs(cor->eta())<2.4)
            {
              corjetmatch_24etacut=true;
              if(tgen.DeltaR(tcor)<dR_24etamucut)
              {
                dR_24etamucut=tgen.DeltaR(tcor);
                ptr_24etamucut=tcor.Pt()/tgen.Pt();
              }
              if(gen->pt()>20 && cor->pt()>20)
              {
                corjetmatch_20pt24etacut=true;
              }
              if(gen->pt()>30 && cor->pt()>30)
              {
                corjetmatch_30pt24etacut=true;
                if(tgen.DeltaR(tcor)<dR_30pt24etamucut)
                {
                  dR_30pt24etamucut=tgen.DeltaR(tcor);
                  ptr_30pt24etamucut=tcor.Pt()/tgen.Pt();
                }
              }
            }
          }
        }
        if(corjetmatch)
        {
          hgenjetpt_ptmucut_corjetmatch->Fill(gen->pt());
                                  hgenjeteta_ptmucut_corjetmatch->Fill(gen->eta());
        }
        if(corjetmatch_24etacut)
        {
          hgenjetpt_pt24etamucut_corjetmatch->Fill(gen->pt());
                                  hgenjeteta_pt24etamucut_corjetmatch->Fill(gen->eta());
          if(tgen.Pt()<100)
          {
            int ii = int(tgen.Pt()/2);
            Tptr_24etamucut[ii]+=ptr_24etamucut;
            Tptr2_24etamucut[ii]+=pow(ptr_24etamucut,2);
            ncorgenjetpt_24etamucut[ii]++;
          }
          if(fabs(tgen.Eta())<2.4)
          {
            int jj = int((tgen.Eta()+2.4)*50/4.8);
            Teta_24etamucut[jj]+=ptr_24etamucut;
            Teta2_24etamucut[jj]+=pow(ptr_24etamucut,2);
            ncorgenjeteta_24etamucut[jj]++;
          }
        }
        if(corjetmatch_20pt24etacut)
        {
          hgenjetpt_20pt24etamucut_corjetmatch->Fill(gen->pt());
                                  hgenjeteta_20pt24etamucut_corjetmatch->Fill(gen->eta());
        }
        if(corjetmatch_30pt24etacut)
        {
          hgenjetpt_30pt24etamucut_corjetmatch->Fill(gen->pt());
                                  hgenjeteta_30pt24etamucut_corjetmatch->Fill(gen->eta());
          if(tgen.Pt()<100)
          {
            int ii = int(tgen.Pt()/2);
            Tptr_30pt24etamucut[ii]+=ptr_30pt24etamucut;
            Tptr2_30pt24etamucut[ii]+=pow(ptr_30pt24etamucut,2);
            ncorgenjetpt_30pt24etamucut[ii]++;
          }
          if(fabs(tgen.Eta())<2.4)
          {
            int jj = int((tgen.Eta()+2.4)*50/4.8);
            Teta_30pt24etamucut[jj]+=ptr_30pt24etamucut;
            Teta2_30pt24etamucut[jj]+=pow(ptr_30pt24etamucut,2);
            ncorgenjeteta_30pt24etamucut[jj]++;
          }
        }
      }
    }
if(debug) cout<<"64"<<endl;
    if(gen->pt()>15 && fabs(gen->eta())<2.4)
    {
      float dR_genmcmu=100;
      int mcmuidx=10;
      for(int i=0;i<4;i++)
      {
        if(mcmu[i].Pt()>10 && mcmu[i].DeltaR(tgen)<dR_genmcmu)
        {
          dR_genmcmu=mcmu[i].DeltaR(tgen);
          mcmuidx=i;
        }
      }
      hgenmcmuptdR_ptetacut->Fill(gen->pt()/mcmu[mcmuidx].Pt(),tgen.DeltaR(mcmu[mcmuidx]));
    }
    if(gen->pt()>30 && fabs(gen->eta())<2.4)
    {
      float dR_genmcmu=100;
      int mcmuidx=10;
      for(int i=0;i<4;i++)
      {
        if(mcmu[i].Pt()>10 && mcmu[i].DeltaR(tgen)<dR_genmcmu)
        {
          dR_genmcmu=mcmu[i].DeltaR(tgen);
          mcmuidx=i;
        }
      }
      hgenmcmuptdR_30ptetacut->Fill(gen->pt()/mcmu[mcmuidx].Pt(),tgen.DeltaR(mcmu[mcmuidx]));
    }
  }
if(debug) cout<<"65"<<endl;
  hngenjet_ptcut->Fill(ngenjet_ptcut);
  hngenjet_ptmucut->Fill(ngenjet_ptmucut);

  for( VertexCollection::const_iterator vt = VtCollection->begin(); vt != VtCollection->end();++vt )
  {
//    cout<<vt->position()<<"1 ";
//    cout<<"recoVertex : "<<vt->x()<<", "<<vt->y()<<", "<<vt->z()<<endl;
  }
//---------------------------tracks(tr000, line000)-----------------------------//
  std::vector<ptvsindex> trptsorted, trptsortedP, trptsortedM;
  std::vector<ptvsindex> trptsortedi, trptsortedPi, trptsortedMi;
  std::vector<ptvsindex> trptsortedi2, trptsortedPi2, trptsortedMi2;
  std::vector<ptvsindex> trptsortedi2j, trptsortedPi2j, trptsortedMi2j;
  std::vector<ptvsindex> trptsortedi2jv, trptsortedPi2jv, trptsortedMi2jv;
  std::vector<ptvsindex> trptsortedi2jvc, trptsortedPi2jvc, trptsortedMi2jvc;
  int tridx=0, tridx2=0, tridx_ptetacut=0, tridx_ptetaisolcut=0, tridx_ptetaisol2cut=0, tridx_ptetaisol2jcut=0;
  int tridx_isol2jv=0, tridx_isol2jvc=0;
//cout<<"track.P : "<<tT(TrCollection).P()<<endl;
//cout<<"nVt : "<<VtCollection->size()<<endl;
//cout<<"nVt3 : "<<VtCollection3->size()<<endl;
//cout<<"nVt2 : "<<VtCollection2->size()<<endl;
//if(vtxon) cout<<"reoVertex.D : "<<tV(VtCollection).Px()<<", "<<tV(VtCollection).Py()<<", "<<tV(VtCollection).Pz()<<", "<<tV(VtCollection).P()<<endl;
  for (TrackCollection::const_iterator tr=TrCollection->begin(); tr!=TrCollection->end(); tr++, tridx++)
  {
if(debug && tridx==0) cout<<"66"<<endl;
//cout<<"track pt : "<<tr->pt()<<", track phi : "<<tr->phi()<<", track vx : "<<tr->vx()<<", track vy : "<<tr->vy()<<", track vz : "<<tr->vz()<<endl;
    htrackpt->Fill(tr->pt()); 
    htracketa->Fill(tr->eta());
    htrackz->Fill(tr->dz());
    htrackz2->Fill( tr->dxy()+x0*sin(tr->phi())-y0*cos(tr->phi()) );
//cout<<tr->dz()<<", "<<( tr->vz() - (tr->vx()*tr->px()+tr->vy()*tr->py())/tr->pt() * (tr->pz()/tr->pt()) )<<", "<<( (tr->vz()-z0) - ( (tr->vx()-x0)*tr->px()+(tr->vy()-y0)*tr->py())/tr->pt() * ((tr->pz()-z0)/tr->pt()) )<<endl;
//cout<<tr->dz()<<", "<<( tr->vz() - (tr->vx()*tr->px()+tr->vy()*tr->py())/tr->pt() * (tr->pz()/tr->pt()) )<<endl;
    htrackxy->Fill(tr->dxy());
    htracksz->Fill(tr->dsz());
    htrackX->Fill(tr->vx()); htrackY->Fill(tr->vy()); htrackZ->Fill(tr->vz());
    htrackXY2D->Fill(tr->vx(),tr->vy());
//    cout<<"dt : "<<dt(&(*tr),&(*tr)).P()<<endl;
if(debug && tridx==0) cout<<"66.1"<<endl;
//    if(int(events)%10==0)
    {
//      double d1=length(&(*tr),VtCollection), l1=length_cor(&(*tr),VtCollection,x0,y0);
      double D1=100, L1=100;
      if(vtxon) D1=length(&(*tr),VtCollection);
      if(vtxon) L1=length_cor(&(*tr),VtCollection,x0,y0);
      if(D1<100) htrvtdxy->Fill(D1);
      if(L1<100) htrvtdxy_cor->Fill(L1);
      double dz=100, dz2=100, dz3=100;
      if(vtxon) dz=lengthz(&(*tr),VtCollection), dz2=lengthz(&(*tr),VtCollection2), dz3=lengthz(&(*tr),VtCollection3);
      if(dz<100) htrvtdz->Fill(dz); if(dz2<100) htrvtdz2->Fill(dz2); if(dz3<100) htrvtdz3->Fill(dz3);
    }
if(debug && tridx==0) cout<<"67"<<endl;
    if(vtxon) htrvtd->Fill(rvd(&(*tr),VtCollection).P());
    if(vtxon) htripvtd->Fill(D(&(*tr),VtCollection).P());
//    htrvtxy->Fill(tr->vertex().x(),tr->vertex().y());

    for (TrackCollection::const_iterator tr2=TrCollection->begin(); tr2!=TrCollection->end(); tr2++, tridx2++)
    {
      if(tridx2>tridx && events<=50)
      {
        htrackdz->Fill(tr->dz()-tr2->dz());
        htrackdxy->Fill(tr->dxy()-tr2->dxy());
        htrackdsz->Fill(tr->dsz()-tr2->dsz());
        htrackdX->Fill(tr->vx()-tr2->vx());
        htrackdY->Fill(tr->vy()-tr2->vy());
        htrackdZ->Fill(tr->vz()-tr2->vz());
      }
    }
if(debug && tridx==0) cout<<"68"<<endl;
    if(tr->pt()>10 && fabs(tr->eta())<2.4)
    {
      TLorentzVector ttr(tr->px(),tr->py(),tr->pz(),tr->p());
      htrackpt_ptetacut->Fill(tr->pt()); 
      htracketa_ptetacut->Fill(tr->eta());
      htrackz_ptetacut->Fill(tr->dz());
      trptsorted.push_back(make_pair(tr->pt(), tridx));
      if(tr->charge()==1 ) trptsortedP.push_back(make_pair(tr->pt(), tridx));
      if(tr->charge()==-1) trptsortedM.push_back(make_pair(tr->pt(), tridx));
      tridx_ptetacut++;

//      htrackisolPt_ptetacut01->Fill(isolPt(ttr,TrCollection));
//      htrackisolPt_ptetacut1->Fill(isolPt(ttr,TrCollection));
      htrackisolPt_ptetacut10->Fill(isolPt(ttr,TrCollection));
      htrackisolPt_ptetacut100->Fill(isolPt(ttr,TrCollection));
//      htrackisolPt_ptetacut1000->Fill(isolPt(ttr,TrCollection));
      htrackisolPt2_ptetacut10->Fill(isolPt2(ttr,TrCollection,MuCollection));
      htrackisolPt2_ptetacut100->Fill(isolPt2(ttr,TrCollection,MuCollection));
if(debug && tridx==0) cout<<"69"<<endl;
      if(isolPt(ttr,TrCollection)<isolPtcut)
      {
        htrackpt_ptetaisolcut->Fill(tr->pt());
        htracketa_ptetaisolcut->Fill(tr->eta());
        trptsortedi.push_back(make_pair(tr->pt(), tridx));
        if(tr->charge()==1 ) trptsortedPi.push_back(make_pair(tr->pt(), tridx));
        if(tr->charge()==-1) trptsortedMi.push_back(make_pair(tr->pt(), tridx));
        tridx_ptetaisolcut++;
      }
      if(isolPt2(ttr,TrCollection,MuCollection)<isolPtcut)
      {
        htrackpt_ptetaisol2cut->Fill(tr->pt());
        htracketa_ptetaisol2cut->Fill(tr->eta());
        htrackz_ptetaisol2cut->Fill(tr->dz());
        trptsortedi2.push_back(make_pair(tr->pt(), tridx));
        if(tr->charge()==1 ) trptsortedPi2.push_back(make_pair(tr->pt(), tridx));
        if(tr->charge()==-1) trptsortedMi2.push_back(make_pair(tr->pt(), tridx));
        tridx_ptetaisol2cut++;
if(debug && tridx==0) cout<<"70"<<endl;
        if(jetisol(ttr,corJets))
        {
if(debug && tridx==0) cout<<"80"<<endl;
          htrackpt_ptetaisol2jcut->Fill(tr->pt());
          htracketa_ptetaisol2jcut->Fill(tr->eta());
          trptsortedi2j.push_back(make_pair(tr->pt(), tridx));
          if(tr->charge()==1 ) trptsortedPi2j.push_back(make_pair(tr->pt(), tridx));
          if(tr->charge()==-1) trptsortedMi2j.push_back(make_pair(tr->pt(), tridx));
          tridx_ptetaisol2jcut++;

          htrackz_ptetaisol2jcut->Fill(tr->dz());
          htrackz2_ptetaisol2jcut->Fill( tr->dxy()+x0*sin(tr->phi())-y0*cos(tr->phi()) );
          htrackxy_ptetaisol2jcut->Fill(tr->dxy());
          htracksz_ptetaisol2jcut->Fill(tr->dsz());
          htrackX_ptetaisol2jcut->Fill(tr->vx());
          htrackY_ptetaisol2jcut->Fill(tr->vy());
          htrackZ_ptetaisol2jcut->Fill(tr->vz());
//    htrackX_beamP->Fill(tr->vx()+x0); htrackY_beamP->Fill(tr->vy()+y0); htrackZ_beamP->Fill(tr->vz()+z0);
//    htrackX_beamM->Fill(tr->vx()+x0); htrackY_beamM->Fill(tr->vy()+y0); htrackZ_beamM->Fill(tr->vz()+z0);
          htrackXY2D_ptetaisol2jcut->Fill(tr->vx(),tr->vy());
          htrackXY2D_ptetaisol2jcut_beamP->Fill(tr->vx()+x0,tr->vy()+y0);
          htrackXY2D_ptetaisol2jcut_beamM->Fill(tr->vx()-x0,tr->vy()-y0);

//          bool trmcmu=false;
//          for(int i=0;i<4;i++) if(mcmu[i].Pt()>10 && fabs(mcmu[i].Eta())<2.4 && mcmu[i].DeltaR(ttr)<0.01 && ((i<=1&&tr->charge()==1)||(i>=2&&tr->charge()==-1))) trmcmu=true;
//cout<<matmcmu(&(*mcmu),&(*tr))<<endl;
          bool matching_mcmu=matmcmu(&(*mcmu),&(*tr));
          double d1=100, l1=100;
          if(vtxon) d1=length(&(*tr),VtCollection), l1=length_cor(&(*tr),VtCollection,x0,y0);
          if(d1<100) 
          {
            htrvtdxy_isol2j->Fill(d1); htrvtdxy_isol2j_in->Fill(d1);
            tridx_isol2jv++;
            if(d1>0.02) htrpt_i2junv->Fill(tr->pt());
//            if(trmcmu) {htrvtdxy_isol2jmc->Fill(d1); htrvtdxy_isol2jmc_in->Fill(d1);}
            if(matching_mcmu) {htrvtdxy_isol2jmc->Fill(d1); htrvtdxy_isol2jmc_in->Fill(d1);}
          }
          if(l1<100) {htrvtdxy_isol2j_cor->Fill(l1); tridx_isol2jvc++;}
          if(d1<0.015 && tr->charge()==1 ) trptsortedPi2jv.push_back(make_pair(tr->pt(), tridx));
          if(d1<0.015 && tr->charge()==-1) trptsortedMi2jv.push_back(make_pair(tr->pt(), tridx));
          if(l1<0.015 && tr->charge()==1 ) trptsortedPi2jvc.push_back(make_pair(tr->pt(), tridx));
          if(l1<0.015 && tr->charge()==-1) trptsortedMi2jvc.push_back(make_pair(tr->pt(), tridx));
          double dz=100, dz2=100, dz3=100;
          if(vtxon) dz=lengthz(&(*tr),VtCollection), dz2=lengthz(&(*tr),VtCollection2), dz3=lengthz(&(*tr),VtCollection3);
          if(dz<100) htrvtdz_isol2j->Fill(dz);
          if(dz2<100) htrvtdz2_isol2j->Fill(dz2);
          if(dz3<100) htrvtdz3_isol2j->Fill(dz3);

          TLorentzVector trvd; if(vtxon) trvd=rvd(&(*tr),VtCollection);
          TLorentzVector tD; if(vtxon) tD=D(&(*tr),VtCollection);
          htrvtd_i2j->Fill(trvd.P());
          htrvtdz_i2j->Fill(trvd.Pz());
          htrvtdxy_i2j->Fill(trvd.Pt());
          htripvtd_i2j->Fill(tD.P());
          htripvtdz_i2j->Fill(tD.Pz());
          htripvtdxy_i2j->Fill(tD.Pt());
          if(matching_mcmu)
          {
            htrvtd_i2jmc->Fill(trvd.P());
            htrvtdz_i2jmc->Fill(trvd.Pz());
            htrvtdxy_i2jmc->Fill(trvd.Pt());
            if(dz<100) htrvtdz_isol2jmc->Fill(dz);
            htripvtd_i2jmc->Fill(tD.P());
            htripvtdz_i2jmc->Fill(tD.Pz());
            htripvtdxy_i2jmc->Fill(tD.Pt());
          }
        }
        else
        {
          TLorentzVector tsta2;
          double dR=0.01;
          for (TrackCollection::const_iterator sta=StaCollection->begin(); sta!=StaCollection->end(); sta++)
          {
            TLorentzVector tsta(sta->px(),sta->py(),sta->pz(),sta->p());
            if(sta->pt()>10 && fabs(sta->eta())<2.4 && ttr.DeltaR(tsta)<dR)
            {
              dR=ttr.DeltaR(tsta);
              tsta2.SetPxPyPzE(sta->px(),sta->py(),sta->pz(),sta->p());
            } 
          }
          if(dR<0.01) htrstapt_ptetaisol2notjcut->Fill(tr->pt(),tsta2.Pt());
        }
if(debug && tridx==0) cout<<"90"<<endl;
      }
      double dR_trgmu=0.1, gmupt_trmatch=0;
      for (TrackCollection::const_iterator gmu=MuCollection->begin(); gmu!=MuCollection->end(); gmu++)
      {
        if(gmu->pt()>10 && fabs(gmu->eta())<2.4)
        {
          TLorentzVector tgmu(gmu->px(),gmu->py(),gmu->pz(),gmu->p());
          if(ttr.DeltaR(tgmu)<dR_trgmu)
          {
            dR_trgmu=ttr.DeltaR(tgmu);
            gmupt_trmatch=gmu->pt();
          }
        }
      }
      if(dR_trgmu<0.01) htrgmupt_ptetacut->Fill(tr->pt(),gmupt_trmatch);

      TLorentzVector tsta2;
      double dR=0.01;
      for (TrackCollection::const_iterator sta=StaCollection->begin(); sta!=StaCollection->end(); sta++)
      {
        TLorentzVector tsta(sta->px(),sta->py(),sta->pz(),sta->p());
        if(sta->pt()>10 && fabs(sta->eta())<2.4 && ttr.DeltaR(tsta)<dR)
        {
          dR=ttr.DeltaR(tsta);
          tsta2.SetPxPyPzE(sta->px(),sta->py(),sta->pz(),sta->p());
        } 
      }
      if(dR<0.01) htrstapt_ptetacut->Fill(tr->pt(),tsta2.Pt());
    }
  }
if(debug) cout<<"100"<<endl;
  hntr_ptetacut->Fill(tridx_ptetacut);
  hntr_ptetaisolcut->Fill(tridx_ptetaisolcut);
  hntr_ptetaisol2cut->Fill(tridx_ptetaisol2cut);
if(debug) cout<<"101"<<endl;
  hntr_ptetaisol2jcut->Fill(tridx_ptetaisol2jcut);
  if(tridx_ptetaisol2jcut>=1) hntr_i2j->Fill(TrCollection->size());
  if(tridx_ptetaisol2jcut==3) 
  {
    hntr_i2j3mu->Fill(TrCollection->size());
    for( VertexCollection::const_iterator vt = VtCollection->begin(); vt != VtCollection->end();++vt )
    {
      if(!(vt->x()==x0 && vt->y()==y0 && vt->z()==z0))
      {
        hvtchi2_BS3mu->Fill(vt->normalizedChi2());
      }
    }
  }
  hntr_isol2j_VS_vt->Fill(tridx_ptetaisol2jcut, VtCollection->size());
//  hntr_isol2j_VS_vt_mat->Fill(tridx_ptetaisol2jcut, tridx_isol2jv);
//  hntr_isol2j_VS_vt_matc->Fill(tridx_ptetaisol2jcut, tridx_isol2jvc);
if(debug) cout<<"102"<<endl;
  sort(trptsorted.begin(), trptsorted.end(), greater<ptvsindex>());
  sort(trptsortedP.begin(), trptsortedP.end(), greater<ptvsindex>());
  sort(trptsortedM.begin(), trptsortedM.end(), greater<ptvsindex>());
  sort(trptsortedi.begin(), trptsortedi.end(), greater<ptvsindex>());
  sort(trptsortedPi.begin(), trptsortedPi.end(), greater<ptvsindex>());
  sort(trptsortedMi.begin(), trptsortedMi.end(), greater<ptvsindex>());
  sort(trptsortedi2.begin(), trptsortedi2.end(), greater<ptvsindex>());
  sort(trptsortedPi2.begin(), trptsortedPi2.end(), greater<ptvsindex>());
  sort(trptsortedMi2.begin(), trptsortedMi2.end(), greater<ptvsindex>());
  sort(trptsortedi2j.begin(), trptsortedi2j.end(), greater<ptvsindex>());
  sort(trptsortedPi2j.begin(), trptsortedPi2j.end(), greater<ptvsindex>());
  sort(trptsortedMi2j.begin(), trptsortedMi2j.end(), greater<ptvsindex>());

  TLorentzVector ttrP[100], ttrM[100], ttrPi[100], ttrMi[100], ttr4[4], ttr4i[4];
  TLorentzVector ttrPi2[100], ttrMi2[100];
  for (unsigned int i=0; i<trptsortedP.size(); i++)
  {
    const Track *trP = & (*TrCollection)[trptsortedP[i].second];
    ttrP[i].SetPxPyPzE(trP->px(), trP->py(), trP->pz(), trP->p());
    if(i<2) ttr4[i].SetPxPyPzE(trP->px(), trP->py(), trP->pz(), trP->p());
  }
  for (unsigned int i=0; i<trptsortedM.size(); i++)
  {
    const Track *trM = & (*TrCollection)[trptsortedM[i].second];
    ttrM[i].SetPxPyPzE(trM->px(), trM->py(), trM->pz(), trM->p());
    if(i<2) ttr4[i+2].SetPxPyPzE(trM->px(), trM->py(), trM->pz(), trM->p());
  }
  if(trptsortedP.size()>=2) hHPPmass_tr->Fill((ttrP[0]+ttrP[1]).M());
  if(trptsortedM.size()>=2) hHMMmass_tr->Fill((ttrM[0]+ttrM[1]).M());
  if(trptsortedP.size()>=2 && trptsortedM.size()>=2)
  {
    hHPPHMMmass_tr->Fill((ttrP[0]+ttrP[1]+ttrM[0]+ttrM[1]).M());
    h2trdeltaphi->Fill(fabs(deltaphi((ttrP[0]+ttrP[1]).Phi(),(ttrM[0]+ttrM[1]).Phi())));
  }

  for (unsigned int i=0; i<trptsortedPi.size(); i++)
  {
    const Track *trPi = & (*TrCollection)[trptsortedPi[i].second];
    ttrPi[i].SetPxPyPzE(trPi->px(), trPi->py(), trPi->pz(), trPi->p());
    if(i<2) ttr4i[i].SetPxPyPzE(trPi->px(), trPi->py(), trPi->pz(), trPi->p());
  }
  for (unsigned int i=0; i<trptsortedMi.size(); i++)
  {
    const Track *trMi = & (*TrCollection)[trptsortedMi[i].second];
    ttrMi[i].SetPxPyPzE(trMi->px(), trMi->py(), trMi->pz(), trMi->p());
    if(i<2) ttr4i[i+2].SetPxPyPzE(trMi->px(), trMi->py(), trMi->pz(), trMi->p());
  }
  if(trptsortedPi.size()>=2) hHPPmass_trisol->Fill((ttrPi[0]+ttrPi[1]).M());
  if(trptsortedMi.size()>=2) hHMMmass_trisol->Fill((ttrMi[0]+ttrMi[1]).M());
  if(trptsortedPi.size()>=2 && trptsortedMi.size()>=2)
  {
    hHPPHMMmass_trisol->Fill((ttrPi[0]+ttrPi[1]+ttrMi[0]+ttrMi[1]).M());
    h2trdeltaphi_isol->Fill(fabs(deltaphi((ttrPi[0]+ttrPi[1]).Phi(),(ttrMi[0]+ttrMi[1]).Phi())));
  }

  for (unsigned int i=0; i<trptsortedPi2.size(); i++)
  {
    const Track *trPi2 = & (*TrCollection)[trptsortedPi2[i].second];
    ttrPi2[i].SetPxPyPzE(trPi2->px(), trPi2->py(), trPi2->pz(), trPi2->p());
//    if(i<2) ttr4i[i].SetPxPyPzE(trPi->px(), trPi->py(), trPi->pz(), trPi->p());
  }
  for (unsigned int i=0; i<trptsortedMi2.size(); i++)
  {
    const Track *trMi2 = & (*TrCollection)[trptsortedMi2[i].second];
    ttrMi2[i].SetPxPyPzE(trMi2->px(), trMi2->py(), trMi2->pz(), trMi2->p());
//    if(i<2) ttr4i[i+2].SetPxPyPzE(trMi->px(), trMi->py(), trMi->pz(), trMi->p());
  }
  if(trptsortedPi2.size()>=2) hHPPmass_trisol2->Fill((ttrPi2[0]+ttrPi2[1]).M());
  if(trptsortedMi2.size()>=2) hHMMmass_trisol2->Fill((ttrMi2[0]+ttrMi2[1]).M());
  if(trptsortedPi2.size()>=2 && trptsortedMi2.size()>=2)
  {
    hHPPHMMmass_trisol2->Fill((ttrPi2[0]+ttrPi2[1]+ttrMi2[0]+ttrMi2[1]).M());
    h2trdeltaphi_isol2->Fill(fabs(deltaphi((ttrPi2[0]+ttrPi2[1]).Phi(),(ttrMi2[0]+ttrMi2[1]).Phi())));
    if(fabs(deltaphi((ttrPi2[0]+ttrPi2[1]).Phi(),(ttrMi2[0]+ttrMi2[1]).Phi()))>dphicut_4mu);
    {
      hHPPmass_trisol2dphi->Fill((ttrPi2[0]+ttrPi2[1]).M());
      hHMMmass_trisol2dphi->Fill((ttrMi2[0]+ttrMi2[1]).M());
      hHPPHMMmass_trisol2dphi->Fill((ttrPi2[0]+ttrPi2[1]+ttrMi2[0]+ttrMi2[1]).M());
    }
  }

if(debug) cout<<"120"<<endl;
  TLorentzVector ttrPi2j[100], ttrMi2j[100];
  double trPi2jz[2]={0,0}, trMi2jz[2]={0,0}, trPi2jZ[2]={0,0}, trMi2jZ[2]={0,0};
//  double trPi2jv[4]={0,}, trMi2jv[4]={0,}, trPi2jvc[4]={0,}, trMi2jvc[4]={0,};
//  bool tracki2jv=true, tracki2jvc=true, trackPi2jv=true, trackPi2jvc=true, trackMi2jv=true, trackMi2jvc=true;
  const Track *tPi2j[2], *tMi2j[2];
  for (unsigned int i=0; i<trptsortedPi2j.size(); i++)
  {
    const Track *trPi2j = & (*TrCollection)[trptsortedPi2j[i].second];
    ttrPi2j[i].SetPxPyPzE(trPi2j->px(), trPi2j->py(), trPi2j->pz(), trPi2j->p());
    double d1=100;
    if(trptsortedPi2j.size()>=2 && trptsortedMi2j.size()>=2)
    {
//      if(matmcmu(&(*mcmu),&(*trPi2j))) d1=length(&(*trPi2j),VtCollection);
      if(vtxon) d1=length(&(*trPi2j),VtCollection);
      if(d1<100)
      {
        htrvtdxy_isol2j4mu->Fill(d1); htrvtdxy_isol2j4mu_in->Fill(d1);
        if(matmcmu(&(*mcmu),&(*trPi2j))) {htrvtdxy_isol2j4mumc->Fill(d1); htrvtdxy_isol2j4mumc_in->Fill(d1);}
      }
    }
    if(i<2)
    {
      trPi2jz[i]=trPi2j->dz(); trPi2jZ[i]=trPi2j->vz();
      tPi2j[i] = & (*TrCollection)[trptsortedPi2j[i].second];
    }
  }
  for (unsigned int i=0; i<trptsortedMi2j.size(); i++)
  {
    const Track *trMi2j = & (*TrCollection)[trptsortedMi2j[i].second];
    ttrMi2j[i].SetPxPyPzE(trMi2j->px(), trMi2j->py(), trMi2j->pz(), trMi2j->p());
    double d1=100;
    if(trptsortedPi2j.size()>=2 && trptsortedMi2j.size()>=2)
    {
//      if(matmcmu(&(*mcmu),&(*trMi2j))) d1=length(&(*trMi2j),VtCollection);
      if(vtxon) d1=length(&(*trMi2j),VtCollection);
      if(d1<100)
      {
        htrvtdxy_isol2j4mu->Fill(d1); htrvtdxy_isol2j4mu_in->Fill(d1);
        if(matmcmu(&(*mcmu),&(*trMi2j))) {htrvtdxy_isol2j4mumc->Fill(d1); htrvtdxy_isol2j4mumc_in->Fill(d1);}
      }
    }
    if(i<2)
    {
      trMi2jz[i]=trMi2j->dz(); trMi2jZ[i]=trMi2j->vz();
      tMi2j[i] = & (*TrCollection)[trptsortedMi2j[i].second];
    }
  }
  if(trptsortedPi2j.size()>=2)
  {
    hHPPmass_trisol2j->Fill((ttrPi2j[0]+ttrPi2j[1]).M());
    if(fabs(trPi2jz[0]-trPi2jz[1])<dzcut_)
    {
      hHPPmass_trisol2jdz->Fill((ttrPi2j[0]+ttrPi2j[1]).M());
      hHPPmass_trisol2jdz_ex->Fill((ttrPi2j[0]+ttrPi2j[1]).M());
    }
//    h2trdz_isol2j->Fill(fabs(trPi2jz[0]-trPi2jz[1]));
    h2trdz_isol2j->Fill(trPi2jz[0]-trPi2jz[1]); h2trdz_isol2j_ex->Fill(trPi2jz[0]-trPi2jz[1]);
//    h2trdZ_isol2j->Fill(trPi2jZ[0]-trPi2jZ[1]); h2trdZ_isol2j_ex->Fill(trPi2jZ[0]-trPi2jZ[1]);
  }
  if(trptsortedMi2j.size()>=2)
  {
    hHMMmass_trisol2j->Fill((ttrMi2j[0]+ttrMi2j[1]).M());
    if(fabs(trMi2jz[0]-trMi2jz[1])<dzcut_)
    {
      hHMMmass_trisol2jdz->Fill((ttrMi2j[0]+ttrMi2j[1]).M());
      hHMMmass_trisol2jdz_ex->Fill((ttrMi2j[0]+ttrMi2j[1]).M());
    }
//    h2trdz_isol2j->Fill(fabs(trMi2jz[0]-trMi2jz[1]));
    h2trdz_isol2j->Fill(trMi2jz[0]-trMi2jz[1]); h2trdz_isol2j_ex->Fill(trMi2jz[0]-trMi2jz[1]);
//    h2trdZ_isol2j->Fill(trMi2jz[0]-trMi2jZ[1]); h2trdz_isol2j_ex->Fill(trMi2jZ[0]-trMi2jZ[1]);
  }
  if(trptsortedPi2j.size()>=2 && trptsortedMi2j.size()>=2)
  {
    hHPPHMMmass_trisol2j->Fill((ttrPi2j[0]+ttrPi2j[1]+ttrMi2j[0]+ttrMi2j[1]).M());
//----------ZZ4mu----------
    double PMmass[4]={0,}, Z=91.1876, PMdphi[4]={0,};
    int k=0;
    for(int i=0;i<2;i++) for(int j=0;j<2;j++)
    {
      PMmass[k]=(ttrPi2j[i]+ttrMi2j[j]).M();
      hPMmassall_trisol2j->Fill(PMmass[k]);
      PMdphi[k]=fabs(deltaphi(ttrPi2j[i].Phi(),ttrMi2j[j].Phi()));
      hPMdphi_trisol2j->Fill(PMdphi[k]);
      k++;
    }
    if( (PMmass[0]+PMmass[3])>(PMmass[1]+PMmass[2]) ) {hPMmassbig_trisol2j->Fill(PMmass[0]); hPMmassbig_trisol2j->Fill(PMmass[3]);}
    else {hPMmassbig_trisol2j->Fill(PMmass[1]); hPMmassbig_trisol2j->Fill(PMmass[2]);}
    if( (fabs(Z-PMmass[0])+fabs(Z-PMmass[3]))<(fabs(Z-PMmass[1])+fabs(Z-PMmass[2])) ) {hPMmassZ_trisol2j->Fill(PMmass[0]); hPMmassZ_trisol2j->Fill(PMmass[3]);}
    else {hPMmassZ_trisol2j->Fill(PMmass[1]); hPMmassZ_trisol2j->Fill(PMmass[2]);}
    if( (PMdphi[0]+PMdphi[3])>(PMdphi[1]+PMdphi[2]) ) {hPMmassdphi_trisol2j->Fill(PMmass[0]); hPMmassdphi_trisol2j->Fill(PMmass[3]);}
    else {hPMmassdphi_trisol2j->Fill(PMmass[1]); hPMmassdphi_trisol2j->Fill(PMmass[2]);}

    htStSd_i2j4mu->Fill(dt(&(*tPi2j[0]),&(*tPi2j[1])).P());
    htStSdz_i2j4mu->Fill(dt(&(*tPi2j[0]),&(*tPi2j[1])).Pz());
    htStSdxy_i2j4mu->Fill(dt(&(*tPi2j[0]),&(*tPi2j[1])).Pt());
    htStSd_i2j4mu->Fill(dt(&(*tMi2j[0]),&(*tMi2j[1])).P());
    htStSdz_i2j4mu->Fill(dt(&(*tMi2j[0]),&(*tMi2j[1])).Pz());
    htStSdxy_i2j4mu->Fill(dt(&(*tMi2j[0]),&(*tMi2j[1])).Pt());

    htPtMd_i2j4mu->Fill(dt(&(*tPi2j[0]),&(*tMi2j[1])).P());
    htPtMdz_i2j4mu->Fill(dt(&(*tPi2j[0]),&(*tMi2j[1])).Pz());
    htPtMdxy_i2j4mu->Fill(dt(&(*tPi2j[0]),&(*tMi2j[1])).Pt());
    htPtMd_i2j4mu->Fill(dt(&(*tMi2j[0]),&(*tPi2j[1])).P());
    htPtMdz_i2j4mu->Fill(dt(&(*tMi2j[0]),&(*tPi2j[1])).Pz());
    htPtMdxy_i2j4mu->Fill(dt(&(*tMi2j[0]),&(*tPi2j[1])).Pt());

    TLorentzVector tcp1P=cp1(&(*tPi2j[0]),&(*tPi2j[1]));
    TLorentzVector tcp2P=cp2(&(*tPi2j[0]),&(*tPi2j[1]));
    TLorentzVector tcp1M=cp1(&(*tMi2j[0]),&(*tMi2j[1]));
    TLorentzVector tcp2M=cp2(&(*tMi2j[0]),&(*tMi2j[1]));

    hmcmuxcpx->Fill(mcmu1->vx(),tcp1P.Px()); hmcmuxcpx->Fill(mcmu1->vx(),tcp2P.Px()); hmcmuxcpx->Fill(mcmu1->vx(),tcp1M.Px()); hmcmuxcpx->Fill(mcmu1->vx(),tcp2M.Px());
    hmcmuycpy->Fill(mcmu1->vy(),tcp1P.Py()); hmcmuycpy->Fill(mcmu1->vy(),tcp2P.Py()); hmcmuycpy->Fill(mcmu1->vy(),tcp1M.Py()); hmcmuycpy->Fill(mcmu1->vy(),tcp2M.Py());
    hmcmuzcpz->Fill(mcmu1->vz(),tcp1P.Pz()); hmcmuzcpz->Fill(mcmu1->vz(),tcp2P.Pz()); hmcmuzcpz->Fill(mcmu1->vz(),tcp1M.Pz()); hmcmuzcpz->Fill(mcmu1->vz(),tcp2M.Pz());
    for( VertexCollection::const_iterator vt = VtCollection->begin(); vt != VtCollection->end();++vt )
    {
      hvtx->Fill(vt->x()); hvty->Fill(vt->y()); hvtz->Fill(vt->z()); hvtxy->Fill(vt->x(),vt->y()); //hmcmurecovtxy->Fill(vt->x()-x0,vt->y()-y0);
      if(!(vt->x()==x0 && vt->y()==y0 && vt->z()==z0))
      {
        hvtxcpx_BS->Fill(vt->x(),tcp1P.Px()); hvtxcpx_BS->Fill(vt->x(),tcp2P.Px()); hvtxcpx_BS->Fill(vt->x(),tcp1M.Px()); hvtxcpx_BS->Fill(vt->x(),tcp2M.Px());
        hvtycpy_BS->Fill(vt->y(),tcp1P.Py()); hvtycpy_BS->Fill(vt->y(),tcp2P.Py()); hvtycpy_BS->Fill(vt->y(),tcp1M.Py()); hvtycpy_BS->Fill(vt->y(),tcp2M.Py());
        hvtzcpz_BS->Fill(vt->z(),tcp1P.Pz()); hvtzcpz_BS->Fill(vt->z(),tcp2P.Pz()); hvtzcpz_BS->Fill(vt->z(),tcp1M.Pz()); hvtzcpz_BS->Fill(vt->z(),tcp2M.Pz());
        hvtchi2_BS4mu->Fill(vt->normalizedChi2());
      }
    }
    hntr_i2j4mu->Fill(TrCollection->size());
  }
//------------------------------

  for (unsigned int i=0; i<trptsortedi2j.size(); i++)
  {
    const Track *tri2j = & (*TrCollection)[trptsortedi2j[i].second];
    TLorentzVector ttri2j; if(vtxon) ttri2j=ip(&(*tri2j),VtCollection);
    for (unsigned int n=0; n<trptsortedi2j.size(); n++)
    {
      const Track *tri2j2 = & (*TrCollection)[trptsortedi2j[n].second];
      TLorentzVector ttri2j2; if(vtxon) ttri2j2=ip(&(*tri2j2),VtCollection);
//      TLorentzVector tip=ip(&(*tri2j),VtCollection)-ip(&(*tri2j2),VtCollection);
      TLorentzVector trrd=rrd(&(*tri2j),&(*tri2j2));
      TLorentzVector tip; if(vtxon) tip=ttri2j-ttri2j2;
      TLorentzVector tdt=dt(&(*tri2j),&(*tri2j2));
//      if(i<4 && n<4 && n>i) htrackdz_ptetaisol2jcut->Fill(fabs(tri2j->dz()-tri2j2->dz()));
      if(i<4 && n<4 && n>i)
      {
        htrackdz_ptetaisol2jcut->Fill(tri2j->dz()-tri2j2->dz());
        htrackdxy_ptetaisol2jcut->Fill(tri2j->dxy()-tri2j2->dxy());
        htrackdsz_ptetaisol2jcut->Fill(tri2j->dsz()-tri2j2->dsz());
        htrackdX_ptetaisol2jcut->Fill(tri2j->vx()-tri2j2->vx());
        htrackdY_ptetaisol2jcut->Fill(tri2j->vy()-tri2j2->vy());
        htrackdZ_ptetaisol2jcut->Fill(tri2j->vz()-tri2j2->vz());
        htrackdZ_ptetaisol2jcut_ex->Fill(tri2j->vz()-tri2j2->vz());

//        htrtripd_i2j->Fill((ttri2j-ttri2j2).P());
//        htrtripdz_i2j->Fill((ttri2j-ttri2j2).Pz());
        htrtrd_i2j->Fill(trrd.P());
        htrtrdz_i2j->Fill(trrd.Pz());
        htrtrdxy_i2j->Fill(trrd.Pt());
        htrtripd_i2j->Fill(tip.P());
        htrtripdz_i2j->Fill(tip.Pz());
        htrtripdxy_i2j->Fill(tip.Pt());
        htrtrsd_i2j->Fill(tdt.P());
        htrtrsdz_i2j->Fill(tdt.Pz());
        htrtrsdxy_i2j->Fill(tdt.Pt());
        if(matmcmu(&(*mcmu),&(*tri2j)) && matmcmu(&(*mcmu),&(*tri2j2)))
        {
//          htrtripd_i2jmc->Fill((ttri2j-ttri2j2).P());
//          htrtripdz_i2jmc->Fill((ttri2j-ttri2j2).Pz());
//          htrtripdxy_i2jmc->Fill((ttri2j-ttri2j2).Pt());
          htrtrd_i2jmc->Fill(trrd.P());
          htrtrdz_i2jmc->Fill(trrd.Pz());
          htrtrdxy_i2jmc->Fill(trrd.Pt());
          htrtripd_i2jmc->Fill(tip.P());
          htrtripdz_i2jmc->Fill(tip.Pz());
          htrtripdxy_i2jmc->Fill(tip.Pt());
          htrtrsd_i2jmc->Fill(tdt.P());
          htrtrsdz_i2jmc->Fill(tdt.Pz());
          htrtrsdxy_i2jmc->Fill(tdt.Pt());
        }
      }
    }
  }
if(debug) cout<<"130"<<endl;
  TLorentzVector ttrPi2jv[100], ttrMi2jv[100];
  bool tracki2jv=true, trackPi2jv=true, trackMi2jv=true;
  for (unsigned int i=0; i<trptsortedPi2jv.size(); i++)
  {
    const Track *trPi2jv = & (*TrCollection)[trptsortedPi2jv[i].second];
    ttrPi2jv[i].SetPxPyPzE(trPi2jv->px(), trPi2jv->py(), trPi2jv->pz(), trPi2jv->p());
    double d1=100;
    if(vtxon) d1=length(&(*trPi2jv),VtCollection);
    if(d1>0.015) {tracki2jv=false; trackPi2jv=false;}
  }
  for (unsigned int i=0; i<trptsortedMi2jv.size(); i++)
  {
    const Track *trMi2jv = & (*TrCollection)[trptsortedMi2jv[i].second];
    ttrMi2jv[i].SetPxPyPzE(trMi2jv->px(), trMi2jv->py(), trMi2jv->pz(), trMi2jv->p());
    double d1=100;
    if(vtxon) d1=length(&(*trMi2jv),VtCollection);
    if(d1>0.015) {tracki2jv=false; trackMi2jv=false;}
  }
  if(trptsortedPi2jv.size()>=2)
  {
    hHPPmass_trisol2jv->Fill((ttrPi2jv[0]+ttrPi2jv[1]).M());
    if(trackPi2jv) hHPPmass_trisol2jvs->Fill((ttrPi2jv[0]+ttrPi2jv[1]).M());
  }
  if(trptsortedMi2jv.size()>=2)
  {
    hHMMmass_trisol2jv->Fill((ttrMi2jv[0]+ttrMi2jv[1]).M());
    if(trackMi2jv) hHMMmass_trisol2jvs->Fill((ttrMi2jv[0]+ttrMi2jv[1]).M());
  }
  if(trptsortedPi2jv.size()>=2 && trptsortedMi2jv.size()>=2)
  {
    hHPPHMMmass_trisol2jv->Fill((ttrPi2jv[0]+ttrPi2jv[1]+ttrMi2jv[0]+ttrMi2jv[1]).M());
    if(tracki2jv) hHPPHMMmass_trisol2jvs->Fill((ttrPi2jv[0]+ttrPi2jv[1]+ttrMi2jv[0]+ttrMi2jv[1]).M());
  }
if(debug) cout<<"131"<<endl;
  TLorentzVector ttrPi2jvc[100], ttrMi2jvc[100];
  bool tracki2jvc=true, trackPi2jvc=true, trackMi2jvc=true;
  for (unsigned int i=0; i<trptsortedPi2jvc.size(); i++)
  {
    const Track *trPi2jvc = & (*TrCollection)[trptsortedPi2jvc[i].second];
    ttrPi2jvc[i].SetPxPyPzE(trPi2jvc->px(), trPi2jvc->py(), trPi2jvc->pz(), trPi2jvc->p());
    double l1=100;
    if(vtxon) l1=length_cor(&(*trPi2jvc),VtCollection,x0,y0);
    if(l1>0.015) {tracki2jvc=false; trackPi2jvc=false;}
  }
  for (unsigned int i=0; i<trptsortedMi2jvc.size(); i++)
  {
    const Track *trMi2jvc = & (*TrCollection)[trptsortedMi2jvc[i].second];
    ttrMi2jvc[i].SetPxPyPzE(trMi2jvc->px(), trMi2jvc->py(), trMi2jvc->pz(), trMi2jvc->p());
    double l1=100;
    if(vtxon)  l1=length_cor(&(*trMi2jvc),VtCollection,x0,y0);
    if(l1>0.015) {tracki2jvc=false; trackMi2jvc=false;}
  }
  if(trptsortedPi2jvc.size()>=2)
  {
    hHPPmass_trisol2jvc->Fill((ttrPi2jvc[0]+ttrPi2jvc[1]).M());
    if(trackPi2jvc) hHPPmass_trisol2jvcs->Fill((ttrPi2jvc[0]+ttrPi2jvc[1]).M());
  }
  if(trptsortedMi2jvc.size()>=2)
  {
    hHMMmass_trisol2jvc->Fill((ttrMi2jvc[0]+ttrMi2jvc[1]).M());
    if(trackMi2jvc) hHMMmass_trisol2jvcs->Fill((ttrMi2jvc[0]+ttrMi2jvc[1]).M());
  }
  if(trptsortedPi2jvc.size()>=2 && trptsortedMi2jvc.size()>=2)
  {
    hHPPHMMmass_trisol2jvc->Fill((ttrPi2jvc[0]+ttrPi2jvc[1]+ttrMi2jvc[0]+ttrMi2jvc[1]).M());
    if(tracki2jvc) hHPPHMMmass_trisol2jvcs->Fill((ttrPi2jvc[0]+ttrPi2jvc[1]+ttrMi2jvc[0]+ttrMi2jvc[1]).M());
  }
if(debug) cout<<"132"<<endl;

//---------------------------global muons(gmu000)-----------------------------//
  std::vector<ptvsindex> muonsptsorted, muonsptsortedP, muonsptsortedM;
  std::vector<ptvsindex> muonsptsortedi, muonsptsortedPi, muonsptsortedMi;
  std::vector<ptvsindex> muonsptsortedi2, muonsptsortedPi2, muonsptsortedMi2;
  std::vector<ptvsindex> muonsptsortedi2j, muonsptsortedPi2j, muonsptsortedMi2j;
  int muidx=0, muidx2=0, muidx_ptetacut=0, muidx_ptetaisolcut=0, muidx_ptetaisol2cut=0, muidx_ptetaisol2jcut=0;
  for (TrackCollection::const_iterator gmu=MuCollection->begin(); gmu!=MuCollection->end(); gmu++, muidx++)
  {
    hgmupt->Fill(gmu->pt()); 
    hgmueta->Fill(gmu->eta());
    hgmuz->Fill(gmu->dz());
    for (TrackCollection::const_iterator gmu2=MuCollection->begin(); gmu2!=MuCollection->end(); gmu2++, muidx2++)
    {
//      if(muidx2>muidx) hgmudz->Fill(fabs(gmu->dz()-gmu2->dz()));
      if(muidx2>muidx) hgmudz->Fill(gmu->dz()-gmu2->dz());
    }
    if(gmu->pt()>10 && fabs(gmu->eta())<2.4)
    {
      TLorentzVector tgmu(gmu->px(),gmu->py(),gmu->pz(),gmu->p());
      hgmupt_ptetacut->Fill(gmu->pt()); 
      hgmueta_ptetacut->Fill(gmu->eta());
      hgmuz_ptetacut->Fill(gmu->dz());
      muonsptsorted.push_back(make_pair(gmu->pt(), muidx));
      if(gmu->charge()==1 ) muonsptsortedP.push_back(make_pair(gmu->pt(), muidx));
      if(gmu->charge()==-1) muonsptsortedM.push_back(make_pair(gmu->pt(), muidx));
      muidx_ptetacut++;

//      hgmuisolPt_ptetacut01->Fill(isolPt_gmu(tgmu,TrCollection));
//      hgmuisolPt_ptetacut1->Fill(isolPt_gmu(tgmu,TrCollection));
      hgmuisolPt_ptetacut10->Fill(isolPt_gmu(tgmu,TrCollection));
      hgmuisolPt_ptetacut100->Fill(isolPt_gmu(tgmu,TrCollection));
//      hgmuisolPt_ptetacut1000->Fill(isolPt_gmu(tgmu,TrCollection));
      hgmuisolPt2_ptetacut10->Fill(isolPt2_gmu(tgmu,TrCollection,MuCollection));
      hgmuisolPt2_ptetacut100->Fill(isolPt2_gmu(tgmu,TrCollection,MuCollection));
      if(isolPt_gmu(tgmu,TrCollection)<isolPtcut)
      {
        hgmupt_ptetaisolcut->Fill(gmu->pt());
        hgmueta_ptetaisolcut->Fill(gmu->eta());
        muonsptsortedi.push_back(make_pair(gmu->pt(), muidx));
        if(gmu->charge()==1 ) muonsptsortedPi.push_back(make_pair(gmu->pt(), muidx));
        if(gmu->charge()==-1) muonsptsortedMi.push_back(make_pair(gmu->pt(), muidx));
        muidx_ptetaisolcut++;
      }
      if(isolPt2_gmu(tgmu,TrCollection,MuCollection)<isolPtcut)
      {
        hgmupt_ptetaisol2cut->Fill(gmu->pt());
        hgmueta_ptetaisol2cut->Fill(gmu->eta());
        hgmuz_ptetaisol2cut->Fill(gmu->dz());
        muonsptsortedi2.push_back(make_pair(gmu->pt(), muidx));
        if(gmu->charge()==1 ) muonsptsortedPi2.push_back(make_pair(gmu->pt(), muidx));
        if(gmu->charge()==-1) muonsptsortedMi2.push_back(make_pair(gmu->pt(), muidx));
        muidx_ptetaisol2cut++;
if(debug && muidx==0) cout<<"140"<<endl;
        if(jetisol(tgmu,corJets))
        {
          hgmupt_ptetaisol2jcut->Fill(gmu->pt());
          hgmueta_ptetaisol2jcut->Fill(gmu->eta());
          muonsptsortedi2j.push_back(make_pair(gmu->pt(), muidx));
          if(gmu->charge()==1 ) muonsptsortedPi2j.push_back(make_pair(gmu->pt(), muidx));
          if(gmu->charge()==-1) muonsptsortedMi2j.push_back(make_pair(gmu->pt(), muidx));
          muidx_ptetaisol2jcut++;
          hgmuz_ptetaisol2jcut->Fill(gmu->dz());
if(debug && muidx==0) cout<<"150"<<endl;
        }
        else
        {
          TLorentzVector tsta2;
          double dR=0.01;
          for (TrackCollection::const_iterator sta=StaCollection->begin(); sta!=StaCollection->end(); sta++)
          {
            TLorentzVector tsta(sta->px(),sta->py(),sta->pz(),sta->p());
            if(sta->pt()>10 && fabs(sta->eta())<2.4 && tgmu.DeltaR(tsta)<dR)
            {
              dR=tgmu.DeltaR(tsta);
              tsta2.SetPxPyPzE(sta->px(),sta->py(),sta->pz(),sta->p());
            } 
          }
          if(dR<0.01) hgmustapt_ptetaisol2notjcut->Fill(gmu->pt(),tsta2.Pt());
        }
      }
//      int ntrackmatch=0;
//      for (TrackCollection::const_iterator tr=TrCollection->begin(); tr!=TrCollection->end(); tr++)
//      {
//        TLorentzVector ttr(tr->px(),tr->py(),tr->pz(),tr->p());
//        if(ttr.DeltaR(tgmu)<0.01 && tr->pt()>10 && fabs(tr->eta())<2.4) ntrackmatch++;
//      }
//      if(ntrackmatch>=2) cout<<"# of matched track with a global muon : "<<ntrackmatch<<"(pt>10, |eta|<2.4)"<<endl;
      TLorentzVector tsta2;
      double dR=0.01;
      for (TrackCollection::const_iterator sta=StaCollection->begin(); sta!=StaCollection->end(); sta++)
      {
        TLorentzVector tsta(sta->px(),sta->py(),sta->pz(),sta->p());
        if(sta->pt()>10 && fabs(sta->eta())<2.4 && tgmu.DeltaR(tsta)<dR)
        {
          dR=tgmu.DeltaR(tsta);
          tsta2.SetPxPyPzE(sta->px(),sta->py(),sta->pz(),sta->p());
        } 
      }
      if(dR<0.01) hgmustapt_ptetacut->Fill(gmu->pt(),tsta2.Pt());
    }
  }
if(debug) cout<<"160"<<endl;
  hngmu_ptetacut->Fill(muidx_ptetacut);
  hngmu_ptetaisolcut->Fill(muidx_ptetaisolcut);
  hngmu_ptetaisol2cut->Fill(muidx_ptetaisol2cut);
  hngmu_ptetaisol2jcut->Fill(muidx_ptetaisol2jcut);
  sort(muonsptsorted.begin(), muonsptsorted.end(), greater<ptvsindex>());
  sort(muonsptsortedP.begin(), muonsptsortedP.end(), greater<ptvsindex>());
  sort(muonsptsortedM.begin(), muonsptsortedM.end(), greater<ptvsindex>());
  sort(muonsptsortedi.begin(), muonsptsortedi.end(), greater<ptvsindex>());
  sort(muonsptsortedPi.begin(), muonsptsortedPi.end(), greater<ptvsindex>());
  sort(muonsptsortedMi.begin(), muonsptsortedMi.end(), greater<ptvsindex>());
  sort(muonsptsortedi2.begin(), muonsptsortedi2.end(), greater<ptvsindex>());
  sort(muonsptsortedPi2.begin(), muonsptsortedPi2.end(), greater<ptvsindex>());
  sort(muonsptsortedMi2.begin(), muonsptsortedMi2.end(), greater<ptvsindex>());
  sort(muonsptsortedi2j.begin(), muonsptsortedi2j.end(), greater<ptvsindex>());
  sort(muonsptsortedPi2j.begin(), muonsptsortedPi2j.end(), greater<ptvsindex>());
  sort(muonsptsortedMi2j.begin(), muonsptsortedMi2j.end(), greater<ptvsindex>());
if(debug) cout<<"170"<<endl;

//  int nmuP=muonsptsortedP.size(), nmuM=muonsptsortedM.size();
//  TLorentzVector tgmuP[nmuP], tgmuM[nmuM], tgmuPi[nmuP], tgmuMi[nmuM];
  TLorentzVector tgmuP[100], tgmuM[100], tgmuPi[100], tgmuMi[100], tgmu4[4], tgmu4i[4];
  TLorentzVector tgmuPi2[100], tgmuMi2[100], tgmuPi2j[100], tgmuMi2j[100];
  for (unsigned int i=0; i<muonsptsortedP.size(); i++)
  {
    const Track *gmuP = & (*MuCollection)[muonsptsortedP[i].second];
    tgmuP[i].SetPxPyPzE(gmuP->px(), gmuP->py(), gmuP->pz(), gmuP->p());
    if(i<2) tgmu4[i].SetPxPyPzE(gmuP->px(), gmuP->py(), gmuP->pz(), gmuP->p());
  }
  for (unsigned int i=0; i<muonsptsortedM.size(); i++)
  {
    const Track *gmuM = & (*MuCollection)[muonsptsortedM[i].second];
    tgmuM[i].SetPxPyPzE(gmuM->px(), gmuM->py(), gmuM->pz(), gmuM->p());
    if(i<2) tgmu4[i+2].SetPxPyPzE(gmuM->px(), gmuM->py(), gmuM->pz(), gmuM->p());
  }
  if(muonsptsortedP.size()>=2) hHPPmass->Fill((tgmuP[0]+tgmuP[1]).M());
  if(muonsptsortedM.size()>=2) hHMMmass->Fill((tgmuM[0]+tgmuM[1]).M());
  if(muonsptsortedP.size()>=2 && muonsptsortedM.size()>=2)
  {
    hHPPHMMmass->Fill((tgmuP[0]+tgmuP[1]+tgmuM[0]+tgmuM[1]).M());
    h2gmudeltaphi->Fill(fabs(deltaphi((tgmuP[0]+tgmuP[1]).Phi(),(tgmuM[0]+tgmuM[1]).Phi())));
  }

  for (unsigned int i=0; i<muonsptsortedPi.size(); i++)
  {
    const Track *gmuPi = & (*MuCollection)[muonsptsortedPi[i].second];
    tgmuPi[i].SetPxPyPzE(gmuPi->px(), gmuPi->py(), gmuPi->pz(), gmuPi->p());
    if(i<2) tgmu4i[i].SetPxPyPzE(gmuPi->px(), gmuPi->py(), gmuPi->pz(), gmuPi->p());
  }
  for (unsigned int i=0; i<muonsptsortedMi.size(); i++)
  {
    const Track *gmuMi = & (*MuCollection)[muonsptsortedMi[i].second];
    tgmuMi[i].SetPxPyPzE(gmuMi->px(), gmuMi->py(), gmuMi->pz(), gmuMi->p());
    if(i<2) tgmu4i[i+2].SetPxPyPzE(gmuMi->px(), gmuMi->py(), gmuMi->pz(), gmuMi->p());
  }
  if(muonsptsortedPi.size()>=2) hHPPmass_isol->Fill((tgmuPi[0]+tgmuPi[1]).M());
  if(muonsptsortedMi.size()>=2) hHMMmass_isol->Fill((tgmuMi[0]+tgmuMi[1]).M());
  if(muonsptsortedPi.size()>=2 && muonsptsortedMi.size()>=2)
  {
    hHPPHMMmass_isol->Fill((tgmuPi[0]+tgmuPi[1]+tgmuMi[0]+tgmuMi[1]).M());
    h2gmudeltaphi_isol->Fill(fabs(deltaphi((tgmuPi[0]+tgmuPi[1]).Phi(),(tgmuMi[0]+tgmuMi[1]).Phi())));
  }

  for (unsigned int i=0; i<muonsptsortedPi2.size(); i++)
  {
    const Track *gmuPi2 = & (*MuCollection)[muonsptsortedPi2[i].second];
    tgmuPi2[i].SetPxPyPzE(gmuPi2->px(), gmuPi2->py(), gmuPi2->pz(), gmuPi2->p());
//    if(i<2) tgmu4i[i].SetPxPyPzE(gmuPi->px(), gmuPi->py(), gmuPi->pz(), gmuPi->p());
  }
  for (unsigned int i=0; i<muonsptsortedMi2.size(); i++)
  {
    const Track *gmuMi2 = & (*MuCollection)[muonsptsortedMi2[i].second];
    tgmuMi2[i].SetPxPyPzE(gmuMi2->px(), gmuMi2->py(), gmuMi2->pz(), gmuMi2->p());
//    if(i<2) tgmu4i[i+2].SetPxPyPzE(gmuMi->px(), gmuMi->py(), gmuMi->pz(), gmuMi->p());
  }
  if(muonsptsortedPi2.size()>=2) hHPPmass_isol2->Fill((tgmuPi2[0]+tgmuPi2[1]).M());
  if(muonsptsortedMi2.size()>=2) hHMMmass_isol2->Fill((tgmuMi2[0]+tgmuMi2[1]).M());
  if(muonsptsortedPi2.size()>=2 && muonsptsortedMi2.size()>=2)
  {
    hHPPHMMmass_isol2->Fill((tgmuPi2[0]+tgmuPi2[1]+tgmuMi2[0]+tgmuMi2[1]).M());
    h2gmudeltaphi_isol2->Fill(fabs(deltaphi((tgmuPi2[0]+tgmuPi2[1]).Phi(),(tgmuMi2[0]+tgmuMi2[1]).Phi())));
    if(fabs(deltaphi((tgmuPi2[0]+tgmuPi2[1]).Phi(),(tgmuMi2[0]+tgmuMi2[1]).Phi()))>dphicut_4mu);
    {
      hHPPmass_isol2dphi->Fill((tgmuPi2[0]+tgmuPi2[1]).M());
      hHMMmass_isol2dphi->Fill((tgmuMi2[0]+tgmuMi2[1]).M());
      hHPPHMMmass_isol2dphi->Fill((tgmuPi2[0]+tgmuPi2[1]+tgmuMi2[0]+tgmuMi2[1]).M());
    }
  }
  for (unsigned int i=0; i<muonsptsortedi2.size(); i++)
  {
    const Track *gmui2 = & (*MuCollection)[muonsptsortedi2[i].second];
    TLorentzVector tgmui2(gmui2->px(), gmui2->py(), gmui2->pz(), gmui2->p());
    double dR_trgmui2=0.1, trpt_gmumatch=0;
    for (unsigned int n=0; n<trptsortedi2.size(); n++)
    {
      const Track *tri2 = & (*TrCollection)[trptsortedi2[n].second];
      TLorentzVector ttri2(tri2->px(), tri2->py(), tri2->pz(), tri2->p());
      if(tgmui2.DeltaR(ttri2)<dR_trgmui2)
      {
        dR_trgmui2=tgmui2.DeltaR(ttri2);
        trpt_gmumatch=tri2->pt();
      }
    }
    if(dR_trgmui2<0.01) htrgmupt_ptetaisol2cut->Fill(trpt_gmumatch,gmui2->pt());
  }

if(debug) cout<<"180"<<endl;
  double gmuPi2jz[2]={0,0}, gmuMi2jz[2]={0,0};
  for (unsigned int i=0; i<muonsptsortedPi2j.size(); i++)
  {
    const Track *gmuPi2j = & (*MuCollection)[muonsptsortedPi2j[i].second];
    tgmuPi2j[i].SetPxPyPzE(gmuPi2j->px(), gmuPi2j->py(), gmuPi2j->pz(), gmuPi2j->p());
    if(i<2) gmuPi2jz[i]=gmuPi2j->dz();
//    if(i<2) tgmu4i[i].SetPxPyPzE(gmuPi->px(), gmuPi->py(), gmuPi->pz(), gmuPi->p());
  }
  for (unsigned int i=0; i<muonsptsortedMi2j.size(); i++)
  {
    const Track *gmuMi2j = & (*MuCollection)[muonsptsortedMi2j[i].second];
    tgmuMi2j[i].SetPxPyPzE(gmuMi2j->px(), gmuMi2j->py(), gmuMi2j->pz(), gmuMi2j->p());
    if(i<2) gmuMi2jz[i]=gmuMi2j->dz();
//    if(i<2) tgmu4i[i+2].SetPxPyPzE(gmuMi->px(), gmuMi->py(), gmuMi->pz(), gmuMi->p());
  }
  if(muonsptsortedPi2j.size()>=2)
  {
    hHPPmass_isol2j->Fill((tgmuPi2j[0]+tgmuPi2j[1]).M());
    if(fabs(gmuPi2jz[0]-gmuPi2jz[1])<dzcut_)
    {
      hHPPmass_isol2jdz->Fill((tgmuPi2j[0]+tgmuPi2j[1]).M());
      hHPPmass_isol2jdz_ex->Fill((tgmuPi2j[0]+tgmuPi2j[1]).M());
    }
//    h2gmudz_isol2j->Fill(fabs(gmuPi2jz[0]-gmuPi2jz[1]));
    h2gmudz_isol2j->Fill(gmuPi2jz[0]-gmuPi2jz[1]); h2gmudz_isol2j_ex->Fill(gmuPi2jz[0]-gmuPi2jz[1]);
  }
  if(muonsptsortedMi2j.size()>=2)
  {
    hHMMmass_isol2j->Fill((tgmuMi2j[0]+tgmuMi2j[1]).M());
    if(fabs(gmuMi2jz[0]-gmuMi2jz[1])<dzcut_)
    {
      hHMMmass_isol2jdz->Fill((tgmuMi2j[0]+tgmuMi2j[1]).M());
      hHMMmass_isol2jdz_ex->Fill((tgmuMi2j[0]+tgmuMi2j[1]).M());
    }
//    h2gmudz_isol2j->Fill(fabs(gmuMi2jz[0]-gmuMi2jz[1]));
    h2gmudz_isol2j->Fill(gmuMi2jz[0]-gmuMi2jz[1]); h2gmudz_isol2j_ex->Fill(gmuMi2jz[0]-gmuMi2jz[1]);
  }
  if(muonsptsortedPi2j.size()>=2 && muonsptsortedMi2j.size()>=2)
  {
    hHPPHMMmass_isol2j->Fill((tgmuPi2j[0]+tgmuPi2j[1]+tgmuMi2j[0]+tgmuMi2j[1]).M());
  }
  for (unsigned int i=0; i<muonsptsortedi2j.size(); i++)
  {
    const Track *gmui2j = & (*MuCollection)[muonsptsortedi2j[i].second];
    TLorentzVector tgmui2j(gmui2j->px(), gmui2j->py(), gmui2j->pz(), gmui2j->p());
    double dR_trgmui2j=0.1, trpt_gmumatchj=0;
    for (unsigned int n=0; n<trptsortedi2j.size(); n++)
    {
      const Track *tri2j = & (*TrCollection)[trptsortedi2j[n].second];
      TLorentzVector ttri2j(tri2j->px(), tri2j->py(), tri2j->pz(), tri2j->p());
      if(tgmui2j.DeltaR(ttri2j)<dR_trgmui2j)
      {
        dR_trgmui2j=tgmui2j.DeltaR(ttri2j);
        trpt_gmumatchj=tri2j->pt();
      }
    }
    if(dR_trgmui2j<0.01) htrgmupt_ptetaisol2jcut->Fill(trpt_gmumatchj,gmui2j->pt());
    for (unsigned int n=0; n<muonsptsortedi2j.size(); n++)
    {
      const Track *gmui2j2 = & (*MuCollection)[muonsptsortedi2j[n].second];
//      if(i<4 && n<4 && n>i) hgmudz_ptetaisol2jcut->Fill(fabs(gmui2j->dz()-gmui2j2->dz()));
      if(i<4 && n<4 && n>i) hgmudz_ptetaisol2jcut->Fill(gmui2j->dz()-gmui2j2->dz());
    }
  }

if(debug) cout<<"190"<<endl;

//---------------------------matching particles(mat000)-----------------------------//
  float dR_match_tr[4]={0.1,0.1,0.1,0.1}, dR_matchi_tr[4]={0.1,0.1,0.1,0.1}, dR_matchi2_tr[4]={0.1,0.1,0.1,0.1};
  float dR_matchi2j_tr[4]={0.1,0.1,0.1,0.1};
  TLorentzVector ttr2[4], ttr2i[4], ttr2i2[4], ttr2i2j[4];
  float tr2i2jz[4]={0};
  float dR_match_gmu[4]={0.1,0.1,0.1,0.1}, dR_matchi_gmu[4]={0.1,0.1,0.1,0.1}, dR_matchi2_gmu[4]={0.1,0.1,0.1,0.1};
  float dR_matchi2j_gmu[4]={0.1,0.1,0.1,0.1};
  float dR_match_gmu3[4]={0.1,0.1,0.1,0.1};
  TLorentzVector tgmu2[4], tgmu2i[4], tgmu3[4], tgmu2i2[4], tgmu2i2j[4]; //matched global muon with MC muon from Higgs
  float gmu2i2jz[4]={0};

  for(int i=0;i<4;i++)
  {
    if(mcmu[i].Pt()>10 && fabs(mcmu[i].Eta())<2.4)
    {
      for (unsigned int n=0; n<trptsorted.size(); n++)
      {
        const Track *tr2 = & (*TrCollection)[trptsorted[n].second];
        TLorentzVector Ttr(tr2->px(), tr2->py(), tr2->pz(), tr2->p());
//        if(mcmu[i].DeltaR(tgmu)<dR_match_gmu[i] && fabs(mcmu[i].Pt()-tgmu.Pt())/mcmu[i].Pt()<0.2)
        if(mcmu[i].DeltaR(Ttr)<dR_match_tr[i])
        {
          if( (i<=1 && tr2->charge()==1)||(i>=2 && tr2->charge()==-1) )
          {
            dR_match_tr[i]=mcmu[i].DeltaR(Ttr);
            ttr2[i].SetPxPyPzE(tr2->px(), tr2->py(), tr2->pz(), tr2->p());
          }
        }
      }
      for (unsigned int n=0; n<trptsortedi.size(); n++)
      {
        const Track *tr2i = & (*TrCollection)[trptsortedi[n].second];
        TLorentzVector Ttri(tr2i->px(), tr2i->py(), tr2i->pz(), tr2i->p());
        if(mcmu[i].DeltaR(Ttri)<dR_matchi_tr[i] && mcmui[i])
        {
          if( (i<=1 && tr2i->charge()==1)||(i>=2 && tr2i->charge()==-1) )
          {
            dR_matchi_tr[i]=mcmu[i].DeltaR(Ttri);
            ttr2i[i].SetPxPyPzE(tr2i->px(), tr2i->py(), tr2i->pz(), tr2i->p());
          }
        }
      }
      if(dR_match_tr[i]<0.01)
      {
        htrackpt_mcmumatch->Fill(ttr2[i].Pt());
        htracketa_mcmumatch->Fill(ttr2[i].Eta());
      }
      if(dR_matchi_tr[i]<0.01)
      {
        htrackpt_mcmumatchisol->Fill(ttr2i[i].Pt());
        htracketa_mcmumatchisol->Fill(ttr2i[i].Eta());
      }

if(debug && i==0) cout<<"200"<<endl;
      int ngmumatch[4]={0,};
      for (unsigned int n=0; n<muonsptsorted.size(); n++)
      {
        const Track *gmu2 = & (*MuCollection)[muonsptsorted[n].second];
        TLorentzVector Tgmu(gmu2->px(), gmu2->py(), gmu2->pz(), gmu2->p());
        if(mcmu[i].DeltaR(Tgmu)<dR_match_gmu[i])
        {
          if( (i<=1 && gmu2->charge()==1)||(i>=2 && gmu2->charge()==-1) )
          {
            dR_match_gmu3[i]=dR_match_gmu[i];
            tgmu3[i].SetPxPyPzE(tgmu2[i].Px(), tgmu2[i].Py(), tgmu2[i].Pz(), tgmu2[i].P());
            dR_match_gmu[i]=mcmu[i].DeltaR(Tgmu);
            tgmu2[i].SetPxPyPzE(gmu2->px(), gmu2->py(), gmu2->pz(), gmu2->p());
          }
        }
        if(mcmu[i].DeltaR(Tgmu)<0.01) ngmumatch[i]++;
      }
      if(ngmumatch[i]>=2) cout<<"# of matched global muon with a MC signal muon : "<<ngmumatch[i]<<endl;
      for (unsigned int n=0; n<muonsptsortedi.size(); n++)
      {
        const Track *gmu2i = & (*MuCollection)[muonsptsortedi[n].second];
        TLorentzVector Tgmui(gmu2i->px(), gmu2i->py(), gmu2i->pz(), gmu2i->p());
        if(mcmu[i].DeltaR(Tgmui)<dR_matchi_gmu[i] && mcmui[i])
        {
          if( (i<=1 && gmu2i->charge()==1)||(i>=2 && gmu2i->charge()==-1) )
          {
            dR_matchi_gmu[i]=mcmu[i].DeltaR(Tgmui);
            tgmu2i[i].SetPxPyPzE(gmu2i->px(), gmu2i->py(), gmu2i->pz(), gmu2i->p());
          }
        }
      }
      if(dR_match_gmu[i]<0.01)
      {
        hgmumcmu_ptratio->Fill(mcmu[i].Pt(),fabs(tgmu2[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt());
        hgmumcmu_ptratio_ex->Fill(mcmu[i].Pt(),fabs(tgmu2[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt());
        hmcmupt_gmumatch->Fill(mcmu[i].Pt());
        hmcmueta_gmumatch->Fill(mcmu[i].Eta());
        hgmupt_mcmumatch->Fill(tgmu2[i].Pt());
        hgmueta_mcmumatch->Fill(tgmu2[i].Eta());
        if(dR_match_gmu3[i]<0.1)
        {
          hgmumcmu_ptratio_2nd->Fill(mcmu[i].Pt(),fabs(tgmu3[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt());
          hgmumcmu_ptratio_2nd_ex->Fill(mcmu[i].Pt(),fabs(tgmu3[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt());
          if(fabs(tgmu2[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt()>0.5)
          {
          hgmumcmu_ptratio_2ndtrue->Fill(mcmu[i].Pt(),fabs(tgmu3[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt());
          hgmumcmu_ptratio_2ndtrue_ex->Fill(mcmu[i].Pt(),fabs(tgmu3[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt());
          }
        }
      }
      if(dR_matchi_gmu[i]<0.01)
      {
        hgmumcmu_ptratioisol->Fill(mcmu[i].Pt(),fabs(tgmu2i[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt());
        hgmumcmu_ptratioisol_ex->Fill(mcmu[i].Pt(),fabs(tgmu2i[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt());
        hmcmupt_gmumatchisol->Fill(mcmu[i].Pt());
        hmcmueta_gmumatchisol->Fill(mcmu[i].Eta());
        hgmupt_mcmumatchisol->Fill(tgmu2i[i].Pt());
        hgmueta_mcmumatchisol->Fill(tgmu2i[i].Eta());
      }
if(debug && i==0) cout<<"210"<<endl;
//--------------------new isolation --------------------//
      for (unsigned int n=0; n<trptsortedi2.size(); n++)
      {
        const Track *tr2i2 = & (*TrCollection)[trptsortedi2[n].second];
        TLorentzVector Ttri2(tr2i2->px(), tr2i2->py(), tr2i2->pz(), tr2i2->p());
        if(mcmu[i].DeltaR(Ttri2)<dR_matchi2_tr[i] && mcmui2[i])
        {
          if( (i<=1 && tr2i2->charge()==1)||(i>=2 && tr2i2->charge()==-1) )
          {
            dR_matchi2_tr[i]=mcmu[i].DeltaR(Ttri2);
            ttr2i2[i].SetPxPyPzE(tr2i2->px(), tr2i2->py(), tr2i2->pz(), tr2i2->p());
          }
        }
      }
      if(dR_matchi2_tr[i]<0.01)
      {
        htrackpt_mcmumatchisol2->Fill(ttr2i2[i].Pt());
        htracketa_mcmumatchisol2->Fill(ttr2i2[i].Eta());
      }
      for (unsigned int n=0; n<muonsptsortedi2.size(); n++)
      {
        const Track *gmu2i2 = & (*MuCollection)[muonsptsortedi2[n].second];
        TLorentzVector Tgmui2(gmu2i2->px(), gmu2i2->py(), gmu2i2->pz(), gmu2i2->p());
        if(mcmu[i].DeltaR(Tgmui2)<dR_matchi2_gmu[i] && mcmui2[i])
        {
          if( (i<=1 && gmu2i2->charge()==1)||(i>=2 && gmu2i2->charge()==-1) )
          {
            dR_matchi2_gmu[i]=mcmu[i].DeltaR(Tgmui2);
            tgmu2i2[i].SetPxPyPzE(gmu2i2->px(), gmu2i2->py(), gmu2i2->pz(), gmu2i2->p());
          }
        }
      }
      if(dR_matchi2_gmu[i]<0.01)
      {
        hgmumcmu_ptratioisol2->Fill(mcmu[i].Pt(),fabs(tgmu2i2[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt());
        hgmumcmu_ptratioisol2_ex->Fill(mcmu[i].Pt(),fabs(tgmu2i2[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt());
        hmcmupt_gmumatchisol2->Fill(mcmu[i].Pt());
        hmcmueta_gmumatchisol2->Fill(mcmu[i].Eta());
        hgmupt_mcmumatchisol2->Fill(tgmu2i2[i].Pt());
        hgmueta_mcmumatchisol2->Fill(tgmu2i2[i].Eta());
      }
if(debug && i==0) cout<<"220"<<endl;
//--------------------jet isolation --------------------//
      for (unsigned int n=0; n<trptsortedi2j.size(); n++)
      {
        const Track *tr2i2j = & (*TrCollection)[trptsortedi2j[n].second];
        TLorentzVector Ttri2j(tr2i2j->px(), tr2i2j->py(), tr2i2j->pz(), tr2i2j->p());
        if(mcmu[i].DeltaR(Ttri2j)<dR_matchi2j_tr[i] && mcmui2j[i])
        {
          if( (i<=1 && tr2i2j->charge()==1)||(i>=2 && tr2i2j->charge()==-1) )
          {
            dR_matchi2j_tr[i]=mcmu[i].DeltaR(Ttri2j);
            ttr2i2j[i].SetPxPyPzE(tr2i2j->px(), tr2i2j->py(), tr2i2j->pz(), tr2i2j->p());
            tr2i2jz[i]=tr2i2j->dz();
          }
        }
      }
      if(dR_matchi2j_tr[i]<0.01)
      {
        htrackpt_mcmumatchisol2j->Fill(ttr2i2j[i].Pt());
        htracketa_mcmumatchisol2j->Fill(ttr2i2j[i].Eta());
      }
      for (unsigned int n=0; n<muonsptsortedi2j.size(); n++)
      {
        const Track *gmu2i2j = & (*MuCollection)[muonsptsortedi2j[n].second];
        TLorentzVector Tgmui2j(gmu2i2j->px(), gmu2i2j->py(), gmu2i2j->pz(), gmu2i2j->p());
        if(mcmu[i].DeltaR(Tgmui2j)<dR_matchi2j_gmu[i] && mcmui2j[i])
        {
          if( (i<=1 && gmu2i2j->charge()==1)||(i>=2 && gmu2i2j->charge()==-1) )
          {
            dR_matchi2j_gmu[i]=mcmu[i].DeltaR(Tgmui2j);
            tgmu2i2j[i].SetPxPyPzE(gmu2i2j->px(), gmu2i2j->py(), gmu2i2j->pz(), gmu2i2j->p());
            gmu2i2jz[i]=gmu2i2j->dz();
          }
        }
      }
      if(dR_matchi2j_gmu[i]<0.01)
      {
        hgmumcmu_ptratioisol2j->Fill(mcmu[i].Pt(),fabs(tgmu2i2j[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt());
        hgmumcmu_ptratioisol2j_ex->Fill(mcmu[i].Pt(),fabs(tgmu2i2j[i].Pt()-mcmu[i].Pt())/mcmu[i].Pt());
//        hmcmupt_gmumatchisol2->Fill(mcmu[i].Pt());
//        hmcmueta_gmumatchisol2->Fill(mcmu[i].Eta());
        hgmupt_mcmumatchisol2j->Fill(tgmu2i2[i].Pt());
        hgmueta_mcmumatchisol2j->Fill(tgmu2i2[i].Eta());
      }
    }
  }
if(debug) cout<<"230"<<endl;

  if(dR_match_tr[0]<0.01&&dR_match_tr[1]<0.01) hHPPmass_trmcmumatch->Fill((ttr2[0]+ttr2[1]).M());
  if(dR_match_tr[2]<0.01&&dR_match_tr[3]<0.01) hHMMmass_trmcmumatch->Fill((ttr2[2]+ttr2[3]).M());
  if(dR_match_tr[0]<0.01&&dR_match_tr[1]<0.01&&dR_match_tr[2]<0.01&&dR_match_tr[3]<0.01) hHPPHMMmass_trmcmumatch->Fill((ttr2[0]+ttr2[1]+ttr2[2]+ttr2[3]).M());
  if(dR_match_gmu[0]<0.01&&dR_match_gmu[1]<0.01) hHPPmass_mcmumatch->Fill((tgmu2[0]+tgmu2[1]).M());
  if(dR_match_gmu[2]<0.01&&dR_match_gmu[3]<0.01) hHMMmass_mcmumatch->Fill((tgmu2[2]+tgmu2[3]).M());
  if(dR_match_gmu[0]<0.01&&dR_match_gmu[1]<0.01&&dR_match_gmu[2]<0.01&&dR_match_gmu[3]<0.01) hHPPHMMmass_mcmumatch->Fill((tgmu2[0]+tgmu2[1]+tgmu2[2]+tgmu2[3]).M());
  if(dR_matchi_tr[0]<0.01&&dR_matchi_tr[1]<0.01) hHPPmass_trmcmumatchisol->Fill((ttr2i[0]+ttr2i[1]).M());
  if(dR_matchi_tr[2]<0.01&&dR_matchi_tr[3]<0.01) hHMMmass_trmcmumatchisol->Fill((ttr2i[2]+ttr2i[3]).M());
  if(dR_matchi_tr[0]<0.01&&dR_matchi_tr[1]<0.01&&dR_matchi_tr[2]<0.01&&dR_matchi_tr[3]<0.01) hHPPHMMmass_trmcmumatchisol->Fill((ttr2i[0]+ttr2i[1]+ttr2i[2]+ttr2i[3]).M());
  if(dR_matchi_gmu[0]<0.01&&dR_matchi_gmu[1]<0.01) hHPPmass_mcmumatchisol->Fill((tgmu2i[0]+tgmu2i[1]).M());
  if(dR_matchi_gmu[2]<0.01&&dR_matchi_gmu[3]<0.01) hHMMmass_mcmumatchisol->Fill((tgmu2i[2]+tgmu2i[3]).M());
  if(dR_matchi_gmu[0]<0.01&&dR_matchi_gmu[1]<0.01&&dR_matchi_gmu[2]<0.01&&dR_matchi_gmu[3]<0.01) hHPPHMMmass_mcmumatchisol->Fill((tgmu2i[0]+tgmu2i[1]+tgmu2i[2]+tgmu2i[3]).M());
  if(dR_matchi2_tr[0]<0.01&&dR_matchi2_tr[1]<0.01) hHPPmass_trmcmumatchisol2->Fill((ttr2i2[0]+ttr2i2[1]).M());
  if(dR_matchi2_tr[2]<0.01&&dR_matchi2_tr[3]<0.01) hHMMmass_trmcmumatchisol2->Fill((ttr2i2[2]+ttr2i2[3]).M());
  if(dR_matchi2_tr[0]<0.01&&dR_matchi2_tr[1]<0.01&&dR_matchi2_tr[2]<0.01&&dR_matchi2_tr[3]<0.01)
  {
    hHPPHMMmass_trmcmumatchisol2->Fill((ttr2i2[0]+ttr2i2[1]+ttr2i2[2]+ttr2i2[3]).M());
    hHPPmass_trmcmumatchisol24mu->Fill((ttr2i2[0]+ttr2i2[1]).M());
    hHMMmass_trmcmumatchisol24mu->Fill((ttr2i2[2]+ttr2i2[3]).M());
    if(fabs(deltaphi((mcmu[0]+mcmu[1]).Phi(),(mcmu[2]+mcmu[3]).Phi()))>dphicut_4mu)
    if(fabs(deltaphi((ttr2i2[0]+ttr2i2[1]).Phi(),(ttr2i2[2]+ttr2i2[3]).Phi()))>dphicut_4mu)
    {
      hHPPmass_trmcmumatchisol2dphi->Fill((ttr2i2[0]+ttr2i2[1]).M());
      hHMMmass_trmcmumatchisol2dphi->Fill((ttr2i2[2]+ttr2i2[3]).M());
      hHPPHMMmass_trmcmumatchisol2dphi->Fill((ttr2i2[0]+ttr2i2[1]+ttr2i2[2]+ttr2i2[3]).M());
    }
    h2tr2mcmudphi_isol2->Fill(fabs(deltaphi((ttr2i2[0]+ttr2i2[1]).Phi(),(ttr2i2[2]+ttr2i2[3]).Phi())),fabs(deltaphi((mcmu[0]+mcmu[1]).Phi(),(mcmu[2]+mcmu[3]).Phi())));
  }
  if(dR_matchi2_gmu[0]<0.01&&dR_matchi2_gmu[1]<0.01) hHPPmass_mcmumatchisol2->Fill((tgmu2i2[0]+tgmu2i2[1]).M());
  if(dR_matchi2_gmu[2]<0.01&&dR_matchi2_gmu[3]<0.01) hHMMmass_mcmumatchisol2->Fill((tgmu2i2[2]+tgmu2i2[3]).M());
  if(dR_matchi2_gmu[0]<0.01&&dR_matchi2_gmu[1]<0.01&&dR_matchi2_gmu[2]<0.01&&dR_matchi2_gmu[3]<0.01)
  {
    hHPPHMMmass_mcmumatchisol2->Fill((tgmu2i2[0]+tgmu2i2[1]+tgmu2i2[2]+tgmu2i2[3]).M());
    hHPPmass_mcmumatchisol24mu->Fill((tgmu2i2[0]+tgmu2i2[1]).M());
    hHMMmass_mcmumatchisol24mu->Fill((tgmu2i2[2]+tgmu2i2[3]).M());
    if(fabs(deltaphi((mcmu[0]+mcmu[1]).Phi(),(mcmu[2]+mcmu[3]).Phi()))>dphicut_4mu)
    if(fabs(deltaphi((tgmu2i2[0]+tgmu2i2[1]).Phi(),(tgmu2i2[2]+tgmu2i2[3]).Phi()))>dphicut_4mu)
    {
      hHPPmass_mcmumatchisol2dphi->Fill((tgmu2i2[0]+tgmu2i2[1]).M());
      hHMMmass_mcmumatchisol2dphi->Fill((tgmu2i2[2]+tgmu2i2[3]).M());
      hHPPHMMmass_mcmumatchisol2dphi->Fill((tgmu2i2[0]+tgmu2i2[1]+tgmu2i2[2]+tgmu2i2[3]).M());
    }
    h2gmu2mcmudphi_isol2->Fill(fabs(deltaphi((tgmu2i2[0]+tgmu2i2[1]).Phi(),(tgmu2i2[2]+tgmu2i2[3]).Phi())),fabs(deltaphi((mcmu[0]+mcmu[1]).Phi(),(mcmu[2]+mcmu[3]).Phi())));
  }
  if(dR_matchi2j_tr[0]<0.01&&dR_matchi2j_tr[1]<0.01)
  {
    hHPPmass_trmcmumatchisol2j->Fill((ttr2i2j[0]+ttr2i2j[1]).M());
    if(fabs(tr2i2jz[0]-tr2i2jz[1])<dzcut_)
    {
      hHPPmass_trmcmumatchisol2jdz->Fill((ttr2i2j[0]+ttr2i2j[1]).M());
      hHPPmass_trmcmumatchisol2jdz_ex->Fill((ttr2i2j[0]+ttr2i2j[1]).M());
    }
  }
  if(dR_matchi2j_tr[2]<0.01&&dR_matchi2j_tr[3]<0.01)
  {
    hHMMmass_trmcmumatchisol2j->Fill((ttr2i2j[2]+ttr2i2j[3]).M());
    if(fabs(tr2i2jz[2]-tr2i2jz[3])<dzcut_)
    {
      hHMMmass_trmcmumatchisol2jdz->Fill((ttr2i2j[2]+ttr2i2j[3]).M());
      hHMMmass_trmcmumatchisol2jdz_ex->Fill((ttr2i2j[2]+ttr2i2j[3]).M());
    }
  }
  if(dR_matchi2j_tr[0]<0.01&&dR_matchi2j_tr[1]<0.01&&dR_matchi2j_tr[2]<0.01&&dR_matchi2j_tr[3]<0.01)
  {
    hHPPHMMmass_trmcmumatchisol2j->Fill((ttr2i2j[0]+ttr2i2j[1]+ttr2i2j[2]+ttr2i2j[3]).M());
  }
  if(dR_matchi2j_gmu[0]<0.01&&dR_matchi2j_gmu[1]<0.01)
  {
    hHPPmass_mcmumatchisol2j->Fill((tgmu2i2j[0]+tgmu2i2j[1]).M());
    if(fabs(gmu2i2jz[0]-gmu2i2jz[1])<dzcut_)
    {
      hHPPmass_mcmumatchisol2jdz->Fill((tgmu2i2j[0]+tgmu2i2j[1]).M());
      hHPPmass_mcmumatchisol2jdz_ex->Fill((tgmu2i2j[0]+tgmu2i2j[1]).M());
    }
  }
  if(dR_matchi2j_gmu[2]<0.01&&dR_matchi2j_gmu[3]<0.01)
  {
    hHMMmass_mcmumatchisol2j->Fill((tgmu2i2j[2]+tgmu2i2j[3]).M());
    if(fabs(gmu2i2jz[2]-gmu2i2jz[3])<dzcut_)
    {
      hHMMmass_mcmumatchisol2jdz->Fill((tgmu2i2j[2]+tgmu2i2j[3]).M());
      hHMMmass_mcmumatchisol2jdz_ex->Fill((tgmu2i2j[2]+tgmu2i2j[3]).M());
    }
  }
  if(dR_matchi2j_gmu[0]<0.01&&dR_matchi2j_gmu[1]<0.01&&dR_matchi2j_gmu[2]<0.01&&dR_matchi2j_gmu[3]<0.01)
  {
    hHPPHMMmass_mcmumatchisol2j->Fill((tgmu2i2j[0]+tgmu2i2j[1]+tgmu2i2j[2]+tgmu2i2j[3]).M());
  }
if(debug) cout<<"240"<<endl;
//---------------------------MET-----------------------------//
  Handle<CaloMETCollection> cmetJets;
  Handle<GenMETCollection> gmetJets;
  Handle<GenMETCollection> gmetJets2;

  iEvent.getByLabel("met",cmetJets);
  iEvent.getByLabel("genMet",gmetJets);
  iEvent.getByLabel("genMetNoNuBSM",gmetJets2);
//cout<<"----------"<<cmetJets->sumet()<<"----------"<<endl;
//cout<<"----------"<<cmetJets->maxEtInEmTowers()<<"----------"<<endl;
//cout<<"----------"<<gmetJets->pt()<<"----------"<<endl;
//cout<<"----------"<<gmetJets->emEnergy()<<"----------"<<endl;
  double calmet=0, calmet3=0, genmet=0, genmet2=0, genmet3=0;
  for( CaloMETCollection::const_iterator cmet = cmetJets->begin(); cmet != cmetJets->end();++cmet )
  {
    calmet=cmet->et();
    calmet3=cmet->et();
    for (unsigned int i=0; i<muonsptsortedi2j.size(); i++)
    {
      const Track *gmu3i2j = & (*MuCollection)[muonsptsortedi2j[i].second];
      calmet3+=gmu3i2j->pt()*0.02;
    }
  }
  for( GenMETCollection::const_iterator gmet = gmetJets->begin(); gmet != gmetJets->end();++gmet )
  {
    genmet=gmet->et();
    genmet3=gmet->et();
    for(int i=0;i<4;i++) if(fabs(mcmu[i].Eta())<5.0) genmet3-=mcmu[i].Et()*0.01;
//cout<<gmet->sumet_()<<endl;
//cout<<gmet->emEnergy()<<endl;
//cout<<gmet->hadEnergy()<<endl;
//cout<<gmet->invisibleEnergy()<<endl;
//cout<<gmet->auxiliaryEnergy()<<endl;
  }
  for( GenMETCollection::const_iterator gmet2 = gmetJets2->begin(); gmet2 != gmetJets2->end();++gmet2 )
  {
    genmet2=gmet2->et();
  }
  hcalmet->Fill(calmet);
  hcalmet01mu->Fill(calmet3);
  hgenmet->Fill(genmet);
  hgenmetnonu->Fill(genmet2);
  hgenmet01mu->Fill(genmet3);
  hcalgenmet->Fill(calmet,genmet);
  hcalgenmetnonu->Fill(calmet,genmet2);
  hcalgenmet01mu->Fill(calmet,genmet3);
  hcal01mugenmet->Fill(calmet3,genmet);
  if(muidx_ptetaisol2jcut==0) hcalmet0gmu->Fill(calmet);
  if(muidx_ptetaisol2jcut==0) hgenmet0gmu->Fill(genmet);
  if(muidx_ptetaisol2jcut==0) hcalgenmet0gmu->Fill(calmet,genmet);
  if(muidx_ptetaisol2jcut==1) hcalmet1gmu->Fill(calmet);
  if(muidx_ptetaisol2jcut==1) hgenmet1gmu->Fill(genmet);
  if(muidx_ptetaisol2jcut==1) hcalgenmet1gmu->Fill(calmet,genmet);
  if(muidx_ptetaisol2jcut==2) hcalmet2gmu->Fill(calmet);
  if(muidx_ptetaisol2jcut==2) hgenmet2gmu->Fill(genmet);
  if(muidx_ptetaisol2jcut==2) hcalgenmet2gmu->Fill(calmet,genmet);
  if(muidx_ptetaisol2jcut==3) hcalmet3gmu->Fill(calmet);
  if(muidx_ptetaisol2jcut==3) hgenmet3gmu->Fill(genmet);
  if(muidx_ptetaisol2jcut==3) hcalgenmet3gmu->Fill(calmet,genmet);
  if(muidx_ptetaisol2jcut==4) hcalmet4gmu->Fill(calmet);
  if(muidx_ptetaisol2jcut==4) hgenmet4gmu->Fill(genmet);
  if(muidx_ptetaisol2jcut==4) hcalgenmet4gmu->Fill(calmet,genmet);

if(debug) cout<<"250"<<endl;
//---------------------------Vertex(line000, vt000)-----------------------------//

//  const reco::VertexCollection tC = *(VtCollection.product());
//  cout << "Reconstructed "<< tC.size() << " vertices (Events : "<<events<<")" << endl ;

//  if (tC.size() >0)
//  {
//    cout<<" PARAMS "<<tC.front().position()<< endl;
//    cout<<" COV "<<tC.front().covariance()<< endl;
//    cout <<"error  " <<tC.front().covariance(2,2)<< endl<< endl;
//  }

  hnvt->Fill(VtCollection->size());
  hnvt2->Fill(VtCollection2->size());
  hnvt3->Fill(VtCollection3->size());
  for( VertexCollection::const_iterator vt = VtCollection->begin(); vt != VtCollection->end();++vt )
  {
//    cout<<vt->position()<<"1 ";
    hvtx->Fill(vt->x()); hvty->Fill(vt->y()); hvtz->Fill(vt->z()); hvtxy->Fill(vt->x(),vt->y()); //hmcmurecovtxy->Fill(vt->x()-x0,vt->y()-y0);
//    if(mcmu1!=0 && fabs(mcmu1->vz())>12) cout<<vt->x()<<", "<<vt->y()<<", "<<vt->z()<<" // ";
//    int i=0; if(vt->x()==x0 && vt->y()==y0 && vt->z()==z0) cout<<VtCollection->size()<<", "<<mcmu1->vx()<<", "<<mcmu1->vy()<<", "<<mcmu1->vz()<<", "<<i<<endl; i++;
    if(vt->x()==x0 && vt->y()==y0 && vt->z()==z0) {hmcmux_BS->Fill(mcmu1->vx()); hmcmuy_BS->Fill(mcmu1->vy()); hmcmuz_BS->Fill(mcmu1->vz()); hmcmuxy_BS->Fill(mcmu1->vx(),mcmu1->vy());}
//    if(VtCollection2->size()==0) cout<<vt->x()<<", "<<vt->y()<<", "<<vt->z()<<endl;
    if(!(vt->x()==x0 && vt->y()==y0 && vt->z()==z0)) hvtz_BS->Fill(vt->z());
//    if(!(vt->x()==x0 && vt->y()==y0 && vt->z()==z0)) cout<<mcmu1->vz()<<", "<<vt->z()<<endl;
    if(vt->x()==x0 && vt->y()==y0 && vt->z()==z0) hntr_BS->Fill(TrCollection->size());
//    if(vt->x()==x0 && vt->y()==y0 && vt->z()==z0) cout<<TrCollection->size()<<", "<<vt->tracksSize()<<endl;
    hvtntr->Fill(vt->tracksSize());
//    if(vt->x()==x0 && vt->y()==y0 && vt->z()==z0) cout<<"ntr_i2j(recoVertex=BeamSpot) "<<tridx_ptetaisol2jcut<<",   ntr "<<TrCollection->size()<<endl;
    if(vt->x()==x0 && vt->y()==y0 && vt->z()==z0) hntr_i2jBS->Fill(tridx_ptetaisol2jcut);
    if(!(vt->x()==x0 && vt->y()==y0 && vt->z()==z0)) hvtchi2_BS->Fill(vt->normalizedChi2());
  }
//cout<<endl;
  for( VertexCollection::const_iterator vt = VtCollection2->begin(); vt != VtCollection2->end();++vt )
  {
//    cout<<vt->position()<<"2 ";
    hvt2x->Fill(vt->x()); hvt2y->Fill(vt->y()); hvt2z->Fill(vt->z()); hvt2xy->Fill(vt->x(),vt->y()); //hmcmurecovt2xy->Fill(vt->x(),vt->y());
//    if(mcmu1!=0 && fabs(mcmu1->vz())>12) cout<<vt->x()<<", "<<vt->y()<<", "<<vt->z()<<" // ";
//    if(vt->x()==x0 && vt->y()==y0 && vt->z()==z0) cout<<"vt2 : "<<mcmu1->vx()<<", "<<mcmu1->vy()<<", "<<mcmu1->vz()<<" // ";
    if(vt->x()==x0 && vt->y()==y0 && vt->z()==z0) {hmcmu2x_BS->Fill(mcmu1->vx()); hmcmu2y_BS->Fill(mcmu1->vy()); hmcmu2z_BS->Fill(mcmu1->vz()); hmcmu2xy_BS->Fill(mcmu1->vx(),mcmu1->vy());}
    if(!(vt->x()==x0 && vt->y()==y0 && vt->z()==z0)) hvt2z_BS->Fill(vt->z());
    hvt2ntr->Fill(vt->tracksSize());
    if(!(vt->x()==x0 && vt->y()==y0 && vt->z()==z0)) hvt2chi2_BS->Fill(vt->normalizedChi2());
  }
//cout<<endl;
  for( VertexCollection::const_iterator vt = VtCollection3->begin(); vt != VtCollection3->end();++vt )
  {
//    cout<<vt->position()<<"3 ";
    hvt3x->Fill(vt->x()); hvt3y->Fill(vt->y()); hvt3z->Fill(vt->z()); hvt3xy->Fill(vt->x(),vt->y()); //hmcmurecovt3xy->Fill(vt->x(),vt->y());
//    if(mcmu1!=0 && fabs(mcmu1->vz())>12) cout<<vt->x()<<", "<<vt->y()<<", "<<vt->z()<<endl;
//    if(vt->x()==x0 && vt->y()==y0 && vt->z()==z0) cout<<"vt3 : "<<mcmu1->vx()<<", "<<mcmu1->vy()<<", "<<mcmu1->vz()<<endl;
    if(vt->x()==x0 && vt->y()==y0 && vt->z()==z0) {hmcmu3x_BS->Fill(mcmu1->vx()); hmcmu3y_BS->Fill(mcmu1->vy()); hmcmu3z_BS->Fill(mcmu1->vz()); hmcmu3xy_BS->Fill(mcmu1->vx(),mcmu1->vy());}
    if(!(vt->x()==x0 && vt->y()==y0 && vt->z()==z0)) hvt3z_BS->Fill(vt->z());
    hvt3ntr->Fill(vt->tracksSize());
    if(!(vt->x()==x0 && vt->y()==y0 && vt->z()==z0)) hvt3chi2_BS->Fill(vt->normalizedChi2());
  }

if(debug) cout<<"280"<<endl;
//  if(int(events)%100==0) cout<<"# of zero vertex = "<<nzerovt<<", # of zero vertex2 = "<<nzerovt2<<", # of zero vertex3 = "<<nzerovt3<<endl;
//  if(int(events)%100==0) cout<<"# of zero vertexz = "<<nzerovtz<<", # of zero vertex2z = "<<nzerovt2z<<", # of zero vertex3z = "<<nzerovt3z<<endl;
//cout<<endl<<endl;


//---------------------------vertex(vt001)-----------------------------//
//  ESHandle<TransientTrackBuilder> theB;
//  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
//  vector<TransientTrack> TTr = (*theB).build(TrCollection);

//  vector<TransientTrack> theTransientTrack = theTTBuilder->build(TrCollection);
//cout<<theTransientTrack.vertex().z()<<endl;//fail
//  AdaptiveVertexReconstructor finder;
//  const AdaptiveVertexReconstructor *finder; //pass
//  VertexReconstructor finder;
//  vector<TransientVertex> vertices = finder.vertices ( vector<TransientTrack> TTr );
//  vector<TransientVertex> vertices = finder->vertices (TTr); //pass
//cout<<vertices.state()<<endl; //fail
//  for (TransientTrack::const_iterator ttr=TTr->begin(); ttr!=TTr->end(); ttr++) //fail
//    vector<TransientVertex> vertices = finder.vertices (ttr);
//  for (TransientVertex::const_iterator vt=vertices->begin(); vt!=vertices->end(); vt++) //fail

//  for (TrackCollection::const_iterator tr=TrCollection->begin(); tr!=TrCollection->end(); tr++)
//  {
//  }
//delete(mcmu1);

}


// ------------ method called once each job just before starting event loop  ------------
void HppMuMuAnalyzer::beginJob(const edm::EventSetup&)
//HppMuMuAnalyzer::beginJob(const edm::EventSetup& es)
{
//  edm::ESHandle<TransientTrackBuilder> builder;
//  es.get<TransientTrackRecord>().get("TransientTrackBuilder",builder);
//  theTTBuilder = builder.product();
  cout<<"Event start"<<endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void HppMuMuAnalyzer::endJob() 
{
  cout<<"end job"<<endl<<"Event "<<events<<endl<<"# of background MC muon = "<<nbackmcmuon<<endl<<"# of background MC muon (pt>10, |eta|<2.4) = "<<nbackmcmuon_ptetacut<<endl;

  double m_pt=0, m_eta=0, s_pt=0, s_eta=0;
  for(int i=0;i<50;i++) 
  {
    if(ncorgenjetpt_24etamucut[i]>0)
    {
      m_pt=Tptr_24etamucut[i]/ncorgenjetpt_24etamucut[i];
      s_pt=sqrt(Tptr2_24etamucut[i]/ncorgenjetpt_24etamucut[i]-pow(m_pt,2));
      hjetresponseVSpt_24etamucut->SetBinContent(i+1,m_pt);
      hjetresolutionVSpt_24etamucut->SetBinContent(i+1,s_pt);
    }
    if(ncorgenjeteta_24etamucut[i]>0)
    {
      m_eta=Teta_24etamucut[i]/ncorgenjeteta_24etamucut[i];
      s_eta=sqrt(Teta2_24etamucut[i]/ncorgenjeteta_24etamucut[i]-pow(m_eta,2));
      hjetresponseVSeta_24etamucut->SetBinContent(i+1,m_eta);
      hjetresolutionVSeta_24etamucut->SetBinContent(i+1,s_eta);
    }
  }
  double m_pt2=0, m_eta2=0, s_pt2=0, s_eta2=0;
  for(int i=0;i<50;i++) 
  {
    if(ncorgenjetpt_30pt24etamucut[i]>0)
    {
      m_pt2=Tptr_30pt24etamucut[i]/ncorgenjetpt_30pt24etamucut[i];
      s_pt2=sqrt(Tptr2_30pt24etamucut[i]/ncorgenjetpt_30pt24etamucut[i]-pow(m_pt2,2));
      hjetresponseVSpt_30pt24etamucut->SetBinContent(i+1,m_pt2);
      hjetresolutionVSpt_30pt24etamucut->SetBinContent(i+1,s_pt2);
    }
    if(ncorgenjeteta_30pt24etamucut[i]>0)
    {
      m_eta2=Teta_30pt24etamucut[i]/ncorgenjeteta_30pt24etamucut[i];
      s_eta2=sqrt(Teta2_30pt24etamucut[i]/ncorgenjeteta_30pt24etamucut[i]-pow(m_eta2,2));
      hjetresponseVSeta_30pt24etamucut->SetBinContent(i+1,m_eta2);
      hjetresolutionVSeta_30pt24etamucut->SetBinContent(i+1,s_eta2);
    }
  }
}

/* vim:set ts=2 sts=2 sw=2: */
