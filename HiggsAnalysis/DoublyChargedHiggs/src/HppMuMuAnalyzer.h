#ifndef HiggsAnalysis_DoublyChargedHiggs_HppCutAnalyzer_H
#define HiggsAnalysis_DoublyChargedHiggs_HppCutAnalyzer_H

// -*- C++ -*-
//
// Package:    DoublyChargedHiggs
// Class:      HppMuMuAnalyzer
//
/**\class HppMuMuAnalyzer HiggsAnalysis/DoublyChargedHiggs/src/HppMuMuAnalyzer.h

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

// system include files
#include <memory>
#include <map>
#include <vector>
#include <string>

#include "TH1F.h"
#include "TH2F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class HppMuMuAnalyzer : public edm::EDAnalyzer 
{
public:
	explicit HppMuMuAnalyzer(const edm::ParameterSet&);
	~HppMuMuAnalyzer();

private:
	virtual void beginJob(const edm::EventSetup&) ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	edm::InputTag Tracks_;
	double dzcut_, min, max, minvtx, maxvtx, minvty, maxvty, minmcvtx, maxmcvtx, minmcvty, maxmcvty;
	double minvt2x, maxvt2x, minvt2y, maxvt2y, minvt3x, maxvt3x, minvt3y, maxvt3y;
	double vtx_cor, vty_cor;
//    std::string JetCorrectionService_;
	std::string CaloJetAlgorithm_;
//    const TransientTrackBuilder * theTTBuilder;
//    std::string GenJetAlgorithm_;

	// ----------member data ---------------------------
	// ----------jet-------------------------------------
	TH1F *hncalojet, *hncalojet_ptcut;
	TH1F *hcalojetpt, *hcalojetpt_ptcut, *hcalojetpt_ptcut_ex;
	TH1F *hcalojeteta, *hcalojeteta_ptcut;

	TH1F *hncorjet, *hncorjet_ptcut, *hncorjet_ptmucut;
	TH1F *hcorjetpt, *hcorjetpt_ptcut, *hcorjetpt_ptcut_ex, *hcorjetpt_ptmucut, *hcorjetpt_ptmucut_ex;
	TH1F *hcorjeteta, *hcorjeteta_ptcut, *hcorjeteta_ptmucut;
	TH2F *hcormcmuptdR_ptetacut, *hcormcmuptdR_30ptetacut;

	TH1F *hngenjet, *hngenjet_ptcut, *hngenjet_ptmucut;
	TH1F *hgenjetpt, *hgenjetpt_ptcut, *hgenjetpt_ptcut_ex, *hgenjetpt_ptmucut, *hgenjetpt_ptmucut_ex;
	TH1F *hgenjetpt_pt24etamucut, *hgenjetpt_20pt24etamucut, *hgenjetpt_30pt24etamucut;
	TH1F *hgenjetpt_ptmucut_calojetmatch, *hgenjetpt_ptmucut_corjetmatch, *hgenjetpt_pt24etamucut_corjetmatch, *hgenjetpt_20pt24etamucut_corjetmatch, *hgenjetpt_30pt24etamucut_corjetmatch;
	TH1F *hgenjeteta, *hgenjeteta_ptcut, *hgenjeteta_ptmucut;
	TH1F *hgenjeteta_pt24etamucut, *hgenjeteta_20pt24etamucut, *hgenjeteta_30pt24etamucut;
	TH1F *hgenjeteta_ptmucut_calojetmatch, *hgenjeteta_ptmucut_corjetmatch, *hgenjeteta_pt24etamucut_corjetmatch, *hgenjeteta_20pt24etamucut_corjetmatch, *hgenjeteta_30pt24etamucut_corjetmatch;
	TH2F *hgenmcmuptdR_ptetacut, *hgenmcmuptdR_30ptetacut;

	TH1F *hjetresponseVSpt_24etamucut, *hjetresponseVSeta_24etamucut;
	TH1F *hjetresolutionVSpt_24etamucut, *hjetresolutionVSeta_24etamucut;
	TH1F *hjetresponseVSpt_30pt24etamucut, *hjetresponseVSeta_30pt24etamucut;
	TH1F *hjetresolutionVSpt_30pt24etamucut, *hjetresolutionVSeta_30pt24etamucut;

	TH1F *htrjetpt, *htrjetpt_mucut;
	TH2F *htrjetcorjetpt, *htrjetcorjetpt_mucut;
	// ----------track-------------------------------------
	TH1F *hntr, *hntr_ptetacut, *hntr_ptetaisolcut, *hntr_ptetaisol2cut, *hntr_ptetaisol2jcut, *hntr_i2j, *hntr_i2j3mu, *hntr_i2j4mu, *hntr_BS, *hntr_i2jBS;
	TH1F *htrackpt, *htrackpt_ptetacut, *htrackpt_ptetaisolcut, *htrackpt_ptetaisol2cut, *htrackpt_ptetaisol2jcut;
	TH1F *htrackpt_mcmumatch, *htrackpt_mcmumatchisol, *htrackpt_mcmumatchisol2, *htrackpt_mcmumatchisol2j;
	TH1F *htracketa, *htracketa_ptetacut, *htracketa_ptetaisolcut, *htracketa_ptetaisol2cut, *htracketa_ptetaisol2jcut;
	TH1F *htracketa_mcmumatch, *htracketa_mcmumatchisol, *htracketa_mcmumatchisol2, *htracketa_mcmumatchisol2j;
	TH1F *htrackisolPt_ptetacut01, *htrackisolPt_ptetacut1, *htrackisolPt_ptetacut10, *htrackisolPt_ptetacut100, *htrackisolPt_ptetacut1000;
	TH1F *htrackisolPt2_ptetacut10, *htrackisolPt2_ptetacut100;
	TH1F *h2trdeltaphi, *h2trdeltaphi_isol, *h2trdeltaphi_isol2;
	TH2F *h2tr2mcmudphi_isol2;
	TH1F *htrpt_i2junv;

	TH1F *htrackz, *htrackz_ptetacut, *htrackz_ptetaisol2cut, *htrackz_ptetaisol2jcut;
	TH1F *htrackz2, *htrackz2_ptetaisol2jcut;
	TH1F *htrackxy, *htrackxy_ptetaisol2jcut;
	TH1F *htracksz, *htracksz_ptetaisol2jcut;
	TH1F *htrackdz, *htrackdz_ptetaisol2jcut;
	TH1F *htrackdxy, *htrackdxy_ptetaisol2jcut;
	TH1F *htrackdsz, *htrackdsz_ptetaisol2jcut;
	TH1F *h2trdz_isol2j, *h2trdz_isol2j_ex;
	TH1F *htrackX, *htrackX_ptetaisol2jcut, *htrackdX, *htrackdX_ptetaisol2jcut;
	TH1F *htrackY, *htrackY_ptetaisol2jcut, *htrackdY, *htrackdY_ptetaisol2jcut;
	TH1F *htrackZ, *htrackZ_ptetaisol2jcut, *htrackdZ, *htrackdZ_ptetaisol2jcut, *htrackdZ_ptetaisol2jcut_ex;
	TH2F *htrackXY2D, *htrackXY2D_ptetaisol2jcut, *htrackXY2D_ptetaisol2jcut_beamP, *htrackXY2D_ptetaisol2jcut_beamM;
	TH1F *h2trdZ_isol2j, *h2trdZ_isol2j_ex;
	TH2F *hbeamXY2D;
	TH1F *hbeamz, *hbeamsigmaz;

	TH1F *hHPPmass_tr, *hHMMmass_tr, *hHPPHMMmass_tr;
	TH1F *hHPPmass_trisol, *hHMMmass_trisol, *hHPPHMMmass_trisol;
	TH1F *hHPPmass_trisol2, *hHMMmass_trisol2, *hHPPHMMmass_trisol2;
	TH1F *hHPPmass_trisol2j, *hHMMmass_trisol2j, *hHPPHMMmass_trisol2j;
	TH1F *hHPPmass_trisol2dphi, *hHMMmass_trisol2dphi, *hHPPHMMmass_trisol2dphi;
	TH1F *hHPPmass_trisol2jdz, *hHMMmass_trisol2jdz;
	TH1F *hHPPmass_trisol2jdz_ex, *hHMMmass_trisol2jdz_ex;
	TH1F *hHPPmass_trisol2jv, *hHMMmass_trisol2jv, *hHPPHMMmass_trisol2jv;
	TH1F *hHPPmass_trisol2jvc, *hHMMmass_trisol2jvc, *hHPPHMMmass_trisol2jvc;
	TH1F *hHPPmass_trisol2jvs, *hHMMmass_trisol2jvs, *hHPPHMMmass_trisol2jvs;
	TH1F *hHPPmass_trisol2jvcs, *hHMMmass_trisol2jvcs, *hHPPHMMmass_trisol2jvcs;

	TH1F *hHPPmass_trmcmumatch, *hHMMmass_trmcmumatch, *hHPPHMMmass_trmcmumatch;
	TH1F *hHPPmass_trmcmumatchisol, *hHMMmass_trmcmumatchisol, *hHPPHMMmass_trmcmumatchisol;
	TH1F *hHPPmass_trmcmumatchisol2, *hHMMmass_trmcmumatchisol2, *hHPPHMMmass_trmcmumatchisol2;
	TH1F *hHPPmass_trmcmumatchisol2j, *hHMMmass_trmcmumatchisol2j, *hHPPHMMmass_trmcmumatchisol2j;
	TH1F *hHPPmass_trmcmumatchisol24mu, *hHMMmass_trmcmumatchisol24mu;
	TH1F *hHPPmass_trmcmumatchisol2dphi, *hHMMmass_trmcmumatchisol2dphi, *hHPPHMMmass_trmcmumatchisol2dphi;
	TH1F *hHPPmass_trmcmumatchisol2jdz, *hHMMmass_trmcmumatchisol2jdz;
	TH1F *hHPPmass_trmcmumatchisol2jdz_ex, *hHMMmass_trmcmumatchisol2jdz_ex;

	TH2F *htrgmupt_ptetacut, *htrgmupt_ptetaisol2cut, *htrgmupt_ptetaisol2jcut;
	TH2F *htrstapt_ptetacut, *htrstapt_ptetaisol2notjcut;

	// ----------global muon-------------------------------------
	TH1F *hngmu, *hngmu_ptetacut, *hngmu_ptetaisolcut, *hngmu_ptetaisol2cut, *hngmu_ptetaisol2jcut;
	TH1F *hgmupt, *hgmupt_ptetacut, *hgmupt_ptetaisolcut, *hgmupt_ptetaisol2cut, *hgmupt_ptetaisol2jcut;
	TH1F *hgmupt_mcmumatch, *hgmupt_mcmumatchisol, *hgmupt_mcmumatchisol2, *hgmupt_mcmumatchisol2j;
	TH1F *hgmueta, *hgmueta_ptetacut, *hgmueta_ptetaisolcut, *hgmueta_ptetaisol2cut, *hgmueta_ptetaisol2jcut;
	TH1F *hgmueta_mcmumatch, *hgmueta_mcmumatchisol, *hgmueta_mcmumatchisol2, *hgmueta_mcmumatchisol2j;
	TH1F *hgmuisolPt_ptetacut01, *hgmuisolPt_ptetacut1, *hgmuisolPt_ptetacut10, *hgmuisolPt_ptetacut100, *hgmuisolPt_ptetacut1000;
	TH1F *hgmuisolPt2_ptetacut10, *hgmuisolPt2_ptetacut100;
	TH1F *h2gmudeltaphi, *h2gmudeltaphi_isol, *h2gmudeltaphi_isol2;
	TH2F *h2gmu2mcmudphi_isol2;

	TH1F *hgmuz, *hgmuz_ptetacut, *hgmuz_ptetaisol2cut, *hgmuz_ptetaisol2jcut;
	TH1F *hgmudz, *hgmudz_ptetaisol2jcut;
	TH1F *h2gmudz_isol2j, *h2gmudz_isol2j_ex;

	TH1F *hHPPmass, *hHMMmass, *hHPPHMMmass;
	TH1F *hHPPmass_isol, *hHMMmass_isol, *hHPPHMMmass_isol;
	TH1F *hHPPmass_isol2, *hHMMmass_isol2, *hHPPHMMmass_isol2;
	TH1F *hHPPmass_isol2j, *hHMMmass_isol2j, *hHPPHMMmass_isol2j;
	TH1F *hHPPmass_isol2dphi, *hHMMmass_isol2dphi, *hHPPHMMmass_isol2dphi;
	TH1F *hHPPmass_isol2jdz, *hHMMmass_isol2jdz;
	TH1F *hHPPmass_isol2jdz_ex, *hHMMmass_isol2jdz_ex;

	TH1F *hHPPmass_mcmumatch, *hHMMmass_mcmumatch, *hHPPHMMmass_mcmumatch;
	TH1F *hHPPmass_mcmumatchisol, *hHMMmass_mcmumatchisol, *hHPPHMMmass_mcmumatchisol;
	TH1F *hHPPmass_mcmumatchisol2, *hHMMmass_mcmumatchisol2, *hHPPHMMmass_mcmumatchisol2;
	TH1F *hHPPmass_mcmumatchisol2j, *hHMMmass_mcmumatchisol2j, *hHPPHMMmass_mcmumatchisol2j;
	TH1F *hHPPmass_mcmumatchisol24mu, *hHMMmass_mcmumatchisol24mu;
	TH1F *hHPPmass_mcmumatchisol2dphi, *hHMMmass_mcmumatchisol2dphi, *hHPPHMMmass_mcmumatchisol2dphi;
	TH1F *hHPPmass_mcmumatchisol2jdz, *hHMMmass_mcmumatchisol2jdz;
	TH1F *hHPPmass_mcmumatchisol2jdz_ex, *hHMMmass_mcmumatchisol2jdz_ex;

	TH2F *hgmustapt_ptetacut, *hgmustapt_ptetaisol2notjcut;
	// ----------MC-------------------------------------
	TH1F *hnmcmu_ptetacut, *hnmcmu_ptetaisolcut, *hnmcmu_ptetaisol2cut, *hnmcmu_ptetaisol2jcut;
	TH1F *hmcmupt, *hmcmupt_ptetacut, *hmcmupt_ptetaisolcut, *hmcmupt_ptetaisol2cut, *hmcmupt_ptetaisol2jcut;
	TH1F *hmcmupt_gmumatch, *hmcmupt_gmumatchisol, *hmcmupt_gmumatchisol2, *hmcmupt_trmatch_gmuunmatch;
	TH1F *hmcmueta, *hmcmueta_ptetacut, *hmcmueta_ptetaisolcut, *hmcmueta_ptetaisol2cut, *hmcmueta_ptetaisol2jcut;
	TH1F *hmcmueta_gmumatch, *hmcmueta_gmumatchisol, *hmcmueta_gmumatchisol2, *hmcmueta_trmatch_gmuunmatch;
	TH1F *hmcisolPt_ptetacut01, *hmcisolPt_ptetacut1, *hmcisolPt_ptetacut10, *hmcisolPt_ptetacut100, *hmcisolPt_ptetacut1000;
	TH1F *hmcisolPt2_ptetacut10, *hmcisolPt2_ptetacut100;
	TH1F *hmcmuisolPt_ptetacut01, *hmcmuisolPt_ptetacut1, *hmcmuisolPt_ptetacut10, *hmcmuisolPt_ptetacut100, *hmcmuisolPt_ptetacut1000;
	TH1F *hmcmuisolPt2_ptetacut10, *hmcmuisolPt2_ptetacut100;
	TH1F *hmcmuisolPt_mcmudR03, *hmcmuisolPt2_mcmudR03;
	TH1F *hmcmudeltaphi, *hmcmudeltaphi2, *hmcHdeltaphi, *h2mcmudeltaphi, *h2mcmudeltaphi_isol, *h2mcmudeltaphi_isol2;
	TH1F *hmcmudR, *hmcmudR_ex, *hmcmudR_allmu;
	TH1F *hmcmudR_upperisolPt_10, *hmcmudR_upperisolPt_20, *hmcmudR_upperisolPt2_10, *hmcmudR_upperisolPt2_20;
	TH2F *hmcmuxy, *hmcmuxy_upper1k, *hmcmuxy_10times1k;

	TH1F *hHPPmass_mc, *hHMMmass_mc, *hHPPHMMmass_mc;
	TH1F *hHPPmass_mcisol, *hHMMmass_mcisol, *hHPPHMMmass_mcisol;
	TH1F *hHPPmass_mcisol2, *hHMMmass_mcisol2, *hHPPHMMmass_mcisol2;
	TH1F *hHPPmass_mcisol2j, *hHMMmass_mcisol2j, *hHPPHMMmass_mcisol2j;
	TH1F *hHPPmass_mcisol2dphi, *hHMMmass_mcisol2dphi, *hHPPHMMmass_mcisol2dphi;

	TH2F *hgmumcmu_ptratio, *hgmumcmu_ptratioisol, *hgmumcmu_ptratioisol2, *hgmumcmu_ptratioisol2j;
	TH2F *hgmumcmu_ptratio_ex, *hgmumcmu_ptratioisol_ex, *hgmumcmu_ptratioisol2_ex, *hgmumcmu_ptratioisol2j_ex;
	TH2F *hgmumcmu_ptratio_2nd, *hgmumcmu_ptratio_2nd_ex, *hgmumcmu_ptratio_2ndtrue, *hgmumcmu_ptratio_2ndtrue_ex;

	// ----------MET-------------------------------------
	TH1F *hcalmet, *hcalmet0gmu, *hcalmet1gmu, *hcalmet2gmu, *hcalmet3gmu, *hcalmet4gmu;
	TH1F *hcalmet01mu;
	TH1F *hgenmet, *hgenmet0gmu, *hgenmet1gmu, *hgenmet2gmu, *hgenmet3gmu, *hgenmet4gmu;
	TH1F *hgenmetnonu, *hgenmet01mu;
	TH2F *hcalgenmet, *hcalgenmet0gmu, *hcalgenmet1gmu, *hcalgenmet2gmu, *hcalgenmet3gmu, *hcalgenmet4gmu;
	TH2F *hcalgenmetnonu, *hcalgenmet01mu, *hcal01mugenmet;

	// ----------ZZ->4mu-------------------------------------
	TH1F *hPMdphi_trisol2j;
	TH1F *hPMmassall_trisol2j, *hPMmassbig_trisol2j, *hPMmassZ_trisol2j, *hPMmassdphi_trisol2j;

	// ----------Vertex-------------------------------------
	TH1F *hnvt, *hvtx, *hvty, *hvtz, *hvtz_BS;
	TH1F *hnvt2, *hvt2x, *hvt2y, *hvt2z, *hvt2z_BS;
	TH1F *hnvt3, *hvt3x, *hvt3y, *hvt3z, *hvt3z_BS;
	TH2F *hvtxy, *hvt2xy, *hvt3xy;
	TH1F *hvtntr, *hvt2ntr, *hvt3ntr;
	TH1F *hmcmuz;
	TH1F *hmcmux_BS, *hmcmuy_BS, *hmcmuz_BS;
	TH1F *hmcmu2x_BS, *hmcmu2y_BS, *hmcmu2z_BS;
	TH1F *hmcmu3x_BS, *hmcmu3y_BS, *hmcmu3z_BS;
	TH2F *hmcmuxy_BS, *hmcmu2xy_BS, *hmcmu3xy_BS;
	TH1F *htrvtdxy, *htrvtdxy_isol2j, *htrvtdxy_cor, *htrvtdxy_isol2j_cor, *htrvtdxy_isol2j_in;
	TH1F *htrvtdxy_isol2jmc, *htrvtdxy_isol2jmc_in;
	TH1F *htrvtdxy_isol2j4mu, *htrvtdxy_isol2j4mu_in, *htrvtdxy_isol2j4mumc, *htrvtdxy_isol2j4mumc_in;
	TH1F *htrvtdz, *htrvtdz_isol2j, *htrvtdz_isol2jmc, *htrvtdz2, *htrvtdz2_isol2j, *htrvtdz3, *htrvtdz3_isol2j;
	TH1F *htrvtd, *htrtrd, *htripvtd, *htrtripd, *htrtrsd;
	TH2F *hntr_isol2j_VS_vt;

	TH1F *htrvtd_i2j, *htrvtdz_i2j, *htrvtdxy_i2j, *htrvtd_i2jmc, *htrvtdz_i2jmc, *htrvtdxy_i2jmc;
	TH1F *htripvtd_i2j, *htripvtdz_i2j, *htripvtdxy_i2j, *htripvtd_i2jmc, *htripvtdz_i2jmc, *htripvtdxy_i2jmc;
	TH1F *htrtrd_i2j, *htrtrdz_i2j, *htrtrdxy_i2j, *htrtrd_i2jmc, *htrtrdz_i2jmc, *htrtrdxy_i2jmc;
	TH1F *htrtripd_i2j, *htrtripdz_i2j, *htrtripdxy_i2j, *htrtripd_i2jmc, *htrtripdz_i2jmc, *htrtripdxy_i2jmc;
	TH1F *htrtrsd_i2j, *htrtrsdz_i2j, *htrtrsdxy_i2j, *htrtrsd_i2jmc, *htrtrsdz_i2jmc, *htrtrsdxy_i2jmc;

	TH1F *htStSd_i2j4mu, *htStSdz_i2j4mu, *htStSdxy_i2j4mu;
	TH1F *htPtMd_i2j4mu, *htPtMdz_i2j4mu, *htPtMdxy_i2j4mu;

	TH1F *hnmcvt;
	TH2F *hmcvtxy10;
	//	TH2F *htrvtxy, *htrvtxy_isol2j;
	TH1F *hmcmurecovtd, *hmcmurecovt2d, *hmcmurecovt3d, *hmcmurecovt2d_ex, *hmcmurecovt3d_ex;
	TH1F *hmcmurecovtcd, *hmcmurecovtc2d, *hmcmurecovtc3d;
	TH2F *hmcmuxvtx, *hmcmuyvty, *hmcmuzvtz;
	TH2F *hmcmuxvtx_BS, *hmcmuyvty_BS, *hmcmuzvtz_BS;
	TH2F *hmcmuxcpx, *hmcmuycpy, *hmcmuzcpz;
	TH2F *hvtxcpx_BS, *hvtycpy_BS, *hvtzcpz_BS;
	//	TH2F *hmcmurecovtxy, *hmcmurecovt2xy, *hmcmurecovt3xy;

	TH1F *hvtchi2_BS, *hvtchi2_BS3mu, *hvtchi2_BS4mu, *hvt2chi2_BS, *hvt3chi2_BS;

	TH1F *hevents;
};

#endif

/* vim:set ts=2 sts=2 sw=2: */
