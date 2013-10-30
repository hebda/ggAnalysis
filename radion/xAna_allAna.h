#ifndef xAnaSummer11_excl_h
#define xAnaSummer11_excl_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <cmath>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include "WeightManager.h"
#include "ConfigReader.h"
#include "EnergyCorrections.h"
#include "MinitreeOut.h"


#include "VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "VertexAnalysis/interface/HggVertexFromConversions.h"
#include "VertexAnalysis/interface/PhotonInfo.h"
#include "VertexAnalysis/interface/ReweightMC.h"
#include "VertexAnalysis/interface/VertexAlgoParameters.h"

#include "JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETObjects/interface/JetResolution.h"
#include "JetMETObjects/interface/Linkdef.h"
#include "JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "JetMETObjects/interface/SimpleJetCorrector.h"

#include <fstream>

const float vcicSTmva[17][4] = {
  { 3.8,     2.2,     1.77,     1.29},    // baseline rel combIso (good vtx)
  {11.7,     3.4,     3.90,     1.84},    // baseline rel combIso (bad vtx)
  { 3.5,     2.2,     2.30,     1.45},    // presel trkIso (good vtx)
  { 0.014,    0.014,  0.034,    0.034},   // presel sigma_ie_ie 
  { 0.082,   0.075,   0.075,    0.075},   // presel H/E       
  { 0.94,    0.36,    0.94,     0.32},    // baseline R9
  { 1,       0.062,   0.97,     0.97},    // baseline dR to trk
  { 1.5,     1.5,     1.5,      1.5},     // baseline pixel     
  { 50,      4,       50,       4},       // presel EtCorrEcalIso
  { 50,      4,       50,       4},       // presel EtCorrHcalIso
  { 50,      4,       50,       4},       // presel EtCorrTrkIso
  { 3,       3,       3,        3},       // presel PUCorrHcalEcal
  { 2.8,     2.8,     2.8,      2.8},     // presel AbsTrkIsoCIC
  { 4,       4,       4,        4},       // presel HollowConeTrkIsoDr03
  { 0.0106,  0.0097,  0.028,    0.027},   // baseline sigma_ie_ie
  { 0.082,   0.062,   0.065,    0.048},   // baseline H/E     
  { 4,       4,       4,        4}        // presel psChargedIso02
};



// ============= CiC4PF ==================

const float vcicST[7][4] = {
  { 6,     4.7,     5.6,     3.6},    //                             
  {10,     6.5,     5.6,     4.4},    //                                                 
  { 3.8,     2.5,    3.1,   2.2},    // Charged PF iso                                                          
  { 0.0108,  0.0102,  0.028,    0.028},   // sigma_ie_ie                                                                  
  { 0.124,   0.092,   0.142,    0.063},   // H/E       
  { 0.94,    0.28,    0.94,     0.24},    // R9
  //  { 1,       0.062,   0.97,     0.97},    // dR to trk
  { 1,       1,     1,     1},    // dR to trk
};



//// delta_Sigma/E
//// May10 
//float Eres[4] = {0.89e-2, 1.99e-2, 4.09e-2, 2.46e-2};
// // 600 /pb July1 email from Paramatti
//float Eres[4] = {1.49e-2, 1.77e-2, 4.65e-2, 3.72e-2};
////  Wieghted Average for May10 and PR
//float Eres[4] = {1.37e-2, 1.81e-2, 4.53e-2, 3.47e-2};

/*
//  May10 and PromptReco July6
//
// ////////////////////////////
// Energy corrected and fitted over (May10 + PromptReco)
float Eres[4] = {1.41e-2, 2.00e-2, 4.74e-2, 3.61e-2};

//// (E_data - E_MC)/E
// 160404 - 163869 May10 ReReco
float Eoffset_0[4] = { 0.47e-2, -0.25e-2, -0.58e-2,  0.10e-2};

// 165071 - 165970 (139.6 pb)  PR
float Eoffset_1[4] = { 0.07e-2, -0.49e-2, -2.49e-2, -0.62e-2};

// 165971 - 166502 (187.9 pb) PR
float Eoffset_2[4] = {-0.03e-2, -0.67e-2, -3.76e-2, -1.33e-2};

// 166503 - 166861 (191.6 pb) PR
float Eoffset_3[4] = {-0.11e-2, -0.63e-2, -4.50e-2, -1.78e-2};

// 166862 - 167784 (252.8 pb) PR
float Eoffset_4[4] = {-0.14e-2, -0.74e-2, -5.61e-2, -2.73e-2};
*/

//  May10 and Re Reco  by Paolo Meridiani
//
// //////////////////////////////////////

/*
float Eres[4] = {0.91e-2, 1.65e-2, 3.23e-2, 2.47e-2};

float Eoffset_0[4] = { 0.21e-2, -0.59e-2,  0.02e-2, -0.15e-2};
float Eoffset_1[4] = { 0.06e-2, -0.55e-2, -0.50e-2, -0.38e-2};
float Eoffset_2[4] = { 0.05e-2, -0.61e-2, -0.34e-2, -0.71e-2};
float Eoffset_3[4] = { 0.04e-2, -0.51e-2,  0.25e-2, -0.38e-2};
float Eoffset_4[4] = {-0.01e-2, -0.65e-2,  0.39e-2, -0.31e-2};
*/


//  

//==== LP ======================
//float Eres[4] = {1.0e-2, 1.69e-2, 3.03e-2, 2.96e-2};
//float Eoffset[4] = { 0.02e-2, -0.03e-2,  0.0e-2, 0.02e-2};
//==============================





// new corrections after the energy regression 
// EB data 
float YMcorr[10] = {1 - 0.830/91.19, // 0
		    1 - 0.550/91.19, // 1
		    1 - 0.495/91.19, // 2
		    1 - 0.560/91.19, // 3
		    1 - 0.463/91.19, // 4
		    1 - 0.457/91.19, // 5
		    1 - 0.393/91.19, // 6
		    1 - 0.460/91.19, // 7
		    1 - 0.133/91.19, // 8
		    1 - 0.219/91.19}; // 9



//--------------------------------------

//If a SC is in 4th module: its energy should be bumped up by this factor: 1 + 1.14/91.19
// scale applied in .C 

float YMsmearExtra = 1;  // this is for basket 4  not used!!!


float YMcorrMC = 1 - 0.055/91.19;
float YMsmear = 0.68 / 91.19;


//--------------------------------------

//EE data  --> no scale corrections needed
float YMEEcorr = 1;

//EE MC ???
float YMEEcorrMC = 1;  

//float YMEEcorrMC = 1 - 1.62/91.19; // mc scale is imperfect

// extra smearing for |eta|>2  (YMEEsmear[1]) to be applied on top of YMEEsmear[0]
float YMEEsmear[2] = {1.6 / 91.19,
                      0.5 / 91.19};



/*
//==== LP Systematics ======================


float sysEres[4]={0.23e-2,0.42e-2,0.51e-2,0.38e-2};
float Eres[4] = {1.0e-2, 1.69e-2, 3.03e-2, 2.96e-2};

float EresPlus[4]={0,0,0,0};
float EresMinus[4]={0,0,0,0};


//======= LP ==================
float Eoffset[4] = { 0.02e-2, 0.03e-2,  0.0e-2, 0.02e-2};

float sysEoffset[4] = { 0.09e-2, 0.44e-2,  0.29e-2, 0.41e-2};


//float Eoffset_jul5[4] = {0.094, -0.06, 1.641, 0.150};
//float Eoffset_aug5[4] = {-0.0846, -0.1263, 2.57, 0.974};
//float Eoffset_pr[4] = {-0.225, -0.266, 3.54, 2.03};

float Eoffset_jul5[4] = {0.03e-2, -0.54e-2, 0.03e-2, -0.24e-2};
float Eoffset_aug5[4] = {-0.2e-2, -0.81e-2, 1.09e-2, 0.43e-2};
float Eoffset_pr[4] = {-0.35e-2, -0.97e-2, 2.15e-2, 1.58e-2};

//==============================



//==== EPS =====================
//float Eres[4] = {0.93e-2, 1.72e-2, 3.15e-2, 2.69e-2};
//float Eoffset[4] = { 0.04e-2, -0.52e-2,  0.22e-2, -0.13e-2};
//==============================


float Eoffset_0[4] = { 0.21e-2, -0.59e-2,  0.02e-2, -0.15e-2};
float Eoffset_1[4] = { 0.06e-2, -0.55e-2, -0.50e-2, -0.38e-2};
float Eoffset_2[4] = { 0.05e-2, -0.61e-2, -0.34e-2, -0.71e-2};
float Eoffset_3[4] = { 0.04e-2, -0.51e-2,  0.25e-2, -0.38e-2};
float Eoffset_4[4] = {-0.01e-2, -0.65e-2,  0.39e-2, -0.31e-2};

*/

/*
*/

const Int_t kMaxeleESEffSigmaRR = 1;
const Int_t kMaxeleESE1 = 1;
const Int_t kMaxeleESE3 = 1;
const Int_t kMaxeleESE5 = 1;
const Int_t kMaxeleESE7 = 1;
const Int_t kMaxeleESE11 = 1;
const Int_t kMaxeleESE21 = 1;
const Int_t kMaxphoESEffSigmaRR = 1;
const Int_t kMaxphoESE1 = 1;
const Int_t kMaxphoESE3 = 1;
const Int_t kMaxphoESE5 = 1;
const Int_t kMaxphoESE7 = 1;
const Int_t kMaxphoESE11 = 1;
const Int_t kMaxphoESE21 = 1;
const Int_t kMaxnPFPho = 1000;
const Int_t kMaxisIsoPFPho = 1;
const Int_t kMaxPFphoE = 1;
const Int_t kMaxPFphoEt = 1;
const Int_t kMaxPFphoEta = 1;
const Int_t kMaxPFphoPhi= 1;
const Int_t kMaxPFphoIsConv = 1;
const Int_t kMaxPFphoTrkIsoHollowDR04 = 1;
const Int_t kMaxPFphoEcalIsoDR04 = 1;
const Int_t kMaxPFphoHcalIsoDR04 = 1;
const Int_t kMaxPFphoHoverE = 1;
const Int_t kMaxnPFElec = 50;
const Int_t kMaxPFelePt = 1;
const Int_t kMaxPFeleEta = 1;
const Int_t kMaxPFelePhi = 1;
const Int_t kMaxPFeleEn = 1;
const Int_t kMaxPFeleCharge = 1;
const Int_t kMaxnPFMu = 50;
const Int_t kMaxPFmuCharge = 1;
const Int_t kMaxPFmuPhi = 1;
const Int_t kMaxPFmuEta = 1;
const Int_t kMaxPFmuPt = 1;
const Int_t kMaxPFmuPz = 1;
const Int_t kMaxPFmuChambers = 1;
const Int_t kMaxPFmuD0 = 1;
const Int_t kMaxPFmuDz = 1;
const Int_t kMaxPFmuChi2NDF = 1;
const Int_t kMaxPFmuNumberOfValidTrkHits = 1;
const Int_t kMaxPFmuNumberOfValidPixelHits = 1;
const Int_t kMaxPFmuNumberOfValidMuonHits = 1;
const Int_t kMaxnPFchad = 1;
const Int_t kMaxPFhad_charge = 1;
const Int_t kMaxPFchad_Eta = 1;
const Int_t kMaxPFchad_Phi = 1;
const Int_t kMaxPFchad_Pt = 1;
const Int_t kMaxPFchad_P = 1;
const Int_t kMaxPFchad_E = 1;
const Int_t kMaxPFchad_Ecal = 1;
const Int_t kMaxPFchad_Hcal = 1;
const Int_t kMaxnPFnhad = 350;
const Int_t kMaxPFnhad_Eta = 1;
const Int_t kMaxPFnhad_Phi = 1;
const Int_t kMaxPFnhad_Pt = 1;
const Int_t kMaxPFnhad_P = 1;
const Int_t kMaxPFnhad_E = 1;
const Int_t kMaxPFnhad_Hcal = 1;
const Int_t maxP = 500;
const Int_t maxTrks = 2000;

class xAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   
   // TH1F *hMCpu;
   
   float phoID_HoE,phoID_covIEtaIEta,
     phoID_tIso1abs,phoID_tIso2abs,phoID_tIso3abs,
     phoID_R9,phoID_absIsoEcal,
     phoID_absIsoHcal,phoID_NVertexes,phoID_ScEta,phoID_EtaWidth,phoID_PhiWidth;
   float phoID_covIEtaIPhi, phoID_S4, phoID_rho, 
     phoID_ESEffSigmaRR, phoID_pfPhoIso03, phoID_pfChIso03, phoID_pfChIso03worst ;

   float Ddipho_masserr,Ddipho_masserrwrongvtx,Ddipho_vtxprob;
   float Ddipho_piT1,Ddipho_piT2,Ddipho_eta1,Ddipho_eta2,Ddipho_cosdPhi,Ddipho_id1,Ddipho_id2;
   
   float myVBFLeadJPt,myVBFSubJPt,myVBFdEta,myVBF_Mjj,myVBFZep,myVBF_Mgg;
   float myVBFdPhi,myVBFDiPhoPtOverM,myVBFLeadPhoPtOverM,myVBFSubPhoPtOverM;

   TMVA::Reader *phoID_2011[2];
   TMVA::Reader *phoID_2012[2];
   TMVA::Reader *phoID_mva[2]; 

   TMVA::Reader *DiscriDiPho_2011;
   TMVA::Reader *DiscriDiPho_2012;
   TMVA::Reader *Ddipho_mva;

   TMVA::Reader *DiscriVBF;
   bool DiscriVBF_UseDiPhoPt;
   bool DiscriVBF_UsePhoPt;
   bool DiscriVBF_useMvaSel;
   string DiscriVBF_Method;
   vector<float> DiscriVBF_cat;
   vector<float> DiscriVBF_sig_cat   ; /// not used for now
   vector<float> DiscriVBF_dipho_cat ; /// not used for now
   vector<float> DiscriVBF_gluglu_cat; /// not used for now

   int doJetRegression;//0 for no regression, 1 for regression with 2 BDTs
   int doControlSample;//0 for normal analysis, 1 for control sample (invert 1 photon ID)
   TMVA::Reader *jetRegres;

   float sumTrackPtInCone(int,int,float,float,float,float,float,float);
   float worstSumTrackPtInCone(int,int&,float,float,float,float,float,float);

   inline void SetConfigReader( const ConfigReader &config );
   inline void SetConfigReader( string configName );

   // Declaration of leaf types
   Int_t           run;
   Long64_t         event;
   Int_t           orbit;
   Int_t           bx;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nHLT;
   Int_t           HLT[2309];   //[nHLT]
   Int_t           HLTIndex[70];
   Float_t         bspotPos[3];
   Int_t           nVtx;
   Float_t         vtx[230][3];   //[nVtx]
   Int_t           vtxNTrk[230];   //[nVtx]
   Int_t           vtxNDF[230];   //[nVtx]
   Float_t         vtxD0[230];   //[nVtx]
   Int_t           IsVtxGood;
   Int_t           nVtxBS;
   float vertProb;
   float vertMVA;
   float vert_ptbal;
   float vert_ptasym;
   float vert_logsumpt2;
   float vert_limPullToConv;
   float vert_nConv;
   Float_t         vtxbs[230][3];   //[nVtxBS]
   Float_t         vtxbsPtMod[230];   //[nVtxBS]
   Float_t         vtxbsSumPt2[230];   //[nVtxBS]
   //vector<vector<float> > *vtxbsTkPxVec;
   //vector<vector<float> > *vtxbsTkPyVec;
   std::vector<std::vector<int> > *vtxbsTkIndex;
   std::vector<std::vector<float> > *vtxbsTkWeight;
   //int             vtxbsTkIndex[10000][1000]
   //float           vtxbsTkWeight[10000][1000]
   Int_t           nTrk;
   Float_t         trkP[10660][3];   //[nTrk]
   Float_t         trkVtx[10660][3];   //[nTrk]
   Float_t         trkd0[10660];   //[nTrk]
   Float_t         trkd0Err[10660];   //[nTrk]
   Float_t         trkdz[10660];   //[nTrk]
   Float_t         trkdzErr[10660];   //[nTrk]
   Float_t         trkPtErr[10660];   //[nTrk]
   Int_t           trkQuality[10660];   //[nTrk]
   Int_t           nGoodTrk;
   Int_t           IsTracksGood;
   Float_t         pdf[7];
   Float_t         pthat;
   Float_t         processID;
   Int_t           nMC;
   Int_t           mcPID[230];   //[nMC]
   Float_t         mcVtx[230][3];   //[nMC]
   Float_t         mcPt[230];   //[nMC]
   Float_t         mcMass[230];   //[nMC]
   Float_t         mcEta[230];   //[nMC]
   Float_t         mcPhi[230];   //[nMC]
   Float_t         mcE[230];   //[nMC]
   Float_t         mcEt[230];   //[nMC]
   Int_t           mcGMomPID[230];   //[nMC]
   Int_t           mcMomPID[230];   //[nMC]
   Float_t         mcMomPt[230];   //[nMC]
   Float_t         mcMomMass[230];   //[nMC]
   Float_t         mcMomEta[230];   //[nMC]
   Float_t         mcMomPhi[230];   //[nMC]
   Int_t           mcIndex[230];   //[nMC]
   Int_t           mcDecayType[230];   //[nMC]
   Float_t         mcCalIsoDR03[230];   //[nMC]
   Float_t         mcTrkIsoDR03[230];   //[nMC]
   Float_t         mcCalIsoDR04[230];   //[nMC]
   Float_t         mcTrkIsoDR04[230];   //[nMC]
   Float_t         genMET;
   Float_t         genMETx;
   Float_t         genMETy;
   Float_t         genMETPhi;
   Int_t           nPUInfo;
   Int_t           nPU[3];   //[nPUInfo]
   Int_t           puBX[3];   //[nPUInfo]
   Int_t           puTrue[3];    //[nPUInfo]
   Float_t         MET;
   Float_t         METx;
   Float_t         METy;
   Float_t         METPhi;
   Float_t         METsumEt;
   /* Float_t         tcMET; */
   /* Float_t         tcMETx; */
   /* Float_t         tcMETy; */
   /* Float_t         tcMETPhi; */
   /* Float_t         tcMETsumEt; */
   /* Float_t         tcMETmEtSig; */
   /* Float_t         tcMETSig; */
   Float_t         pfMET;
   Float_t         pfMETx;
   Float_t         pfMETy;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Float_t         recoPfMET;
   Float_t         recoPfMETx;
   Float_t         recoPfMETy;
   Float_t         recoPfMETPhi;
   Float_t         recoPfMETsumEt;
   Float_t         recoPfMETmEtSig;
   Float_t         recoPfMETSig;
   Int_t           nEle;
   Int_t           eleTrg[1000][13];   //[nEle]
   Int_t           eleID[1000][12];   //[nEle]
   Int_t           eleClass[1000];   //[nEle]
   Int_t           eleCharge[1000];   //[nEle]
   Float_t         eleEn[1000];   //[nEle]
   Float_t         eleSCRawEn[1000];   //[nEle]
   Float_t         eleSCEn[1000];   //[nEle]
   Float_t         elePt[1000];   //[nEle]
   Float_t         eleEta[1000];   //[nEle]
   Float_t         eleEtaVtx[1000][100];   //[nEle]
   Float_t         elePhi[1000];   //[nEle]
   Float_t         elePhiVtx[1000][100];   //[nEle]
   Float_t         eleSCEta[1000];   //[nEle]
   Float_t         eleSCPhi[1000];   //[nEle]
   Float_t         eleSCEtaWidth[1000];   //[nEle]
   Float_t         eleSCPhiWidth[1000];   //[nEle]
   Float_t         eleVtx[1000][3];   //[nEle]
   Float_t         eleHoverE[1000];   //[nEle]
   Float_t         eleEoverP[1000];   //[nEle]
   Float_t         elePin[1000];   //[nEle]
   Float_t         elePout[1000];   //[nEle]
   Float_t         eleBrem[1000];   //[nEle]
   Float_t         eledEtaAtVtx[1000];   //[nEle]
   Float_t         eledPhiAtVtx[1000];   //[nEle]
   Float_t         eleSigmaIEtaIEta[1000];   //[nEle]
   Float_t         eleSigmaIEtaIPhi[1000];   //[nEle]
   Float_t         eleSigmaIPhiIPhi[1000];   //[nEle]
   Float_t         eleD0[1000];   //[nEle]
   Float_t         eleDz[1000];   //[nEle]
   Float_t         eleD0Vtx[1000][100];   //[nEle]
   Float_t         eleDzVtx[1000][100];   //[nEle]
   Float_t         eleEmax[1000];   //[nEle]
   Float_t         eleE3x3[1000];   //[nEle]
   Float_t         eleE5x5[1000];   //[nEle]
   Float_t         eleE2x5Right[1000];   //[nEle]
   Float_t         eleE2x5Left[1000];   //[nEle]
   Float_t         eleE2x5Top[1000];   //[nEle]
   Float_t         eleE2x5Bottom[1000];   //[nEle]
   Float_t         eleSeedTime[1000];   //[nEle]
   Int_t           eleRecoFlag[1000];   //[nEle]
   Int_t           eleGenIndex[1000];   //[nEle]
   Int_t           eleGenGMomPID[1000];   //[nEle]
   Int_t           eleGenMomPID[1000];   //[nEle]
   Float_t         eleGenMomPt[1000];   //[nEle]
   Float_t         eleIsoTrkDR03[1000];   //[nEle]
   Float_t         eleIsoEcalDR03[1000];   //[nEle]
   Float_t         eleIsoHcalDR03[1000];   //[nEle]
   Float_t         eleIsoHcalSolidDR03[1000];   //[nEle]
   Float_t         eleIsoTrkDR04[1000];   //[nEle]
   Float_t         eleIsoEcalDR04[1000];   //[nEle]
   Float_t         eleIsoHcalDR04[1000];   //[nEle]
   Float_t         eleIsoHcalSolidDR04[1000];   //[nEle]
   Int_t           eleMissHits[1000];   //[nEle]
   Float_t         eleConvDist[1000];   //[nEle]
   Float_t         eleConvDcot[1000];   //[nEle]
   Float_t         elePFChIso03[1000];
   Float_t         elePFPhoIso03[1000];
   Float_t         elePFNeuIso03[1000];
   Float_t         eleIDMVANonTrig[1000];
   Float_t         eleEcalEn[1000];
   Int_t           eleConvVtxFit[1000];
   Float_t         eleESEffSigmaRR[1000];   //[nEle]
   Float_t         eleESE1[1000][2];   //[nEle]
   Float_t         eleESE3[1000][2];   //[nEle]
   Float_t         eleESE5[1000][2];   //[nEle]
   Float_t         eleESE7[1000][2];   //[nEle]
   Float_t         eleESE11[1000][2];   //[nEle]
   Float_t         eleESE21[1000][2];   //[nEle]
   Int_t           nPho;
   Int_t           phoTrg[1000][8];   //[nPho]
   Bool_t          phoIsPhoton[1000];   //[nPho]
   Float_t         phoSCPos[1000][3];   //[nPho]
   Float_t         phoCaloPos[1000][3];   //[nPho]
   Float_t         phoE[1000];   //[nPho]
   Float_t         phoRegrE[1000];   //[nPho]
   Float_t         phoRegrErr[1000];   //[nPho]
   Float_t         phoRegrSmear[1000];   //[nPho]
   Float_t         phoEt[1000];   //[nPho]
   Float_t         phoEta[1000];   //[nPho]
   Float_t         phoVtx[1000][3];   //[nPho]
   Float_t         phoPhi[1000];   //[nPho]
   Float_t         phoEtVtx[1000][100];   //[nPho]
   Float_t         phoEtaVtx[1000][100];   //[nPho]
   Float_t         phoPhiVtx[1000][100];   //[nPho]
   Float_t         phoR9[1000];   //[nPho]
   Float_t         phoCetaCorrE[1000];   //[nPho]
   Float_t         phoCetaCorrEt[1000];   //[nPho]
   Float_t         phoBremCorrE[1000];   //[nPho]
   Float_t         phoBremCorrEt[1000];   //[nPho]
   Float_t         phoFullCorrE[1000];   //[nPho]
   Float_t         phoFullCorrEt[1000];   //[nPho]
   Float_t         phoTrkIsoSolidDR03[1000];   //[nPho]
   Float_t         phoTrkIsoHollowDR03[1000];   //[nPho]
   Float_t         phoEcalIsoDR03[1000];   //[nPho]
   Float_t         phoHcalIsoDR03[1000];   //[nPho]
   Float_t         phoHcalIsoSolidDR03[1000];   //[nPho]
   Float_t         phoTrkIsoSolidDR04[1000];   //[nPho]
   Float_t         phoTrkIsoHollowDR04[1000];   //[nPho]
   Float_t         phoTrkIsoHollowNoDzDR03[1000];   //[nPho]
   Float_t         phoTrkIsoHollowNoDzDR04[1000];   //[nPho]
   Float_t         phoCiCTrkIsoDR03[1000][100];   //[nPho]
   Float_t         phoCiCTrkIsoDR04[1000][100];   //[nPho]
   Float_t         phoCiCTrkIsoDR03Vtx[1000][100];   //[nPho]
   Float_t         phoCiCTrkIsoDR04Vtx[1000][100];   //[nPho]
   Float_t         phoCiCdRtoTrk[1000];   //[nPho]
   Float_t         phoEcalIsoDR04[1000];   //[nPho]
   Float_t         phoHcalIsoDR04[1000];   //[nPho]
   Float_t         phoHcalIsoSolidDR04[1000];   //[nPho]
   Float_t         phoHoverE[1000];   //[nPho]
   Float_t         phoSigmaIEtaIEta[1000];   //[nPho]
   Float_t         phoSigmaIEtaIPhi[1000];   //[nPho]
   Float_t         phoSigmaIPhiIPhi[1000];   //[nPho]
   Int_t           phoEleVeto[1000];   //[nPho]   ----  2012

   Float_t  phoCiCPF4phopfIso03[1000];
   Float_t  phoCiCPF4phopfIso04[1000];
   Float_t  phoCiCPF4chgpfIso02[1000][100];
   Float_t  phoCiCPF4chgpfIso03[1000][100];
   Float_t  phoCiCPF4chgpfIso04[1000][100];
   
   Float_t  phoCiCPF4chgpfIso03AppVeto[1000][100];
   Float_t  phoCiCPF4chgpfIso04AppVeto[1000][100];
   Float_t  phoCiCPF4phopfIso03AppVeto[1000];
   Float_t  phoCiCPF4phopfIso04AppVeto[1000];
   
   Float_t  phoCiCPF4phopfIso03AppVetoClean[1000];
   Float_t  phoCiCPF4phopfIso04AppVetoClean[1000];
   
   Float_t  phoCiCPF4chgpfIso03Mod[1000][100];
   Float_t  phoCiCPF4chgpfIso04Mod[1000][100];

   Float_t         phoEmax[1000];   //[nPho]
   Float_t         phoS4ratio[1000];   //[nPho]
   Float_t         phoE2x2[1000];   //[nPho]
   Float_t         phoE3x3[1000];   //[nPho]
   Float_t         phoE5x5[1000];   //[nPho]
   Float_t         phoE2x5Right[1000];   //[nPho]
   Float_t         phoE2x5Left[1000];   //[nPho]
   Float_t         phoE2x5Top[1000];   //[nPho]
   Float_t         phoE2x5Bottom[1000];   //[nPho]
   Float_t         phoSeedTime[1000];   //[nPho]
   Int_t           phoSeedDetId1[1000];   //[nPho]
   Int_t           phoSeedDetId2[1000];   //[nPho]
   Int_t           phoRecoFlag[1000];   //[nPho]
   Int_t           phoPos[1000];   //[nPho]
   Int_t           phoGenIndex[1000];   //[nPho]
   Int_t           phoGenGMomPID[1000];   //[nPho]
   Int_t           phoGenMomPID[1000];   //[nPho]
   Float_t         phoGenMomPt[1000];   //[nPho]
   Float_t         phoSCE[1000];   //[nPho]
   Float_t         phoSCRawE[1000];   //[nPho]
   Float_t         phoSCEt[1000];   //[nPho]
   Float_t         phoSCEta[1000];   //[nPho]
   Float_t         phoSCPhi[1000];   //[nPho]
   Float_t         phoSCEtaWidth[1000];   //[nPho]
   Float_t         phoSCPhiWidth[1000];   //[nPho]
   Float_t         phoSCBrem[1000];   //[nPho]
   Int_t           phoOverlap[1000];   //[nPho]
   Int_t           phohasPixelSeed[1000];   //[nPho]
   Int_t           phoIsConv[1000];   //[nPho]
   Int_t           phoNConv[1000];   //[nPho]
  
   Float_t         phoConvInvMass[1000];   //[nPho]
   Float_t         phoConvCotTheta[1000];   //[nPho]
   Float_t         phoConvEoverP[1000];   //[nPho]
   Float_t         phoConvZofPVfromTrks[1000];   //[nPho]
   Float_t         phoConvMinDist[1000];   //[nPho]
   Float_t         phoConvdPhiAtVtx[1000];   //[nPho]
   Float_t         phoConvdPhiAtCalo[1000];   //[nPho]
   Float_t         phoConvdEtaAtCalo[1000];   //[nPho]
   Float_t         phoConvTrkd0[1000][2];   //[nPho]
   Float_t         phoConvTrkPin[1000][2];   //[nPho]
   Float_t         phoConvTrkPout[1000][2];   //[nPho]
   Float_t         phoConvTrkdz[1000][2];   //[nPho]
   Float_t         phoConvTrkdzErr[1000][2];   //[nPho]
   Float_t         phoConvChi2[1000];   //[nPho]
   Float_t         phoConvChi2Prob[1000];   //[nPho]
   Int_t           phoConvNTrks[1000];   //[nPho]
   Float_t         phoConvCharge[1000][2];   //[nPho]
   Float_t         phoConvValidVtx[1000];   //[nPho]
   Float_t         phoConvLikeLihood[1000];   //[nPho]
   Float_t         phoConvP4[1000][4];   //[nPho]
   Float_t         phoConvVtx[1000][3];   //[nPho]
   Float_t         phoConvVtxErr[1000][3];   //[nPho]
   Float_t         phoConvPairMomentum[1000][3];   //[nPho]
   Float_t         phoConvRefittedMomentum[1000][3];   //[nPho]
   Float_t         phoESEffSigmaRR[1000][3];   //[nPho]
   Float_t         phoESE1[1000][2];   //[nPho]
   Float_t         phoESE3[1000][2];   //[nPho]
   Float_t         phoESE5[1000][2];   //[nPho]
   Float_t         phoESE7[1000][2];   //[nPho]
   Float_t         phoESE11[1000][2];   //[nPho]
   Float_t         phoESE21[1000][2];   //[nPho]
   Int_t           nMu;
   Int_t           muTrg[230][6];   //[nMu]
   Float_t         muEta[230];   //[nMu]
   Float_t         muPhi[230];   //[nMu]
   Int_t           muCharge[230];   //[nMu]
   Float_t         muPt[230];   //[nMu]
   Float_t         muPz[230];   //[nMu]
   Int_t           muGenIndex[230];   //[nMu]
   Float_t         muIsoTrk[230];   //[nMu]
   Float_t         muIsoCalo[230];   //[nMu]
   Float_t         muIsoEcal[230];   //[nMu]
   Float_t         muIsoHcal[230];   //[nMu]
   Float_t         muChi2NDF[230];   //[nMu]
   Float_t         muEmVeto[230];   //[nMu]
   Float_t         muHadVeto[230];   //[nMu]
   Int_t           muType[230];   //[nMu]
   Bool_t          muID[230][6];   //[nMu]
   Float_t         muD0[230];   //[nMu]
   Float_t         muDz[230];   //[nMu]
   Float_t         muD0Vtx[230][100];   //[nMu]
   Float_t         muDzVtx[230][100];   //[nMu]
   Float_t         muVtx[230][3];     //[nMu]
   Int_t           muNumberOfValidTrkHits[230];   //[nMu]
   Int_t           muNumberOfValidPixelHits[230];   //[nMu]
   Int_t           muNumberOfValidMuonHits[230];   //[nMu]
   Int_t           muStations[230];   //[nMu]
   Int_t           muChambers[230];   //[nMu]

   Int_t muNumberOfValidTrkLayers[230];

   Float_t  muPFIsoR04_CH[230];
   Float_t  muPFIsoR04_NH[230];
   Float_t  muPFIsoR04_Pho[230];
   Float_t  muPFIsoR04_PU[230];
   Float_t  muPFIsoR04_CPart[230];
   Float_t  muPFIsoR04_NHHT[230];
   Float_t  muPFIsoR04_PhoHT[230];

   Float_t  muPFIsoR03_CH[230];
   Float_t  muPFIsoR03_NH[230];
   Float_t  muPFIsoR03_Pho[230];
   Float_t  muPFIsoR03_PU[230];
   Float_t  muPFIsoR03_CPart[230];
   Float_t  muPFIsoR03_NHHT[230];
   Float_t  muPFIsoR03_PhoHT[230];


   Int_t           nPFPho_;
   Int_t           isIsoPFPho_[kMaxnPFPho];   //[nPFPho_]
   Int_t           PFphoE_[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFphoEt_[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFphoEta_[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFphoPhi_[kMaxnPFPho];   //[nPFPho_]
   Int_t           PFphoIsConv_[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFphoTrkIsoHollowDR04_[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFphoEcalIsoDR04_[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFphoHcalIsoDR04_[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFphoHoverE_[kMaxnPFPho];   //[nPFPho_]
   Int_t           nPFElec_;
   Float_t         PFelePt_[kMaxnPFElec];   //[nPFElec_]
   Float_t         PFeleEta_[kMaxnPFElec];   //[nPFElec_]
   Float_t         PFelePhi_[kMaxnPFElec];   //[nPFElec_]
   Float_t         PFeleEn_[kMaxnPFElec];   //[nPFElec_]
   Int_t           PFeleCharge[kMaxnPFElec];   //[nPFElec_]
   Int_t           nPFMu_;
   Int_t           PFmuCharge_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuPhi_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuEta_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuPt_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuPz_[kMaxnPFMu];   //[nPFMu_]
   Int_t           PFmuChambers_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuD0_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuDz_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuChi2NDF_[kMaxnPFMu];   //[nPFMu_]
   Int_t           PFmuNumberOfValidTrkHits_[kMaxnPFMu];   //[nPFMu_]
   Int_t           PFmuNumberOfValidPixelHits_[kMaxnPFMu];   //[nPFMu_]
   Int_t           PFmuNumberOfValidMuonHits_[kMaxnPFMu];   //[nPFMu_]
   Int_t           nPFchad_;
   Int_t           PFhad_charge_[kMaxnPFchad];   //[nPFchad_]
   Float_t         PFchad_Eta_[kMaxnPFchad];   //[nPFchad_]
   Float_t         PFchad_Phi_[kMaxnPFchad];   //[nPFchad_]
   Float_t         PFchad_Pt_[kMaxnPFchad];   //[nPFchad_]
   Float_t         PFchad_P_[kMaxnPFchad];   //[nPFchad_]
   Float_t         PFchad_E_[kMaxnPFchad];   //[nPFchad_]
   Float_t         PFchad_Ecal_[kMaxnPFchad];   //[nPFchad_]
   Float_t         PFchad_Hcal_[kMaxnPFchad];   //[nPFchad_]
   Int_t           nPFnhad_;
   Float_t         PFnhad_Eta_[kMaxnPFnhad];   //[nPFnhad_]
   Float_t         PFnhad_Phi_[kMaxnPFnhad];   //[nPFnhad_]
   Float_t         PFnhad_Pt_[kMaxnPFnhad];   //[nPFnhad_]
   Float_t         PFnhad_P_[kMaxnPFnhad];   //[nPFnhad_]
   Float_t         PFnhad_E_[kMaxnPFnhad];   //[nPFnhad_]
   Float_t         PFnhad_Hcal_[kMaxnPFnhad];   //[nPFnhad_]
   Float_t         rho2012;
   Float_t         rho25;
   Float_t         rho44;
   Float_t         rho25_muPFiso;
   Int_t           nJet;
   Int_t           jetTrg[470][14];   //[nJet]
   Float_t         jetEn[470];   //[nJet]
   Float_t         jetPt[470];   //[nJet]
   Float_t         jetEta[470];   //[nJet]
   Float_t         jetPhi[470];   //[nJet]
   Float_t         jetEt[470];   //[nJet]
   Float_t         jetRawPt[470];   //[nJet]
   Float_t         jetRawPtSmeared[470];   //[nJet]
   Float_t         jetRawEn[470];   //[nJet]
   Float_t         jetRawEnSmeared[470];   //[nJet]
   Float_t         jetArea[470];   //[nJet]
   Float_t         jetCHF[470];   //[nJet]
   Float_t         jetNHF[470];   //[nJet]
   Float_t         jetCEF[470];   //[nJet]
   Float_t         jetNEF[470];   //[nJet]
   Float_t         jetMEF[470];   //[nJet]
   Int_t           jetNCH[470];   //[nJet]
   Float_t         jetHFHAE[470];   //[nJet]
   Float_t         jetHFEME[470];   //[nJet]
   Int_t           jetNConstituents[470];   //[nJet]
   Float_t         jetLeadTrackPt[470];   //[nJet]
   Float_t         jetSoftLeptPt[470];   //[nJet]
   Float_t         jetSoftLeptPtRel[470];   //[nJet]
   Float_t         jetSoftLeptdR[470];   //[nJet]
   Float_t         jetVtx3deL[470];   //[nJet]
   Float_t         jetVtxMass[470];   //[nJet]
   Float_t         jetTrackCountHiEffBJetTags[470];   //[nJet]
   Float_t         jetTrackCountHiPurBJetTags[470];   //[nJet]
   Float_t         jetSimpleSVHiEffBJetTags[470];   //[nJet]
   Float_t         jetSimpleSVHiPurBJetTags[470];   //[nJet]
   Float_t         jetCombinedSecondaryVtxBJetTags[470];   //[nJet]
   Float_t         jetCombinedSecondaryVtxMVABJetTags[470];   //[nJet]
   Float_t         jetJetProbabilityBJetTags[470];   //[nJet]
   Float_t         jetJetBProbabilityBJetTags[470];   //[nJet]
   Float_t         jetTrackCountingHighPurBJetTags[470];   //[nJet]
   Float_t         jetJECUnc[470];   //[nJet]
   Int_t           jetPartonID[470];   //[nJet]
   Int_t           jetGenJetIndex[470];   //[nJet]
   Float_t         jetGenJetEn[470];   //[nJet]
   Float_t         jetGenJetPt[470];   //[nJet]
   Float_t         jetGenJetEta[470];   //[nJet]
   Float_t         jetGenJetPhi[470];   //[nJet]
   Int_t           jetGenPartonID[470];   //[nJet]
   Float_t         jetGenEn[470];   //[nJet]
   Float_t         jetGenPt[470];   //[nJet]
   Float_t         jetGenEta[470];   //[nJet]
   Float_t         jetGenPhi[470];   //[nJet]
   Float_t         jetVtxPt[470];
   Float_t         jetVtx3dL[470];
   Float_t         jetDPhiMETJet[470];
   //vars for regression reader:
   Float_t jetRegPt, jetRegEta, jetRegEMFrac, jetRegHadFrac, jetRegNConstituents, jetRegVtx3dL, jetRegMET, jetRegDPhiMETJet, jetRegRho25;

   /// Variables for PU Jet ID
   Float_t         jetMVAs[500][4]; // [4] -> [3] for new ntuples before May30
   Int_t           jetWPLevels[500][4]; // [4] -> [3] for new ntuples before May30
   Float_t         jetMVAsExt[500][4][100]; // [4] -> [3] for new ntuples before May30
   Int_t           jetWPLevelsExt[500][4][100]; // [4] -> [3] for new ntuples before May30
   // std::vector<float * > jetMVAs_;
   // std::vector<int * >   jetWPLevels_;
   // std::vector<std::vector<std::vector<float> > * > jetMVAsExt_;
   // std::vector<std::vector<std::vector<int> > * > jetWPLevelsExt_;
   Float_t         jetDRMean[500];
   Float_t         jetDR2Mean[500];
   Float_t         jetDZ[500];
   Float_t         jetFrac01[500];
   Float_t         jetFrac02[500];
   Float_t         jetFrac03[500];
   Float_t         jetFrac04[500];
   Float_t         jetFrac05[500];
   Float_t         jetFrac06[500];
   Float_t         jetFrac07[500];
   Float_t         jetBeta[500];
   Float_t         jetBetaStarCMG[500];
   Float_t         jetBetaStarClassic[500];
   Float_t         jetBetaExt[500][100];
   Float_t         jetBetaStarCMGExt[500][100];
   Float_t         jetBetaStarClassicExt[500][100];
   Float_t         jetNNeutrals[500];
   Float_t         jetNCharged[500];
   Bool_t          jetPFLooseId[500];
   /// Variables for Low Pt Jets
   Int_t           nLowPtJet;
   Float_t         jetLowPtEn[500];
   Float_t         jetLowPtPt[500];
   Float_t         jetLowPtEta[500];
   Float_t         jetLowPtPhi[500];
   Float_t         jetLowPtCharge[500];
   Float_t         jetLowPtEt[500];
   Float_t         jetLowPtRawPt[500];
   Float_t         jetLowPtRawEn[500];
   Float_t         jetLowPtArea[500];
   Float_t         jetLowPtGenJetEn[500];
   Float_t         jetLowPtGenJetPt[500];
   Float_t         jetLowPtGenJetEta[500];
   Float_t         jetLowPtGenJetPhi[500];
   Int_t           jetLowPtPartonID[500];
   Int_t           jetLowPtGenPartonID[500];
   Float_t         jetLowPtGenEn[500];
   Float_t         jetLowPtGenPt[500];
   Float_t         jetLowPtGenEta[500];
   Float_t         jetLowPtGenPhi[500];
   /// Zee
   Int_t           nZee;
   Float_t         ZeeMass[1500];   //[nZee]
   Float_t         ZeePt[1500];   //[nZee]
   Float_t         ZeeEta[1500];   //[nZee]
   Float_t         ZeePhi[1500];   //[nZee]
   Int_t           ZeeLeg1Index[1500];   //[nZee]
   Int_t           ZeeLeg2Index[1500];   //[nZee]
   Int_t           nZmumu;
   Float_t         ZmumuMass[2];   //[nZmumu]
   Float_t         ZmumuPt[2];   //[nZmumu]
   Float_t         ZmumuEta[2];   //[nZmumu]
   Float_t         ZmumuPhi[2];   //[nZmumu]
   Int_t           ZmumuLeg1Index[2];   //[nZmumu]
   Int_t           ZmumuLeg2Index[2];   //[nZmumu]
   Int_t           nWenu;
   Float_t         WenuMassTPfMET[6];   //[nWenu]
   Float_t         WenuEtPfMET[6];   //[nWenu]
   Float_t         WenuACopPfMET[6];   //[nWenu]
   Int_t           WenuEleIndex[6];   //[nWenu]
   Int_t           nWmunu;
   Float_t         WmunuMassTPfMET[3];   //[nWmunu]
   Float_t         WmunuEtPfMET[3];   //[nWmunu]
   Float_t         WmunuACopPfMET[3];   //[nWmunu]
   Int_t           WmunuMuIndex[3];   //[nWmunu]
   Int_t           nConv;
   Float_t         convP4[1000][4];   //[nConv]
   Float_t         convVtx[1000][3];   //[nConv]
   Float_t         convVtxErr[1000][3];   //[nConv]
   Float_t         convPairMomentum[1000][3];   //[nConv]
   Float_t         convRefittedMomentum[1000][3];   //[nConv]
   Int_t           convNTracks[1000];   //[nConv]
   Float_t         convPairInvMass[1000];   //[nConv]
   Float_t         convPairCotThetaSep[1000];   //[nConv]
   Float_t         convEoverP[1000];   //[nConv]
   Float_t         convDistOfMinApproach[1000];   //[nConv]
   Float_t         convDPhiTrksAtVtx[1000];   //[nConv]
   Float_t         convDPhiTrksAtEcal[1000];   //[nConv]
   Float_t         convDEtaTrksAtEcal[1000];   //[nConv]
   Float_t         convDxy[1000];   //[nConv]
   Float_t         convDz[1000];   //[nConv]
   Float_t         convLxy[1000];   //[nConv]
   Float_t         convLz[1000];   //[nConv]
   Float_t         convZofPrimVtxFromTrks[1000];   //[nConv]
   Int_t           convNHitsBeforeVtx[1000][2];   //[nConv]
   Int_t           convNSharedHits[1000];   //[nConv]
   Int_t           convValidVtx[1000];   //[nConv]
   Float_t         convMVALikelihood[1000];   //[nConv]
   Float_t         convChi2[1000];   //[nConv]
   Float_t         convChi2Probability[1000];   //[nConv]
   Float_t         convTk1Dz[1000];   //[nConv]
   Float_t         convTk2Dz[1000];   //[nConv]
   Float_t         convTk1DzErr[1000];   //[nConv]
   Float_t         convTk2DzErr[1000];   //[nConv]
   Int_t           convCh1Ch2[1000];   //[nConv]
   Float_t         convTk1D0[1000];   //[nConv]
   Float_t         convTk1Pout[1000];   //[nConv]
   Float_t         convTk1Pin[1000];   //[nConv]
   Float_t         convTk2D0[1000];   //[nConv]
   Float_t         convTk2Pout[1000];   //[nConv]
   Float_t         convTk2Pin[1000];   //[nConv]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_HLT;   //!
   TBranch        *b_HLTIndex;   //!
   TBranch        *b_bspotPos;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vtxNTrk;   //!
   TBranch        *b_vtxNDF;   //!
   TBranch        *b_vtxD0;   //!
   TBranch        *b_IsVtxGood;   //!
   TBranch        *b_nVtxBS;   //!
   TBranch        *b_vtxbs;   //!
   TBranch        *b_vtxbsPtMod;   //!
   TBranch        *b_vtxbsSumPt2;   //!
   TBranch        *b_vtxbsTkPxVec;   //!
   TBranch        *b_vtxbsTkPyVec;   //!
   TBranch        *b_vtxbsTkIndex;   //!
   TBranch        *b_vtxbsTkWeight;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkP;   //!
   TBranch        *b_trkVtx;   //!
   TBranch        *b_trkd0;   //!
   TBranch        *b_trkd0Err;   //!
   TBranch        *b_trkdz;   //!
   TBranch        *b_trkdzErr;   //!
   TBranch        *b_trkPtErr;   //!
   TBranch        *b_trkQuality;   //!
   TBranch        *b_nGoodTrk;   //!
   TBranch        *b_IsTracksGood;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcIndex;   //!
   TBranch        *b_mcDecayType;   //!
   TBranch        *b_mcCalIsoDR03;   //!
   TBranch        *b_mcTrkIsoDR03;   //!
   TBranch        *b_mcCalIsoDR04;   //!
   TBranch        *b_mcTrkIsoDR04;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMETx;   //!
   TBranch        *b_genMETy;   //!
   TBranch        *b_genMETPhi;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METx;   //!
   TBranch        *b_METy;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_METsumEt;   //!
   /* TBranch        *b_tcMET;   //! */
   /* TBranch        *b_tcMETx;   //! */
   /* TBranch        *b_tcMETy;   //! */
   /* TBranch        *b_tcMETPhi;   //! */
   /* TBranch        *b_tcMETsumEt;   //! */
   /* TBranch        *b_tcMETmEtSig;   //! */
   /* TBranch        *b_tcMETSig;   //! */
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETx;   //!
   TBranch        *b_pfMETy;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_recoPfMET;   //!
   TBranch        *b_recoPfMETx;   //!
   TBranch        *b_recoPfMETy;   //!
   TBranch        *b_recoPfMETPhi;   //!
   TBranch        *b_recoPfMETsumEt;   //!
   TBranch        *b_recoPfMETmEtSig;   //!
   TBranch        *b_recoPfMETSig;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleTrg;   //!
   TBranch        *b_eleID;   //!
   TBranch        *b_eleClass;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleEtaVtx;   //!
   TBranch        *b_elePhiVtx;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleVtx;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_elePin;   //!
   TBranch        *b_elePout;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleSigmaIEtaIPhi;   //!
   TBranch        *b_eleSigmaIPhiIPhi;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleD0Vtx;   //!
   TBranch        *b_eleDzVtx;   //!
   TBranch        *b_eleEmax;   //!
   TBranch        *b_eleE3x3;   //!
   TBranch        *b_eleE5x5;   //!
   TBranch        *b_eleE2x5Right;   //!
   TBranch        *b_eleE2x5Left;   //!
   TBranch        *b_eleE2x5Top;   //!
   TBranch        *b_eleE2x5Bottom;   //!
   TBranch        *b_eleSeedTime;   //!
   TBranch        *b_eleRecoFlag;   //!
   TBranch        *b_eleGenIndex;   //!
   TBranch        *b_eleGenGMomPID;   //!
   TBranch        *b_eleGenMomPID;   //!
   TBranch        *b_eleGenMomPt;   //!
   TBranch        *b_eleIsoTrkDR03;   //!
   TBranch        *b_eleIsoEcalDR03;   //!
   TBranch        *b_eleIsoHcalDR03;   //!
   TBranch        *b_eleIsoHcalSolidDR03;   //!
   TBranch        *b_eleIsoTrkDR04;   //!
   TBranch        *b_eleIsoEcalDR04;   //!
   TBranch        *b_eleIsoHcalDR04;   //!
   TBranch        *b_eleIsoHcalSolidDR04;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleConvDist;   //!
   TBranch        *b_eleConvDcot;   //!
   TBranch        *b_elePFChIso03;
   TBranch        *b_elePFPhoIso03;
   TBranch        *b_elePFNeuIso03;
   TBranch        *b_eleIDMVANonTrig;
   TBranch        *b_eleEcalEn;
   TBranch        *b_eleConvVtxFit;
   TBranch        *b_eleESEffSigmaRR;   //!
   TBranch        *b_eleESE1;   //!
   TBranch        *b_eleESE3;   //!
   TBranch        *b_eleESE5;   //!
   TBranch        *b_eleESE7;   //!
   TBranch        *b_eleESE11;   //!
   TBranch        *b_eleESE21;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoTrg;   //!
   TBranch        *b_phoIsPhoton;   //!
   TBranch        *b_phoSCPos;   //!
   TBranch        *b_phoCaloPos;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoRegrE;   //!
   TBranch        *b_phoRegrErr;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoVtx;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEtVtx;   //!
   TBranch        *b_phoEtaVtx;   //!
   TBranch        *b_phoPhiVtx;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoCetaCorrE;   //!
   TBranch        *b_phoCetaCorrEt;   //!
   TBranch        *b_phoBremCorrE;   //!
   TBranch        *b_phoBremCorrEt;   //!
   TBranch        *b_phoFullCorrE;   //!
   TBranch        *b_phoFullCorrEt;   //!
   TBranch        *b_phoTrkIsoSolidDR03;   //!
   TBranch        *b_phoTrkIsoHollowDR03;   //!
   TBranch        *b_phoEcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoSolidDR03;   //!
   TBranch        *b_phoTrkIsoSolidDR04;   //!
   TBranch        *b_phoTrkIsoHollowDR04;   //!
   TBranch        *b_phoTrkIsoHollowNoDzDR03;   //!
   TBranch        *b_phoTrkIsoHollowNoDzDR04;   //!
   TBranch        *b_phoCiCTrkIsoDR03;   //!
   TBranch        *b_phoCiCTrkIsoDR04;   //!
   //TBranch        *b_phoCiCTrkIsoDR03Vtx;   //!
   //TBranch        *b_phoCiCTrkIsoDR04Vtx;   //!


  



   TBranch        *b_phoCiCPF4phopfIso03;
   TBranch        *b_phoCiCPF4phopfIso04;
   TBranch        *b_phoCiCPF4chgpfIso02;
   TBranch        *b_phoCiCPF4chgpfIso03;
   TBranch        *b_phoCiCPF4chgpfIso04;
   TBranch        *b_phoCiCPF4chgpfIso03AppVeto;
   TBranch        *b_phoCiCPF4chgpfIso04AppVeto;
   TBranch        *b_phoCiCPF4phopfIso03AppVeto;
   TBranch        *b_phoCiCPF4phopfIso04AppVeto;
   TBranch        *b_phoCiCPF4phopfIso03AppVetoClean;
   TBranch        *b_phoCiCPF4phopfIso04AppVetoClean;
   TBranch        *b_phoCiCPF4chgpfIso03Mod;
   TBranch        *b_phoCiCPF4chgpfIso04Mod;
   //  TBranch        *b_;
   // TBranch        *b_;
   // TBranch        *b_;
   // TBranch        *b_;

   TBranch        *b_phoEleVeto;   // ----  2012
   TBranch        *b_phoCiCdRtoTrk;   //!
   TBranch        *b_phoEcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoSolidDR04;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoSigmaIEtaIPhi;   //!
   TBranch        *b_phoSigmaIPhiIPhi;   //!
   TBranch        *b_phoEmax;   //!
   TBranch        *b_phoE2x2;   //!
   TBranch        *b_phoE3x3;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoE2x5Right;   //!
   TBranch        *b_phoE2x5Left;   //!
   TBranch        *b_phoE2x5Top;   //!
   TBranch        *b_phoE2x5Bottom;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedDetId1;   //!
   TBranch        *b_phoSeedDetId2;   //!
   TBranch        *b_phoRecoFlag;   //!
   TBranch        *b_phoPos;   //!
   TBranch        *b_phoGenIndex;   //!
   TBranch        *b_phoGenGMomPID;   //!
   TBranch        *b_phoGenMomPID;   //!
   TBranch        *b_phoGenMomPt;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoSCEt;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phoOverlap;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoIsConv;   //!
   TBranch        *b_phoNConv;   //!
   TBranch        *b_phoConvInvMass;   //!
   TBranch        *b_phoConvCotTheta;   //!
   TBranch        *b_phoConvEoverP;   //!
   TBranch        *b_phoConvZofPVfromTrks;   //!
   TBranch        *b_phoConvMinDist;   //!
   TBranch        *b_phoConvdPhiAtVtx;   //!
   TBranch        *b_phoConvdPhiAtCalo;   //!
   TBranch        *b_phoConvdEtaAtCalo;   //!
   TBranch        *b_phoConvTrkd0;   //!
   TBranch        *b_phoConvTrkPin;   //!
   TBranch        *b_phoConvTrkPout;   //!
   TBranch        *b_phoConvTrkdz;   //!
   TBranch        *b_phoConvTrkdzErr;   //!
   TBranch        *b_phoConvChi2;   //!
   TBranch        *b_phoConvChi2Prob;   //!
   TBranch        *b_phoConvNTrks;   //!
   TBranch        *b_phoConvCharge;   //!
   TBranch        *b_phoConvValidVtx;   //!
   TBranch        *b_phoConvLikeLihood;   //!
   TBranch        *b_phoConvP4;   //!
   TBranch        *b_phoConvVtx;   //!
   TBranch        *b_phoConvVtxErr;   //!
   TBranch        *b_phoConvPairMomentum;   //!
   TBranch        *b_phoConvRefittedMomentum;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_phoESE1;   //!
   TBranch        *b_phoESE3;   //!
   TBranch        *b_phoESE5;   //!
   TBranch        *b_phoESE7;   //!
   TBranch        *b_phoESE11;   //!
   TBranch        *b_phoESE21;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muTrg;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muPz;   //!
   TBranch        *b_muGenIndex;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muIsoCalo;   //!
   TBranch        *b_muIsoEcal;   //!
   TBranch        *b_muIsoHcal;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muEmVeto;   //!
   TBranch        *b_muHadVeto;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muID;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muVtx;   //!
   TBranch        *b_muD0Vtx;   //!
   TBranch        *b_muDzVtx;   //!   
   TBranch        *b_muNumberOfValidTrkHits;   //!
   TBranch        *b_muNumberOfValidPixelHits;   //!
   TBranch        *b_muNumberOfValidMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muChambers;   //!

  //------------ 2012 --------------------
   TBranch        *b_muNumberOfValidTrkLayers;
   TBranch        *b_muPFIsoR04_CH;
   TBranch        *b_muPFIsoR04_NH;
   TBranch        *b_muPFIsoR04_Pho;
   TBranch        *b_muPFIsoR04_PU;
   TBranch        *b_muPFIsoR04_CPart;
   TBranch        *b_muPFIsoR04_NHHT;
   TBranch        *b_muPFIsoR04_PhoHT;
   
   TBranch        *b_muPFIsoR03_CH;
   TBranch        *b_muPFIsoR03_NH;
   TBranch        *b_muPFIsoR03_Pho;
   TBranch        *b_muPFIsoR03_PU;
   TBranch        *b_muPFIsoR03_CPart;
   TBranch        *b_muPFIsoR03_NHHT;
   TBranch        *b_muPFIsoR03_PhoHT;
   //------------------------


   TBranch        *b_nPFPho_;   //!
   TBranch        *b_isIsoPFPho_;   //!
   TBranch        *b_PFphoE_;   //!
   TBranch        *b_PFphoEt_;   //!
   TBranch        *b_PFphoEta_;   //!
   TBranch        *b_PFphoPhi_;   //!
   TBranch        *b_PFphoIsConv_;   //!
   TBranch        *b_PFphoTrkIsoHollowDR04_;   //!
   TBranch        *b_PFphoEcalIsoDR04_;   //!
   TBranch        *b_PFphoHcalIsoDR04_;   //!
   TBranch        *b_PFphoHoverE_;   //!
   TBranch        *b_nPFElec_;   //!
   TBranch        *b_PFelePt_;   //!
   TBranch        *b_PFeleEta_;   //!
   TBranch        *b_PFelePhi_;   //!
   TBranch        *b_PFeleEn_;   //!
   TBranch        *b_PFeleCharge;   //!
   TBranch        *b_nPFMu_;   //!
   TBranch        *b_PFmuCharge_;   //!
   TBranch        *b_PFmuPhi_;   //!
   TBranch        *b_PFmuEta_;   //!
   TBranch        *b_PFmuPt_;   //!
   TBranch        *b_PFmuPz_;   //!
   TBranch        *b_PFmuChambers_;   //!
   TBranch        *b_PFmuD0_;   //!
   TBranch        *b_PFmuDz_;   //!
   TBranch        *b_PFmuChi2NDF_;   //!
   TBranch        *b_PFmuNumberOfValidTrkHits_;   //!
   TBranch        *b_PFmuNumberOfValidPixelHits_;   //!
   TBranch        *b_PFmuNumberOfValidMuonHits_;   //!
   TBranch        *b_nPFchad_;   //!
   TBranch        *b_PFhad_charge_;   //!
   TBranch        *b_PFchad_Eta_;   //!
   TBranch        *b_PFchad_Phi_;   //!
   TBranch        *b_PFchad_Pt_;   //!
   TBranch        *b_PFchad_P_;   //!
   TBranch        *b_PFchad_E_;   //!
   TBranch        *b_PFchad_Ecal_;   //!
   TBranch        *b_PFchad_Hcal_;   //!
   TBranch        *b_nPFnhad_;   //!
   TBranch        *b_PFnhad_Eta_;   //!
   TBranch        *b_PFnhad_Phi_;   //!
   TBranch        *b_PFnhad_Pt_;   //!
   TBranch        *b_PFnhad_P_;   //!
   TBranch        *b_PFnhad_E_;   //!
   TBranch        *b_PFnhad_Hcal_;   //!
   TBranch        *b_rho2012;   //!
   TBranch        *b_rho25;   //!
   TBranch        *b_rho44;   //!
   TBranch        *b_rho25_muPFiso;  //  2012
   TBranch        *b_nJet;   //!
   TBranch        *b_jetTrg;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetEt;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetArea;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetMEF;   //!
   TBranch        *b_jetNCH;   //!
   TBranch        *b_jetHFHAE;   //!
   TBranch        *b_jetHFEME;   //!
   TBranch        *b_jetNConstituents;   //!
   TBranch        *b_jetLeadTrackPt;   //!
   TBranch        *b_jetSoftLeptPt;   //!
   TBranch        *b_jetSoftLeptPtRel;   //!
   TBranch        *b_jetSoftLeptdR;   //!
   TBranch        *b_jetVtx3deL;   //!
   TBranch        *b_jetVtxMass;   //!
   TBranch        *b_jetTrackCountHiEffBJetTags;   //!
   TBranch        *b_jetTrackCountHiPurBJetTags;   //!
   TBranch        *b_jetSimpleSVHiEffBJetTags;   //!
   TBranch        *b_jetSimpleSVHiPurBJetTags;   //!
   TBranch        *b_jetCombinedSecondaryVtxBJetTags;   //!
   TBranch        *b_jetCombinedSecondaryVtxMVABJetTags;   //!
   TBranch        *b_jetJetProbabilityBJetTags;   //!
   TBranch        *b_jetJetBProbabilityBJetTags;   //!
   TBranch        *b_jetTrackCountingHighPurBJetTags;   //!
   TBranch        *b_jetJECUnc;   //!
   TBranch        *b_jetPartonID;   //!
   TBranch        *b_jetGenJetIndex;   //!
   TBranch        *b_jetGenJetEn;   //!
   TBranch        *b_jetGenJetPt;   //!
   TBranch        *b_jetGenJetEta;   //!
   TBranch        *b_jetGenJetPhi;   //!
   TBranch        *b_jetGenPartonID;   //!
   TBranch        *b_jetGenEn;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetGenEta;   //!
   TBranch        *b_jetGenPhi;   //!
   TBranch        *b_jetVtxPt;
   TBranch        *b_jetVtx3dL;
   TBranch        *b_jetDPhiMETJet;
   TBranch        *b_jetMVAs;
   TBranch        *b_jetWPLevels;
   TBranch        *b_jetMVAsExt;
   TBranch        *b_jetWPLevelsExt;
   TBranch        *b_jetDRMean;
   TBranch        *b_jetDR2Mean;
   TBranch        *b_jetDZ;
   TBranch        *b_jetFrac01;
   TBranch        *b_jetFrac02;
   TBranch        *b_jetFrac03;
   TBranch        *b_jetFrac04;
   TBranch        *b_jetFrac05;
   TBranch        *b_jetFrac06;
   TBranch        *b_jetFrac07;
   TBranch        *b_jetBeta;
   TBranch        *b_jetBetaStarCMG;
   TBranch        *b_jetBetaStarClassic;
   TBranch        *b_jetBetaExt;
   TBranch        *b_jetBetaStarCMGExt;
   TBranch        *b_jetBetaStarClassicExt;
   TBranch        *b_jetNNeutrals;
   TBranch        *b_jetNCharged;
   TBranch        *b_jetPFLooseId;
   TBranch        *b_nLowPtJet;
   TBranch        *b_jetLowPtEn;
   TBranch        *b_jetLowPtPt;
   TBranch        *b_jetLowPtEta;
   TBranch        *b_jetLowPtPhi;
   TBranch        *b_jetLowPtCharge;
   TBranch        *b_jetLowPtEt;
   TBranch        *b_jetLowPtRawPt;
   TBranch        *b_jetLowPtRawEn;
   TBranch        *b_jetLowPtArea;
   TBranch        *b_jetLowPtGenJetEn;
   TBranch        *b_jetLowPtGenJetPt;
   TBranch        *b_jetLowPtGenJetEta;
   TBranch        *b_jetLowPtGenJetPhi;
   TBranch        *b_jetLowPtPartonID;
   TBranch        *b_jetLowPtGenPartonID;
   TBranch        *b_jetLowPtGenEn;
   TBranch        *b_jetLowPtGenPt;
   TBranch        *b_jetLowPtGenEta;
   TBranch        *b_jetLowPtGenPhi;
   TBranch        *b_nZee;   //!
   TBranch        *b_ZeeMass;   //!
   TBranch        *b_ZeePt;   //!
   TBranch        *b_ZeeEta;   //!
   TBranch        *b_ZeePhi;   //!
   TBranch        *b_ZeeLeg1Index;   //!
   TBranch        *b_ZeeLeg2Index;   //!
   TBranch        *b_nZmumu;   //!
   TBranch        *b_ZmumuMass;   //!
   TBranch        *b_ZmumuPt;   //!
   TBranch        *b_ZmumuEta;   //!
   TBranch        *b_ZmumuPhi;   //!
   TBranch        *b_ZmumuLeg1Index;   //!
   TBranch        *b_ZmumuLeg2Index;   //!
   TBranch        *b_nWenu;   //!
   TBranch        *b_WenuMassTPfMET;   //!
   TBranch        *b_WenuEtPfMET;   //!
   TBranch        *b_WenuACopPfMET;   //!
   TBranch        *b_WenuEleIndex;   //!
   TBranch        *b_nWmunu;   //!
   TBranch        *b_WmunuMassTPfMET;   //!
   TBranch        *b_WmunuEtPfMET;   //!
   TBranch        *b_WmunuACopPfMET;   //!
   TBranch        *b_WmunuMuIndex;   //!
   TBranch        *b_nConv;   //!
   TBranch        *b_convP4;   //!
   TBranch        *b_convVtx;   //!
   TBranch        *b_convVtxErr;   //!
   TBranch        *b_convPairMomentum;   //!
   TBranch        *b_convRefittedMomentum;   //!
   TBranch        *b_convNTracks;   //!
   TBranch        *b_convPairInvMass;   //!
   TBranch        *b_convPairCotThetaSep;   //!
   TBranch        *b_convEoverP;   //!
   TBranch        *b_convDistOfMinApproach;   //!
   TBranch        *b_convDPhiTrksAtVtx;   //!
   TBranch        *b_convDPhiTrksAtEcal;   //!
   TBranch        *b_convDEtaTrksAtEcal;   //!
   TBranch        *b_convDxy;   //!
   TBranch        *b_convDz;   //!
   TBranch        *b_convLxy;   //!
   TBranch        *b_convLz;   //!
   TBranch        *b_convZofPrimVtxFromTrks;   //!
   TBranch        *b_convNHitsBeforeVtx;   //!
   TBranch        *b_convNSharedHits;   //!
   TBranch        *b_convValidVtx;   //!
   TBranch        *b_convMVALikelihood;   //!
   TBranch        *b_convChi2;   //!
   TBranch        *b_convChi2Probability;   //!
   TBranch        *b_convTk1Dz;   //!
   TBranch        *b_convTk2Dz;   //!
   TBranch        *b_convTk1DzErr;   //!
   TBranch        *b_convTk2DzErr;   //!
   TBranch        *b_convCh1Ch2;   //!
   TBranch        *b_convTk1D0;   //!
   TBranch        *b_convTk1Pout;   //!
   TBranch        *b_convTk1Pin;   //!
   TBranch        *b_convTk2D0;   //!
   TBranch        *b_convTk2Pout;   //!
   TBranch        *b_convTk2Pin;   //!


//   xAna(Int_t mode = 0);
   xAna(TTree *tree=0);
   virtual ~xAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   // virtual void     HggTreeWriteLoop(const char* filename, Int_t mode = 0, bool correctVertex=true, bool correctEnergy=true, bool setRho0=false);

   /// FC: add some function for reweighting
   void SetWeightManager( WeightManager *weight_manager);
   /// FC: Trigger selection (file TriggerSelection.h)
   bool PassTriggerSelection(void);

   /// FC: MVA function
   float PhoID_MVA( int i, int ivtx ) ;
   float DiPhoID_MVA( const TLorentzVector &glead, const TLorentzVector &gtrail, const TLorentzVector &hcand, 
		      const MassResolution & massResoCalc, float vtxProb,
		      const float &idlead, const float &idtrail );
   float JetRegression( int iJet );

   void ReweightMC_phoIdVar( int i );
   void ReweightMC_phoIdVarBaseline( int i );

   /// FC: electron and conversion matching
   bool isMatchedToAPromptElec( int ipho, float &dRclosestElec );
   bool isAGoodConv( int ipho );
   bool isInGAP_EB( int ipho);
   
   void Setup_MVA(void);
   
   /// FC: met selection
   bool MetTagSelection( const TLorentzVector &g1, const TLorentzVector &g2, const vector<int> &jindex );

   /// FC: MC truth functions
   bool isZgamma(void);
   void fillMCtruthInfo_HiggsSignal(void);

   
   virtual int     HggTreeWriteLoop(const char* filename, int ijob, bool correctVertex=true, 
				     bool correctEnergy=true, bool setRho0=false
				     );
   //virtual void     ZggTreeWriteLoop(const char* filename, Int_t mode = 0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual Int_t photonID(int ipho, int ivtx, float corEt_ipho = -1, bool invEleVeto = false );
   virtual Int_t cicPhotonID(int ipho, int ivtx, float corEt_ipho = -1 );
   virtual Int_t cic4PFPhotonID(int ipho, int ivtx, float corEt_ipho = -1, bool invEleVeto = false );
   virtual Int_t mvaPhotonID(int ipho, int ivtx, float corEt_ipho = -1, bool invEleVeto = false );
   virtual void  fillPhotonVariablesToMiniTree(int ipho, int ivtx, int index );

//   virtual Int_t    romePhotonID(int ipho);
   virtual Int_t    phocat(int i);
   virtual Int_t    phocatCIC(int i);
   virtual Int_t    cat4(float r1, float r2, float eta1, float eta2);
   virtual Double_t deltaPhi(Double_t phi1, Double_t phi2);
   virtual Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
   virtual void     InitHists();

   virtual Int_t   matchPhotonToConversion( int lpho);
   virtual Float_t   etaTransformation(  float EtaParticle , float vertex);
   virtual PhotonInfo  fillPhotonInfo(  int iPho1, bool useAllConvs);
   //virtual int    getSelectedVertex(int iPho1, int iPho2, bool useAllConvs=true);
   //  int    getSelectedVertex(int iPho1, int iPho2, bool useAllConvs=true);
   vector<int>    getSelectedVertex(int iPho1, int iPho2, bool useAllConvs=true, bool useMva=true);

   //virtual Int_t photonCat(int i) ;
   //virtual bool preselectPhoton(int i) ;
   virtual bool preselectPhoton(int i, float Et) ;
   virtual TVector3 getCorPhotonTVector3(int i, int selVtx);
 


   virtual Double_t getHqTWeight(Float_t HiggsMCMass, Float_t HiggsMCPt);
   
   virtual bool selectTwoPhotons( int &i, int &j, int selVtx, 
				  TLorentzVector &g1, TLorentzVector &g2, bool elecVeto[2] );
   
   virtual Int_t*   selectElectrons2011(Int_t iPho1, Int_t iPho2, Double_t wei, Int_t vtxInd);
   virtual Int_t*   selectMuons2011(Int_t iPho1, Int_t iPho2, Double_t wei, Int_t vtxInd);
   virtual vector<int> selectElectrons(Int_t iPho1, Int_t iPho2, Double_t wei, Int_t vtxInd);
   virtual vector<int> selectMuons(Int_t iPho1, Int_t iPho2, Double_t wei, Int_t vtxInd);
   virtual vector<int> selectElectronsHCP2012(double wei, int &iLepVtx );
   virtual vector<int> selectMuonsHCP2012(    double wei, int &iLepVtx);

   vector<int>   selectJets(Int_t iPho1, Int_t iPho2, Int_t nVtxJetID, Double_t wei, Int_t selVtx);
   vector<int>   selectJetsJEC(Int_t iPho1, Int_t iPho2, Int_t nVtxJetID, Double_t wei, Int_t selVtx);
   vector<int>   selectJetsRadion(Int_t nVtxJetID, Int_t selVtx);
   void dijetSelection( const TLorentzVector &corPho1,const TLorentzVector &corPho2,
			const vector<int> &goodJetsIndex,
			const Double_t &wei, int ivtx,
			int &vbftag, int &hstratag, int &catjet );
   double correctJetPt(int i);

   void           applyJEC4SoftJets();
   void           METSmearCorrection(Int_t iPho1, Int_t iPho2);
   double         getMETJetSmearFactor(Float_t metJetEta);
   double         ErrEt( double Et, double Eta);
   void           METPhiShiftCorrection();
   void           METScaleCorrection(Int_t iPho1, Int_t iPho2);
   double         getMETJetScaleFactor(Float_t metJetEta);
   //----------------------
   

   std::ofstream _xcheckTextFile;
   bool DoCorrectVertex_;
   bool DoCorrectEnergy_;
   bool DoSetRho0_;
   bool DoDebugEvent;
   Int_t mode_;

   float lumi_,xsec_;
   int nEvts_;

   WeightManager *_weight_manager;
   
   HggVertexAnalyzer        *_vtxAna;
   HggVertexFromConversions *_vtxAnaFromConversions;
   VertexAlgoParameters     _vtxAlgoParams;
   TMVA::Reader *_perVtxReader;
   std::string _perVtxMvaWeights;
   std::string _perVtxMvaMethod;
   TMVA::Reader *_perEvtReader;
   std::string _perEvtMvaWeights;
   std::string _perEvtMvaMethod;

   ConfigReader *_config;
   MiniTree *_minitree;

   //ofstream outTextFile;
   TFile outFile;
   TFile *fout;
   TH1F *hTotEvents;
   TH1D *data_npu_estimated;

   TH1F *hITPU_wei;
   TH1F *hITPU;

   TH1F *hITPUTrue_wei;
   TH1F *hITPUTrue;

   TH1F *h_pho1Et;
   TH1F *h_pho2Et;
   TH1F *h_pho3Et;

   TH1F *hPhoCorrEt1;
   TH1F *hPhoCorrEt2;

   TH1F *hPhoCorrEt1_WH;
   TH1F *hPhoCorrEt2_WH;
 
   TH1F *hPhoCorrEt1_ZH;
   TH1F *hPhoCorrEt2_ZH;
   
   //------------------ lepton tag----------
   
   TH1F *hDMZmin;
   TH1F *hDPhi_gele_g;

   // TH2F *hDPhi_geleg_gele;


   TH2F *hDPhi_geleg_gele;
   
   TH1F *hcat4;
   TH1F *hggMass_eleVeto;
   TH1F *hggMass_eleVeto_jetVeto;
   
   TH1F *hDRtoTrk;
   TH1F *hggMass_muTag;
   TH1F *hggMass_eleTag;
   TH1F *hggMassWH_muTag;
   TH1F *hggMassWH_eleTag;
   TH1F *hggMassZH_muTag;
   TH1F *hggMassZH_eleTag;
   TH1F *hggMass_eletag_Zveto;
   TH1F *hggMassWH_eleTag_Zveto;
   TH1F *hggMassZH_eleTag_Zveto;
   TH1F *hDR_Pho1Ele;
   TH1F *hDR_Pho2Ele;
   TH1F *hDR_SCPho1Ele;
   TH1F *hDR_SCPho2Ele;
   TH1F *hDRmin_PhoEle;
   TH1F *hDRmin_SCPhoEle;
   TH1F *hDRmin_PhoMu;

   TH1F *hAllElePt;
   TH1F *hAllSCEleEta;
   TH1F *hEleCombRelIso_EB;
   TH1F *hEleCombRelIso_EE;
   TH1F *hSelEle1Pt;
   TH1F *hSelSCEle1Eta; 
   TH1F *hSelEle2Pt; 
   TH1F *hSelSCEle2Eta;

   TH1F *helePFRelIso_BID_B;
   TH1F *helePFRelIso_AID_B;
   TH1F *helePFRelIso_BID_EC;
   TH1F *helePFRelIso_AID_EC;

   TH1F *hmuPFRelIso_BID;
   TH1F *hmuPFRelIso_AID;

   TH1F *hMuCombRelIso;
   TH1F *hAllMuPt;
   TH1F *hAllMuEta;
  

   TH1F *hSelMu1Pt;
   TH1F *hSelMu1Eta; 
   TH1F *hSelMu2Pt; 
   TH1F *hSelMu2Eta;
  
   TH1F *h_pfMET;
   TH1F *h_tcMET;
   TH1F *h_CaloMET;

   TH1F *h_pfMET_lepTag;
   TH1F *h_tcMET_lepTag;
   TH1F *h_CaloMET_lepTag;

   TH1F *hMassEleGamGam;
   TH1F *hPtEleGamGam;

   TH1F *hMassMuGamGam;
   TH1F *hPtMuGamGam;


   TH1F *hMassEleEleGam1;
   TH1F *hMassEleEleGam2;
   TH1F *hPtEleEleGam1;
   TH1F *hPtEleEleGam2;
   TH1F *hMassEleGam1;
   TH1F *hMassEleGam2;
   TH1F *hPtEleGam1;
   TH1F *hPtEleGam2;

   TH1F *hMassEleFakePho;
   TH1F *hMassEleRealPho;



  //------------------two jet tag----------


   TH1F *hJetNHF;
   TH1F *hJetNEF;
   TH1F *hJetNConst;

   TH1F *hJetCHF;
   TH1F *hJetNCH;
   TH1F *hJetCEF;

   TH1F *hDRAllJetg1;
   TH1F *hDRAllJetg2;
   TH1F *hAllJetPt;
   TH1F *hAllJetEta;

   TH1F *hDRjet1g1;
   TH1F *hDRjet1g2;
   TH1F *hDRjet2g1;
   TH1F *hDRjet2g2;

   TH1F *hDPhijet1g1;
   TH1F *hDPhijet1g2;
   TH1F *hDPhijet2g1;
   TH1F *hDPhijet2g2;

   TH1F *hDPhijet1jet2;


   TH1F *hggMass_twoJetTag;
   TH1F *hggMass;
   TH1F *hggPt_twoJetTag;

   TH1F *hMassJetJet;   
   TH1F *hJet1Pt;
   TH1F *hJet2Pt;
   TH1F *hJet1Eta;
   TH1F *hJet2Eta;
   TH1F *hDeltaEtaJets;
   // TH2F *hDeltaEtaJetsVsMgg;
   //   TH2F *hMassJetJetVsMgg;
   TH1F *hZeppenfeld;
   TH1F *hDeltaPhiDiJetDiPho;
  
   TH1F *hZeppenfeld_bc;
   TH1F *hDeltaEtaJets_bc;
   TH1F *hJetEtaProd_bc;
   TH1F *hMassJetJet_bc;   
   TH1F *hDeltaPhiDiJetDiPho_bc;

   //TH1F *hMgg_twoJetTag_etaprod;
   TH1F *hMgg_twoJetTag_jet1Pt;
   TH1F *hMgg_twoJetTag_deltaEta;
   TH1F *hMgg_twoJetTag_zep;
   TH1F *hMgg_twoJetTag_Mjets;
   TH1F *hMgg_twoJetTag_dPhi;
  
   TH1F *hJet1Pt_HSTRA_bc;
   TH1F *hDeltaEtaJets_HSTRA_bc;
   TH1F *hZeppenfeld_HSTRA_bc;
   TH1F *hMassJetJet_HSTRA_bc;

   TH1F *hMgg_HSTRA_twoJetTag;
   TH1F *hMgg_HSTRA_twoJetTag_jet1Pt;
   TH1F *hMgg_HSTRA_twoJetTag_deltaEta;
   TH1F *hMgg_HSTRA_twoJetTag_zep;
   TH1F *hMgg_HSTRA_twoJetTag_Mjets;


   //----------------------------------------

};

/* #endif */

/* #ifdef xAna_cxx */


void xAna::SetConfigReader( const ConfigReader &config ) {
  if( _config != 0 ) delete _config;
  _config = new ConfigReader(config.configName());
  
}
void xAna::SetConfigReader( string configName ) {
  if( _config != 0 ) delete _config;
  _config = new ConfigReader(configName);
}

xAna::xAna(TTree *tree) {
  
  for( int i = 0 ; i < 2 ; i++ ) {
    phoID_2011[i] = 0;
    phoID_2012[i] = 0;
  }
  _config = 0;
  Init(tree);
  
}


Double_t xAna::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {

  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);

  return sqrt(dEta*dEta+dPhi*dPhi);
}

Double_t xAna::deltaPhi(Double_t phi1, Double_t phi2) {

  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();

  return dPhi;
}


Double_t xAna::getHqTWeight(Float_t HiggsMCMass, Float_t HiggsMCPt) {
  // HiggsMCMass=120;
  Double_t wei = 1;
  Double_t x = (Double_t) HiggsMCPt;
  

  Double_t par[6];
  if (HiggsMCMass == 90) {
    par[0] = 1.5512; par[1] = -0.0415444; par[2] = 0.00110316; par[3] = -1.51299e-05; par[4] = 7.51189e-08; par[5] = 0.574011;
  } else if (HiggsMCMass == 95) {
    par[0] = 1.50275; par[1] = -0.0342925; par[2] = 0.000796335; par[3] = -9.8689e-06; par[4] = 4.48063e-08; par[5] = 0.570917;
  } else if (HiggsMCMass == 100) {
    par[0] = 1.45047; par[1] = -0.0284948; par[2] = 0.000624926; par[3] = -7.66006e-06; par[4] = 3.42959e-08; par[5] = 0.574828;
  } else if (HiggsMCMass == 105) {
    par[0] = 1.29301; par[1] = -0.018701; par[2] = 0.000466271; par[3] = -6.21643e-06; par[4] = 2.84922e-08; par[5] = 0.659637;
  } else if (HiggsMCMass == 110) {
    par[0] = 1.41094; par[1] = -0.0216905; par[2] = 0.000379204; par[3] = -4.22212e-06; par[4] = 1.79075e-08; par[5] = 0.569545;
  } else if (HiggsMCMass == 115) {
    par[0] = 1.38536; par[1] = -0.0196079; par[2] = 0.000340012; par[3] = -3.85773e-06; par[4] = 1.64619e-08; par[5] = 0.569461;
  } else if (HiggsMCMass == 120) {
    par[0] = 1.40129; par[1] = -0.0206659; par[2] = 0.000339413; par[3] = -3.33799e-06; par[4] = 1.23947e-08; par[5] = 0.575397;
  } else if (HiggsMCMass == 121) {
    par[0] = 1.41674; par[1] = -0.0216618; par[2] = 0.000382855; par[3] = -4.06627e-06; par[4] = 1.59015e-08; par[5] = 0.563047;
  } else if (HiggsMCMass == 123) {
    par[0] = 1.38972; par[1] = -0.0189709; par[2] = 0.000287065; par[3] = -2.7141e-06; par[4] = 9.80252e-09; par[5] = 0.575012;
  } else if (HiggsMCMass == 125) {
    par[0] = 1.36709; par[1] = -0.0164974; par[2] = 0.000211968; par[3] = -1.82729e-06; par[4] = 6.32798e-09; par[5] = 0.588449;
  } else if (HiggsMCMass == 130) {
    par[0] = 1.40301; par[1] = -0.0204618; par[2] = 0.00034626; par[3] = -3.42286e-06; par[4] = 1.24308e-08; par[5] = 0.57312;
  } else if (HiggsMCMass == 135) {
    par[0] = 1.36409; par[1] = -0.0175548; par[2] = 0.000286312; par[3] = -2.86484e-06; par[4] = 1.05347e-08; par[5] = 0.560269;
  } else if (HiggsMCMass == 140) {
    par[0] = 1.32729; par[1] = -0.0133144; par[2] = 0.000156537; par[3] = -1.31663e-06; par[4] = 4.32718e-09; par[5] = 0.590839;
  } else if (HiggsMCMass == 150) {
    par[0] = 1.3531; par[1] = -0.0152357; par[2] = 0.000211272; par[3] = -1.84944e-06; par[4] = 6.00008e-09; par[5] = 0.582344;
  } else if (HiggsMCMass == 155) {
    par[0] = 1.33302; par[1] = -0.0132161; par[2] = 0.000154738; par[3] = -1.20301e-06; par[4] = 3.58784e-09; par[5] = 0.569991;
  } else if (HiggsMCMass == 160) {
    par[0] = 1.30238; par[1] = -0.0113894; par[2] = 0.000130351; par[3] = -1.05532e-06; par[4] = 3.16211e-09; par[5] = 0.598551;
  } else {
    return wei; 
  }

  wei = par[0] + par[1] * x + par[2] * pow(x, 2) + par[3] * pow(x, 3) + par[4] * pow(x, 4);
  if (x > HiggsMCMass) wei = par[5]; 

  return wei;
}

void xAna::SetWeightManager(WeightManager * weight_manager) {

  _weight_manager = weight_manager;
}


xAna::~xAna()
{
  ////  if( _config ) delete _config;
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t xAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t xAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void xAna::InitHists() {

  /* float xpt[11] = {10, 15, 20, 30, 40, 50, 60, 80, 100, 140, 200}; */
  /* float xvt[8]  = { 0,  2,  4,  6,  8, 10, 12,  20}; */

  hTotEvents =   new TH1F("hTotEvents",   "hTotEvents",3,0,3);
  //hITPU_wei =  new TH1F("hITPU_wei",   "hITPU_wei",60,0,60);
  // hITPU =  new TH1F("hITPU",   "hITPU",60,0,60);

  
  hITPU_wei =  new TH1F("hITPU_wei",   "hITPU_wei",100,0,100);
  hITPU =  new TH1F("hITPU",   "hITPU",100,0,100);

  hITPUTrue_wei =  new TH1F("hITPUTrue_wei",   "hITPUTrue_wei",100,0,100);
  hITPUTrue =  new TH1F("hITPUTrue",   "hITPUTrue",100,0,100);

  h_pho1Et =  new TH1F("h_pho1Et",   "h_pho1Et",200,0,200);
  h_pho2Et =  new TH1F("h_pho2Et",   "h_pho2Et",200,0,200);
  h_pho3Et =  new TH1F("h_pho3Et",   "h_pho3Et",200,0,200);
  hPhoCorrEt1 =  new TH1F("hPhoCorrEt1",   "hPhoCorrEt1",200,0,200);
  hPhoCorrEt2 =  new TH1F("hPhoCorrEt2",   "hPhoCorrEt2",200,0,200);
  hPhoCorrEt1_WH =  new TH1F("hPhoCorrEt1_WH",   "hPhoCorrEt1_WH",200,0,200);
  hPhoCorrEt2_WH =  new TH1F("hPhoCorrEt2_WH",   "hPhoCorrEt2_WH",200,0,200);
  hPhoCorrEt1_ZH =  new TH1F("hPhoCorrEt1_ZH",   "hPhoCorrEt1_ZH",200,0,200);
  hPhoCorrEt2_ZH =  new TH1F("hPhoCorrEt2_ZH",   "hPhoCorrEt2_ZH",200,0,200);

  /*
  hmc_pt_eb    = new TH1F("hmc_pt_eb",   "Pt EB",  10, xpt);
  hmc_pt_ee    = new TH1F("hmc_pt_ee",   "Pt EE",  10, xpt);
  hmc_eta      = new TH1F("hmc_eta",     "Eta",    20, -3, 3);
  hmc_npu_eb   = new TH1F("hmc_npu_eb",  "nPU EB",  7, xvt);
  hmc_npu_ee   = new TH1F("hmc_npu_ee",  "nPU EE",  7, xvt);
  hmc_nvtx_eb  = new TH1F("hmc_nvtx_eb", "nVtx EB", 7, xvt);
  hmc_nvtx_ee  = new TH1F("hmc_nvtx_ee", "nVtx EE", 7, xvt);

  hdata_pt_eb    = new TH1F("hdata_pt_eb",   "Pt EB",  10, xpt);
  hdata_pt_ee    = new TH1F("hdata_pt_ee",   "Pt EE",  10, xpt);
  hdata_eta      = new TH1F("hdata_eta",     "Eta",    20, -3, 3);
  hdata_npu_eb   = new TH1F("hdata_npu_eb",  "nPU EB",  7, xvt);
  hdata_npu_ee   = new TH1F("hdata_npu_ee",  "nPU EE",  7, xvt);
  hdata_nvtx_eb  = new TH1F("hdata_nvtx_eb", "nVtx EB", 7, xvt);
  hdata_nvtx_ee  = new TH1F("hdata_nvtx_ee", "nVtx EE", 7, xvt);

  hreco_bad_pt_eb    = new TH1F("hreco_bad_pt_eb",   "Pt EB",  10, xpt);
  hreco_bad_pt_ee    = new TH1F("hreco_bad_pt_ee",   "Pt EE",  10, xpt);
  hreco_bad_eta      = new TH1F("hreco_bad_eta",     "Eta",    20, -3, 3);
  hreco_bad_npu_eb   = new TH1F("hreco_bad_npu_eb",  "nPU EB",  7, xvt);
  hreco_bad_npu_ee   = new TH1F("hreco_bad_npu_ee",  "nPU EE",  7, xvt);
  hreco_bad_nvtx_eb  = new TH1F("hreco_bad_nvtx_eb", "nVtx EB", 7, xvt);
  hreco_bad_nvtx_ee  = new TH1F("hreco_bad_nvtx_ee", "nVtx EE", 7, xvt);

  hreco_pass_pt_eb    = new TH1F("hreco_pass_pt_eb",   "Pt EB",  10, xpt);
  hreco_pass_pt_ee    = new TH1F("hreco_pass_pt_ee",   "Pt EE",  10, xpt);
  hreco_pass_eta      = new TH1F("hreco_pass_eta",     "Eta",    20, -3, 3);
  hreco_pass_npu_eb   = new TH1F("hreco_pass_npu_eb",  "nPU EB",  7, xvt);
  hreco_pass_npu_ee   = new TH1F("hreco_pass_npu_ee",  "nPU EE",  7, xvt);
  hreco_pass_nvtx_eb  = new TH1F("hreco_pass_nvtx_eb", "nVtx EB", 7, xvt);
  hreco_pass_nvtx_ee  = new TH1F("hreco_pass_nvtx_ee", "nVtx EE", 7, xvt);
  */


  //=========================== lepton tag==============================
  hDPhi_gele_g = new TH1F("hDPhi_gele_g", "hDPhi_gele_g", 70, -3.5, 3.5);

  hDPhi_geleg_gele = new TH2F("hDPhi_geleg_gele", "hDPhi_geleg_gele", 70, -3.5, 3.5, 70,-3.5, 3.5);
  
  hcat4 = new TH1F("hcat4","hcat4",4,0,4);
  hDMZmin = new TH1F("hDMZmin","hDMZmin",500,0,50);
  hDRtoTrk = new TH1F("hDRtoTrk","hDRtoTrk", 300,0,30);
  hMassEleGamGam   = new TH1F("hMassEleGamGam", "M_{#gamma#gamma e}", 500, 0, 500);
  hMassMuGamGam   = new TH1F("hMassMuGamGam", "M_{#gamma#gamma #mu}", 500, 0, 500);
  hPtEleGamGam   = new TH1F("hPtEleGamGam", "p_{T}(#gamma#gamma e)", 500, 0, 500);
  hPtMuGamGam   = new TH1F("hPtMuGamGam", "p_{t}(#gamma#gamma #mu)", 500, 0, 500);


  hggMass_eleVeto = new TH1F("hggMass_eleVeto", "M_{#gamma#gamma}", 1000, 0, 500);
  hggMass_eleVeto_jetVeto = new TH1F("hggMass_eleVeto_jetVeto", "M_{#gamma#gamma}", 1000, 0, 500);

  hggMass_muTag     = new TH1F("hggMass_muTag", "M_{#gamma#gamma}", 1000, 0, 500);
  hggMass_eleTag     = new TH1F("hggMass_eleTag", "M_{#gamma#gamma}", 1000, 0, 500);
  hggMass_eletag_Zveto = new TH1F("hggMass_eleTag_Zveto", "M_{#gamma#gamma} Zveto", 1000, 0, 500);


  hggMassWH_eleTag_Zveto = new TH1F("hggMassWH_eleTag_Zveto", "M_{#gamma#gamma} Zveto", 1000, 0, 500);
  hggMassZH_eleTag_Zveto = new TH1F("hggMassZH_eleTag_Zveto", "M_{#gamma#gamma} Zveto", 1000, 0, 500);

  hggMassWH_muTag     = new TH1F("hggMassWH_muTag", "M_{#gamma#gamma}", 1000, 0,500);
  hggMassWH_eleTag     = new TH1F("hggMassWH_eleTag", "M_{#gamma#gamma}", 1000, 0, 500);
  hggMassZH_muTag     = new TH1F("hggMassZH_muTag", "M_{#gamma#gamma}", 1000, 0, 500);
  hggMassZH_eleTag     = new TH1F("hggMassZH_eleTag", "M_{#gamma#gamma}", 1000, 0, 500);
  hDR_Pho1Ele     = new TH1F("hDR_Pho1Ele", "#DeltaR(#gamma_{1},e)", 60, 0,6);
  hDR_Pho2Ele     = new TH1F("hDR_Pho2Ele", "#DeltaR(#gamma_{2},e)", 60, 0, 6);
  hDR_SCPho1Ele     = new TH1F("hDR_SCPho1Ele", "#DeltaR(#gamma_{1}(SC),e(SC))", 60, 0,6);
  hDR_SCPho2Ele     = new TH1F("hDR_SCPho2Ele", "#DeltaR(#gamma_{2}(SC),e(SC))", 60, 0, 6);
  hDRmin_PhoMu     = new TH1F("hDRmin_PhoMu", "#DeltaR_{min}(#gamma,#mu)", 60, 0, 6);
  hDRmin_PhoEle     = new TH1F("hDRmin_PhoEle", "#DeltaR_{min}(#gamma,e) respect to teh selected vertex", 60, 0, 6);
  hDRmin_SCPhoEle     = new TH1F("hDRmin_SCPhoEle", "#DeltaR_{min}(#gamma(SC),e(SC))", 60, 0, 6);
  hAllElePt  = new TH1F("hAllElePt", "p_{T}(e)", 200, 0, 200);
  hAllSCEleEta  = new TH1F("hAllSCEleEta", "#eta(e)", 60, -3, 3);
  hEleCombRelIso_EB =  new TH1F("hEleCombRelIso_EB", "CombRelIso(e) EB", 200, 0, 0.2);
  hEleCombRelIso_EE =  new TH1F("hEleCombRelIso_EE", "CombRelIso(e) EE", 200, 0, 0.2);

 /*  hPFRelIso_BID_B =  new TH1F("hPFRelIso_BID_B", "CombRelPFIso(e) EB", 200, 0, 0.2); */
/*   hPFRelIso_AID_B =  new TH1F("hPFRelIso_AID_B", "CombRelPFIso(e) EB", 200, 0, 0.2); */
/*   hPFRelIso_BID_EC =  new TH1F("hPFRelIso_BID_EC", "CombRelPFIso(e) EC", 200, 0, 0.2); */
/*   hPFRelIso_AID_EC =  new TH1F("hPFRelIso_AID_EC", "CombRelPFIso(e) EC", 200, 0, 0.2); */


  helePFRelIso_BID_B =  new TH1F("helePFRelIso_BID_B", "CombRelPFIso(e) EB", 200, 0, 0.2);
  helePFRelIso_AID_B =  new TH1F("helePFRelIso_AID_B", "CombRelPFIso(e) EB", 200, 0, 0.2);
  helePFRelIso_BID_EC =  new TH1F("helePFRelIso_BID_EC", "CombRelPFIso(e) EC", 200, 0, 0.2);
  helePFRelIso_AID_EC =  new TH1F("helePFRelIso_AID_EC", "CombRelPFIso(e) EC", 200, 0, 0.2);
  
  hmuPFRelIso_BID =  new TH1F("hmuPFRelIso_BID", "CombRelPFIso(#mu) EC", 200, 0, 0.2);
  hmuPFRelIso_AID =  new TH1F("hmuPFRelIso_AID", "CombRelPFIso(#mu) EC", 200, 0, 0.2);


  hSelEle1Pt  = new TH1F("hSelEle1Pt", "p_{T}(e)", 200, 0, 200);
  hSelSCEle1Eta  = new TH1F("hSelSCEle1Eta", "#eta(e)", 60, -3, 3);
  hSelEle2Pt  = new TH1F("hSelEle2Pt", "p_{T}(e)", 200, 0, 200);
  hSelSCEle2Eta  = new TH1F("hSelSCEle2Eta", "#eta(e)", 60, -3, 3);

  hAllMuPt  = new TH1F("hAllMuPt", "p_{T}(#mu)", 200, 0, 200);
  hAllMuEta  = new TH1F("hAllMuEta", "#eta(#mu)", 60, -3, 3);
  hMuCombRelIso =  new TH1F("hMuCombRelIso", "CombRelIso(#mu)" , 200, 0, 2);

  hSelMu1Pt  = new TH1F("hSelMu1Pt", "p_{T}(#mu)", 200, 0, 200);
  hSelMu1Eta  = new TH1F("hSelMu1Eta", "#eta(#mu)", 60, -3, 3);
  hSelMu2Pt  = new TH1F("hSelMu2Pt", "p_{T}(#mu)", 200, 0, 200);
  hSelMu2Eta  = new TH1F("hSelMu2Eta", "#eta(#mu)", 60, -3, 3);


  hMassEleEleGam1    = new TH1F("hMassEleEleGam1", "M(ee#gamma_{1})", 500, 0, 500);
  hMassEleEleGam2    = new TH1F("hMassEleEleGam2", "M(ee#gamma_{2})", 1000, 0, 500);
  hMassEleGam1    = new TH1F("hMassEleGam1", "M(e#gamma_{1})", 500, 0, 500);
  hMassEleGam2    = new TH1F("hMassEleGam2", "M(e#gamma_{2})", 500, 0, 500);


  hMassEleFakePho = new TH1F("hMassEleFakePho", "M(e gamma_{fake}) (GeV/c^{2})", 500, 0, 500)  ;
  hMassEleRealPho = new TH1F("hMassEleRealPho", "M(e gamma_{real}) (GeV/c^{2})", 500, 0, 500);

  hPtEleEleGam1    = new TH1F("hPtEleEleGam1", "p_{T}(ee#gamma_{1})", 500, 0, 500);
  hPtEleEleGam2    = new TH1F("hPtEleEleGam2", "p_{T}(ee#gamma_{2})", 1000, 0, 500);
  hPtEleGam1    = new TH1F("hPtEleGam1", "p_{T}(e#gamma_{1})", 500, 0, 500);
  hPtEleGam2    = new TH1F("hPtEleGam2", "p_{T}(e#gamma_{2})", 500, 0, 500);






  h_pfMET_lepTag = new TH1F("h_pfMET_lepTag", "PFMET lepTag", 200, 0, 200);
  h_tcMET_lepTag = new TH1F("h_tcMET_lepTag", "TCMET lepTag", 200, 0, 200);
  h_CaloMET_lepTag = new TH1F("h_CaloMET_lepTag", "Calo MET lepTag", 200, 0, 200);


  h_pfMET = new TH1F("h_pfMET", "PFMET ", 200, 0, 200);
  h_tcMET = new TH1F("h_tcMET", "TCMET ", 200, 0, 200);
  h_CaloMET = new TH1F("h_CaloMET", "Calo MET ", 200, 0, 200);
  //================================================
  hggMass     = new TH1F("hggMass", "M_{#gamma#gamma}", 1000, 0, 500);
  //==================================================

  hJetNHF =  new TH1F("hJetNHF", "hJetNHF", 100, 0, 1);
  hJetNEF =  new TH1F("hJetNEF", "hJetNEF", 100, 0, 1);
  hJetNConst =  new TH1F("hJetNConst", "hJetNConst", 10, 0, 10);

  hJetCHF =  new TH1F("hJetCHF", "hJetCHF", 100, 0, 1);
  hJetNCH =  new TH1F("hJetNCH", "hJetNCH", 100, 0, 10);
  hJetCEF =  new TH1F("hJetCEF", "hJetCEF", 100, 0, 1);



  hDRAllJetg1 = new TH1F("hDRAllJetg1", "hDRAllJetg1", 70, 0, 7);
  hDRAllJetg2 = new TH1F("hDRAllJetg2", "hDRAllJetg2", 70, 0, 7);
  hAllJetPt = new TH1F("hAllJetPt", "hAllJetPt", 200, 0, 200);
  hAllJetEta = new TH1F("hAllJetEta", "hAllJetEta", 100, -5, 5);

  hDRjet1g1   = new TH1F("hDRjet1g1", "hDRjet1g1", 70, 0, 7);
  hDRjet1g2   = new TH1F("hDRjet1g2", "hDRjet1g2", 70, 0, 7);
  hDRjet2g1   = new TH1F("hDRjet2g1", "hDRjet2g1", 70, 0, 7);
  hDRjet2g2   = new TH1F("hDRjet2g2", "hDRjet2g2", 70, 0, 7);

  hDPhijet1g1   = new TH1F("hDPhijet1g1", "hDPhijet1g1", 70, -3.5, 3.5);
  hDPhijet1g2   = new TH1F("hDPhijet1g2", "hDPhijet1g2", 70, -3.5, 3.5);
  hDPhijet2g1   = new TH1F("hDPhijet2g1", "hDPhijet2g1",70, -3.5, 3.5);
  hDPhijet2g2   = new TH1F("hDPhijet2g2", "hDPhijet2g2", 70, -3.5, 3.5);
  hDPhijet1jet2   = new TH1F("hDPhijet1jet2", "hDPhijet1jet2", 70, -3.5, 3.5);

  hggMass_twoJetTag     = new TH1F("hggMass_twoJetTag", "M_{#gamma#gamma}", 1000, 0, 1000);
  hggPt_twoJetTag    = new TH1F("hggPt_twoJetTag", "Pt_{#gamma#gamma}", 100, 0, 1000);
  hJet1Pt    = new TH1F("hJet1Pt", "p_{T}(jet1)", 300, 0, 300);
  hJet2Pt    = new TH1F("hJet2Pt", "p_{T}(jet2)", 300, 0, 300);
 
  hJet1Eta    = new TH1F("hJet1Eta", "#eta (jet1)", 100, -5, 5);
  hJet2Eta    = new TH1F("hJet2Eta", "#eta (jet2)", 100, -5, 5);


  hDeltaEtaJets    = new TH1F("hDeltaEtaJets", "#Delta #eta (jet1,jet2)", 80, 0, 8);
  hDeltaEtaJets_bc    = new TH1F("hDeltaEtaJets_bc", "#Delta #eta (jet1,jet2) bc", 80, 0, 8);

 
  hJetEtaProd_bc   = new TH1F("hJetEtaProd_bc", "hJetEtaProd_bc", 3, -1,2 );

  //hDeltaEtaJetsVsMgg    = new TH2F("hDeltaEtaJetsVsMgg", "#Delta #eta (jet1,jet2) vs M(#gamma #gamma)", 1000, 0, 1000, 80, 0, 8);
  //hMassJetJetVsMgg  = new TH2F("hMassJetJetVsMgg", " M(jet1,jet2) vs M(#gamma #gamma)", 1000, 0, 1000, 200, 0, 2000);

  hZeppenfeld = new TH1F("hZeppenfeld", "hZeppenfeld", 200, -10, 10);
  hZeppenfeld_bc = new TH1F("hZeppenfeld_bc", "hZeppenfeld bc", 200, -10, 10);

  hMassJetJet = new TH1F("hMassJetJet", "M(jet_{1},#gamma_{1}) [GeV]", 200, 0, 2000); 
  hMassJetJet_bc = new TH1F("hMassJetJet_bc", "M(jet_{1},#gamma_{1}) [GeV]  bc", 200, 0, 2000); 

  hDeltaPhiDiJetDiPho   = new TH1F("hDeltaPhiDiJetDiPho", "dPhi 2jet 2pho", 70, -3.5, 3.5);
  hDeltaPhiDiJetDiPho_bc = new TH1F("hDeltaPhiDiJetDiPho_bc", "dPhi 2jet 2pho bc", 70, -3.5, 3.5);



  //hMgg_twoJetTag_etaprod   = new TH1F("hMgg_twoJetTag_etaprod",  "M_{#gamma#gamma} #eta_{1} #times #eta_{2}<0", 1000, 0, 1000);
  hMgg_twoJetTag_jet1Pt    = new TH1F("hMgg_twoJetTag_jet1Pt",   "M_{#gamma#gamma} jet1 pt > 30",               1000, 0, 1000);
  hMgg_twoJetTag_deltaEta  = new TH1F("hMgg_twoJetTag_deltaEta", "M_{#gamma#gamma} #Delta #eta > 3.5",          1000, 0, 1000);
  hMgg_twoJetTag_zep       = new TH1F("hMgg_twoJetTag_zep",      "M_{#gamma#gamma} Zeppenfeld < 2.5",           1000, 0, 1000);
  hMgg_twoJetTag_Mjets     = new TH1F("hMgg_twoJetTag_Mjets",    "M_{#gamma#gamma} M(jet1,jet2) > 350",         1000, 0, 1000);
  hMgg_twoJetTag_dPhi      = new TH1F("hMgg_twoJetTag_dPhi",     "M_{#gamma#gamma} dPhi(jj,gg) > 2.6",          1000, 0, 1000);


  hJet1Pt_HSTRA_bc               = new TH1F("hJet1Pt_HSTRA_bc", "p_{T}(jet1) bc", 300, 0, 300);
  hDeltaEtaJets_HSTRA_bc         = new TH1F("hDeltaEtaJets_HSTRA_bc", "#Delta #eta (jet1,jet2) bc", 80, 0, 8);
  hZeppenfeld_HSTRA_bc           = new TH1F("hZeppenfeld_HSTRA_bc",          "hZeppenfeld bc",                            200, -10, 10);
  hMassJetJet_HSTRA_bc           = new TH1F("hMassJetJet_HSTRA_bc",          "M(jet_{1},#gamma_{1}) [GeV]  bc",           200, 0, 2000); 

  hMgg_HSTRA_twoJetTag           = new TH1F("hMgg_HSTRA_twoJetTag",          "M_{#gamma#gamma} pho1 pt > 60",            1000, 0, 1000);
  hMgg_HSTRA_twoJetTag_jet1Pt    = new TH1F("hMgg_HSTRA_twoJetTag_jet1Pt",   "M_{#gamma#gamma} jet1 pt > 30",            1000, 0, 1000);
  hMgg_HSTRA_twoJetTag_deltaEta  = new TH1F("hMgg_HSTRA_twoJetTag_deltaEta", "M_{#gamma#gamma} #Delta #eta < 2.5",       1000, 0, 1000);
  hMgg_HSTRA_twoJetTag_zep       = new TH1F("hMgg_HSTRA_twoJetTag_zep",      "M_{#gamma#gamma} Zeppenfeld < 1.5",        1000, 0, 1000);
  hMgg_HSTRA_twoJetTag_Mjets     = new TH1F("hMgg_HSTRA_twoJetTag_Mjets",    "M_{#gamma#gamma} 55 < M(jet1,jet2) < 115", 1000, 0, 1000);


}

void xAna::Init(TTree *tree) {
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   //vtxbsTkPxVec = 0;
   //vtxbsTkPyVec = 0;
   DiscriVBF_UseDiPhoPt = true ;
   DiscriVBF_UsePhoPt   = true ;
   DiscriVBF_useMvaSel  = false;
   vtxbsTkIndex = 0;
   vtxbsTkWeight = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("HLT", HLT, &b_HLT);
   fChain->SetBranchAddress("HLTIndex", HLTIndex, &b_HLTIndex);
   fChain->SetBranchAddress("bspotPos", bspotPos, &b_bspotPos);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("vtx", vtx, &b_vtx);
   fChain->SetBranchAddress("vtxNTrk", vtxNTrk, &b_vtxNTrk);
   fChain->SetBranchAddress("vtxNDF", vtxNDF, &b_vtxNDF);
   fChain->SetBranchAddress("vtxD0", vtxD0, &b_vtxD0);
   fChain->SetBranchAddress("IsVtxGood", &IsVtxGood, &b_IsVtxGood);
   fChain->SetBranchAddress("nVtxBS", &nVtxBS, &b_nVtxBS);
   fChain->SetBranchAddress("vtxbs", vtxbs, &b_vtxbs);
   fChain->SetBranchAddress("vtxbsPtMod", vtxbsPtMod, &b_vtxbsPtMod);
   fChain->SetBranchAddress("vtxbsSumPt2", vtxbsSumPt2, &b_vtxbsSumPt2);
   //fChain->SetBranchAddress("vtxbsTkPxVec", &vtxbsTkPxVec, &b_vtxbsTkPxVec);
   //fChain->SetBranchAddress("vtxbsTkPyVec", &vtxbsTkPyVec, &b_vtxbsTkPyVec);
   fChain->SetBranchAddress("vtxbsTkIndex", &vtxbsTkIndex, &b_vtxbsTkIndex);
   fChain->SetBranchAddress("vtxbsTkWeight", &vtxbsTkWeight, &b_vtxbsTkWeight);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("trkP", trkP, &b_trkP);
   fChain->SetBranchAddress("trkVtx", trkVtx, &b_trkVtx);
   fChain->SetBranchAddress("trkd0", trkd0, &b_trkd0);
   fChain->SetBranchAddress("trkd0Err", trkd0Err, &b_trkd0Err);
   fChain->SetBranchAddress("trkdz", trkdz, &b_trkdz);
   fChain->SetBranchAddress("trkdzErr", trkdzErr, &b_trkdzErr);
   fChain->SetBranchAddress("trkPtErr", trkPtErr, &b_trkPtErr);
   fChain->SetBranchAddress("trkQuality", trkQuality, &b_trkQuality);
   fChain->SetBranchAddress("nGoodTrk", &nGoodTrk, &b_nGoodTrk);
   fChain->SetBranchAddress("IsTracksGood", &IsTracksGood, &b_IsTracksGood);
   fChain->SetBranchAddress("pdf", pdf, &b_pdf);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("processID", &processID, &b_processID);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcVtx", mcVtx, &b_mcVtx);
   fChain->SetBranchAddress("mcPt", mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMass", mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcEta", mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcGMomPID", mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcMomPID", mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcMomPt", mcMomPt, &b_mcMomPt);
   fChain->SetBranchAddress("mcMomMass", mcMomMass, &b_mcMomMass);
   fChain->SetBranchAddress("mcMomEta", mcMomEta, &b_mcMomEta);
   fChain->SetBranchAddress("mcMomPhi", mcMomPhi, &b_mcMomPhi);
   fChain->SetBranchAddress("mcIndex", mcIndex, &b_mcIndex);
   fChain->SetBranchAddress("mcDecayType", mcDecayType, &b_mcDecayType);
   fChain->SetBranchAddress("mcCalIsoDR03", mcCalIsoDR03, &b_mcCalIsoDR03);
   fChain->SetBranchAddress("mcTrkIsoDR03", mcTrkIsoDR03, &b_mcTrkIsoDR03);
   fChain->SetBranchAddress("mcCalIsoDR04", mcCalIsoDR04, &b_mcCalIsoDR04);
   fChain->SetBranchAddress("mcTrkIsoDR04", mcTrkIsoDR04, &b_mcTrkIsoDR04);
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   /* fChain->SetBranchAddress("genMETx", &genMETx, &b_genMETx); */
   /* fChain->SetBranchAddress("genMETy", &genMETy, &b_genMETy); */
   fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", puTrue, &b_puTrue);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   /* fChain->SetBranchAddress("METx", &METx, &b_METx); */
   /* fChain->SetBranchAddress("METy", &METy, &b_METy); */
   fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
   fChain->SetBranchAddress("METsumEt", &METsumEt, &b_METsumEt);
   /* fChain->SetBranchAddress("tcMET", &tcMET, &b_tcMET); */
   /* fChain->SetBranchAddress("tcMETx", &tcMETx, &b_tcMETx); */
   /* fChain->SetBranchAddress("tcMETy", &tcMETy, &b_tcMETy); */
   /* fChain->SetBranchAddress("tcMETPhi", &tcMETPhi, &b_tcMETPhi); */
   /* fChain->SetBranchAddress("tcMETsumEt", &tcMETsumEt, &b_tcMETsumEt); */
   /* fChain->SetBranchAddress("tcMETmEtSig", &tcMETmEtSig, &b_tcMETmEtSig); */
   /* fChain->SetBranchAddress("tcMETSig", &tcMETSig, &b_tcMETSig); */
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   /* fChain->SetBranchAddress("pfMETx", &pfMETx, &b_pfMETx); */
   /* fChain->SetBranchAddress("pfMETy", &pfMETy, &b_pfMETy); */
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("recoPfMET", &recoPfMET, &b_recoPfMET);
   /* fChain->SetBranchAddress("recoPfMETx", &recoPfMETx, &b_recoPfMETx); */
   /* fChain->SetBranchAddress("recoPfMETy", &recoPfMETy, &b_recoPfMETy); */
   fChain->SetBranchAddress("recoPfMETPhi", &recoPfMETPhi, &b_recoPfMETPhi);
   fChain->SetBranchAddress("recoPfMETsumEt", &recoPfMETsumEt, &b_recoPfMETsumEt);
   fChain->SetBranchAddress("recoPfMETmEtSig", &recoPfMETmEtSig, &b_recoPfMETmEtSig);
   fChain->SetBranchAddress("recoPfMETSig", &recoPfMETSig, &b_recoPfMETSig);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleTrg", eleTrg, &b_eleTrg);
   fChain->SetBranchAddress("eleID", eleID, &b_eleID);
   fChain->SetBranchAddress("eleClass", eleClass, &b_eleClass);
   fChain->SetBranchAddress("eleCharge", eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleEn", eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleSCRawEn", eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEn", eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("elePt", elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleEtaVtx", eleEtaVtx, &b_eleEtaVtx);
   fChain->SetBranchAddress("elePhiVtx", elePhiVtx, &b_elePhiVtx);
   fChain->SetBranchAddress("eleSCEta", eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCEtaWidth", eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleVtx", eleVtx, &b_eleVtx);
   fChain->SetBranchAddress("eleHoverE", eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("elePin", elePin, &b_elePin);
   fChain->SetBranchAddress("elePout", elePout, &b_elePout);
   fChain->SetBranchAddress("eleBrem", eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSigmaIEtaIPhi", eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
   fChain->SetBranchAddress("eleSigmaIPhiIPhi", eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
   fChain->SetBranchAddress("eleD0", eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleD0Vtx", eleD0Vtx, &b_eleD0Vtx);
   fChain->SetBranchAddress("eleDzVtx", eleDzVtx, &b_eleDzVtx);
   fChain->SetBranchAddress("eleEmax", eleEmax, &b_eleEmax);
   fChain->SetBranchAddress("eleE3x3", eleE3x3, &b_eleE3x3);
   fChain->SetBranchAddress("eleE5x5", eleE5x5, &b_eleE5x5);
   fChain->SetBranchAddress("eleE2x5Right", eleE2x5Right, &b_eleE2x5Right);
   fChain->SetBranchAddress("eleE2x5Left", eleE2x5Left, &b_eleE2x5Left);
   fChain->SetBranchAddress("eleE2x5Top", eleE2x5Top, &b_eleE2x5Top);
   fChain->SetBranchAddress("eleE2x5Bottom", eleE2x5Bottom, &b_eleE2x5Bottom);
   fChain->SetBranchAddress("eleSeedTime", eleSeedTime, &b_eleSeedTime);
   fChain->SetBranchAddress("eleRecoFlag", eleRecoFlag, &b_eleRecoFlag);
   fChain->SetBranchAddress("eleGenIndex", eleGenIndex, &b_eleGenIndex);
   fChain->SetBranchAddress("eleGenGMomPID", eleGenGMomPID, &b_eleGenGMomPID);
   fChain->SetBranchAddress("eleGenMomPID", eleGenMomPID, &b_eleGenMomPID);
   fChain->SetBranchAddress("eleGenMomPt", eleGenMomPt, &b_eleGenMomPt);
   fChain->SetBranchAddress("eleIsoTrkDR03", eleIsoTrkDR03, &b_eleIsoTrkDR03);
   fChain->SetBranchAddress("eleIsoEcalDR03", eleIsoEcalDR03, &b_eleIsoEcalDR03);
   fChain->SetBranchAddress("eleIsoHcalDR03", eleIsoHcalDR03, &b_eleIsoHcalDR03);
   fChain->SetBranchAddress("eleIsoHcalSolidDR03", eleIsoHcalSolidDR03, &b_eleIsoHcalSolidDR03);
   fChain->SetBranchAddress("eleIsoTrkDR04", eleIsoTrkDR04, &b_eleIsoTrkDR04);
   fChain->SetBranchAddress("eleIsoEcalDR04", eleIsoEcalDR04, &b_eleIsoEcalDR04);
   fChain->SetBranchAddress("eleIsoHcalDR04", eleIsoHcalDR04, &b_eleIsoHcalDR04);
   fChain->SetBranchAddress("eleIsoHcalSolidDR04", eleIsoHcalSolidDR04, &b_eleIsoHcalSolidDR04);
   fChain->SetBranchAddress("eleMissHits", eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleConvDist", eleConvDist, &b_eleConvDist);
   fChain->SetBranchAddress("eleConvDcot", eleConvDcot, &b_eleConvDcot);
   fChain->SetBranchAddress("elePFChIso03", elePFChIso03, &b_elePFChIso03);
   fChain->SetBranchAddress("elePFPhoIso03", elePFPhoIso03, &b_elePFPhoIso03);
   fChain->SetBranchAddress("elePFNeuIso03", elePFNeuIso03, &b_elePFNeuIso03);
   fChain->SetBranchAddress("eleIDMVANonTrig", eleIDMVANonTrig, &b_eleIDMVANonTrig);

   fChain->SetBranchAddress("eleEcalEn", eleEcalEn, &b_eleEcalEn);
   fChain->SetBranchAddress("eleConvVtxFit", eleConvVtxFit, &b_eleConvVtxFit);
   fChain->SetBranchAddress("eleESEffSigmaRR", eleESEffSigmaRR, &b_eleESEffSigmaRR);
   /* fChain->SetBranchAddress("eleESE1", eleESE1, &b_eleESE1); */
   /* fChain->SetBranchAddress("eleESE3", eleESE3, &b_eleESE3); */
   /* fChain->SetBranchAddress("eleESE5", eleESE5, &b_eleESE5); */
   /* fChain->SetBranchAddress("eleESE7", eleESE7, &b_eleESE7); */
   /* fChain->SetBranchAddress("eleESE11", eleESE11, &b_eleESE11); */
   /* fChain->SetBranchAddress("eleESE21", eleESE21, &b_eleESE21); */
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoTrg", phoTrg, &b_phoTrg);
   fChain->SetBranchAddress("phoIsPhoton", phoIsPhoton, &b_phoIsPhoton);
   fChain->SetBranchAddress("phoSCPos", phoSCPos, &b_phoSCPos);
   fChain->SetBranchAddress("phoCaloPos", phoCaloPos, &b_phoCaloPos);
   fChain->SetBranchAddress("phoE", phoE, &b_phoE);
   fChain->SetBranchAddress("phoRegrE", phoRegrE, &b_phoRegrE);
   fChain->SetBranchAddress("phoRegrEerr", phoRegrErr, &b_phoRegrErr);
   fChain->SetBranchAddress("phoEt", phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoVtx", phoVtx, &b_phoVtx);
   fChain->SetBranchAddress("phoPhi", phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoEtVtx", phoEtVtx, &b_phoEtVtx);
   fChain->SetBranchAddress("phoEtaVtx", phoEtaVtx, &b_phoEtaVtx);
   fChain->SetBranchAddress("phoPhiVtx", phoPhiVtx, &b_phoPhiVtx);
   fChain->SetBranchAddress("phoR9", phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoCetaCorrE", phoCetaCorrE, &b_phoCetaCorrE);
   fChain->SetBranchAddress("phoCetaCorrEt", phoCetaCorrEt, &b_phoCetaCorrEt);
   fChain->SetBranchAddress("phoBremCorrE", phoBremCorrE, &b_phoBremCorrE);
   fChain->SetBranchAddress("phoBremCorrEt", phoBremCorrEt, &b_phoBremCorrEt);
   fChain->SetBranchAddress("phoFullCorrE", phoFullCorrE, &b_phoFullCorrE);
   fChain->SetBranchAddress("phoFullCorrEt", phoFullCorrEt, &b_phoFullCorrEt);
   fChain->SetBranchAddress("phoTrkIsoSolidDR03", phoTrkIsoSolidDR03, &b_phoTrkIsoSolidDR03);
   fChain->SetBranchAddress("phoTrkIsoHollowDR03", phoTrkIsoHollowDR03, &b_phoTrkIsoHollowDR03);
   fChain->SetBranchAddress("phoEcalIsoDR03", phoEcalIsoDR03, &b_phoEcalIsoDR03);
   fChain->SetBranchAddress("phoHcalIsoDR03", phoHcalIsoDR03, &b_phoHcalIsoDR03);
   fChain->SetBranchAddress("phoHcalIsoSolidDR03", phoHcalIsoSolidDR03, &b_phoHcalIsoSolidDR03);
   fChain->SetBranchAddress("phoTrkIsoSolidDR04", phoTrkIsoSolidDR04, &b_phoTrkIsoSolidDR04);
   fChain->SetBranchAddress("phoTrkIsoHollowDR04", phoTrkIsoHollowDR04, &b_phoTrkIsoHollowDR04);
   fChain->SetBranchAddress("phoTrkIsoHollowNoDzDR03", phoTrkIsoHollowNoDzDR03, &b_phoTrkIsoHollowNoDzDR03);
   fChain->SetBranchAddress("phoTrkIsoHollowNoDzDR04", phoTrkIsoHollowNoDzDR04, &b_phoTrkIsoHollowNoDzDR04);
   fChain->SetBranchAddress("phoCiCTrkIsoDR03", phoCiCTrkIsoDR03, &b_phoCiCTrkIsoDR03);
   fChain->SetBranchAddress("phoCiCTrkIsoDR04", phoCiCTrkIsoDR04, &b_phoCiCTrkIsoDR04);
   //fChain->SetBranchAddress("phoCiCTrkIsoDR03Vtx", phoCiCTrkIsoDR03Vtx, &b_phoCiCTrkIsoDR03Vtx);
   //fChain->SetBranchAddress("phoCiCTrkIsoDR04Vtx", phoCiCTrkIsoDR04Vtx, &b_phoCiCTrkIsoDR04Vtx);
   fChain->SetBranchAddress("phoCiCdRtoTrk", phoCiCdRtoTrk, &b_phoCiCdRtoTrk);
   fChain->SetBranchAddress("phoEcalIsoDR04", phoEcalIsoDR04, &b_phoEcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoDR04", phoHcalIsoDR04, &b_phoHcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoSolidDR04", phoHcalIsoSolidDR04, &b_phoHcalIsoSolidDR04);

   fChain->SetBranchAddress("phoCiCPF4phopfIso03", phoCiCPF4phopfIso03, &b_phoCiCPF4phopfIso03);
   fChain->SetBranchAddress("phoCiCPF4phopfIso04", phoCiCPF4phopfIso04, &b_phoCiCPF4phopfIso04);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso02", phoCiCPF4chgpfIso02, &b_phoCiCPF4chgpfIso02);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso03", phoCiCPF4chgpfIso03, &b_phoCiCPF4chgpfIso03);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso04", phoCiCPF4chgpfIso04, &b_phoCiCPF4chgpfIso04);

   fChain->SetBranchAddress("phoCiCPF4phopfIso03AppVeto", phoCiCPF4phopfIso03AppVeto, &b_phoCiCPF4phopfIso03AppVeto);
   fChain->SetBranchAddress("phoCiCPF4phopfIso04AppVeto", phoCiCPF4phopfIso04AppVeto, &b_phoCiCPF4phopfIso04AppVeto);

   fChain->SetBranchAddress("phoCiCPF4chgpfIso03AppVeto", phoCiCPF4chgpfIso03AppVeto, &b_phoCiCPF4chgpfIso03AppVeto);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso04AppVeto", phoCiCPF4chgpfIso04AppVeto, &b_phoCiCPF4chgpfIso04AppVeto);

   fChain->SetBranchAddress("phoCiCPF4phopfIso03AppVetoClean", phoCiCPF4phopfIso03AppVetoClean, &b_phoCiCPF4phopfIso03AppVetoClean);
   fChain->SetBranchAddress("phoCiCPF4phopfIso04AppVetoClean", phoCiCPF4phopfIso04AppVetoClean, &b_phoCiCPF4phopfIso04AppVetoClean);

   fChain->SetBranchAddress("phoCiCPF4chgpfIso03Mod", phoCiCPF4chgpfIso03Mod, &b_phoCiCPF4chgpfIso03Mod);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso04Mod", phoCiCPF4chgpfIso04Mod, &b_phoCiCPF4chgpfIso04Mod);

   fChain->SetBranchAddress("phoEleVeto",phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoHoverE", phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoSigmaIEtaIPhi", phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
   fChain->SetBranchAddress("phoSigmaIPhiIPhi", phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
   fChain->SetBranchAddress("phoEmax", phoEmax, &b_phoEmax);
   fChain->SetBranchAddress("phoE2x2", phoE2x2, &b_phoE2x2);
   fChain->SetBranchAddress("phoE3x3", phoE3x3, &b_phoE3x3);
   fChain->SetBranchAddress("phoE5x5", phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoE2x5Right", phoE2x5Right, &b_phoE2x5Right);
   fChain->SetBranchAddress("phoE2x5Left", phoE2x5Left, &b_phoE2x5Left);
   fChain->SetBranchAddress("phoE2x5Top", phoE2x5Top, &b_phoE2x5Top);
   fChain->SetBranchAddress("phoE2x5Bottom", phoE2x5Bottom, &b_phoE2x5Bottom);
   fChain->SetBranchAddress("phoSeedTime", phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedDetId1", phoSeedDetId1, &b_phoSeedDetId1);
   fChain->SetBranchAddress("phoSeedDetId2", phoSeedDetId2, &b_phoSeedDetId2);
   fChain->SetBranchAddress("phoRecoFlag", phoRecoFlag, &b_phoRecoFlag);
   fChain->SetBranchAddress("phoPos", phoPos, &b_phoPos);
   fChain->SetBranchAddress("phoGenIndex", phoGenIndex, &b_phoGenIndex);
   fChain->SetBranchAddress("phoGenGMomPID", phoGenGMomPID, &b_phoGenGMomPID);
   fChain->SetBranchAddress("phoGenMomPID", phoGenMomPID, &b_phoGenMomPID);
   fChain->SetBranchAddress("phoGenMomPt", phoGenMomPt, &b_phoGenMomPt);
   fChain->SetBranchAddress("phoSCE", phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoSCEt", phoSCEt, &b_phoSCEt);
   fChain->SetBranchAddress("phoSCEta", phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phoOverlap", phoOverlap, &b_phoOverlap);
   fChain->SetBranchAddress("phohasPixelSeed", phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoIsConv", phoIsConv, &b_phoIsConv);
   fChain->SetBranchAddress("phoNConv", phoNConv, &b_phoNConv);
   fChain->SetBranchAddress("phoConvInvMass", phoConvInvMass, &b_phoConvInvMass);
   fChain->SetBranchAddress("phoConvCotTheta", phoConvCotTheta, &b_phoConvCotTheta);
   fChain->SetBranchAddress("phoConvEoverP", phoConvEoverP, &b_phoConvEoverP);
   fChain->SetBranchAddress("phoConvZofPVfromTrks", phoConvZofPVfromTrks, &b_phoConvZofPVfromTrks);
   fChain->SetBranchAddress("phoConvMinDist", phoConvMinDist, &b_phoConvMinDist);
   fChain->SetBranchAddress("phoConvdPhiAtVtx", phoConvdPhiAtVtx, &b_phoConvdPhiAtVtx);
   fChain->SetBranchAddress("phoConvdPhiAtCalo", phoConvdPhiAtCalo, &b_phoConvdPhiAtCalo);
   fChain->SetBranchAddress("phoConvdEtaAtCalo", phoConvdEtaAtCalo, &b_phoConvdEtaAtCalo);
   fChain->SetBranchAddress("phoConvTrkd0", phoConvTrkd0, &b_phoConvTrkd0);
   fChain->SetBranchAddress("phoConvTrkPin", phoConvTrkPin, &b_phoConvTrkPin);
   fChain->SetBranchAddress("phoConvTrkPout", phoConvTrkPout, &b_phoConvTrkPout);
   fChain->SetBranchAddress("phoConvTrkdz", phoConvTrkdz, &b_phoConvTrkdz);
   fChain->SetBranchAddress("phoConvTrkdzErr", phoConvTrkdzErr, &b_phoConvTrkdzErr);
   fChain->SetBranchAddress("phoConvChi2", phoConvChi2, &b_phoConvChi2);
   fChain->SetBranchAddress("phoConvChi2Prob", phoConvChi2Prob, &b_phoConvChi2Prob);
   fChain->SetBranchAddress("phoConvNTrks", phoConvNTrks, &b_phoConvNTrks);
   fChain->SetBranchAddress("phoConvCharge", phoConvCharge, &b_phoConvCharge);
   fChain->SetBranchAddress("phoConvValidVtx", phoConvValidVtx, &b_phoConvValidVtx);
   fChain->SetBranchAddress("phoConvLikeLihood", phoConvLikeLihood, &b_phoConvLikeLihood);
   fChain->SetBranchAddress("phoConvP4", phoConvP4, &b_phoConvP4);
   fChain->SetBranchAddress("phoConvVtx", phoConvVtx, &b_phoConvVtx);
   fChain->SetBranchAddress("phoConvVtxErr", phoConvVtxErr, &b_phoConvVtxErr);
   fChain->SetBranchAddress("phoConvPairMomentum", phoConvPairMomentum, &b_phoConvPairMomentum);
   fChain->SetBranchAddress("phoConvRefittedMomentum", phoConvRefittedMomentum, &b_phoConvRefittedMomentum);
   fChain->SetBranchAddress("phoESEffSigmaRR", phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("phoESE1", phoESE1, &b_phoESE1);
   fChain->SetBranchAddress("phoESE3", phoESE3, &b_phoESE3);
   fChain->SetBranchAddress("phoESE5", phoESE5, &b_phoESE5);
   fChain->SetBranchAddress("phoESE7", phoESE7, &b_phoESE7);
   fChain->SetBranchAddress("phoESE11", phoESE11, &b_phoESE11);
   fChain->SetBranchAddress("phoESE21", phoESE21, &b_phoESE21);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muTrg", muTrg, &b_muTrg);
   fChain->SetBranchAddress("muEta", muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", muCharge, &b_muCharge);
   fChain->SetBranchAddress("muPt", muPt, &b_muPt);
   fChain->SetBranchAddress("muPz", muPz, &b_muPz);
   fChain->SetBranchAddress("muGenIndex", muGenIndex, &b_muGenIndex);
   fChain->SetBranchAddress("muIsoTrk", muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muIsoCalo", muIsoCalo, &b_muIsoCalo);
   fChain->SetBranchAddress("muIsoEcal", muIsoEcal, &b_muIsoEcal);
   fChain->SetBranchAddress("muIsoHcal", muIsoHcal, &b_muIsoHcal);
   fChain->SetBranchAddress("muChi2NDF", muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muEmVeto", muEmVeto, &b_muEmVeto);
   fChain->SetBranchAddress("muHadVeto", muHadVeto, &b_muHadVeto);
   fChain->SetBranchAddress("muType", muType, &b_muType);
   fChain->SetBranchAddress("muID", muID, &b_muID);
   fChain->SetBranchAddress("muD0", muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", muDz, &b_muDz);
   fChain->SetBranchAddress("muVtx", muVtx, &b_muVtx);
   fChain->SetBranchAddress("muD0Vtx", muD0Vtx, &b_muD0Vtx);
   fChain->SetBranchAddress("muDzVtx", muDzVtx, &b_muDzVtx);
   fChain->SetBranchAddress("muNumberOfValidTrkHits", muNumberOfValidTrkHits, &b_muNumberOfValidTrkHits);
   fChain->SetBranchAddress("muNumberOfValidPixelHits", muNumberOfValidPixelHits, &b_muNumberOfValidPixelHits);
   fChain->SetBranchAddress("muNumberOfValidMuonHits", muNumberOfValidMuonHits, &b_muNumberOfValidMuonHits);
   fChain->SetBranchAddress("muStations", muStations, &b_muStations);
   fChain->SetBranchAddress("muChambers", muChambers, &b_muChambers);

   //-------2012----------
   fChain->SetBranchAddress("muNumberOfValidTrkLayers", muNumberOfValidTrkLayers,&b_muNumberOfValidTrkLayers);
   fChain->SetBranchAddress("muPFIsoR04_CH",muPFIsoR04_CH,&b_muPFIsoR04_CH);
   fChain->SetBranchAddress("muPFIsoR04_NH",muPFIsoR04_NH,&b_muPFIsoR04_NH);
   fChain->SetBranchAddress("muPFIsoR04_Pho",muPFIsoR04_Pho,&b_muPFIsoR04_Pho);
   fChain->SetBranchAddress("muPFIsoR04_PU",muPFIsoR04_PU,&b_muPFIsoR04_PU);
   fChain->SetBranchAddress("muPFIsoR04_NHHT",muPFIsoR04_NHHT,&b_muPFIsoR04_NHHT);
   fChain->SetBranchAddress("muPFIsoR04_PhoHT",muPFIsoR04_PhoHT,&b_muPFIsoR04_PhoHT);
   fChain->SetBranchAddress("muPFIsoR03_CH",muPFIsoR03_CH,&b_muPFIsoR03_CH);
   fChain->SetBranchAddress("muPFIsoR03_NH",muPFIsoR03_NH,&b_muPFIsoR03_NH);
   fChain->SetBranchAddress("muPFIsoR03_Pho",muPFIsoR03_Pho,&b_muPFIsoR03_Pho);
   fChain->SetBranchAddress("muPFIsoR03_PU",muPFIsoR03_PU,&b_muPFIsoR03_PU);
   fChain->SetBranchAddress("muPFIsoR03_NHHT",muPFIsoR03_NHHT,&b_muPFIsoR03_NHHT);
   fChain->SetBranchAddress("muPFIsoR03_PhoHT",muPFIsoR03_PhoHT,&b_muPFIsoR03_PhoHT);


   fChain->SetBranchAddress("nPFPho_", &nPFPho_, &b_nPFPho_);
   fChain->SetBranchAddress("isIsoPFPho_", isIsoPFPho_, &b_isIsoPFPho_);
   fChain->SetBranchAddress("PFphoE_", PFphoE_, &b_PFphoE_);
   fChain->SetBranchAddress("PFphoEt_", PFphoEt_, &b_PFphoEt_);
   fChain->SetBranchAddress("PFphoEta_", PFphoEta_, &b_PFphoEta_);
   fChain->SetBranchAddress("PFphoPhi_", PFphoPhi_, &b_PFphoPhi_);
   fChain->SetBranchAddress("PFphoIsConv_", PFphoIsConv_, &b_PFphoIsConv_);
   fChain->SetBranchAddress("PFphoTrkIsoHollowDR04_", PFphoTrkIsoHollowDR04_, &b_PFphoTrkIsoHollowDR04_);
   fChain->SetBranchAddress("PFphoEcalIsoDR04_", PFphoEcalIsoDR04_, &b_PFphoEcalIsoDR04_);
   fChain->SetBranchAddress("PFphoHcalIsoDR04_", PFphoHcalIsoDR04_, &b_PFphoHcalIsoDR04_);
   fChain->SetBranchAddress("PFphoHoverE_", PFphoHoverE_, &b_PFphoHoverE_);
   fChain->SetBranchAddress("nPFElec_", &nPFElec_, &b_nPFElec_);
   fChain->SetBranchAddress("PFelePt_", PFelePt_, &b_PFelePt_);
   fChain->SetBranchAddress("PFeleEta_", PFeleEta_, &b_PFeleEta_);
   fChain->SetBranchAddress("PFelePhi_", PFelePhi_, &b_PFelePhi_);
   fChain->SetBranchAddress("PFeleEn_", PFeleEn_, &b_PFeleEn_);
   fChain->SetBranchAddress("PFeleCharge", PFeleCharge, &b_PFeleCharge);
   fChain->SetBranchAddress("nPFMu_", &nPFMu_, &b_nPFMu_);
   fChain->SetBranchAddress("PFmuCharge_", PFmuCharge_, &b_PFmuCharge_);
   fChain->SetBranchAddress("PFmuPhi_", PFmuPhi_, &b_PFmuPhi_);
   fChain->SetBranchAddress("PFmuEta_", PFmuEta_, &b_PFmuEta_);
   fChain->SetBranchAddress("PFmuPt_", PFmuPt_, &b_PFmuPt_);
   fChain->SetBranchAddress("PFmuPz_", PFmuPz_, &b_PFmuPz_);
   fChain->SetBranchAddress("PFmuChambers_", PFmuChambers_, &b_PFmuChambers_);
   fChain->SetBranchAddress("PFmuD0_", PFmuD0_, &b_PFmuD0_);
   fChain->SetBranchAddress("PFmuDz_", PFmuDz_, &b_PFmuDz_);
   fChain->SetBranchAddress("PFmuChi2NDF_", PFmuChi2NDF_, &b_PFmuChi2NDF_);
   fChain->SetBranchAddress("PFmuNumberOfValidTrkHits_", PFmuNumberOfValidTrkHits_, &b_PFmuNumberOfValidTrkHits_);
   fChain->SetBranchAddress("PFmuNumberOfValidPixelHits_", PFmuNumberOfValidPixelHits_, &b_PFmuNumberOfValidPixelHits_);
   fChain->SetBranchAddress("PFmuNumberOfValidMuonHits_", PFmuNumberOfValidMuonHits_, &b_PFmuNumberOfValidMuonHits_);
   fChain->SetBranchAddress("nPFchad_", &nPFchad_, &b_nPFchad_);
   fChain->SetBranchAddress("PFhad_charge_", &PFhad_charge_, &b_PFhad_charge_);
   fChain->SetBranchAddress("PFchad_Eta_", &PFchad_Eta_, &b_PFchad_Eta_);
   fChain->SetBranchAddress("PFchad_Phi_", &PFchad_Phi_, &b_PFchad_Phi_);
   fChain->SetBranchAddress("PFchad_Pt_", &PFchad_Pt_, &b_PFchad_Pt_);
   fChain->SetBranchAddress("PFchad_P_", &PFchad_P_, &b_PFchad_P_);
   fChain->SetBranchAddress("PFchad_E_", &PFchad_E_, &b_PFchad_E_);
   fChain->SetBranchAddress("PFchad_Ecal_", &PFchad_Ecal_, &b_PFchad_Ecal_);
   fChain->SetBranchAddress("PFchad_Hcal_", &PFchad_Hcal_, &b_PFchad_Hcal_);
   fChain->SetBranchAddress("nPFnhad_", &nPFnhad_, &b_nPFnhad_);
   fChain->SetBranchAddress("PFnhad_Eta_", PFnhad_Eta_, &b_PFnhad_Eta_);
   fChain->SetBranchAddress("PFnhad_Phi_", PFnhad_Phi_, &b_PFnhad_Phi_);
   fChain->SetBranchAddress("PFnhad_Pt_", PFnhad_Pt_, &b_PFnhad_Pt_);
   fChain->SetBranchAddress("PFnhad_P_", PFnhad_P_, &b_PFnhad_P_);
   fChain->SetBranchAddress("PFnhad_E_", PFnhad_E_, &b_PFnhad_E_);
   fChain->SetBranchAddress("PFnhad_Hcal_", PFnhad_Hcal_, &b_PFnhad_Hcal_);
   fChain->SetBranchAddress("rho2012", &rho2012, &b_rho2012);
   fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
   fChain->SetBranchAddress("rho44", &rho44, &b_rho44);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetTrg", jetTrg, &b_jetTrg);
   fChain->SetBranchAddress("jetEn", jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetEt", jetEt, &b_jetEt);
   fChain->SetBranchAddress("jetRawPt", jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetArea", jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetCHF", jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetMEF", jetMEF, &b_jetMEF);
   fChain->SetBranchAddress("jetNCH", jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetHFHAE", jetHFHAE, &b_jetHFHAE);
   fChain->SetBranchAddress("jetHFEME", jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConstituents", jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("jetLeadTrackPt", jetLeadTrackPt, &b_jetLeadTrackPt);
   fChain->SetBranchAddress("jetSoftLeptPt", jetSoftLeptPt, &b_jetSoftLeptPt);
   fChain->SetBranchAddress("jetSoftLeptPtRel", jetSoftLeptPtRel, &b_jetSoftLeptPtRel);
   fChain->SetBranchAddress("jetSoftLeptdR", jetSoftLeptdR, &b_jetSoftLeptdR);
   fChain->SetBranchAddress("jetVtx3deL", jetVtx3deL, &b_jetVtx3deL);
   fChain->SetBranchAddress("jetVtxMass", jetVtxMass, &b_jetVtxMass);
   fChain->SetBranchAddress("jetTrackCountHiEffBJetTags", jetTrackCountHiEffBJetTags, &b_jetTrackCountHiEffBJetTags);
   fChain->SetBranchAddress("jetTrackCountHiPurBJetTags", jetTrackCountHiPurBJetTags, &b_jetTrackCountHiPurBJetTags);
   fChain->SetBranchAddress("jetSimpleSVHiEffBJetTags", jetSimpleSVHiEffBJetTags, &b_jetSimpleSVHiEffBJetTags);
   fChain->SetBranchAddress("jetSimpleSVHiPurBJetTags", jetSimpleSVHiPurBJetTags, &b_jetSimpleSVHiPurBJetTags);
   fChain->SetBranchAddress("jetCombinedSecondaryVtxBJetTags", jetCombinedSecondaryVtxBJetTags, &b_jetCombinedSecondaryVtxBJetTags);
   fChain->SetBranchAddress("jetCombinedSecondaryVtxMVABJetTags", jetCombinedSecondaryVtxMVABJetTags, &b_jetCombinedSecondaryVtxMVABJetTags);
   fChain->SetBranchAddress("jetJetProbabilityBJetTags", jetJetProbabilityBJetTags, &b_jetJetProbabilityBJetTags);
   fChain->SetBranchAddress("jetJetBProbabilityBJetTags", jetJetBProbabilityBJetTags, &b_jetJetBProbabilityBJetTags);
   fChain->SetBranchAddress("jetTrackCountingHighPurBJetTags", jetTrackCountingHighPurBJetTags, &b_jetTrackCountingHighPurBJetTags);
   fChain->SetBranchAddress("jetJECUnc", jetJECUnc, &b_jetJECUnc);
   fChain->SetBranchAddress("jetPartonID", jetPartonID, &b_jetPartonID);
   fChain->SetBranchAddress("jetGenJetIndex", jetGenJetIndex, &b_jetGenJetIndex);
   fChain->SetBranchAddress("jetGenJetEn", jetGenJetEn, &b_jetGenJetEn);
   fChain->SetBranchAddress("jetGenJetPt", jetGenJetPt, &b_jetGenJetPt);
   fChain->SetBranchAddress("jetGenJetEta", jetGenJetEta, &b_jetGenJetEta);
   fChain->SetBranchAddress("jetGenJetPhi", jetGenJetPhi, &b_jetGenJetPhi);
   fChain->SetBranchAddress("jetGenPartonID", jetGenPartonID, &b_jetGenPartonID);
   fChain->SetBranchAddress("jetGenEn", jetGenEn, &b_jetGenEn);
   fChain->SetBranchAddress("jetGenPt", jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenPhi", jetGenPhi, &b_jetGenPhi);
   fChain->SetBranchAddress("jetVtxPt", jetVtxPt, &b_jetVtxPt);
   fChain->SetBranchAddress("jetVtx3dL", jetVtx3dL, &b_jetVtx3dL);
   fChain->SetBranchAddress("jetDPhiMETJet", jetDPhiMETJet, &b_jetDPhiMETJet);

   fChain->SetBranchAddress("jetMVAs", jetMVAs, &b_jetMVAs);
   fChain->SetBranchAddress("jetWPLevels", jetWPLevels, &b_jetWPLevels);
   fChain->SetBranchAddress("jetMVAsExt", jetMVAsExt, &b_jetMVAsExt);
   fChain->SetBranchAddress("jetWPLevelsExt", jetWPLevelsExt, &b_jetWPLevelsExt);
   fChain->SetBranchAddress("jetDRMean", jetDRMean, &b_jetDRMean);
   fChain->SetBranchAddress("jetDR2Mean", jetDR2Mean, &b_jetDR2Mean);
   fChain->SetBranchAddress("jetDZ", jetDZ, &b_jetDZ);
   fChain->SetBranchAddress("jetFrac01", jetFrac01, &b_jetFrac01);
   fChain->SetBranchAddress("jetFrac02", jetFrac02, &b_jetFrac02);
   fChain->SetBranchAddress("jetFrac03", jetFrac03, &b_jetFrac03);
   fChain->SetBranchAddress("jetFrac04", jetFrac04, &b_jetFrac04);
   fChain->SetBranchAddress("jetFrac05", jetFrac05, &b_jetFrac05);
   fChain->SetBranchAddress("jetFrac06", jetFrac06, &b_jetFrac06);
   fChain->SetBranchAddress("jetFrac07", jetFrac07, &b_jetFrac07);
   fChain->SetBranchAddress("jetBeta", jetBeta, &b_jetBeta);
   fChain->SetBranchAddress("jetBetaStarCMG", jetBetaStarCMG, &b_jetBetaStarCMG);
   fChain->SetBranchAddress("jetBetaStarClassic", jetBetaStarClassic, &b_jetBetaStarClassic);
   fChain->SetBranchAddress("jetBetaExt", jetBetaExt, &b_jetBetaExt);
   fChain->SetBranchAddress("jetBetaStarCMGExt", jetBetaStarCMGExt, &b_jetBetaStarCMGExt);
   fChain->SetBranchAddress("jetBetaStarClassicExt", jetBetaStarClassicExt, &b_jetBetaStarClassicExt);
   fChain->SetBranchAddress("jetNNeutrals", jetNNeutrals, &b_jetNNeutrals);
   fChain->SetBranchAddress("jetNCharged", jetNCharged, &b_jetNCharged);
   fChain->SetBranchAddress("jetPFLooseId", jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("nLowPtJet", &nLowPtJet, &b_nLowPtJet);
   fChain->SetBranchAddress("jetLowPtEn", jetLowPtEn, &b_jetLowPtEn);
   fChain->SetBranchAddress("jetLowPtPt", jetLowPtPt, &b_jetLowPtPt);
   fChain->SetBranchAddress("jetLowPtEta", jetLowPtEta, &b_jetLowPtEta);
   fChain->SetBranchAddress("jetLowPtPhi", jetLowPtPhi, &b_jetLowPtPhi);
   fChain->SetBranchAddress("jetLowPtCharge", jetLowPtCharge, &b_jetLowPtCharge);
   fChain->SetBranchAddress("jetLowPtEt", jetLowPtEt, &b_jetLowPtEt);
   fChain->SetBranchAddress("jetLowPtRawPt", jetLowPtRawPt, &b_jetLowPtRawPt);
   fChain->SetBranchAddress("jetLowPtRawEn", jetLowPtRawEn, &b_jetLowPtRawEn);
   fChain->SetBranchAddress("jetLowPtArea", jetLowPtArea, &b_jetLowPtArea);
   fChain->SetBranchAddress("jetLowPtGenJetEn", jetLowPtGenJetEn, &b_jetLowPtGenJetEn);
   fChain->SetBranchAddress("jetLowPtGenJetPt", jetLowPtGenJetPt, &b_jetLowPtGenJetPt);
   fChain->SetBranchAddress("jetLowPtGenJetEta", jetLowPtGenJetEta, &b_jetLowPtGenJetEta);
   fChain->SetBranchAddress("jetLowPtGenJetPhi", jetLowPtGenJetPhi, &b_jetLowPtGenJetPhi);
   fChain->SetBranchAddress("jetLowPtPartonID", jetLowPtPartonID, &b_jetLowPtPartonID);
   fChain->SetBranchAddress("jetLowPtGenPartonID", jetLowPtGenPartonID, &b_jetLowPtGenPartonID);
   fChain->SetBranchAddress("jetLowPtGenEn", jetLowPtGenEn, &b_jetLowPtGenEn);
   fChain->SetBranchAddress("jetLowPtGenPt", jetLowPtGenPt, &b_jetLowPtGenPt);
   fChain->SetBranchAddress("jetLowPtGenEta", jetLowPtGenEta, &b_jetLowPtGenEta);
   fChain->SetBranchAddress("jetLowPtGenPhi", jetLowPtGenPhi, &b_jetLowPtGenPhi);
   /* fChain->SetBranchAddress("nZee", &nZee, &b_nZee); */
   /* fChain->SetBranchAddress("ZeeMass", ZeeMass, &b_ZeeMass); */
   /* fChain->SetBranchAddress("ZeePt", ZeePt, &b_ZeePt); */
   /* fChain->SetBranchAddress("ZeeEta", ZeeEta, &b_ZeeEta); */
   /* fChain->SetBranchAddress("ZeePhi", ZeePhi, &b_ZeePhi); */
   /* fChain->SetBranchAddress("ZeeLeg1Index", ZeeLeg1Index, &b_ZeeLeg1Index); */
   /* fChain->SetBranchAddress("ZeeLeg2Index", ZeeLeg2Index, &b_ZeeLeg2Index); */
   /* fChain->SetBranchAddress("nZmumu", &nZmumu, &b_nZmumu); */
   /* fChain->SetBranchAddress("ZmumuMass", ZmumuMass, &b_ZmumuMass); */
   /* fChain->SetBranchAddress("ZmumuPt", ZmumuPt, &b_ZmumuPt); */
   /* fChain->SetBranchAddress("ZmumuEta", ZmumuEta, &b_ZmumuEta); */
   /* fChain->SetBranchAddress("ZmumuPhi", ZmumuPhi, &b_ZmumuPhi); */
   /* fChain->SetBranchAddress("ZmumuLeg1Index", ZmumuLeg1Index, &b_ZmumuLeg1Index); */
   /* fChain->SetBranchAddress("ZmumuLeg2Index", ZmumuLeg2Index, &b_ZmumuLeg2Index); */
   /* fChain->SetBranchAddress("nWenu", &nWenu, &b_nWenu); */
   /* fChain->SetBranchAddress("WenuMassTPfMET", WenuMassTPfMET, &b_WenuMassTPfMET); */
   /* fChain->SetBranchAddress("WenuEtPfMET", WenuEtPfMET, &b_WenuEtPfMET); */
   /* fChain->SetBranchAddress("WenuACopPfMET", WenuACopPfMET, &b_WenuACopPfMET); */
   /* fChain->SetBranchAddress("WenuEleIndex", WenuEleIndex, &b_WenuEleIndex); */
   /* fChain->SetBranchAddress("nWmunu", &nWmunu, &b_nWmunu); */
   /* fChain->SetBranchAddress("WmunuMassTPfMET", WmunuMassTPfMET, &b_WmunuMassTPfMET); */
   /* fChain->SetBranchAddress("WmunuEtPfMET", WmunuEtPfMET, &b_WmunuEtPfMET); */
   /* fChain->SetBranchAddress("WmunuACopPfMET", WmunuACopPfMET, &b_WmunuACopPfMET); */
   /* fChain->SetBranchAddress("WmunuMuIndex", WmunuMuIndex, &b_WmunuMuIndex); */
   fChain->SetBranchAddress("nConv", &nConv, &b_nConv);
   fChain->SetBranchAddress("convP4", convP4, &b_convP4);
   fChain->SetBranchAddress("convVtx", convVtx, &b_convVtx);
   fChain->SetBranchAddress("convVtxErr", convVtxErr, &b_convVtxErr);
   fChain->SetBranchAddress("convPairMomentum", convPairMomentum, &b_convPairMomentum);
   fChain->SetBranchAddress("convRefittedMomentum", convRefittedMomentum, &b_convRefittedMomentum);
   fChain->SetBranchAddress("convNTracks", convNTracks, &b_convNTracks);
   fChain->SetBranchAddress("convPairInvMass", convPairInvMass, &b_convPairInvMass);
   fChain->SetBranchAddress("convPairCotThetaSep", convPairCotThetaSep, &b_convPairCotThetaSep);
   fChain->SetBranchAddress("convEoverP", convEoverP, &b_convEoverP);
   fChain->SetBranchAddress("convDistOfMinApproach", convDistOfMinApproach, &b_convDistOfMinApproach);
   fChain->SetBranchAddress("convDPhiTrksAtVtx", convDPhiTrksAtVtx, &b_convDPhiTrksAtVtx);
   fChain->SetBranchAddress("convDPhiTrksAtEcal", convDPhiTrksAtEcal, &b_convDPhiTrksAtEcal);
   fChain->SetBranchAddress("convDEtaTrksAtEcal", convDEtaTrksAtEcal, &b_convDEtaTrksAtEcal);
   fChain->SetBranchAddress("convDxy", convDxy, &b_convDxy);
   fChain->SetBranchAddress("convDz", convDz, &b_convDz);
   fChain->SetBranchAddress("convLxy", convLxy, &b_convLxy);
   fChain->SetBranchAddress("convLz", convLz, &b_convLz);
   fChain->SetBranchAddress("convZofPrimVtxFromTrks", convZofPrimVtxFromTrks, &b_convZofPrimVtxFromTrks);
   fChain->SetBranchAddress("convNHitsBeforeVtx", convNHitsBeforeVtx, &b_convNHitsBeforeVtx);
   fChain->SetBranchAddress("convNSharedHits", convNSharedHits, &b_convNSharedHits);
   fChain->SetBranchAddress("convValidVtx", convValidVtx, &b_convValidVtx);
   fChain->SetBranchAddress("convMVALikelihood", convMVALikelihood, &b_convMVALikelihood);
   fChain->SetBranchAddress("convChi2", convChi2, &b_convChi2);
   fChain->SetBranchAddress("convChi2Probability", convChi2Probability, &b_convChi2Probability);
   fChain->SetBranchAddress("convTk1Dz", convTk1Dz, &b_convTk1Dz);
   fChain->SetBranchAddress("convTk2Dz", convTk2Dz, &b_convTk2Dz);
   fChain->SetBranchAddress("convTk1DzErr", convTk1DzErr, &b_convTk1DzErr);
   fChain->SetBranchAddress("convTk2DzErr", convTk2DzErr, &b_convTk2DzErr);
   fChain->SetBranchAddress("convCh1Ch2", convCh1Ch2, &b_convCh1Ch2);
   fChain->SetBranchAddress("convTk1D0", convTk1D0, &b_convTk1D0);
   fChain->SetBranchAddress("convTk1Pout", convTk1Pout, &b_convTk1Pout);
   fChain->SetBranchAddress("convTk1Pin", convTk1Pin, &b_convTk1Pin);
   fChain->SetBranchAddress("convTk2D0", convTk2D0, &b_convTk2D0);
   fChain->SetBranchAddress("convTk2Pout", convTk2Pout, &b_convTk2Pout);
   fChain->SetBranchAddress("convTk2Pin", convTk2Pin, &b_convTk2Pin);
   Notify();

}

Bool_t xAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void xAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t xAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}




// ---------------------------------------------------------------------------------------------------------------------------------------------
// Vertex Analysis
// ---------------------------------------------------------------------------------------------------------------------------------------------
class ggVertexInfo : public VertexInfoAdapter
{
public:
  ggVertexInfo( const vector<float> &vtx_x, const vector<float> &vtx_y, const vector<float> &vtx_z,
		const vector<int> &vtx_ntks,
 		const vector<float> &tk_px, const vector<float> &tk_py, const vector<float> &tk_pz,
		const vector<float> &tk_pterr, const vector<int> &tk_quality,
		vector<vector<unsigned short> > & vtx_tkind,
		vector<vector<float> >    & vtx_tkweight
		);
  
  virtual int nvtx() const    { return _vtx_n; };
  virtual int ntracks() const { return _tk_n; };
  
  virtual bool hasVtxTracks() const { return true; }
  virtual const unsigned short * vtxTracks(int ii) const { return &(_vtx_tkind)[ii][0]; };
  virtual int vtxNTracks(int ii) const { return _vtx_ntks[ii]; };
  virtual const float * vtxTkWeights(int ii) const { return &(_vtx_tkweight)[ii][0]; };
  
  virtual float tkpx(int ii) const { return _tk_px[ii]; };
  virtual float tkpy(int ii) const { return _tk_py[ii]; };
  virtual float tkpz(int ii) const { return _tk_pz[ii]; };
  
  virtual float tkPtErr(int ii) const { return _tk_pterr[ii]; };
  virtual int   tkVtxId(int ii) const { return -1; };
  
  virtual float tkWeight(int ii, int jj) const { return _vtx_tkweight[jj][ii]; };
  
  virtual float vtxx(int ii) const { return _vtx_x[ii]; };
  virtual float vtxy(int ii) const { return _vtx_y[ii]; };
  virtual float vtxz(int ii) const { return _vtx_z[ii]; };
  
  virtual float tkd0(int ii, int jj) const { return 0.; }; // FIXME
  virtual float tkd0Err(int ii, int jj) const { return 1.; };  // FIXME
  
  virtual float tkdz(int ii, int jj) const { return 0.; };  // FIXME
  virtual float tkdzErr(int ii, int jj) const { return 1.; };  // FIXME
  
  virtual bool tkIsHighPurity(int ii) const { return ( _tk_quality[ii] & (1<<2) ) >> 2; };
  
  // virtual ~ggVertexInfo();
  
 private:
  int _vtx_n;
  vector<float> _vtx_x,_vtx_y,_vtx_z;
  vector<int>   _vtx_ntks;
  int _tk_n;
  vector<vector<unsigned short> > _vtx_tkind;
  vector<vector<float> > _vtx_tkweight;
  vector<float> _tk_px, _tk_py, _tk_pz, _tk_pterr;
  vector<int>  _tk_quality;
  
};





/*   ///FC: do not use this anymore

class ggTupleVertexInfo : public VertexInfoAdapter
{
public:

	ggTupleVertexInfo(int nvtx, float * vtxx, float * vtxy, float * vtxz, 
			   int ntracks, float * tkpx, float * tkpy, float * tkpz,
			   float * tkPtErr, int * tkVtxId, float * tkWeight,
			   float * tkd0, float * tkd0Err, float * tkdz, float * tkdzErr,
			  const unsigned short ** vtxTracks, int* vtxNTracks, const float** vtxTkWeights,
			   int  * tkIsHighPurity);
	
	virtual int nvtx() const    { return nvtx_; };
	virtual int ntracks() const { return ntracks_; };

	virtual bool hasVtxTracks()  const { return true; };
	virtual const unsigned short * vtxTracks(int i) const { return vtxTracks_[i]; };
	virtual int vtxNTracks(int i)  const { return vtxNTracks_[i]; };
	virtual const float * vtxTkWeights(int i) const { return vtxTrackWeights_[i]; };

	virtual float tkpx(int ii) const { return tkpx_ != 0 ? tkpx_[ii] : 0.; };
	virtual float tkpy(int ii) const { return tkpx_ != 0 ? tkpy_[ii] : 0.; };
	virtual float tkpz(int ii) const { return tkpx_ != 0 ? tkpz_[ii] : 0.; };
	
	virtual float tkPtErr(int ii) const { return tkPtErr_  != 0 ? tkPtErr_[ii] : 999.; };
	virtual int   tkVtxId(int ii) const { return tkVtxId_  != 0 ? tkVtxId_[ii] : 999; };

	virtual float tkWeight(int ii, int jj) const { return tkWeight_ != 0 ? tkWeight_[ii]*(float)( tkVtxId(ii) == jj) : 0.; };
	
	virtual float vtxx(int ii) const { return vtxx_ != 0 ? vtxx_[ii] : 0.; };
	virtual float vtxy(int ii) const { return vtxy_ != 0 ? vtxy_[ii] : 0.; };
	virtual float vtxz(int ii) const { return vtxz_ != 0 ? vtxz_[ii] : 0.; };

	//virtual float tkd0(int ii, int jj) const { assert(tkVtxId(ii) == jj); return tkd0_ != 0 ? tkd0_[ii] : 0.; };
	//virtual float tkd0Err(int ii, int jj) const { assert(tkVtxId(ii) == jj); return tkd0Err_ != 0 ? tkd0Err_[ii] : 0.; };
	virtual float tkd0(int ii, int jj) const { return 0.;}
	virtual float tkd0Err(int ii, int jj) const { return 1.;}

	//virtual float tkdz(int ii, int jj) const { assert(tkVtxId(ii) == jj); return tkdz_ != 0 ? tkdz_[ii] : 0.; };
	//virtual float tkdzErr(int ii, int jj) const { assert(tkVtxId(ii) == jj); return tkdzErr_ != 0 ? tkdzErr_[ii] : 0.; };
	virtual float tkdz(int ii, int jj) const { assert(tkVtxId(ii) == jj); return 0.; };
	virtual float tkdzErr(int ii, int jj) const { assert(tkVtxId(ii) == jj); return 1.; };

	//virtual bool tkIsHighPurity(int ii) const { return tkIsHighPurity_ != 0 ? tkIsHighPurity_[ii] : 0.; };
	virtual bool tkIsHighPurity(int ii) const { return (tkIsHighPurity_[ii] & (1 << 2) ) >> 2; };

	virtual ~ggTupleVertexInfo(){}
	
private:
	int nvtx_;
	float * vtxx_;
	float * vtxy_;
	float * vtxz_;

	int ntracks_;
	float * tkpx_;
	float * tkpy_;
	float * tkpz_;
	float * tkPtErr_;
	int * tkVtxId_;	
	float * tkWeight_;
	float * tkd0_;
	float * tkd0Err_;
	float * tkdz_;
	float * tkdzErr_;

	const unsigned short ** vtxTracks_;
	int * vtxNTracks_;
	const float ** vtxTrackWeights_;

	int * tkIsHighPurity_;
};


*/





#endif // #ifdef xAna_cxx
