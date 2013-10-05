#ifndef  __minitree_h
#define __minitree_h

#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

class MiniTree {
 public:
  MiniTree( const char * filename, bool read = false );
  ~MiniTree(void);

  /// TFile
  TFile *mtree_file;
  /// main tree with di-photon system variable
  TChain *mtree_mainTree;
  TChain *mtree_miniTree;
  TChain *mtree_mcTrue;

  /// dedicated tree to specific channel analysis
  TChain *mtree_elecTree;
  TChain *mtree_muonTree;
  TChain *mtree_dijetTree;

  bool mtree_fillLepTree;
  bool mtree_fillDijetTree;

  /// radion => HH studies
  TChain *mtree_radion;

  bool addSyncVariables;      

  int mtree_ievt; /// pointer to event in main tree

  void createBranches(void);
  void setBranchesAddresses(void);
  void initEvent(void);
  void fill(void);
  void fillMCtrueOnly(void);

  void end(void);

  ///////////////////////////////////////////////////////////
  /////////////// Declare all variables /////////////////////
  ///////////////////////////////////////////////////////////

  /// MC weights
  float mc_wei;
  float mc_wBSz;
  float mc_wPhoEffi;
  float mc_wTrigEffi;
  float mc_wPU, mc_wVtxId;
  float mc_wHQT;
  float mc_wXsec, mc_wNgen;

  /// truth variables
  float mc_mH, mc_cThetaStar_CS;

  /// event id
  Int_t    mtree_runNum;
  Long64_t mtree_evtNum;
  Int_t    mtree_lumiSec;
  Float_t  weight;
  
  /// general variables
  Float_t mtree_rho , mtree_rho25;
  Float_t mtree_zVtx;
  Int_t   mtree_nVtx;
  Int_t   mtree_ivtx1, mtree_ivtx2, mtree_ivtx3;
  Float_t mtree_mass, mtree_pt, mtree_piT, mtree_massDefVtx;
  Float_t mtree_y;
  Float_t mtree_massResoTot,mtree_massResoAng,mtree_massResoEng;
  Float_t mtree_massResoRightVtx, mtree_massResoRandVtx;
  Float_t mtree_vtxProb, mtree_vtxMva ;
  Float_t mtree_diphoMva;
  Float_t mtree_rawMet, mtree_rawMetPhi;
  Float_t mtree_corMet, mtree_corMetPhi;
  Int_t   mtree_catMva;
  Int_t   mtree_catBase;
  Int_t   mtree_vbfTag;
  Int_t   mtree_vbfCat;
  Int_t   mtree_metTag;
  Int_t   mtree_lepTag;
  Int_t   mtree_lepCat;
  Int_t   mtree_hvjjTag;  
  

  /// for photon specific use array[2]
  /// + pho1, pho2 only for variables going in the fitter (roofit)
  float mtree_eg[2],mtree_ptg[2], mtree_eta[2], mtree_phi[2],mtree_sceta[2];
  float mtree_relResOverE[2], mtree_relSmearing[2];
  int   mtree_isElec[2];
  int   mtree_cat[2];
  
  /// spin specific
  float mtree_cThetaLead_heli, mtree_cThetaTrail_heli,mtree_cThetaStar_CS;
  float mtree_pt1, mtree_pt2, mtree_mvaid1, mtree_mvaid2, mtree_eta1, mtree_eta2, mtree_r91, mtree_r92;
  float mtree_pit1, mtree_pit2;
  float mtree_minR9, mtree_minPhoIdEB, mtree_minPhoIdEE;
  float mtree_minSCEta, mtree_maxSCEta;

  /// phoid variables
  float mtree_mvaid[2], mtree_s4Ratio[2];
  float mtree_phoIso1[2], mtree_phoIso2[2], mtree_phoIso3[2]; 
  float mtree_coviEiE[2], mtree_coviEiP[2], mtree_r9[2], mtree_drTrk[2];
  float mtree_etCorrEcalIso[2], mtree_etCorrHcalIso[2], mtree_etCorrTrkIso[2];
  float mtree_pfChgIso02_selVtx[2];
  float mtree_esEffSigmaRR[2];
  float mtree_HoE[2];
    
  /// photon PF isolation
  float mtree_pfPhoIso03[2], mtree_pfPhoIso04[2];
  float mtree_pfChgIso03_selVtx[2], mtree_pfChgIso04_selVtx[2];
  float mtree_pfChgIso03_badVtx[2], mtree_pfChgIso03_defVtx[2];
  float mtree_pfChgIso04_badVtx[2], mtree_pfChgIso04_defVtx[2];

  
  /// sync only variables (float only var for simpl)
  float convInd[2], convNtk[2], convValidVtx[2];
  float convChi2Prob[2], convPt[2];

  float phoID_covIEtaIPhi[2], phoID_S4[2], phoID_ESEffSigmaRR[2],  phoID_pfPhoIso03[2];
  float phoID_pfChIso03[2], phoID_pfChIso03worst[2],phoID_E2x2[2], phoID_E5x5[2];
      
  float vtxId[3];
  float vtxMva[3], vtxDz[3];
  float vtxPtBal, vtxAsym, vtxLogSumPt2, vtxPullToConv, vtxNConv;


  /// dijet selection variables
  float dijet_dEtaJJ, dijet_zeppenfeld, dijet_mJJ, dijet_dphiGGJJ;
  float dijet_mvaVbf, dijet_mvaVbf_sig, dijet_mvaVbf_dipho, dijet_mvaVbf_gluglu;
  float dijet_ptj1, dijet_ptj2;
  float dijet_ptj[2], dijet_jec[2], dijet_etj[2];
  float dijet_eta[2], dijet_phi[2], dijet_en[2];
  float dijet_betaStar[2], dijet_rms[2];

  /// radion analysis variables
  Long64_t radion_evtNum;
  TLorentzVector *radion_gamma1, *radion_gamma2;
  TClonesArray *radion_jets;
  TArrayF *radion_bJetTags;

  /// vector to do synchronisation printout
  vector<string> sync_iName, sync_lName,sync_fName, sync_mtreeName;
  vector<int*> sync_iVal;
  vector<Long64_t*> sync_lVal;
  vector<float*> sync_fVal;
  //    virtual void loopSync(void);
  void setSynchVariables(void);

};


MiniTree::MiniTree( const char * filenames, bool read ) {
  radion_gamma1 = new TLorentzVector;
  radion_gamma2 = new TLorentzVector;
  radion_jets = new TClonesArray("TLorentzVector");
  radion_bJetTags = new TArrayF;

  radion_gamma1->Class()->IgnoreTObjectStreamer();
  radion_gamma2->Class()->IgnoreTObjectStreamer();
  radion_jets->Class()->IgnoreTObjectStreamer();
  radion_bJetTags->Class()->IgnoreTObjectStreamer();

  addSyncVariables = true;
  
  if( !read ) {
    mtree_ievt = 0;
    mtree_file = new TFile( filenames,"RECREATE", "H->gg input tree for unbinned maximum-likelihood fit");
    
    mtree_mainTree  = (TChain*) new TTree( "HToGG","main minitree for unbinned fit" );
    mtree_elecTree  = (TChain*) new TTree( "HToGG_elevar","electron variables ONLY" );
    mtree_muonTree  = (TChain*) new TTree( "HToGG_muovar","muon variables ONLY"     );
    mtree_dijetTree = (TChain*) new TTree( "HToGG_dijet" ,"dijet variables ONLY"    );
    mtree_mcTrue    = (TChain*) new TTree( "HToGG_mcTrue" ,"dijet variables ONLY"    );
    mtree_radion    = (TChain*) new TTree( "RToHH", "studies of radion => HH decays" );

    createBranches();
    initEvent();
  } else {
    /// reread already existing minitree
    mtree_mainTree  = new TChain( "HToGG" );
    mtree_elecTree  = new TChain( "HToGG_elevar");
    mtree_muonTree  = new TChain( "HToGG_muovar");
    mtree_dijetTree = new TChain( "HToGG_dijet" );
    mtree_mcTrue    = new TChain( "HToGG_mcTrue" );

    /// adding files in filenames may need a loop if file name too complex
    mtree_mainTree ->Add(filenames);
    mtree_elecTree ->Add(filenames);
    mtree_muonTree ->Add(filenames);
    mtree_dijetTree->Add(filenames);
    mtree_mcTrue   ->Add(filenames);

    setBranchesAddresses();
  }
}

void MiniTree::fillMCtrueOnly(void) {
   mtree_mcTrue->Fill();
}

void MiniTree::fill(void) {
  /// roofti doesnot really like arrays
  mtree_pt1  = mtree_ptg[0];
  mtree_pt2  = mtree_ptg[1];
  mtree_pit1 = mtree_ptg[0]/mtree_mass;
  mtree_pit2 = mtree_ptg[1]/mtree_mass;

  //Do not fill useless trees to cut down on file size.
  //mtree_mainTree->Fill();
  //if( mtree_fillLepTree && mtree_lepCat  == 0 ) mtree_elecTree->Fill();
  //if( mtree_fillLepTree && mtree_lepCat  == 1 ) mtree_muonTree->Fill();
  //if( mtree_fillDijetTree                     ) mtree_dijetTree->Fill();
  
  mtree_radion->Fill();

  // increase internal ievt pointer
  mtree_ievt++;
}


void MiniTree::initEvent(void) {

  mtree_fillLepTree = false;
  mtree_fillDijetTree = false;

  /// MC weights
  mc_wei = mc_wBSz = 1;
  mc_wPhoEffi = mc_wTrigEffi = 1;
  mc_wPU   = 1;
  mc_wHQT  = 1; 
  mc_wXsec = 1;
  mc_wNgen = 1;
  mc_wVtxId = 1;
  
  /// MC truth variables
  mc_mH = -1;
  mc_cThetaStar_CS = -999;

  /// event id
  mtree_runNum = mtree_evtNum = mtree_lumiSec = -1;
  
  /// weight section  
  weight = 1;
  
  /// general variables
  mtree_rho = mtree_rho25 = -1;
  mtree_zVtx = -999;
  mtree_nVtx = -1;
  mtree_ivtx1 = mtree_ivtx2 = mtree_ivtx3 = -1;
  mtree_mass  = mtree_pt = mtree_piT = mtree_massDefVtx = -1;
  mtree_massResoTot = mtree_massResoAng = mtree_massResoEng = -1;
  mtree_massResoRightVtx = mtree_massResoRandVtx = -1;
  mtree_y = -999;
  mtree_vtxProb  = -1;
  mtree_diphoMva = -1;
  mtree_rawMet = mtree_rawMetPhi = -1;
  mtree_corMet = mtree_corMetPhi = -1;
  mtree_catMva = -1;
  mtree_catBase=-1;
  /// set exclusive tagging to zero
  mtree_vbfTag = mtree_metTag = mtree_lepTag = 0;
  mtree_hvjjTag = 0;  
  mtree_vbfCat = -1;
  mtree_lepCat = -1;
  
  mtree_cThetaLead_heli = mtree_cThetaTrail_heli = mtree_cThetaStar_CS = -999;

  /// photon specific use array[2]
  for( int i = 0; i < 2; i++ ) {
    mtree_eg[i] = mtree_ptg[i] = mtree_eta[i] = mtree_phi[i] = mtree_sceta[i] = -999;
    mtree_relResOverE[i] = mtree_relSmearing[i] = -999;
    mtree_cat[i]    = -1;
    mtree_isElec[i] = -1;
    
    mtree_pt1 = mtree_pt2 = -999;
    mtree_pit1 = mtree_pit2 = -999;
    mtree_mvaid1 = mtree_mvaid2 = -999;
    mtree_eta1 = mtree_eta2 = -999;
    mtree_r91 = mtree_r92 = -999;
    
    /// phoid variables
    mtree_mvaid[i] = mtree_s4Ratio[i] = -999;
    mtree_phoIso1[i] = mtree_phoIso2[i] = mtree_phoIso3[i] = -999; 
    mtree_coviEiE[i] = mtree_coviEiP[i] = mtree_r9[i] = mtree_drTrk[i] = -999;
    mtree_etCorrEcalIso[i] = mtree_etCorrHcalIso[i] = mtree_etCorrTrkIso[i] = -999;
    mtree_pfChgIso02_selVtx[i] = -999;
    mtree_esEffSigmaRR[i] = -999;
    mtree_HoE[i] = -999;
    
    /// photon PF isolation
    mtree_pfPhoIso03[i] = mtree_pfPhoIso04[i] = -999;
    mtree_pfChgIso03_selVtx[i] = mtree_pfChgIso04_selVtx[i] = -999;
    mtree_pfChgIso03_badVtx[i] = mtree_pfChgIso03_defVtx[i] = -999;
    mtree_pfChgIso04_badVtx[i] = mtree_pfChgIso04_defVtx[i] = -999;
    
    /// dijet selection variables
    dijet_dEtaJJ = dijet_zeppenfeld = dijet_mJJ = dijet_dphiGGJJ = dijet_ptj1 = dijet_ptj2 = -999;
    dijet_ptj[i] = dijet_jec[i] = dijet_etj[i] = -999;
    dijet_eta[i] = dijet_phi[i] = dijet_en[i] = -999;
    dijet_betaStar[i] = dijet_rms[i] = -999;

    /// radion selection variables
    radion_evtNum = -1;
    radion_gamma1->SetXYZT(0, 0, 0, 0);
    radion_gamma2->SetXYZT(0, 0, 0, 0);
    radion_jets->Clear();
    radion_bJetTags->Set(0);

    if( addSyncVariables ) {
      convInd[i] = convNtk[i] = convValidVtx[i] = -999;
      convChi2Prob[i] = convPt[i] = -999;      
      
      phoID_covIEtaIPhi[i] = phoID_S4[i] = phoID_ESEffSigmaRR[i] = phoID_pfPhoIso03[i] = -999;
      phoID_pfChIso03[i]   = phoID_pfChIso03worst[i] = phoID_E2x2[i] = phoID_E5x5[i] = -999;

    }

  }

  if( addSyncVariables ) {
    vtxPtBal =vtxAsym = vtxLogSumPt2 =  vtxPullToConv = -999;
    for( int iv = 0 ; iv < 3; iv ++ ) vtxId[iv] = vtxMva[iv] = vtxDz[iv] = -999;
  }

  dijet_mvaVbf        = -1;
  dijet_mvaVbf_sig    = -1;
  dijet_mvaVbf_dipho  = -1;
  dijet_mvaVbf_gluglu = -1;
    
}

void MiniTree::end(void) {
  mtree_file->Write() ;
  mtree_file->Close();  
}

MiniTree::~MiniTree(void) {
  delete radion_gamma1;
  delete radion_gamma2;
  delete radion_jets;
  delete radion_bJetTags;
}



void MiniTree::createBranches(void) {
  cout << " ------------- MiniTree set output branches ------------" << endl;
  mtree_mcTrue->Branch("mH"           , &mc_mH          , "mH/F");
  mtree_mcTrue->Branch("cThetaStar_CS", &mc_cThetaStar_CS, "cThetaStar_CS/F");
  
  //----------------------- MC weights
  mtree_mainTree->Branch("mH"       , &mc_mH      , "mH/F");
  mtree_mainTree->Branch("cThetaStar_CS_truth", &mc_cThetaStar_CS, "cThetaStar_CS_truth/F");
  mtree_mainTree->Branch("wei"      , &mc_wei      , "wei/F");
  mtree_mainTree->Branch("wPU"      , &mc_wPU      , "wPU/F");
  mtree_mainTree->Branch("wHQT"     , &mc_wHQT     , "wHQT/F");
  mtree_mainTree->Branch("wXsec"    , &mc_wXsec    , "wXsec/F");
  mtree_mainTree->Branch("wNgen"    , &mc_wNgen    , "wNgen/F");
  mtree_mainTree->Branch("wBSz"     , &mc_wBSz     , "wBSz/F");
  mtree_mainTree->Branch("wVtxId"   , &mc_wVtxId   , "wVtxId/F");
  mtree_mainTree->Branch("wPhoEffi" , &mc_wPhoEffi , "wPhoEffi/F");
  mtree_mainTree->Branch("wTrigEffi", &mc_wTrigEffi, "wTrigEffi/F");


  //-----------------------
  mtree_mainTree->Branch("iMainEvt", &mtree_ievt    , "iMainEvt/I" ); ///pointer to main tree
  mtree_mainTree->Branch("runNum"  , &mtree_runNum  , "runNum/I"   );
  mtree_mainTree->Branch("evtNum"  , &mtree_evtNum  , "evtNum/L" );
  mtree_mainTree->Branch("lumiSec" , &mtree_lumiSec , "lumiSec/I"  );
  mtree_mainTree->Branch("rho"     , &mtree_rho     , "rho/F");
  mtree_mainTree->Branch("rho25"   , &mtree_rho25   , "rho25/F");

  mtree_mainTree->Branch("zVtx"    , &mtree_zVtx  , "zVtx/F");
  mtree_mainTree->Branch("nVtx"    , &mtree_nVtx  , "nVtx/I");
  mtree_mainTree->Branch("iVtx1"   , &mtree_ivtx1 , "iVtx1/I");
  mtree_mainTree->Branch("iVtx2"   , &mtree_ivtx2 , "iVtx2/I");
  mtree_mainTree->Branch("iVtx3"   , &mtree_ivtx3 , "iVtx3/I");
  mtree_mainTree->Branch("vtxProb" , &mtree_vtxProb, "vtxProb/F");
  mtree_mainTree->Branch("vtxMva"  , &mtree_vtxMva , "vtxMva/F");
  
  // skeleton
  //  mtree_mainTree->Branch("x", &mtree_x ,"x/F"  );
  //--------------- diphoton event variables
  mtree_mainTree->Branch("mass"      , &mtree_mass       , "mass/F");
  mtree_mainTree->Branch("massDefVtx", &mtree_massDefVtx , "massDefVtx/F");
  mtree_mainTree->Branch("pt"        , &mtree_pt         , "pt/F"   );
  mtree_mainTree->Branch("y"         , &mtree_y          , "y/F"    );
  mtree_mainTree->Branch("piT"       , &mtree_piT        , "piT/F"  );
  mtree_mainTree->Branch("pit1"      , &mtree_pit1       , "pit1/F"  );
  mtree_mainTree->Branch("pit2"      , &mtree_pit2       , "pit2/F"  );
  mtree_mainTree->Branch("minR9"     , &mtree_minR9      , "minR9/F"  );
  mtree_mainTree->Branch("maxSCEta"  , &mtree_maxSCEta   , "maxSCEta/F"  );
  mtree_mainTree->Branch("minSCEta"  , &mtree_minSCEta   , "minSCEta/F"  );
  mtree_mainTree->Branch("minPhoIdEB", &mtree_minPhoIdEB , "minPhoIdEB/F"  );
  mtree_mainTree->Branch("minPhoIdEE", &mtree_minPhoIdEE , "minPhoIdEE/F"  );
  mtree_mainTree->Branch("massResoTot"     , &mtree_massResoTot      , "massResoTot/F"  );
  mtree_mainTree->Branch("massResoAng"     , &mtree_massResoAng      , "massResoAng/F"  );
  mtree_mainTree->Branch("massResoEng"     , &mtree_massResoEng      , "massResoEng/F"  );
  mtree_mainTree->Branch("massResoRightVtx", &mtree_massResoRightVtx , "massResoRightVtx/F" );
  mtree_mainTree->Branch("massResoRandVtx" , &mtree_massResoRandVtx  , "massResoRandVtx/F"  );
  mtree_mainTree->Branch("diphoMva" , &mtree_diphoMva  , "diphoMva/F"  );
  mtree_mainTree->Branch("rawMet"   , &mtree_rawMet    , "rawMet/F"  );
  mtree_mainTree->Branch("corMet"   , &mtree_corMet    , "corMet/F"  );
  mtree_mainTree->Branch("rawMetPhi", &mtree_rawMetPhi , "rawMetPhi/F"  );
  mtree_mainTree->Branch("corMetPhi", &mtree_corMetPhi , "corMetPhi/F"  );

  mtree_mainTree->Branch("catMva" , &mtree_catMva  , "catMva/I"  );
  mtree_mainTree->Branch("catBase", &mtree_catBase , "catBase/I" );
  mtree_mainTree->Branch("lepTag" , &mtree_lepTag  , "lepTag/I"  );
  mtree_mainTree->Branch("vbfTag" , &mtree_vbfTag  , "vbfTag/I"  );
  mtree_mainTree->Branch("metTag" , &mtree_metTag  , "metTag/I"  );
  mtree_mainTree->Branch("hvjjTag", &mtree_hvjjTag , "hvjjTag/I" );
  mtree_mainTree->Branch("lepCat" , &mtree_lepCat  , "lepCat/I"  );
  mtree_mainTree->Branch("vbfCat" , &mtree_vbfCat  , "vbfCat/I"  );
  mtree_mainTree->Branch("cThetaLead_heli" ,  &mtree_cThetaLead_heli , "cThetaLead_heli/F"  );
  mtree_mainTree->Branch("cThetaTrail_heli",  &mtree_cThetaTrail_heli, "cThetaTrail_heli/F"  );
  mtree_mainTree->Branch("cThetaStar_CS"   ,  &mtree_cThetaStar_CS   , "cThetaStar_CS/F"  );
  mtree_mainTree->Branch("dEtaJJ"    , &dijet_dEtaJJ     , "dEtaJJ/F"    );
  mtree_mainTree->Branch("zeppenfeld", &dijet_zeppenfeld , "zeppenfeld/F");
  mtree_mainTree->Branch("mJJ"       , &dijet_mJJ        , "mJJ/F"       );
  mtree_mainTree->Branch("dphiGGJJ"  , &dijet_dphiGGJJ   , "dphiGGJJ/F"  );
  mtree_mainTree->Branch("jet1Pt", &dijet_ptj1 , "jet1Pt/F" );
  mtree_mainTree->Branch("jet2Pt", &dijet_ptj2 , "jet2Pt/F" );


  mtree_mainTree->Branch("eg" , mtree_eg , "eg[2]/F"  );
  mtree_mainTree->Branch("ptg", mtree_ptg, "ptg[2]/F" );
  mtree_mainTree->Branch("eta", mtree_eta, "eta[2]/F" );
  mtree_mainTree->Branch("phi", mtree_phi, "phi[2]/F" );
  mtree_mainTree->Branch("cat", mtree_cat, "cat[2]/I" );
  mtree_mainTree->Branch("sceta", mtree_sceta ,"sceta[2]/F"  );
  mtree_mainTree->Branch("relResOverE", mtree_relResOverE ,"relResOverE[2]/F"  );
  mtree_mainTree->Branch("relSmearing" , mtree_relSmearing  ,"relSmearing[2]/F"  );
  mtree_mainTree->Branch("isElec", mtree_isElec ,"isElec[2]/I"  );

  
/// phoid variables
  mtree_mainTree->Branch("mvaid"  , mtree_mvaid   , "mvaid[2]/F"    );
  mtree_mainTree->Branch("s4Ratio", mtree_s4Ratio , "s4Ratio[2]/F"  );
  mtree_mainTree->Branch("phoIso1", mtree_phoIso1 , "phoIso1[2]/F"  );
  mtree_mainTree->Branch("phoIso2", mtree_phoIso2 , "phoIso2[2]/F"  );
  mtree_mainTree->Branch("phoIso3", mtree_phoIso3 , "phoIso3[2]/F"  );
  mtree_mainTree->Branch("coviEiE", mtree_coviEiE , "coviEiE[2]/F"  );
  mtree_mainTree->Branch("coviEiP", mtree_coviEiP , "coviEiP[2]/F"  );
  mtree_mainTree->Branch("r9"     , mtree_r9      , "r9[2]/F"  );
  mtree_mainTree->Branch("HoE"    , mtree_HoE     , "HoE[2]/F"  );
  mtree_mainTree->Branch("drTrk"  , mtree_drTrk   , "drTrk[2]/F"  );
  mtree_mainTree->Branch("mvaid1" , &mtree_mvaid1 , "mvaid1/F"    );
  mtree_mainTree->Branch("mvaid2" , &mtree_mvaid2 , "mvaid2/F"    );
  mtree_mainTree->Branch("eta1"   , &mtree_eta1   , "eta1/F"    );
  mtree_mainTree->Branch("eta2"   , &mtree_eta2   , "eta2/F"    );
  mtree_mainTree->Branch("r91"    , &mtree_r91    , "r91/F"    );
  mtree_mainTree->Branch("r92"    , &mtree_r92    , "r92/F"    );

  /// sync only variables
  if( addSyncVariables ) {
    cout << "    ******* adding sync variables to minitree " << endl;
    mtree_mainTree->Branch("etCorrEcalIso", mtree_etCorrEcalIso, "etCorrEcalIso[2]/F");
    mtree_mainTree->Branch("etCorrHcalIso", mtree_etCorrHcalIso, "etCorrHcalIso[2]/F");
    mtree_mainTree->Branch("etCorrTrkIso" , mtree_etCorrTrkIso , "etCorrTrkIso[2]/F" );
    mtree_mainTree->Branch("esEffSigmaRR" , mtree_esEffSigmaRR , "esEffSigmaRR[2]/F"  );
    
    /// photon PF isolation
    mtree_mainTree->Branch("pfChgIso02_selVtx", mtree_pfChgIso02_selVtx , "pfChgIso02_selVtx[2]/F"  );
    mtree_mainTree->Branch("pfChgIso03_selVtx", mtree_pfChgIso03_selVtx , "pfChgIso03_selVtx[2]/F"  );
    mtree_mainTree->Branch("pfChgIso03_badVtx", mtree_pfChgIso03_badVtx , "pfChgIso03_badVtx[2]/F"  );
    mtree_mainTree->Branch("pfChgIso03_defVtx", mtree_pfChgIso03_defVtx , "pfChgIso03_defVtx[2]/F"  );
    mtree_mainTree->Branch("pfChgIso04_selVtx", mtree_pfChgIso04_selVtx , "pfChgIso04_selVtx[2]/F"  );
    mtree_mainTree->Branch("pfChgIso04_badVtx", mtree_pfChgIso04_badVtx , "pfChgIso04_badVtx[2]/F"  );
    mtree_mainTree->Branch("pfChgIso04_defVtx", mtree_pfChgIso04_defVtx , "pfChgIso04_defVtx[2]/F"  );
    mtree_mainTree->Branch("pfPhoIso03", mtree_pfPhoIso03 ,"pfPhoIso03[2]/F"  );
    mtree_mainTree->Branch("pfPhoIso04", mtree_pfPhoIso04 ,"pfPhoIso04[2]/F"  );
    
    
    mtree_mainTree->Branch("convInd", convInd, "convInd[2]/F");
    mtree_mainTree->Branch("convNtk", convNtk, "convNtk[2]/F");
    mtree_mainTree->Branch("convPt" , convPt ,"convPt[2]/F" );
    mtree_mainTree->Branch("convValidVtx", convValidVtx,"convValidVtx[2]/F" );
    mtree_mainTree->Branch("convChi2Prob", convChi2Prob, "convChi2Prob[2]/F");
    mtree_mainTree->Branch("vId"  , vtxId   , "vId[3]/F");
    mtree_mainTree->Branch("vDz"  , vtxDz   , "vDz[3]/F");
    mtree_mainTree->Branch("vMva" , vtxMva  , "vMva[3]/F");
    mtree_mainTree->Branch("vAsym"      , &vtxAsym      ,"vAsym/F"      );
    mtree_mainTree->Branch("vPtBal"     , &vtxPtBal     ,"vPtBal/F"     );
    mtree_mainTree->Branch("vLogSumPt2" , &vtxLogSumPt2 ,"vLogSumPt2/F" );
    mtree_mainTree->Branch("vPullToConv", &vtxPullToConv,"vPullToConv/F");
    mtree_mainTree->Branch("vNConv"     , &vtxNConv     ,"vNConv/F"     );
    mtree_mainTree->Branch("phoID_covIEtaIPhi" , phoID_covIEtaIPhi ,"phoID_covIEtaIPhi[2]/F");
    mtree_mainTree->Branch("phoID_S4"          , phoID_S4          ,"phoID_S4[2]/F"         );
    mtree_mainTree->Branch("phoID_ESEffSigmaRR", phoID_ESEffSigmaRR,"phoID_ESEffSigmaRR[2]/F");
    mtree_mainTree->Branch("phoID_pfPhoIso03"  , phoID_pfPhoIso03  ,"phoID_pfPhoIso03[2]/F");
    mtree_mainTree->Branch("phoID_pfChIso03"   , phoID_pfChIso03   ,"phoID_pfChIso03[2]/F");
    mtree_mainTree->Branch("phoID_pfChIso03worst", phoID_pfChIso03worst,"phoID_pfChIso03worst[2]/F");
    mtree_mainTree->Branch("phoID_E2x2", phoID_E2x2,"phoID_E2x2[2]/F");
    mtree_mainTree->Branch("phoID_E5x5", phoID_E5x5,"phoID_E5x5[2]/F");
  }

  /// dijet selection variables
  // skeleton
  //  mtree_dijetTree->Branch("x", &dijet_x ,"x/F"  );
  mtree_dijetTree->Branch("iMainEvt"  , &mtree_ievt       , "iMainEvt/I"  ); ///pointer to main tree
  mtree_dijetTree->Branch("mGG"       , &mtree_mass       , "mgg/F"       ); 
  mtree_dijetTree->Branch("mvaGG"     , &mtree_diphoMva   , "mvaGG/F"     ); 
  mtree_dijetTree->Branch("catBaseGG" , &mtree_catBase    , "catBaseGG/I" ); 
  mtree_dijetTree->Branch("dEtaJJ"    , &dijet_dEtaJJ     , "dEtaJJ/F"    );
  mtree_dijetTree->Branch("zeppenfeld", &dijet_zeppenfeld , "zeppenfeld/F");
  mtree_dijetTree->Branch("mJJ"       , &dijet_mJJ        , "mJJ/F"       );
  mtree_dijetTree->Branch("dphiGGJJ"  , &dijet_dphiGGJJ   , "dphiGGJJ/F"  );
  mtree_dijetTree->Branch("vbfTag"    , &mtree_vbfTag     , "vbfTag/I"    );
  mtree_dijetTree->Branch("vbfCat"    , &mtree_vbfCat     , "vbfCat/I"    );
  mtree_dijetTree->Branch("hvjjTag"   , &mtree_hvjjTag    , "hvjj/I"      );
  mtree_dijetTree->Branch("mvaVbf"       , &dijet_mvaVbf        , "mvaVbf/F"       );
  mtree_dijetTree->Branch("mvaVbf_sig"   , &dijet_mvaVbf_sig    , "mvaVbf_sig/F"   );
  mtree_dijetTree->Branch("mvaVbf_dipho" , &dijet_mvaVbf_dipho  , "mvaVbf_dipho/F" );
  mtree_dijetTree->Branch("mvaVbf_gluglu", &dijet_mvaVbf_gluglu , "mvaVbf_gluglu/F");

  mtree_dijetTree->Branch("jet1Pt", &dijet_ptj1 , "jet1Pt/F" );
  mtree_dijetTree->Branch("jet2Pt", &dijet_ptj2 , "jet2Pt/F" );
  mtree_dijetTree->Branch("etj", dijet_etj , "etj[2]/F" );
  mtree_dijetTree->Branch("ptj", dijet_ptj , "ptj[2]/F" );
  mtree_dijetTree->Branch("jec", dijet_jec , "jec[2]/F" );
  mtree_dijetTree->Branch("en" , dijet_en  , "en[2]/F"  );
  mtree_dijetTree->Branch("eta", dijet_eta , "eta[2]/F" );
  mtree_dijetTree->Branch("phi", dijet_phi , "phi[2]/F" );
  mtree_dijetTree->Branch("rms", dijet_rms , "rms[2]/F" );
  mtree_dijetTree->Branch("betaStar", dijet_betaStar , "betaStar[2]/F"  );

  mtree_dijetTree->Branch("massResoTot", &mtree_massResoTot, "massResoTot/F");
  mtree_dijetTree->Branch("massResoAng", &mtree_massResoAng, "massResoAng/F");
  mtree_dijetTree->Branch("massResoEng", &mtree_massResoEng, "massResoEng/F");
  mtree_dijetTree->Branch("massResoRightVtx", &mtree_massResoRightVtx , "massResoRightVtx/F" );
  mtree_dijetTree->Branch("massResoRandVtx" , &mtree_massResoRandVtx, "massResoRandVtx/F" );
  mtree_dijetTree->Branch("diphoMva" , &mtree_diphoMva, "diphoMva/F" );
  mtree_dijetTree->Branch("wei", &mc_wei, "wei/F");
  mtree_dijetTree->Branch("wPU", &mc_wPU, "wPU/F");
  mtree_dijetTree->Branch("wHQT", &mc_wHQT, "wHQT/F");
  mtree_dijetTree->Branch("wXsec", &mc_wXsec, "wXsec/F");
  mtree_dijetTree->Branch("wNgen", &mc_wNgen, "wNgen/F");
  mtree_dijetTree->Branch("wBSz", &mc_wBSz, "wBSz/F");
  mtree_dijetTree->Branch("wVtxId", &mc_wVtxId, "wVtxId/F");
  mtree_dijetTree->Branch("wPhoEffi" , &mc_wPhoEffi , "wPhoEffi/F");
  mtree_dijetTree->Branch("wTrigEffi", &mc_wTrigEffi, "wTrigEffi/F");

  mtree_dijetTree->Branch("mvaid1" , &mtree_mvaid1 , "mvaid1/F"    );
  mtree_dijetTree->Branch("mvaid2" , &mtree_mvaid2 , "mvaid2/F"    );
  mtree_dijetTree->Branch("eta1"   , &mtree_eta1   , "eta1/F"    );
  mtree_dijetTree->Branch("eta2"   , &mtree_eta2   , "eta2/F"    );
  mtree_dijetTree->Branch("r91"    , &mtree_r91    , "r91/F"    );
  mtree_dijetTree->Branch("r92"    , &mtree_r92    , "r92/F"    );
  mtree_dijetTree->Branch("mvaid"  , mtree_mvaid   , "mvaid[2]/F"    );
  mtree_dijetTree->Branch("s4Ratio", mtree_s4Ratio , "s4Ratio[2]/F"  );
  mtree_dijetTree->Branch("phoIso1", mtree_phoIso1 , "phoIso1[2]/F"  );
  mtree_dijetTree->Branch("phoIso2", mtree_phoIso2 , "phoIso2[2]/F"  );
  mtree_dijetTree->Branch("phoIso3", mtree_phoIso3 , "phoIso3[2]/F"  );
  mtree_dijetTree->Branch("coviEiE", mtree_coviEiE , "coviEiE[2]/F"  );
  mtree_dijetTree->Branch("coviEiP", mtree_coviEiP , "coviEiP[2]/F"  );
  mtree_dijetTree->Branch("r9", mtree_r9 ,"r9[2]/F"  );
  mtree_dijetTree->Branch("HoE", mtree_HoE ,"HoE[2]/F"  );
  mtree_dijetTree->Branch("drTrk", mtree_drTrk ,"drTrk[2]/F"  );



  /// lepton selection variables
  // skeleton
  mtree_elecTree->Branch("iMainEvt"  , &mtree_ievt       , "iMainEvt/I"  ); ///pointer to main tree
  mtree_elecTree->Branch("mGG"       , &mtree_mass       , "mgg/F"       ); 
  mtree_elecTree->Branch("mvaGG"     , &mtree_diphoMva   , "mvaGG/F"     ); 
  mtree_elecTree->Branch("catBaseGG" , &mtree_catBase    , "catBaseGG/I" ); 

  mtree_muonTree->Branch("iMainEvt"  , &mtree_ievt       , "iMainEvt/I"  ); ///pointer to main tree
  mtree_muonTree->Branch("mGG"       , &mtree_mass       , "mgg/F"       ); 
  mtree_muonTree->Branch("mvaGG"     , &mtree_mvaid      , "mvaGG/F"     ); 
  mtree_muonTree->Branch("catBaseGG" , &mtree_catBase    , "catBaseGG/I" ); 

  mtree_radion->Branch("evtNum", &radion_evtNum);
  mtree_radion->Branch("gamma1", &radion_gamma1);
  mtree_radion->Branch("gamma2", &radion_gamma2);
  mtree_radion->Branch("jets", &radion_jets);
  mtree_radion->Branch("bJetTags", &radion_bJetTags);
  mtree_radion->Branch("wei", &mc_wei);
}












//////////////////////////////////////////////////////////////////////////////////
/////////////////////// MiniTree reading specific functions //////////////////////
//////////////////////////////////////////////////////////////////////////////////
void MiniTree::setSynchVariables(void) {
  // type:0	run:190645	lumi:20	event:10264809	r9_1:8.9807e-01	sceta_1:2.1313e+00	hoe_1:3.1121e-02	sigieie_1:2.5966e-02	ecaliso_1:-6.4750e-01	hcal_1:4.6044e-02	trckiso_1:-1.0280e-01	chpfiso_1:0.0000e+00	r9_2:9.3988e-01	sceta_2:1.7397e+00	hoe_2:0.0000e+00	sigieie_2:2.5945e-02	ecaliso_2:-5.8803e-01	hcal_2:1.4802e+00	trckiso_2:-9.7490e-02	chpfiso_2:0.0000e+00	mgg:0.0000e+00	ptH:3.2530e+00	phoid_1:1.4830e-01	phoid_2:1.8015e-01	phoeta_1:2.1477e+00	phoeta_2:1.7567e+00	sigmrv:3.8392e+00	sigmwv:4.4829e+00	pt_1:5.1399e+01	pt_2:4.8745e+01	vtxprob:9.9487e-01	cosdphi:-6.6283e-01	diphoBDT:-4.2062e-01	e_1:2.2312e+02	e_2:1.4541e+02	vbfevt:0	evcat:0

  sync_iName.push_back("run"  ); sync_iVal.push_back( &mtree_runNum ); sync_mtreeName.push_back("runNum");
  sync_iName.push_back("lumi" ); sync_iVal.push_back( &mtree_lumiSec); sync_mtreeName.push_back("lumiSec");
  sync_iName.push_back("evcat"); sync_iVal.push_back( &mtree_catMva ); sync_mtreeName.push_back("catMva");
  sync_iName.push_back("evcatBase"); sync_iVal.push_back( &mtree_catBase ); sync_mtreeName.push_back("catBase");
  sync_lName.push_back("event"); sync_lVal.push_back( &mtree_evtNum ); sync_mtreeName.push_back("evtNum");
  sync_fName.push_back("rho"  ); sync_fVal.push_back( &mtree_rho    ); sync_mtreeName.push_back("rho");

  string phoInd[2] = {"_1","_2"};
  for( int i = 0; i < 2; i++ ) {
    sync_fName.push_back("r9"      + phoInd[i] ); sync_fVal.push_back( &mtree_r9[i] );     if(i==0 ) sync_mtreeName.push_back("r9");
    sync_fName.push_back("sceta"   + phoInd[i] ); sync_fVal.push_back( &mtree_sceta[i] );  if(i==0 ) sync_mtreeName.push_back("sceta");
    sync_fName.push_back("hoe"     + phoInd[i] ); sync_fVal.push_back( &mtree_HoE[i] );     if(i==0 ) sync_mtreeName.push_back("HoE");
    sync_fName.push_back("sigieie" + phoInd[i] ); sync_fVal.push_back( &mtree_coviEiE[i] );  if(i==0 ) sync_mtreeName.push_back("coviEiE");
    sync_fName.push_back("ecaliso" + phoInd[i] ); sync_fVal.push_back( &mtree_etCorrEcalIso[i] );  if(i==0 ) sync_mtreeName.push_back("etCorrEcalIso");
    sync_fName.push_back("hcaliso" + phoInd[i] ); sync_fVal.push_back( &mtree_etCorrHcalIso[i] ); if(i==0 ) sync_mtreeName.push_back("etCorrHcalIso");
    sync_fName.push_back("trckiso" + phoInd[i] ); sync_fVal.push_back( &mtree_etCorrTrkIso[i]  ); if(i==0 ) sync_mtreeName.push_back("etCorrTrkIso");
    sync_fName.push_back("chpfiso2"+ phoInd[i] ); sync_fVal.push_back( &mtree_pfChgIso02_selVtx[i]  ); if(i==0 ) sync_mtreeName.push_back("pfChgIso02_selVtx");
    sync_fName.push_back("phoid"   + phoInd[i] ); sync_fVal.push_back( &mtree_mvaid[i]  ); if(i==0 ) sync_mtreeName.push_back("mvaid");
    sync_fName.push_back("phoeta"  + phoInd[i] ); sync_fVal.push_back( &mtree_eta[i]  ); if(i==0 ) sync_mtreeName.push_back("eta");
    sync_fName.push_back("pt"      + phoInd[i] ); sync_fVal.push_back( &mtree_ptg[i]  ); if(i==0 ) sync_mtreeName.push_back("ptg");
    sync_fName.push_back("e"       + phoInd[i] ); sync_fVal.push_back( &mtree_eg[i]  ); if(i==0 ) sync_mtreeName.push_back("eg");
    sync_fName.push_back("eerr"    + phoInd[i] ); sync_fVal.push_back( &mtree_relResOverE[i]  ); if(i==0 ) sync_mtreeName.push_back("relResOverE");
    sync_fName.push_back("eerrsmeared" + phoInd[i] ); sync_fVal.push_back( &mtree_relSmearing[i]  ); if(i==0 ) sync_mtreeName.push_back("relSmearing");
  }
  string cvInd[2] = {"1","2"};
  for( int i = 0; i < 2; i++ ) {
    sync_fName.push_back("convindex"   + cvInd[i] ); sync_fVal.push_back( &convInd[i] ); if(i==0 ) sync_mtreeName.push_back("convInd");
    sync_fName.push_back("convNtrk"    + cvInd[i] ); sync_fVal.push_back( &convNtk[i] ); if(i==0 ) sync_mtreeName.push_back("convNtk");
    //    sync_fName.push_back("convpt"  +  cvInd[i] ); sync_fVal.push_back( &convPt[i]  ); if(i==0 ) sync_mtreeName.push_back("convPt" );
    sync_fName.push_back("convChiProb" + cvInd[i] ); sync_fVal.push_back( &convChi2Prob[i]  ); if(i==0 ) sync_mtreeName.push_back("convChi2Prob" );
  }
  sync_fName.push_back("mgg"     ); sync_fVal.push_back( &mtree_mass ); sync_mtreeName.push_back("mass");
  sync_fName.push_back("ptH"     ); sync_fVal.push_back( &mtree_pt ); sync_mtreeName.push_back("pt");
  sync_fName.push_back("sigmrv"  ); sync_fVal.push_back( &mtree_massResoRightVtx); sync_mtreeName.push_back("massResoRightVtx");
  sync_fName.push_back("sigmwv"  ); sync_fVal.push_back( &mtree_massResoRandVtx ); sync_mtreeName.push_back("massResoRandVtx");
  sync_fName.push_back("vtxprob" ); sync_fVal.push_back( &mtree_vtxProb ); sync_mtreeName.push_back("vtxProb");
  sync_fName.push_back("diphoBDT" ); sync_fVal.push_back( &mtree_diphoMva ); sync_mtreeName.push_back("diphoMva");

  string vtxInd[3] = {"1","2","3"};
  for( int iv = 0; iv < 3; iv++ ) {
    sync_fName.push_back("vertexId"    + vtxInd[iv] ); sync_fVal.push_back( &vtxId[iv] ); if( iv == 0 ) sync_mtreeName.push_back("vId");
    sync_fName.push_back("vertexMva"   + vtxInd[iv] ); sync_fVal.push_back( &vtxMva[iv]); if( iv == 0 ) sync_mtreeName.push_back("vMva");
    sync_fName.push_back("vertexdeltaz"+ vtxInd[iv] ); sync_fVal.push_back( &vtxDz[iv] ); if( iv == 0 ) sync_mtreeName.push_back("vDz");
  }
  sync_fName.push_back("ptbal"  ); sync_fVal.push_back( &vtxPtBal ); sync_mtreeName.push_back("vPtBal");
  sync_fName.push_back("ptasym" ); sync_fVal.push_back( &vtxAsym  ); sync_mtreeName.push_back("vAsym");
  sync_fName.push_back("logspt2"); sync_fVal.push_back( &vtxLogSumPt2 ); sync_mtreeName.push_back("vLogSumPt2");
  sync_fName.push_back("p2conv" ); sync_fVal.push_back( &vtxPullToConv); sync_mtreeName.push_back("vPullToConv");
  sync_fName.push_back("nconv"  ); sync_fVal.push_back( &vtxNConv    ); sync_mtreeName.push_back("vNConv");

}


void MiniTree::setBranchesAddresses(void) {
  cout << " ------------- MiniTree set input branches ------------" << endl;
    //----------------------- MC true
  mtree_mcTrue->Branch("mH"       , &mc_wei      );
  mtree_mcTrue->Branch("cThetaStarTruth_CS", &mc_cThetaStar_CS); 

  //----------------------- MC weights
  mtree_mainTree->SetBranchAddress("wei"      , &mc_wei      );
  mtree_mainTree->SetBranchAddress("wPU"      , &mc_wPU      );
  mtree_mainTree->SetBranchAddress("wHQT"     , &mc_wHQT     );
  mtree_mainTree->SetBranchAddress("wXsec"    , &mc_wXsec    );
  mtree_mainTree->SetBranchAddress("wNgen"    , &mc_wNgen    );
  mtree_mainTree->SetBranchAddress("wBSz"     , &mc_wBSz     );
  mtree_mainTree->SetBranchAddress("wVtxId"   , &mc_wVtxId   );
  mtree_mainTree->SetBranchAddress("wPhoEffi" , &mc_wPhoEffi );
  mtree_mainTree->SetBranchAddress("wTrigEffi", &mc_wTrigEffi);


  //-----------------------
  mtree_mainTree->SetBranchAddress("iMainEvt", &mtree_ievt    );
  mtree_mainTree->SetBranchAddress("runNum"  , &mtree_runNum  );
  mtree_mainTree->SetBranchAddress("evtNum"  , &mtree_evtNum  );
  mtree_mainTree->SetBranchAddress("lumiSec" , &mtree_lumiSec );
  mtree_mainTree->SetBranchAddress("rho"     , &mtree_rho     );
  mtree_mainTree->SetBranchAddress("rho25"   , &mtree_rho25   );

  mtree_mainTree->SetBranchAddress("zVtx"    , &mtree_zVtx  );
  mtree_mainTree->SetBranchAddress("nVtx"    , &mtree_nVtx  );
  mtree_mainTree->SetBranchAddress("iVtx1"   , &mtree_ivtx1 );
  mtree_mainTree->SetBranchAddress("iVtx2"   , &mtree_ivtx2 );
  mtree_mainTree->SetBranchAddress("iVtx3"   , &mtree_ivtx3 );
  mtree_mainTree->SetBranchAddress("vtxProb" , &mtree_vtxProb);
  mtree_mainTree->SetBranchAddress("vtxMva"  , &mtree_vtxMva );
  
  // skeleton
  //  mtree_mainTree->SetBranchAddress("x", &mtree_x ,"x/F"  );
  //--------------- diphoton event variables
  mtree_mainTree->SetBranchAddress("mass"      , &mtree_mass      );
  mtree_mainTree->SetBranchAddress("massDefVtx", &mtree_massDefVtx);
  mtree_mainTree->SetBranchAddress("pt"        , &mtree_pt        );
  mtree_mainTree->SetBranchAddress("y"         , &mtree_y         );
  mtree_mainTree->SetBranchAddress("piT"       , &mtree_piT       );
  mtree_mainTree->SetBranchAddress("pit1"      , &mtree_pit1      );
  mtree_mainTree->SetBranchAddress("pit2"      , &mtree_pit2      );
  mtree_mainTree->SetBranchAddress("massResoTot"     , &mtree_massResoTot      );
  mtree_mainTree->SetBranchAddress("massResoAng"     , &mtree_massResoAng      );
  mtree_mainTree->SetBranchAddress("massResoEng"     , &mtree_massResoEng      );
  mtree_mainTree->SetBranchAddress("massResoRightVtx", &mtree_massResoRightVtx );
  mtree_mainTree->SetBranchAddress("massResoRandVtx" , &mtree_massResoRandVtx  );
  mtree_mainTree->SetBranchAddress("diphoMva" , &mtree_diphoMva  );
  mtree_mainTree->SetBranchAddress("rawMet"   , &mtree_rawMet    );
  mtree_mainTree->SetBranchAddress("corMet"   , &mtree_corMet    );
  mtree_mainTree->SetBranchAddress("rawMetPhi", &mtree_rawMetPhi );
  mtree_mainTree->SetBranchAddress("corMetPhi", &mtree_corMetPhi );

  mtree_mainTree->SetBranchAddress("catMva" , &mtree_catMva  );
  mtree_mainTree->SetBranchAddress("catBase", &mtree_catBase );
  mtree_mainTree->SetBranchAddress("lepTag" , &mtree_lepTag  );
  mtree_mainTree->SetBranchAddress("vbfTag" , &mtree_vbfTag  );
  mtree_mainTree->SetBranchAddress("metTag" , &mtree_metTag  );
  mtree_mainTree->SetBranchAddress("hvjjTag", &mtree_hvjjTag );
  mtree_mainTree->SetBranchAddress("lepCat" , &mtree_lepCat  );
  mtree_mainTree->SetBranchAddress("vbfCat" , &mtree_vbfCat  );
  mtree_mainTree->SetBranchAddress("cThetaLead_heli" ,  &mtree_cThetaLead_heli );
  mtree_mainTree->SetBranchAddress("cThetaTrail_heli",  &mtree_cThetaTrail_heli);
  mtree_mainTree->SetBranchAddress("cThetaStar_CS"   ,  &mtree_cThetaStar_CS   );
  mtree_mainTree->SetBranchAddress("dEtaJJ"    , &dijet_dEtaJJ     );
  mtree_mainTree->SetBranchAddress("zeppenfeld", &dijet_zeppenfeld );
  mtree_mainTree->SetBranchAddress("mJJ"       , &dijet_mJJ        );
  mtree_mainTree->SetBranchAddress("dphiGGJJ"  , &dijet_dphiGGJJ   );
  mtree_mainTree->SetBranchAddress("jet1Pt"  , &dijet_ptj1   );
  mtree_mainTree->SetBranchAddress("jet2Pt"  , &dijet_ptj2   );

  mtree_mainTree->SetBranchAddress("eg" , mtree_eg );
  mtree_mainTree->SetBranchAddress("ptg", mtree_ptg);
  mtree_mainTree->SetBranchAddress("eta", mtree_eta);
  mtree_mainTree->SetBranchAddress("phi", mtree_phi);
  mtree_mainTree->SetBranchAddress("cat", mtree_cat);
  mtree_mainTree->SetBranchAddress("sceta", mtree_sceta);
  mtree_mainTree->SetBranchAddress("relResOverE", mtree_relResOverE);
  mtree_mainTree->SetBranchAddress("relSmearing" , mtree_relSmearing);
  mtree_mainTree->SetBranchAddress("isElec", mtree_isElec);

  
/// phoid variables
  mtree_mainTree->SetBranchAddress("mvaid"  , mtree_mvaid   );
  mtree_mainTree->SetBranchAddress("s4Ratio", mtree_s4Ratio );
  mtree_mainTree->SetBranchAddress("phoIso1", mtree_phoIso1 );
  mtree_mainTree->SetBranchAddress("phoIso2", mtree_phoIso2 );
  mtree_mainTree->SetBranchAddress("phoIso3", mtree_phoIso3 );
  mtree_mainTree->SetBranchAddress("coviEiE", mtree_coviEiE );
  mtree_mainTree->SetBranchAddress("coviEiP", mtree_coviEiP );
  mtree_mainTree->SetBranchAddress("r9", mtree_r9 );
  mtree_mainTree->SetBranchAddress("HoE", mtree_HoE );
  mtree_mainTree->SetBranchAddress("drTrk", mtree_drTrk );
  mtree_mainTree->SetBranchAddress("etCorrEcalIso", mtree_etCorrEcalIso);
  mtree_mainTree->SetBranchAddress("etCorrHcalIso", mtree_etCorrHcalIso);
  mtree_mainTree->SetBranchAddress("etCorrTrkIso" , mtree_etCorrTrkIso );
  mtree_mainTree->SetBranchAddress("esEffSigmaRR" , mtree_esEffSigmaRR );
  mtree_mainTree->SetBranchAddress("mvaid1"  , &mtree_mvaid1   );
  mtree_mainTree->SetBranchAddress("mvaid2"  , &mtree_mvaid2   );
  mtree_mainTree->SetBranchAddress("eta1"    , &mtree_eta1   );
  mtree_mainTree->SetBranchAddress("eta2"    , &mtree_eta2   );
  mtree_mainTree->SetBranchAddress("r91"     , &mtree_r91   );
  mtree_mainTree->SetBranchAddress("r92"     , &mtree_r92   );

  /// photon PF isolation
  mtree_mainTree->SetBranchAddress("pfChgIso02_selVtx", mtree_pfChgIso02_selVtx );
  mtree_mainTree->SetBranchAddress("pfChgIso03_selVtx", mtree_pfChgIso03_selVtx );
  mtree_mainTree->SetBranchAddress("pfChgIso03_badVtx", mtree_pfChgIso03_badVtx );
  mtree_mainTree->SetBranchAddress("pfChgIso03_defVtx", mtree_pfChgIso03_defVtx );
  mtree_mainTree->SetBranchAddress("pfChgIso04_selVtx", mtree_pfChgIso04_selVtx );
  mtree_mainTree->SetBranchAddress("pfChgIso04_badVtx", mtree_pfChgIso04_badVtx );
  mtree_mainTree->SetBranchAddress("pfChgIso04_defVtx", mtree_pfChgIso04_defVtx );
  mtree_mainTree->SetBranchAddress("pfPhoIso03", mtree_pfPhoIso03 );
  mtree_mainTree->SetBranchAddress("pfPhoIso04", mtree_pfPhoIso04 );

  /// add sync variables
  /// sync only variables
  if( addSyncVariables ) {
    mtree_mainTree->SetBranchAddress("convInd", convInd );
    mtree_mainTree->SetBranchAddress("convNtk", convNtk );
    mtree_mainTree->SetBranchAddress("convPt" , convPt  );
    mtree_mainTree->SetBranchAddress("convValidVtx", convValidVtx );
    mtree_mainTree->SetBranchAddress("convChi2Prob", convChi2Prob );
    mtree_mainTree->SetBranchAddress("vId"  , vtxId    );
    mtree_mainTree->SetBranchAddress("vDz"  , vtxDz    );
    mtree_mainTree->SetBranchAddress("vMva" , vtxMva   );
    mtree_mainTree->SetBranchAddress("vAsym", &vtxAsym );
    mtree_mainTree->SetBranchAddress("vPtBal"     , &vtxPtBal     );
    mtree_mainTree->SetBranchAddress("vLogSumPt2" , &vtxLogSumPt2 );
    mtree_mainTree->SetBranchAddress("vPullToConv", &vtxPullToConv);
    mtree_mainTree->SetBranchAddress("vNConv"     , &vtxNConv     );

    mtree_mainTree->Branch("phoID_covIEtaIPhi"   , phoID_covIEtaIPhi );
    mtree_mainTree->Branch("phoID_S4"            , phoID_S4          );
    mtree_mainTree->Branch("phoID_ESEffSigmaRR"  , phoID_ESEffSigmaRR);
    mtree_mainTree->Branch("phoID_pfPhoIso03"    , phoID_pfPhoIso03  );
    mtree_mainTree->Branch("phoID_pfChIso03"     , phoID_pfChIso03   );
    mtree_mainTree->Branch("phoID_pfChIso03worst", phoID_pfChIso03worst);
    mtree_mainTree->Branch("phoID_E2x2", phoID_E2x2);
    mtree_mainTree->Branch("phoID_E5x5", phoID_E5x5);
  }
    
  /// dijet selection variables
  // skeleton
  //  mtree_dijetTree->SetBranchAddress("x", &dijet_x ,"x/F"  );
  mtree_dijetTree->SetBranchAddress("iMainEvt"  , &mtree_ievt       );
  mtree_dijetTree->SetBranchAddress("dEtaJJ"    , &dijet_dEtaJJ     );
  mtree_dijetTree->SetBranchAddress("zeppenfeld", &dijet_zeppenfeld );
  mtree_dijetTree->SetBranchAddress("mJJ"       , &dijet_mJJ        );
  mtree_dijetTree->SetBranchAddress("dphiGGJJ"  , &dijet_dphiGGJJ   );
  mtree_dijetTree->SetBranchAddress("vbfTag"    , &mtree_vbfTag     );
  mtree_dijetTree->SetBranchAddress("vbfCat"    , &mtree_vbfCat     );
  mtree_dijetTree->SetBranchAddress("hvjjTag"   , &mtree_hvjjTag    );
  mtree_dijetTree->SetBranchAddress("mvaVbf"       , &dijet_mvaVbf        );
  mtree_dijetTree->SetBranchAddress("mvaVbf_sig"   , &dijet_mvaVbf_sig    );
  mtree_dijetTree->SetBranchAddress("mvaVbf_dipho" , &dijet_mvaVbf_dipho  );
  mtree_dijetTree->SetBranchAddress("mvaVbf_gluglu", &dijet_mvaVbf_gluglu );

  mtree_dijetTree->SetBranchAddress("etj", dijet_etj );
  mtree_dijetTree->SetBranchAddress("ptj", dijet_ptj );
  mtree_dijetTree->SetBranchAddress("jec", dijet_jec );
  mtree_dijetTree->SetBranchAddress("en" , dijet_en  );
  mtree_dijetTree->SetBranchAddress("eta", dijet_eta );
  mtree_dijetTree->SetBranchAddress("phi", dijet_phi );
  mtree_dijetTree->SetBranchAddress("rms", dijet_rms );
  mtree_dijetTree->SetBranchAddress("betaStar", dijet_betaStar );
/// phoid variables
  mtree_dijetTree->SetBranchAddress("jet1Pt"  , &dijet_ptj1   );
  mtree_dijetTree->SetBranchAddress("jet2Pt"  , &dijet_ptj2   );
  mtree_dijetTree->SetBranchAddress("mvaid1"  , &mtree_mvaid1   );
  mtree_dijetTree->SetBranchAddress("mvaid2"  , &mtree_mvaid2   );
  mtree_dijetTree->SetBranchAddress("eta1"    , &mtree_eta1   );
  mtree_dijetTree->SetBranchAddress("eta2"    , &mtree_eta2   );
  mtree_dijetTree->SetBranchAddress("r91"     , &mtree_r91   );
  mtree_dijetTree->SetBranchAddress("r92"     , &mtree_r92   );
  mtree_dijetTree->SetBranchAddress("mvaid"  , mtree_mvaid   );
  mtree_dijetTree->SetBranchAddress("s4Ratio", mtree_s4Ratio );
  mtree_dijetTree->SetBranchAddress("phoIso1", mtree_phoIso1 );
  mtree_dijetTree->SetBranchAddress("phoIso2", mtree_phoIso2 );
  mtree_dijetTree->SetBranchAddress("phoIso3", mtree_phoIso3 );
  mtree_dijetTree->SetBranchAddress("coviEiE", mtree_coviEiE );
  mtree_dijetTree->SetBranchAddress("coviEiP", mtree_coviEiP );
  mtree_dijetTree->SetBranchAddress("r9", mtree_r9 );
  mtree_dijetTree->SetBranchAddress("HoE", mtree_HoE );
  mtree_dijetTree->SetBranchAddress("drTrk", mtree_drTrk );

}



#endif
