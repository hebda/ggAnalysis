#define xAna_cxx
#include "UtilsDiPhoton.h"
#include "TriggerSelection.h"
#include "MvaComputation.h"
#include "EnergyCorrections.h"
#include "EnergyScaleReader.h"
#include "JetSelection.h"
#include "METSelection.h"
#include "LeptonSelection.h"
#include "PhotonId.h"
#include "MinitreeOut.h"
#include "MCtruthFunctions.h"
//#include "EventId.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <Math/VectorUtil.h>
#include <TGraphAsymmErrors.h>

#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace std;

#include "TSystem.h" // For memory check
ProcInfo_t pi; // For memory check



int xAna::HggTreeWriteLoop(const char* filename, int ijob, 
					  bool correctVertex, bool correctEnergy, 
					  bool setRho0
					  ) {

  //bool invDRtoTrk = false;


  if( _config == 0 ) {
    cout << " config file was not set properly... bail out" << endl;
    return -1;
  }

  if (fChain == 0) return -1;
  
  Long64_t nentries = fChain->GetEntriesFast();
  //  nentries = 10000;
  cout << "nentries: " << nentries << endl;  
  Long64_t entry_start = ijob    *_config->nEvtsPerJob();
  Long64_t entry_stop  = (ijob+1)*_config->nEvtsPerJob();
  if( _config->nEvtsPerJob() < 0 ) {
    entry_stop = nentries;
  }
  if( entry_stop  > nentries ) entry_stop = nentries;
  cout << "   *** doing entries from: " << entry_start << " -> " << entry_stop << endl;
  if( entry_start > entry_stop ) return -1;

  
  EnergyScaleReader enScaleSkimEOS; /// skim EOS bugged so need to undo the energy scale in the skim
  EnergyScaleReader enScale;
  //  enScaleSkimEOS.setup( "ecalCalibFiles/EnergyScale2012_Lisbon_9fb.txt" );
  enScale.setup( _config->energyScaleFile() );


  Float_t HiggsMCMass =  _weight_manager->getCrossSection()->getHiggsMass();
  Float_t HiggsMCPt   = -1;
  bool isHiggsSignal = false;
  if( HiggsMCMass > 0 ) isHiggsSignal = true;

  mode_ = _weight_manager->getCrossSection()->mode();
  isData = false;
  if(mode_==-1) isData = true;        

//  mode_ = ijob; //  

  DoCorrectVertex_ = correctVertex;
  DoCorrectEnergy_ = correctEnergy;
  DoSetRho0_ = setRho0;

  doJetRegression = _config->getDoJetRegression();
  doControlSample = _config->getDoControlSample();
  
  phoID_2011[0] = new TMVA::Reader("!Color:Silent");
  phoID_2011[1] = new TMVA::Reader("!Color:Silent");
  phoID_2012[0] = new TMVA::Reader("!Color:Silent");
  phoID_2012[1] = new TMVA::Reader("!Color:Silent");
  DiscriDiPho_2011 = new TMVA::Reader("!Color:Silent");
  DiscriDiPho_2012 = new TMVA::Reader("!Color:Silent");
  if(doJetRegression==1) {
    jetRegres0 = new TMVA::Reader("!Color:Silent");
    jetRegres1 = new TMVA::Reader("!Color:Silent");
  }
  else if(doJetRegression==2) jetRegres = new TMVA::Reader("!Color:Silent");
  Setup_MVA();
  
  if( _config->setup() == "ReReco2011" )  for( int i = 0 ; i < 2; i++ ) phoID_mva[i] = phoID_2011[i];
  else                                    for( int i = 0 ; i < 2; i++ ) phoID_mva[i] = phoID_2012[i];

  MassResolution massResoCalc;
  massResoCalc.addSmearingTerm();
  if( _config->setup() == "ReReco2011" )  Ddipho_mva = DiscriDiPho_2011;
  else                                    Ddipho_mva = DiscriDiPho_2012;

  float Ddipho_cat[5]; Ddipho_cat[4] = -1; 
  if( _config->setup() == "ReReco2011" ) { Ddipho_cat[0] = 0.89; Ddipho_cat[1] = 0.72; Ddipho_cat[2] = 0.55; Ddipho_cat[3] = +0.05; }
  else                                   { Ddipho_cat[0] = 0.91; Ddipho_cat[1] = 0.79; Ddipho_cat[2] = 0.49; Ddipho_cat[3] = -0.05; }
  //  else                                   { Ddipho_cat[0] = 0.88; Ddipho_cat[1] = 0.71; Ddipho_cat[2] = 0.50; Ddipho_cat[3] = -0.05; }

  DiscriVBF_UseDiPhoPt = true;
  DiscriVBF_UsePhoPt   = true;
  DiscriVBF_cat.resize(2); DiscriVBF_cat[0] = 0.985; DiscriVBF_cat[1] = 0.93;
  DiscriVBF_useMvaSel  = _config->doVBFmvaCat();


  


  /// depending on the selection veto or not on electrons (can do muele, elemu,eleele)
  bool vetoElec[2] = {true,true};
  if( _config->invertElectronVeto() ) { vetoElec[0] = false; vetoElec[1] = false; }
  if( _config->isEleGamma()         ) { vetoElec[0] = false; vetoElec[1] = true ; }
  if( _config->isGammaEle()         ) { vetoElec[0] = true ; vetoElec[1] = false; }
  cout << " --------- veto electron config -----------" << endl;
  cout << " Leading  Pho VetoElec: " << vetoElec[0] << endl;
  cout << " Trailing Pho VetoElec: " << vetoElec[1] << endl;
  

  DoDebugEvent = true;
  
  bool DoPreselection = true;
  //  bool DoPrint = true;
  
  
  TString VertexFileNamePrefix;
  
  TRandom3 *rnd = new TRandom3();
  rnd->SetSeed(0);
  
  /// output tree and cross check file
   _xcheckTextFile.open(TString(filename)+".xcheck.txt");
   // _xcheckTextFile = cout;

  _minitree = new MiniTree( filename );
  TTree * tSkim = 0;
  if( _config->doSkimming() ) tSkim = (TTree*) fChain->CloneTree(0);

  InitHists();
  _minitree->mc_wXsec = _weight_manager->xSecW();
  _minitree->mc_wNgen = 100000./_weight_manager->getNevts();
  if( isData ) {
    _minitree->mc_wXsec = 1;
    _minitree->mc_wNgen = 1;
  }
    
  Int_t isprompt0 = -1;
  Int_t isprompt1 = -1;
  
  set<Long64_t> syncEvt;

  cout <<" ================ mode  "  << mode_   <<" ===============================  "<<endl;
  /// setupType has to be passed via config file
  //  photonOverSmearing overSmearICHEP("Test52_ichep");
  photonOverSmearing overSmearHCP( "oversmear_hcp2012" );
  photonOverSmearing overSmear(    _config->setup()    );
  int overSmearSyst = _config->getEnergyOverSmearingSyst();

  Long64_t nbytes = 0, nb = 0;
  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////// Start the loop ////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  vector<int>    nEvts;
  vector<string> nCutName;
  nCutName.push_back("TOTAL          :"); nEvts.push_back(0);
  nCutName.push_back("2 gammas       :"); nEvts.push_back(0);
  nCutName.push_back("triggers       :"); nEvts.push_back(0);
  nCutName.push_back("nan weight     :"); nEvts.push_back(0);
  nCutName.push_back("presel kin cuts:"); nEvts.push_back(0);
  nCutName.push_back("pass 2 gam incl:"); nEvts.push_back(0);
  nCutName.push_back("pass all       :"); nEvts.push_back(0);
  nCutName.push_back("   --> pass 2 gam lep: "); nEvts.push_back(0);
  nCutName.push_back("   --> pass 2 gam vbf: "); nEvts.push_back(0);

  for (Long64_t jentry=entry_start; jentry< entry_stop ; ++jentry) {
    Long64_t ientry = LoadTree(jentry);

    if (ientry < -9999) break;
    if( jentry % 10000 == 0) cout<<"processing event "<<jentry<<endl;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    /// reset minitree variables
    _minitree->initEvent();

    /// study mc truth block
    //if( !isData ) {
      //fillMCtruthInfo_HiggsSignal();//This line causes things to break in MC.
      //_minitree->fillMCtrueOnly();
    //}

    /// reco analysis
    unsigned icutlevel = 0;
    nEvts[icutlevel++]++;
    if( nPho < 2 ) continue; 
    nEvts[icutlevel++]++;
    
    /// set synchronisation flag
    if( syncEvt.find(event) != syncEvt.end() ) DoDebugEvent = true;
    else                                       DoDebugEvent = false;
    if( DoDebugEvent ) _xcheckTextFile << "==========================================================" << endl;
    if( DoDebugEvent ) _xcheckTextFile << "================= debugging event: " << event << endl;
    
    /// PU & BSz reweightings
    if( !isData ) {
      int itpu = 1;  /// 0 without OOT PU - 1 with OOT PU
      _minitree->mc_wPU  = _weight_manager->puW( nPU[itpu] );
     // PUwei  = _weight_manager->puWTrue( puTrue[itpu] );      
    }
    hTotEvents->Fill(1,_minitree->mc_wPU);


    
    bool sigWH = false ;
    bool sigZH = false ;    
    int  mc_whzh_type = 0;
    _minitree->mc_wHQT = 1;
    if( isHiggsSignal ) {
      if ( _weight_manager->getCrossSection()->getMCType() == "vh" ) 
      for( Int_t i=0; i < nMC && i <= 1; ++i ) { 
	  if( abs(mcPID[i]) == 24 ) sigWH=true;
	  if( abs(mcPID[i]) == 23 ) sigZH=true;
	}
      if( sigWH ) mc_whzh_type = 1;
      if( sigZH ) mc_whzh_type = 2;
      
      for( Int_t i=0; i<nMC; ++i)
      if ( abs(mcPID[i]) == 25 ) HiggsMCPt   = mcPt[i];

      if( _weight_manager->getCrossSection()->getMCType() == "ggh" && 
	  _config->setup().find( "ReReco2011" ) != string::npos )
	_minitree->mc_wHQT = getHqTWeight(HiggsMCMass, HiggsMCPt);      
    }
    


    
    if ((mode_ == 2 || mode_ ==1 ||  mode_ == 18 || mode_ == 19) && processID==18) continue;      
    
    // Remove double counting in gamma+jets and QCDjets
    
    Int_t mcIFSR_pho = 0;
    Int_t mcPartonic_pho = 0;
    if (mode_ == 1 || mode_ == 2 || mode_ == 3  ||  mode_ == 18 ||  mode_ == 19 ) {
      for (Int_t i=0; i<nMC; ++i) {    
	if (mcPID[i] == 22 && (fabs(mcMomPID[i]) < 6 || mcMomPID[i] == 21)) mcIFSR_pho++;
	if (mcPID[i] == 22 && mcMomPID[i] == 22) mcPartonic_pho++;
      }
    }
    
    // if pythia is used for diphoton.. no IFSR removing from QCD and Gjets!!!!!!      
    if ((mode_==1 || mode_ == 2  ||   mode_ == 18 ) && mcIFSR_pho >= 1 && mcPartonic_pho >=1) continue;
    if ((mode_ == 3 ||  mode_ == 19 )&& mcIFSR_pho == 2) continue;   

 
    bool prompt2= false;
    bool prompt1= false;
    bool prompt0= false;
    vertProb = -1;
    
    if (mode_ == 2  || mode_ == 1    ||   mode_ == 18     ){
      if      ( mcPartonic_pho >= 1 &&  mcIFSR_pho >= 1 ) prompt2 = true;
      else if ( mcPartonic_pho >= 1 &&  mcIFSR_pho == 0 ) prompt1 = true;
      else if ( mcPartonic_pho == 0 &&  mcIFSR_pho == 0 ) prompt0 = true;
    } else if(mode_ == 3 ||  mode_ == 19   ){
      if      ( mcIFSR_pho >= 2 ) prompt2 = true;
      else if ( mcIFSR_pho == 1 ) prompt1 = true;
      else if ( mcIFSR_pho == 0 ) prompt0 = true;
      if( prompt1 ) _minitree->mc_wXsec = 1.3*_weight_manager->xSecW();
    }
    
    if(mode_==1 || mode_==2 || mode_==3 || mode_==18 || mode_==19){
      if(prompt0)isprompt0=1;
      else isprompt0=0;
      if(prompt1)isprompt1=1;
      else isprompt1=0;
    }
    
    if( mode_ == 20 && isZgamma() ) continue;

    /// wei weight is just temporary and may not contain all info.
    float wei = _minitree->mc_wXsec * _minitree->mc_wPU  
      * _minitree->mc_wNgen * _minitree->mc_wHQT;

    if( isData && !PassTriggerSelection() ) continue;      nEvts[icutlevel++]++;
    if( std::isinf( wei ) || std::isnan( wei ) )continue;  nEvts[icutlevel++]++;
    
    
    
    //// ********************* define S4 variable **********************////
    for( int i=0; i<nPho; ++i){
      if( _config->setup() == "ReReco2011" ) phoS4ratio[i] = 1;
      else                                   phoS4ratio[i] = phoE2x2[i] / phoE5x5[i];
    }
    //// ************************************************************* ////

    if( !isData ) {
      //// ************** MC corrections (mostly MC) ******************* //// 
      // 1. energy shifting / smearing
      for( int i=0; i<nPho; ++i)
	if( fabs(phoSCEta[i]) <= 2.5 ) {
	  float smearing = overSmear.randOverSmearing(phoSCEta[i],phoR9[i],isInGAP_EB(i),overSmearSyst);
	  phoRegrE[i] *= (1 + smearing); 
	  phoE[i]     *= (1 + smearing); 
	  /// from MassFactorized in gglobe:   energyCorrectedError[ipho] *=(l.pho_isEB[ipho]) ? 1.07 : 1.045 ;
	  float smearFactor = 1;
	  if( _config->setup() == "ReReco2011" ) smearFactor = fabs(phoSCEta[i]) < 1.45 ? 1.07: 1.045;
	  phoRegrErr[i] *= smearFactor;
	}  
      
      
      // 2. reweighting of photon id variables (R9...)
      for (int i=0; i<nPho; ++i) ReweightMC_phoIdVar(i);
      //// ************************************************************* ////
    }
    
    //// ********** Apply regression energy ************* ////
    float phoStdE[500];
    for( int i=0; i<nPho; ++i)
      if( fabs(phoSCEta[i]) <= 2.5 ) {
	if( isData ){
	  float enCorrSkim = 1;//enScaleSkimEOS.energyScale( phoR9[i], phoSCEta[i], run);
	  float phoEnScale = enScale.energyScale( phoR9[i], phoSCEta[i], run)/enCorrSkim;
	  phoRegrE[i]  *= phoEnScale;
	  phoE[i]      *= phoEnScale;
	}
	phoStdE[i] = phoE[i];
	phoE[i]   = phoRegrE[i];
	
	/// transform calo position abd etaVtx, phiVtx with SC position
	for( int x = 0 ; x < 3; x++ )  phoCaloPos[i][x] = phoSCPos[i][x];
	for( int ivtx = 0 ; ivtx < nVtxBS; ivtx++ ) {
	  TVector3 xxi = getCorPhotonTVector3(i,ivtx);
	  phoEtaVtx[i][ivtx] = xxi.Eta();
	  phoPhiVtx[i][ivtx] = xxi.Phi();
	}
	
	/// additionnal smearing to go to data energy resolution
	phoRegrSmear[i] = phoE[i]*overSmearHCP.meanOverSmearing(phoSCEta[i],phoR9[i],isInGAP_EB(i),0);
	phoEt[i] = phoE[i] / cosh(phoEta[i]);  
      }    
    //// ************************************************* ////
     

    /// lepton selection
    int iElecVtx(-1), iMuonVtx(-1);
    vector<int> elecIndex = selectElectronsHCP2012( wei, iElecVtx );
    vector<int> muonIndex = selectMuonsHCP2012(     wei, iMuonVtx );

    vector<int>             event_vtx; 
    vector<int>             event_ilead ;
    vector<int>             event_itrail;
    vector<TLorentzVector>  event_plead ;
    vector<TLorentzVector>  event_ptrail;
    TLorentzVector leptag;
    
    int lepCat = -1;
    bool exitLoop = false;
    for( int ii = 0     ; ii < nPho ; ++ii ) {
      for( int jj = (ii+1); jj < nPho ; ++jj ) {
	// Preselection 2nd leg
	if (DoPreselection && !preselectPhoton(ii,phoEt[ii])) continue;
	if (DoPreselection && !preselectPhoton(jj,phoEt[jj])) continue;
	
	/// define i, j locally, so when they are inverted this does not mess up the loop
	int i = ii; 
	int j = jj;
	if(phoEt[j] > phoEt[i]){ i = jj; j = ii; }	  
	
	// Select vertex
	int selVtx  = 0;
	TLorentzVector gi,gj;
	double mij(-1);
	
	vector<int> selVtxIncl;      
	selVtxIncl = getSelectedVertex(i,j,true, true);
	selVtx = selVtxIncl[0];
	if( selVtx < 0 ) continue;

	/// check lepton tag
	if( muonIndex.size() > 0 ) {
	  //selVtx = iMuonVtx;	  
	  leptag.SetPtEtaPhiM( muPt[muonIndex[0]], muEta[muonIndex[0]], muPhi[muonIndex[0]],0);
	  if( selectTwoPhotons(i,j,selVtx,gi,gj,vetoElec) ) {
	    lepCat = 0;
	  }
	}
	
	/// check electron tag only if muon tag failed
	if( elecIndex.size() > 0 && lepCat == -1 ) {
	  //selVtx = iElecVtx;
	  leptag.SetPtEtaPhiM( elePt[elecIndex[0]], eleEta[elecIndex[0]], elePhi[elecIndex[0]],0);
	  if( selectTwoPhotons(i,j,selVtx,gi,gj,vetoElec) ) {
	    lepCat = 1;
	    if( fabs( (leptag+gi).M() - 91.19 ) < 10 || fabs( (leptag+gj).M() - 91.19 ) < 10 ) lepCat = -1;
	    /// this is not actually the cut, but it should be but ok (dR(pho,eleTrk) > 1 no dR(pho,anyTrk) > 1 )
	    if(phoCiCdRtoTrk[i] < 1 || phoCiCdRtoTrk[j] <1 ) lepCat = -1;
	  }
	}
	

	if( lepCat >= 0 ) {
	  mij = (gi+gj).M();
          if( _config->analysisType() == "MVA" )
	        if( gi.Pt() / mij < 45./120. || gj.Pt() / mij < 30./120. ) lepCat = -1;
          if( _config->analysisType() == "baselineCiC4PF" )
	        if( gi.Pt() / mij < 45./120. || gj.Pt() < 25. ) lepCat = -1;
	  if( leptag.DeltaR(gi) < 1.0  || leptag.DeltaR(gj) < 1.0  ) lepCat = -1;
	}
	if( lepCat >= 0 ) {
	  cout << " ****** keep leptag event pts[photons] i: " << i << " j: " << j << "   -> event: " << ientry <<  endl;
	  cout << "        leptonCat: " << lepCat << endl;
	  /// if one pair passes lepton tag then no need to go further
	  /// fill in variables
	  event_vtx.resize(   1); event_vtx[0]    = selVtx;
	  event_ilead.resize( 1); event_ilead[0]  = i;
	  event_itrail.resize(1); event_itrail[0] = j;
	  event_plead.resize( 1); event_plead[0]  = gi;
	  event_ptrail.resize(1); event_ptrail[0] = gj;
	  exitLoop = true;
	  break;	
	} else {
	  /// inclusive + VBF + MetTag preselection
	  
	  /// apply kinematic cuts
	  if( !selectTwoPhotons(i,j,selVtx,gi,gj,vetoElec) ) continue;
	  /// drop photon pt cuts from inclusive cuts (add them in categorisation only)
	  mij = (gi+gj).M();	
	  if( ! _config->doSkimming() && 
	      ( gi.Pt()     < _config->pt1Cut() || gi.Pt() < _config->pt2Cut() ||
		gi.Pt()/mij < _config->loosePt1MCut() ||
		gj.Pt()/mij < _config->loosePt2MCut() ) ) continue;
	}
	
	/// here i = lead; j = trail (selectTwoPhotons does that)
	/// fill in variables
	event_vtx.push_back( selVtx );
	event_ilead.push_back(   i  );
	event_itrail.push_back(  j  );
	event_plead.push_back(  gi  );
	event_ptrail.push_back( gj  );   	  
      }
      if( exitLoop ) break;
    }
      
    if(event_ilead.size()==0 || event_itrail.size()==0)continue;
    if(event_vtx.size() == 0 ) continue;  // no pairs selected
    nEvts[icutlevel++]++;

    // Now decide which photon-photon pair in the event to select
    // for lepton tag (size of arrays is 1, so dummy selection)
    unsigned int selectedPair = 0;
    float sumEt = -1;
    for (unsigned int p=0; p < event_vtx.size(); p++) {
      float tempSumEt = event_plead[p].Pt() + event_ptrail[p].Pt();
      if( tempSumEt > sumEt) {
	sumEt = tempSumEt;
      selectedPair = p;
      }
    } 
    

    int ilead  = event_ilead[  selectedPair ];
    int itrail = event_itrail[ selectedPair ];
    int selVtx = event_vtx[    selectedPair ];
    TLorentzVector glead  = event_plead[selectedPair];
    TLorentzVector gtrail = event_ptrail[selectedPair];
    TLorentzVector hcand  = glead + gtrail;

    int CAT4 = cat4(phoR9[ilead], phoR9[itrail], phoSCEta[ilead], phoSCEta[itrail]);
    
    if( glead.Pt()  < _config->pt1Cut() ) continue;
    if( gtrail.Pt() < _config->pt2Cut() ) continue;
    if( hcand.M()   < _config->mggCut() ) continue;
    nEvts[icutlevel++]++;

    _minitree->mtree_runNum = run;
    _minitree->mtree_evtNum = event;
    _minitree->mtree_lumiSec = lumis;

    //    TLorentzVector g1,g2;
    //    g1.SetPtEtaPhiM(phoE[ilead ]/cosh(phoEta[ilead]),phoEta[ilead ], phoPhi[ilead ], 0);
    //    g2.SetPtEtaPhiM(phoE[itrail]/cosh(phoEta[ilead]),phoEta[itrail], phoPhi[itrail], 0);
    //    TLorentzVector lgg = g1 + g2;
    //_minitree->mtree_massDefVtx = lgg.M();
    _minitree->mtree_massDefVtx = hcand.M()/sqrt(phoE[ilead]*phoE[itrail])*sqrt(phoStdE[ilead]*phoStdE[itrail]);

    // calc again vertex to get correct vertex probability (can be different from last one in loop)
    vector<int> sortedVertex = getSelectedVertex( ilead, itrail, true, true);
    if( sortedVertex.size() < 2 ) sortedVertex.push_back(-1);
    if( sortedVertex.size() < 3 ) sortedVertex.push_back(-1);

    if( lepCat >= 0 ) {
      /// lepton tag
      sortedVertex[0] = selVtx;  
      sortedVertex[1] = -1; 
      sortedVertex[2] = -1; 
      //      vertProb = 1;
    } 

    _minitree->mtree_rho   = rho2012;
    _minitree->mtree_rho25 = rho25;    
    _minitree->mtree_zVtx  = vtxbs[selVtx][2];
    _minitree->mtree_nVtx  = nVtxBS;
    _minitree->mtree_ivtx1 = selVtx;
    _minitree->mtree_ivtx2 = sortedVertex[1];
    _minitree->mtree_ivtx3 = sortedVertex[2];
    _minitree->mtree_vtxProb = vertProb;
    _minitree->mtree_vtxMva  = vertMVA;

    _minitree->mtree_mass  = hcand.M();
    _minitree->mtree_pt    = hcand.Pt();
    _minitree->mtree_piT   = hcand.Pt()/hcand.M();
    _minitree->mtree_y     = hcand.Rapidity();

    /// spin variables
    TLorentzVector gtmp1 = glead;  gtmp1.Boost( -hcand.BoostVector() );
    TLorentzVector gtmp2 = gtrail; gtmp2.Boost( -hcand.BoostVector() );
    _minitree->mtree_cThetaLead_heli  = cos( gtmp1.Angle(hcand.BoostVector()) );
    _minitree->mtree_cThetaTrail_heli = cos( gtmp2.Angle(hcand.BoostVector()) );
    _minitree->mtree_cThetaStar_CS    = 2*(glead.E()*gtrail.Pz() - gtrail.E()*glead.Pz())/(hcand.M()*sqrt(hcand.M2()+hcand.Pt()*hcand.Pt()));
    
    /// fill photon id variables in main tree
    _minitree->mtree_minR9      = +999;
    _minitree->mtree_minPhoIdEB = +999;
    _minitree->mtree_minPhoIdEE = +999;
    _minitree->mtree_maxSCEta   =   -1;
    _minitree->mtree_minSCEta   = +999;
    for( int i = 0 ; i < 2; i++ ) {
      int ipho = -1;
      if( i == 0 ) ipho = ilead;
      if( i == 1 ) ipho = itrail;
      fillPhotonVariablesToMiniTree( ipho, selVtx, i );      
      if( _minitree->mtree_r9[i] < _minitree->mtree_minR9 ) _minitree->mtree_minR9 = _minitree->mtree_r9[i];
      if( fabs( _minitree->mtree_sceta[i] ) <  1.5 &&  _minitree->mtree_mvaid[i] <  _minitree->mtree_minPhoIdEB ) _minitree->mtree_minPhoIdEB = _minitree->mtree_mvaid[i];
      if( fabs( _minitree->mtree_sceta[i] ) >= 1.5 &&  _minitree->mtree_mvaid[i] <  _minitree->mtree_minPhoIdEE ) _minitree->mtree_minPhoIdEE = _minitree->mtree_mvaid[i];
      if( fabs( _minitree->mtree_sceta[i] ) > _minitree->mtree_maxSCEta ) _minitree->mtree_maxSCEta =  fabs(_minitree->mtree_sceta[i]);
      if( fabs( _minitree->mtree_sceta[i] ) < _minitree->mtree_minSCEta ) _minitree->mtree_minSCEta =  fabs(_minitree->mtree_sceta[i]);
    }
    
    //------------ compute diphoton mva (add var to minitree inside function) ----------------//
    massResoCalc.setP4CalPosVtxResoSmear( glead,gtrail, 
					  TVector3(phoCaloPos[ilead ][0], phoCaloPos[ilead ][1],phoCaloPos[ilead ][2]),
					  TVector3(phoCaloPos[itrail][0], phoCaloPos[itrail][1],phoCaloPos[itrail][2]),
					  TVector3(vtxbs[selVtx][0], vtxbs[selVtx][1], vtxbs[selVtx][2]),
					  _minitree->mtree_relResOverE, _minitree->mtree_relSmearing );

    
    _minitree->mtree_massResoTot = massResoCalc.relativeMassResolutionFab_total( vertProb );    
    _minitree->mtree_massResoEng = massResoCalc.relativeMassResolutionFab_energy( );    
    _minitree->mtree_massResoAng = massResoCalc.relativeMassResolutionFab_angular();    
   
    float diphotonmva = DiPhoID_MVA( glead, gtrail, hcand, massResoCalc, vertProb, 
				     _minitree->mtree_mvaid[0],_minitree->mtree_mvaid[1] );

    
    // ---- Start categorisation ---- //
    _minitree->mtree_lepTag  = 0;
    _minitree->mtree_metTag  = 0;
    _minitree->mtree_vbfTag  = 0;
    _minitree->mtree_hvjjTag = 0;

    // 1. lepton tag
    if( lepCat >= 0 ) {
      _minitree->mtree_lepTag = 1;
      _minitree->mtree_lepCat = lepCat;
      _minitree->mtree_fillLepTree  = true;
      // if( lepCat == 0 ) fillMuonTagVariables(lepTag,glead,gtrail,elecIndex[0],selVtx);
      // if( lepCat == 1 ) fillElecTagVariables(lepTag,glead,gtrail,muonIndex[0],selVtx);
    }

    // 3. met tag (For one flavor of jet energy regression, MET needs to be corrected first.)
    _minitree->mtree_rawMet    = recoPfMET;
    _minitree->mtree_rawMetPhi = recoPfMETPhi;
    //    if( !isData ) {
      /// bug in data skim, met correction already applied
      //3.1 Soft Jet correction (FC?)    
      applyJEC4SoftJets();
      //3.2 smearing
      if( !isData ) METSmearCorrection(ilead, itrail);
      //3.3 shifting (even in data? but different in data and MC)
      METPhiShiftCorrection();
      //3.4 scaling
      if( isData) METScaleCorrection(ilead, itrail);
      //    }
    _minitree->mtree_corMet    = recoPfMET;
    _minitree->mtree_corMetPhi = recoPfMETPhi;
    
    // 2. dijet tag
    Int_t nVtxJetID = -1;
    if( _config->doPUJetID() ) nVtxJetID = nVtxBS;
    //vector<int> goodJetsIndex = selectJets( ilead, itrail, nVtxJetID, wei, selVtx );
    vector<int> goodJetsIndex = selectJetsJEC(  ilead, itrail, nVtxJetID, wei, selVtx );
    int vbftag(-1),hstratag(-1),catjet(-1);
    dijetSelection( glead, gtrail, goodJetsIndex, wei, selVtx, vbftag, hstratag, catjet);
    _minitree->mtree_vbfTag  = vbftag;
    _minitree->mtree_hvjjTag = hstratag;
    _minitree->mtree_vbfCat  = catjet;      

    // 2x. radion analysis

    // take the very same photon candidates as in the H->GG analysis above
    _minitree->radion_evtNum = event;
    *(_minitree->radion_gamma1) = glead;
    *(_minitree->radion_gamma2) = gtrail;

    vector<int> goodJetsIndexRadion = selectJetsRadion(nVtxJetID, selVtx);
    _minitree->radion_bJetTags->Set(goodJetsIndexRadion.size());

    for (unsigned i = 0; i < goodJetsIndexRadion.size(); i++) {
      int iJet = goodJetsIndexRadion[i];

      //_minitree->radion_bJetTags->AddAt(jetCombinedSecondaryVtxMVABJetTags[iJet], i);
      _minitree->radion_bJetTags->AddAt(jetCombinedSecondaryVtxBJetTags[iJet], i);

      TLorentzVector* jet = new((*(_minitree->radion_jets))[_minitree->radion_jets->GetEntriesFast()]) TLorentzVector();
      jet->SetPtEtaPhiE(jetPt[iJet], jetEta[iJet], jetPhi[iJet], jetEn[iJet]);
    }

    // Continue step 3.
    if (  MetTagSelection(glead,gtrail,goodJetsIndex) && _minitree->mtree_vbfTag == 0 && _minitree->mtree_lepTag == 0   && 
	 _minitree->mtree_corMet > 70.  &&	 
	 fabs( phoSCEta[ilead] ) < 1.45 && fabs( phoSCEta[itrail]) < 1.45 ) _minitree->mtree_metTag = 1;

    //----------- categorisation (for now incl + vbf + lep + met) -----------//
    _minitree->mtree_catBase = CAT4;
    _minitree->mtree_catMva  = -1;
    for( int icat_mva = 0 ; icat_mva < 4; icat_mva++ )
      if( diphotonmva >= Ddipho_cat[icat_mva] ) { _minitree->mtree_catMva = icat_mva; break; }
    if       ( _minitree->mtree_lepTag == 1 ) {
      _minitree->mtree_catBase = _minitree->mtree_lepCat + 6;
      _minitree->mtree_catMva  = _minitree->mtree_lepCat + 6;
    } else if( _minitree->mtree_vbfTag == 1 ) {
      _minitree->mtree_catBase = _minitree->mtree_vbfCat + 4;
      _minitree->mtree_catMva  = _minitree->mtree_vbfCat + 4;
    } else if( _minitree->mtree_metTag == 1 ) {
      _minitree->mtree_catBase = 8;
      _minitree->mtree_catMva  = 8;
    }
    if( diphotonmva < Ddipho_cat[3] )  _minitree->mtree_catMva  = -1;
    
    /// photon pt cut was dropped from the inclusive cuts
    if(  _minitree->mtree_catMva  < 4 && (glead.Pt()/hcand.M() < 1./3. || gtrail.Pt()/hcand.M() < 1./4.) ) 
       _minitree->mtree_catMva  = -1;

    if(  _minitree->mtree_catBase < 4 && (glead.Pt()/hcand.M() < 1./3. || gtrail.Pt()/hcand.M() < 1./4.) ) 
       _minitree->mtree_catBase  = -1;
    
    //------------ MC weighting factors and corrections    
    if( !isData ) {
      /// needs to be recomputed for each event (because reinit var for each event)
      _minitree->mc_wNgen = 100000./_weight_manager->getNevts();
      _minitree->mc_wXsec  = _weight_manager->xSecW(mc_whzh_type);      
      if( ( mode_ == 3 ||  mode_ == 19 ) && isprompt1 == 1 ) _minitree->mc_wXsec  *= 1.3;
      _minitree->mc_wBSz   = _weight_manager->bszW(vtxbs[selVtx][2]-mcVtx[0][2]);
      _minitree->mc_wVtxId = _weight_manager->vtxPtCorrW( hcand.Pt() );
      
      /// photon identification
      float wPhoId[] = { 1., 1.};
      for( int i = 0 ; i < 2; i++ ) {
	wPhoId[i] = 1;
	int index = -1;
	if( i == 0 ) index = ilead;
	if( i == 1 ) index = itrail;
	wPhoId[i] *= _weight_manager->phoIdPresel(phoR9[index],phoSCEta[index]);
	if( _config->analysisType() == "baselineCiC4PF" ) wPhoId[i] *= _weight_manager->phoIdCiC(phoR9[index],phoSCEta[index]);
	if( _config->analysisType() == "MVA"            ) wPhoId[i] *= _weight_manager->phoIdMVA(phoR9[index],phoSCEta[index]);
	if( vetoElec[i] ) wPhoId[i] *= _weight_manager->phoIdEleVeto(phoR9[index],phoSCEta[index]);	      
      }
      _minitree->mc_wPhoEffi   = wPhoId[0]*wPhoId[1];
      
      /// trigger efficiency
      _minitree->mc_wTrigEffi = 0.9968; /// FIX ME

      /// cross section volontary not included in total weight
      _minitree->mc_wei = 
	_minitree->mc_wPU       *
	_minitree->mc_wBSz      *
	_minitree->mc_wHQT      *  /// = 1 but in 2011
	_minitree->mc_wVtxId    *
	_minitree->mc_wPhoEffi  *
	_minitree->mc_wTrigEffi *
	_minitree->mc_wNgen;

      wei =  _minitree->mc_wei;
    }
    nEvts[icutlevel++]++;
    if( _minitree->mtree_lepTag ) nEvts[icutlevel+0]++;
    if( _minitree->mtree_vbfTag ) nEvts[icutlevel+1]++;
    
    //// following crap can be removed when the synchornization will be successfull. 
    if( DoDebugEvent ) {
      _minitree->setSynchVariables();
      _xcheckTextFile << " pho1 pos: " << phoPos[ilead]  << endl;
      _xcheckTextFile << " ieta1 : " << phoSeedDetId1[ilead] << " iphi1: " << phoSeedDetId2[ilead] << endl;
      _xcheckTextFile << " pho2 pos: " << phoPos[itrail] << endl;
      _xcheckTextFile << " ieta1 : " << phoSeedDetId1[itrail] << " iphi1: " << phoSeedDetId2[itrail] << endl;

      _xcheckTextFile << " oversmear 1: " << overSmear.meanOverSmearing(phoSCEta[ilead],phoR9[ilead]  ,isInGAP_EB(ilead),0) << endl;
      _xcheckTextFile << " oversmear 2: " << overSmear.meanOverSmearing(phoSCEta[itrail],phoR9[itrail],isInGAP_EB(itrail),0) << endl;
      _xcheckTextFile << " mass reso (eng only): " << _minitree->mtree_massResoEng << endl;
      _xcheckTextFile << " mass reso (ang only): " << _minitree->mtree_massResoAng << endl;

      for( unsigned i = 0 ; i < _minitree->sync_iName.size(); i++ )
	cout  << _minitree->sync_iName[i] << ":" << *_minitree->sync_iVal[i]  << "  ";
      for( unsigned i = 0 ; i < _minitree->sync_lName.size(); i++ )
	cout << _minitree->sync_lName[i] << ":" << *_minitree->sync_lVal[i]  << "  ";
      for( unsigned i = 0 ; i < _minitree->sync_fName.size(); i++ )
	cout << _minitree->sync_fName[i] << ":" << *_minitree->sync_fVal[i]  << "  ";
      cout << endl;	
      

      cout << "  myVtx = " << selVtx << "; 1=" << sortedVertex[1] << " 2= " << sortedVertex[2] << endl;
      for( int ivtx = 0; ivtx < nVtxBS; ivtx++ ) {
	cout << " etas[ ivtx = " << ivtx << "] = " << phoEtaVtx[ilead][ivtx] << " - " << phoEtaVtx[itrail][ivtx] << "  - zVtx = " << vtxbs[ivtx][2] << endl;

      }
    }

    _minitree->fill();    
    if( _config->doSkimming() && tSkim ) {
      /// undo all modifs before filling the skim
      fChain->GetEntry(jentry);
      tSkim->Fill();
    }
    //---------------      
    
  //   }    
  }// end for entry
  if( _config->doSkimming() ) {
    _minitree->mtree_file->cd();
    TH1F *hEvents = new TH1F("hEvents","hEvents",2,0,2);
    hEvents->SetBinContent(1,_weight_manager->getNevts());
    hEvents->SetBinContent(2,nEvts[nEvts.size()-1]);
    hEvents->Write();
    tSkim->Write();
    _minitree->mtree_file->Close();
  } else  _minitree->end();

  for( unsigned icut = 0 ; icut < nEvts.size(); icut++ )
    cout << nCutName[icut] << nEvts[icut] << endl;
  
  delete rnd;  

  return nEvts[nEvts.size()-1];
}



Int_t xAna::cat4(float r1, float r2, float eta1, float eta2) {

  int cat(-1);
  float rMin(0.), etaMax(10.);

 if (fabs(eta1)>fabs(eta2))
    etaMax = fabs(eta1);
  else
    etaMax = fabs(eta2);

  if (r1<r2)
    rMin = r1;
  else
    rMin = r2;

// cat:  R9_min>0.94, R9_min<0.94, |eta|_max - barrel, |eta|_max - endcap;

  float Rmin(0.94);
  float EtaMax(1.479);

// Make cat:

  if(etaMax < EtaMax) {
    if (rMin > Rmin)
      cat = 0;
    else
      cat = 1;
  } else {
    if (rMin > Rmin)
      cat = 2;
    else
      cat = 3;
  }
  
  return cat;
  
}


Float_t xAna::etaTransformation(  float EtaParticle , float vertexZ)  {


  //---Definitions for ECAL
  const float R_ECAL           = 136.5;
  const float Z_Endcap         = 328.0;
  const float etaBarrelEndcap  = 1.479; 
   
  //---ETA correction

  float Theta = 0.0  ; 
  float ZEcal = R_ECAL*sinh(EtaParticle)+vertexZ;

  if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
  if(Theta<0.0) Theta = Theta+ TMath::Pi() ;
  double ETA = - log(tan(0.5*Theta));
         
  if( fabs(ETA) > etaBarrelEndcap )
    {
      float Zend = Z_Endcap ;
      if(EtaParticle<0.0 )  Zend = -Zend ;
      float Zlen = Zend - vertexZ ;
      float RR = Zlen/sinh(EtaParticle); 
      Theta = atan(RR/Zend);
      if(Theta<0.0) Theta = Theta+ TMath::Pi() ;
      ETA = - log(tan(0.5*Theta));		      
    } 
  //---Return the result
  return ETA;
  //---end
}




bool xAna::preselectPhoton(int i, float Et){

  bool rc = true;

  if (fabs(phoSCEta[i]) > 1.4442 && fabs(phoSCEta[i]) < 1.566) return false; 
  if (fabs(phoSCEta[i]) > 2.5) return false; 

  if (Et < 20) return false; 

  return rc;

}

TVector3 xAna::getCorPhotonTVector3(int i, int selVtx) {
  float vtx_X, vtx_Y, vtx_Z;
  vtx_X = vtxbs[selVtx][0];
  vtx_Y = vtxbs[selVtx][1];
  vtx_Z = vtxbs[selVtx][2];
  TVector3 vtxVec(vtx_X,vtx_Y,vtx_Z);

  // float pho1_X, pho1_Y, pho1_Z;
  // pho1_X = phoCaloPos[i][0];
  // pho1_Y = phoCaloPos[i][1];
  // pho1_Z = phoCaloPos[i][2];
  // TVector3 phoVec(pho1_X,pho1_Y,pho1_Z);
  // // Temporary SC pos for old ggNtuple
  // phoVec.SetPtEtaPhi(sqrt(pow(pho1_X,2)+pow(pho1_Y,2)), phoSCEta[i], phoSCPhi[i]);

  // For new version ggNtuple with SC pos
  float pho1_X = phoSCPos[i][0];
  float pho1_Y = phoSCPos[i][1];
  float pho1_Z = phoSCPos[i][2];
  TVector3 phoVec(pho1_X,pho1_Y,pho1_Z);


  return phoVec - vtxVec;
}



