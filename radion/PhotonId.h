#ifndef __photonid_h__
#define __photonid_h__

Int_t xAna::phocatCIC(int i){
  Int_t r9 = (phoR9[i] > 0.9) ? 0 : 1;
  Int_t eta = (fabs(phoSCEta[i]) < 1.479) ? 0 : 1;

  return r9 + 2*eta;
}

Int_t xAna::phocat(int i) {

  Int_t r9 = (phoR9[i] > 0.94) ? 0 : 1;
  Int_t eta = (fabs(phoSCEta[i]) < 1.479) ? 0 : 1;

  return r9 + 2*eta;
}

Int_t xAna::photonID(int ipho, int ivtx, float corEt_ipho, bool invEleVeto ) {
  
  //  if(invEleVeto){
  //}
  
  //  cout <<"ipho  " << ipho   << "    invEleVeto  "<<  invEleVeto <<endl;
  if( _config->analysisType() == "MVA"      ) return mvaPhotonID(ipho,ivtx,corEt_ipho,invEleVeto);
  if( _config->analysisType() == "baseline" ) return cicPhotonID(ipho,ivtx,corEt_ipho);
  if( _config->analysisType() == "baselineCiC4PF" ) {
    if( !mvaPhotonID(ipho,ivtx,corEt_ipho,invEleVeto     ) ) return 0;
    return cic4PFPhotonID(ipho,ivtx,corEt_ipho,invEleVeto);  
  }
  
  return 0;
}


Int_t xAna::cic4PFPhotonID(int ipho, int ivtx, float corEt_ipho, bool invEleVeto){
    
  Int_t passID = 1;
  Float_t Et = phoEt[ipho]; 
  if(corEt_ipho != -1) {
    Et = corEt_ipho;
  }
 
   
  Et = phoE[ipho] / cosh(phoEtaVtx[ipho][ivtx]);
  float effectiveArea         = 0.09;
  float effectiveAreaWorstVtx = 0.23;
  float pfiso_charged_chosen_vertex = phoCiCPF4chgpfIso03[ipho][ivtx];
  float pfiso_photon03 = phoCiCPF4phopfIso03[ipho];
  float pfiso_photon04 = phoCiCPF4phopfIso04[ipho];

  Int_t badVtx = 0;
  float  pfiso_charged_worst_vertex = -99;
  for( int iv=0; iv<nVtxBS ; iv++) {
    if(phoCiCPF4chgpfIso04[ipho][iv] > pfiso_charged_worst_vertex) {
      pfiso_charged_worst_vertex = phoCiCPF4chgpfIso04[ipho][iv];
      badVtx = iv;
    }
  }
  
  /// replica of the bug in gglobe
  //  float pho_tkiso_badvtx = worstSumTrackPtInCone(ipho,badVtx,0,0.4,0.02,0.0,1.0,0.1);
  
  float pfiso_charged_over_et        = (pfiso_charged_chosen_vertex)*50./Et;  
  float isosum_chosen_vertex_over_et = (pfiso_charged_chosen_vertex + pfiso_photon03 + 2.5 - rho2012*effectiveArea)*50./Et;  
  float isosum_worst_vertex_over_et = 
    (pfiso_charged_worst_vertex + pfiso_photon04 + 2.5 - rho2012*effectiveAreaWorstVtx)*50./( phoE[ipho]/cosh(phoEtaVtx[ipho][badVtx]) );


  Int_t cat = phocat(ipho);
  if( isosum_chosen_vertex_over_et > vcicST[0][cat] ) return 0;
  if( isosum_worst_vertex_over_et  > vcicST[1][cat] ) return 0;
  if( pfiso_charged_over_et  > vcicST[2][cat] ) return 0;
  if( phoSigmaIEtaIEta[ipho] > vcicST[3][cat] ) return 0;
  if( phoHoverE[ipho]        > vcicST[4][cat] ) return 0;
  if( phoR9[ipho]            < vcicST[5][cat] ) return 0;

  if (phoCiCdRtoTrk[ipho] == 0) phoCiCdRtoTrk[ipho] = 99;
  //  hDRtoTrk->Fill(phoCiCdRtoTrk[ipho],wei);
  
  if( _config->noElectronVeto() ) return passID;

  int isElec = phoEleVeto[ipho];
  if(  invEleVeto && isElec == 0 ) return 0;
  if( !invEleVeto && isElec == 1 ) return 0;

  // if(_config->isOnlyLepTag()){
  //   if( !_config->invertElectronVeto() && !invEleVeto_){
  //     //  if( phoEleVeto[ipho]==1 || phoCiCdRtoTrk[ipho]<1  ||  phohasPixelSeed[ipho]==1 )return 0;
  //    if( phoEleVeto[ipho]==1)return 0;
  //   }
  //   if(invEleVeto_ ||  _config->invertElectronVeto()){
  //     // if( phoEleVeto[ipho]==0  &&  phoCiCdRtoTrk[ipho]>1  && phohasPixelSeed[ipho]==0) return 0;
  //     if( phoEleVeto[ipho]==0  ) return 0;
  //   }  
  //}
  return passID;  
}

void xAna::fillPhotonVariablesToMiniTree(int ipho, int ivtx, int index ){ 

  float Et = phoE[ipho]/cosh(phoEtaVtx[ipho][ivtx]);

  _minitree->mtree_r9[index]    = phoR9[ipho];
  _minitree->mtree_sceta[index] = phoSCEta[ipho];
  _minitree->mtree_eg[index]    = phoE[ipho];
  _minitree->mtree_ptg[index]   = phoE[ipho]/cosh(phoEtaVtx[ipho][ivtx]);

  _minitree->mtree_eta[index]   = phoEtaVtx[ipho][ivtx];
  _minitree->mtree_phi[index]   = phoPhiVtx[ipho][ivtx];
  /* /// debugging */
  /* _minitree->mtree_eta[index]   = phoSeedDetId1[ipho]; */
  /* _minitree->mtree_phi[index]   = phoSeedDetId2[ipho];; */
  
  _minitree->mtree_relResOverE[index] = phoRegrErr[ipho]/phoE[ipho]; 
  _minitree->mtree_relSmearing[index] = phoRegrSmear[ipho]/phoE[ipho]; 

  _minitree->mtree_coviEiE[index] = phoSigmaIEtaIEta[ipho];
  _minitree->mtree_HoE[index]     = phoHoverE[ipho];
  _minitree->mtree_drTrk[index]   = phoCiCdRtoTrk[ipho];
  
  _minitree->mtree_cat[index]     = phocatCIC(ipho); // (r9 < 0.90 / r9>0.90 )

  /// preselection photon id cuts
  _minitree->mtree_isElec[index]            = phoEleVeto[ipho];
  _minitree->mtree_etCorrEcalIso[index]     = phoEcalIsoDR03[ipho] - 0.012*Et;
  _minitree->mtree_etCorrHcalIso[index]     = phoHcalIsoDR03[ipho] - 0.005*Et;
  _minitree->mtree_etCorrTrkIso[index]      = phoTrkIsoHollowDR03[ipho] - 0.002*Et;
  _minitree->mtree_pfChgIso02_selVtx[index] = phoCiCPF4chgpfIso02[ipho][ivtx];

  /// mva photon id specific 2012
  _minitree->mtree_coviEiP[index]      = phoSigmaIEtaIPhi[ipho];
  _minitree->mtree_s4Ratio[index]      = phoS4ratio[ipho];
  _minitree->mtree_esEffSigmaRR[index] = phoESEffSigmaRR[ipho][0];

  //// PF isolation dedicated to CiCpf id (also used in photon Id MVA )
  /// compute isolation
  float effectiveArea         = 0.09;
  float effectiveAreaWorstVtx = 0.23;
  float pfiso_charged_chosen_vertex = phoCiCPF4chgpfIso03[ipho][ivtx];
  float pfiso_photon03 = phoCiCPF4phopfIso03[ipho];
  float pfiso_photon04 = phoCiCPF4phopfIso04[ipho]; 
  
  Int_t badVtx03   = 0;
  Int_t badVtx04   = 0;
  float  pfiso_charged03_worst_vertex = -99;
  float  pfiso_charged04_worst_vertex = -99;
  for( int iv=0; iv<nVtxBS ; iv++) {
    // worst vtx for chIso03
    if(phoCiCPF4chgpfIso03[ipho][iv] > pfiso_charged03_worst_vertex) {
      pfiso_charged03_worst_vertex = phoCiCPF4chgpfIso03[ipho][iv];
      badVtx03 = iv;
    }
   
    // worst vtx for chIso04
    if(phoCiCPF4chgpfIso04[ipho][iv] > pfiso_charged04_worst_vertex) {
      pfiso_charged04_worst_vertex = phoCiCPF4chgpfIso04[ipho][iv];
      badVtx04 = iv;
    }
  } 

  float pfiso_charged_over_et        = (pfiso_charged_chosen_vertex)*50./phoE[ipho]/cosh(phoEtaVtx[ipho][ivtx]);
  float isosum_chosen_vertex_over_et = (pfiso_charged_chosen_vertex  + pfiso_photon03 + 2.5 
					- rho2012*effectiveArea)*50./Et;
  float isosum_worst_vertex_over_et  = (pfiso_charged04_worst_vertex + pfiso_photon04 + 2.5 
					- rho2012*effectiveAreaWorstVtx)*50./phoE[ipho]/cosh(phoEtaVtx[ipho][badVtx04]);
  
  _minitree->mtree_phoIso1[index] = isosum_chosen_vertex_over_et; // CombinedPFIso1
  _minitree->mtree_phoIso2[index] = isosum_worst_vertex_over_et;  // CombinedPFIso2
  _minitree->mtree_phoIso3[index] = pfiso_charged_over_et;        // CombinedPFIso3

  /// particle flow photon isolation
  _minitree->mtree_pfPhoIso03[index] = phoCiCPF4phopfIso03[ipho]; /// mvaPhotonID
  _minitree->mtree_pfPhoIso04[index] = phoCiCPF4phopfIso04[ipho]; 
  _minitree->mtree_pfChgIso03_selVtx[index] = phoCiCPF4chgpfIso03[ipho][ivtx];     /// mvaPhotonID
  _minitree->mtree_pfChgIso03_badVtx[index] = phoCiCPF4chgpfIso03[ipho][badVtx03]; /// mvaPhotonID
  _minitree->mtree_pfChgIso03_defVtx[index] = phoCiCPF4chgpfIso03[ipho][0];

  _minitree->mtree_pfChgIso04_selVtx[index] = phoCiCPF4chgpfIso04[ipho][ivtx];
  _minitree->mtree_pfChgIso04_badVtx[index] = phoCiCPF4chgpfIso04[ipho][badVtx04];
  _minitree->mtree_pfChgIso04_defVtx[index] = phoCiCPF4chgpfIso04[ipho][0];

  PhoID_MVA(ipho,ivtx);
  _minitree->mtree_mvaid[index] = -1;
  if( fabs(phoSCEta[ipho]) < 1.45 ) _minitree->mtree_mvaid[index] = phoID_mva[0]->EvaluateMVA("BDT");
  else                              _minitree->mtree_mvaid[index] = phoID_mva[1]->EvaluateMVA("BDT");

  _minitree->mtree_mvaid1 = _minitree->mtree_mvaid[0];
  _minitree->mtree_mvaid2 = _minitree->mtree_mvaid[1];
  _minitree->mtree_eta1   = _minitree->mtree_sceta[0];
  _minitree->mtree_eta2   = _minitree->mtree_sceta[1];
  _minitree->mtree_r91    = _minitree->mtree_r9[0];
  _minitree->mtree_r92    = _minitree->mtree_r9[1];

  ///------ synch only variables.
  if( _minitree->addSyncVariables ) {
    int iconv = matchPhotonToConversion(ipho);
    _minitree->convInd[index] = iconv;
    _minitree->convNtk[index] = iconv>=0 ? convNTracks[iconv] : -1;
    _minitree->convValidVtx[index] = iconv>=0 && convNTracks[iconv] == 2 ? convValidVtx[iconv] : -1;
    _minitree->convChi2Prob[index] = iconv>=0 && convNTracks[iconv] == 2 ? convChi2Probability[iconv] : -1;
    if( iconv >= 0 ) {
      TVector3 qconv;
      if     ( convNTracks[iconv] == 2 ) qconv.SetXYZ( convRefittedMomentum[iconv][0],
						       convRefittedMomentum[iconv][1],
						       convRefittedMomentum[iconv][2] );
      else if( convNTracks[iconv] == 1 ) qconv.SetXYZ( convPairMomentum[iconv][0],
						       convPairMomentum[iconv][1],
						       convPairMomentum[iconv][2] ); 
      _minitree->convPt[index] = qconv.Pt();	
    }
  }
  _minitree->phoID_covIEtaIPhi[index]    =  phoSigmaIEtaIPhi[ipho];
  _minitree->phoID_S4[index]             = phoS4ratio[ipho];
  _minitree->phoID_ESEffSigmaRR[index]   = phoESEffSigmaRR[ipho][0];
  _minitree->phoID_pfPhoIso03[index]     = phoCiCPF4phopfIso03[ipho];
  _minitree->phoID_pfChIso03[index]      = phoCiCPF4chgpfIso03[ipho][ivtx];
  _minitree->phoID_pfChIso03worst[index] = pfiso_charged03_worst_vertex;
  _minitree->phoID_E2x2[index]           = phoE2x2[ipho];
  _minitree->phoID_E5x5[index]           = phoE5x5[ipho];
  
}


Int_t xAna::cicPhotonID(int ipho, int ivtx, float corEt_ipho ) {
  Int_t passID = 1;
  Float_t Et = phoEt[ipho]; 
  if(corEt_ipho != -1) {
    Et = corEt_ipho;
  }
  
  Float_t rhoAbsFactor = 1;
  if(DoSetRho0_) rhoAbsFactor = 0;

  Float_t rhofracbad = 0.52*rhoAbsFactor, rhofrac = 0.17*rhoAbsFactor, combIso = 0*rhoAbsFactor, trkIso = 0*rhoAbsFactor;
  Int_t cat = phocat(ipho);
  
  /*
  // find out the tracker isolation with "worst vertex"
  Float_t maxTrkIso = -99;
  for (int i=0; i<nVtxBS; ++i) {
  if (phoCiCTrkIsoDR04[ipho][i] >  maxTrkIso) {
  maxTrkIso = phoCiCTrkIsoDR04[ipho][i];
  }
  }
  
  // relative combIso (bad vertex)
  combIso = (maxTrkIso + phoEcalIsoDR04[ipho] + phoHcalIsoDR04[ipho] - rho25*rhofracbad) * 50./Et;
  */
  
  // relative combIso (good vertex)
  combIso = (phoCiCTrkIsoDR03[ipho][ivtx] + phoEcalIsoDR03[ipho] + phoHcalIsoDR04[ipho] - rho25*rhofrac) * 50./Et; 
  if (combIso > vcicSTmva[0][cat]) return 0;
  
  
  // find out the tracker isolation with "worst vertex"
  Float_t maxTrkIso = -99;
  Int_t badVtx = 0;
  for (int i=0; i<nVtxBS; ++i) {
    if (phoCiCTrkIsoDR04[ipho][i] >  maxTrkIso) {
      maxTrkIso = phoCiCTrkIsoDR04[ipho][i];
      badVtx = i;
    }
  }
  
  // relative combIso (bad vertex)
  combIso = (maxTrkIso + phoEcalIsoDR04[ipho] + phoHcalIsoDR04[ipho] - rho25*rhofracbad) * 50./(phoE[ipho]/cosh(phoEtaVtx[ipho][badVtx]));
  //combIso = (maxTrkIso + phoEcalIsoDR04[ipho] + phoHcalIsoDR04[ipho] - rho25*rhofracbad) * 50./phoEtVtx[ipho][badVtx];
  //combIso = (maxTrkIso + phoEcalIsoDR04[ipho] + phoHcalIsoDR04[ipho] - rho25*rhofracbad) * 50./Et;
  if (combIso > vcicSTmva[1][cat]) return 0;
  
  // relative trkIso (good vertex)
  trkIso = phoCiCTrkIsoDR03[ipho][ivtx] * 50./Et;
  if (trkIso > vcicSTmva[2][cat]) return 0;
  
  // Sigma_IetaIeta
  if (phoSigmaIEtaIEta[ipho] > vcicSTmva[14][cat]) return 0;
  // H/E
  if (phoHoverE[ipho] > vcicSTmva[15][cat]) return 0;
  // R9
  if (phoR9[ipho] < vcicSTmva[5][cat])  return 0;
  // dR to electron track
  if (phoCiCdRtoTrk[ipho] == 0) phoCiCdRtoTrk[ipho] = 99;
  
  if(  _config->invertElectronVeto() && phoCiCdRtoTrk[ipho] > vcicSTmva[6][cat] ) return 0;
  if( !_config->invertElectronVeto() && phoCiCdRtoTrk[ipho] < vcicSTmva[6][cat] ) return 0;
  
  return passID;
}

Int_t xAna::mvaPhotonID(int ipho, int ivtx, float corEt_ipho, bool invEleVeto ) {


  
  // cout <<" preselection     ipho  "<< ipho <<"  invEleVeto  "<<   invEleVeto <<endl;
  Int_t passID = 1;
  Int_t cat = phocatCIC(ipho);

  Float_t Et = phoEt[ipho]; 
  if(corEt_ipho != -1) {
    Et = corEt_ipho;  
  }

  Et = phoE[ipho] / cosh(phoEtaVtx[ipho][ivtx]);
  float EtCorrEcalIso = phoEcalIsoDR03[ipho]      - 0.012*Et;
  float EtCorrHcalIso = phoHcalIsoDR03[ipho]      - 0.005*Et;
  float EtCorrTrkIso  = phoTrkIsoHollowDR03[ipho] - 0.002*Et;
  float  trkPfIso02   = phoCiCPF4chgpfIso02[ipho][ivtx];
  
  int isElec =  phoEleVeto[ipho];//isMatchedToAPromptElec(ipho, dr );

  if( phoSigmaIEtaIEta[ipho]  >= vcicSTmva[ 3][cat] ) return 0;
  if( phoHoverE[ipho]         >= vcicSTmva[ 4][cat] ) return 0;
  //if( EtCorrEcalIso           >= vcicSTmva[ 8][cat] ) return 0;
  if( EtCorrHcalIso           >= vcicSTmva[ 9][cat] ) return 0;
  if( EtCorrTrkIso            >= vcicSTmva[10][cat] ) return 0;


  if(run==191247  && event== 550253379){
    cout << " invEleVeto:  " <<  invEleVeto << " phoEleVeto[ipho]  "   <<  phoEleVeto[ipho] 
	 << " isElec "       << isElec      << " phoCiCdRtoTrk[ipho] " <<  phoCiCdRtoTrk[ipho] <<endl;   
  }


  /// 2012 preselection only
  if( trkPfIso02            >=  vcicSTmva[16][cat] ) return 0;

  /// 2011 preselection only
  //  float PUCorrHcalEcal       = EtCorrEcalIso + EtCorrHcalIso - 0.17*rho25;
  //  float HollowConeTrkIsoDr03 = phoTrkIsoHollowDR03[ipho];
  /// the 2 are not exactly identical though very very close
  /// the second gives slightly closer results to gglobe
  // float AbsTrkIsoCIC         = sumTrackPtInCone(ipho,ivtx,0,0.3,0.02,0.0,1.0,0.1);
  //  float AbsTrkIsoCIC         = phoCiCTrkIsoDR03[ipho][ivtx];
  //  if( PUCorrHcalEcal         >= vcicSTmva[11][cat] ) return 0;
  //  if( AbsTrkIsoCIC           >= vcicSTmva[12][cat] ) return 0;
  //  if( HollowConeTrkIsoDr03   >= vcicSTmva[13][cat] ) return 0;

  if (phoCiCdRtoTrk[ipho] == 0) phoCiCdRtoTrk[ipho] = 99;
  // if( phoCiCdRtoTrk[ipho] > 0.15 ) phoCiCdRtoTrk[ipho] = 0.15;

  // float dr = -1;
  //  hDRtoTrk->Fill(phoCiCdRtoTrk[ipho],wei);

  if( _config->noElectronVeto() ) return passID;
    
  if(  invEleVeto && isElec == 0 ) return 0;
  if( !invEleVeto && isElec == 1 ) return 0;

  return passID;
}

#endif
