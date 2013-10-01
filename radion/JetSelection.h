#ifndef JetSelection_h__
#define JetSelection_h__

#include "xAna_allAna.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

vector<int> xAna::selectJets(Int_t iPho1, Int_t iPho2, Int_t nvtx, Double_t wei, Int_t selVtx) {
  vector<int> jetIndex; 
  for( int i = 0 ; i < nJet; ++i ) { 
    //  if(i>20)continue;
    hAllJetPt->Fill(  jetPt[i],wei);
    hAllJetEta->Fill(jetEta[i],wei);
    
    float dR_jetg1 = deltaR(phoEta[iPho1], phoPhi[iPho1],jetEta[i],jetPhi[i]);
    float dR_jetg2 = deltaR(phoEta[iPho2], phoPhi[iPho2],jetEta[i],jetPhi[i]);
    //cout <<   " dR_jetg1  "<<  dR_jetg1  <<" dR_jetg2  "<<  dR_jetg2   <<endl;
    
    hDRAllJetg1->Fill(dR_jetg1,wei);
    hDRAllJetg2->Fill(dR_jetg2,wei);
    
    if( jetPt[i] > 25 ){
      
      if( dR_jetg1 < 0.5 || dR_jetg2 < 0.5 ) continue;    
      
      hJetNHF->Fill(jetNHF[i],wei);
      hJetNEF->Fill(jetNEF[i],wei);
      hJetNConst->Fill(jetNConstituents[i],wei);
      
      hJetCHF->Fill(jetCHF[i],wei);
      hJetNCH->Fill(jetNCH[i],wei);
      hJetCEF->Fill(jetCEF[i],wei);

      // 2011
      if( nvtx <= 0 ) {     
        if ( fabs(jetEta[i]) > 4.7 ) continue;
      }
      // PU Jet ID 2012
      else if( nvtx > 0 ) {
        if ( fabs(jetEta[i]) > 4.7 ) continue;
        // if ( (jetWPLevels[i][2] & 4) == 0 ) continue; // [0] simple MVA, [1] full MVA, [2] cutbased, [3] PhilV1; 1(= 1<<0): tight, 2(= 1<<1): medium, 4(= 2<<2): loose /// !! for ntuples after May30
        if (      fabs(jetEta[i]) < 2.5  ) { if (jetBetaStarClassicExt[i][selVtx]/log(nvtx-0.64) > 0.2 || jetDR2Mean[i] > 0.06 ) continue; }
        else if ( fabs(jetEta[i]) < 2.75 ) { if (jetBetaStarClassicExt[i][selVtx]/log(nvtx-0.64) > 0.3 || jetDR2Mean[i] > 0.05 ) continue; }
        else if ( fabs(jetEta[i]) < 3.0  ) { if (                                                         jetDR2Mean[i] > 0.05 ) continue; }
        else                               { if (                                                         jetDR2Mean[i] > 0.055) continue; }
      }

      jetIndex.push_back( i );
      
    }	     
  }   
  return jetIndex;  
}

vector<int> xAna::selectJetsJEC(Int_t iPho1, Int_t iPho2, Int_t nvtx, Double_t wei, Int_t selVtx) {
  vector<int> jetIndex;
  vector<int> jetSortedIndex;

  int jecMode = _weight_manager->getCrossSection()->mode();
  /// ------------------------------- JEC --------------------------------------------------//
  /// step 1
  // Create the JetCorrectorParameter objects, the order does not matter.
  // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
  JetCorrectorParameters *ResJetPar = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_MC_L2L3Residual_AK5PF.txt");
  JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_MC_L3Absolute_AK5PF.txt");
  JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_MC_L2Relative_AK5PF.txt");
  JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_MC_L1FastJet_AK5PF.txt");

  if (jecMode<0) { // ======== Offline DATA JEC Setup ======== //
    ResJetPar = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_DATA_L2L3Residual_AK5PF.txt");
    L3JetPar  = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_DATA_L3Absolute_AK5PF.txt");
    L2JetPar  = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_DATA_L2Relative_AK5PF.txt");
    L1JetPar  = new JetCorrectorParameters("./jecAK5PF/Summer12_V7_DATA_L1FastJet_AK5PF.txt");
  }

  /// step 2
  //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
  vector<JetCorrectorParameters> vPar;
  vPar.push_back(*L1JetPar);
  vPar.push_back(*L2JetPar);
  vPar.push_back(*L3JetPar);
  if (jecMode<0) vPar.push_back(*ResJetPar);

  /// step 3
  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(vPar);

  for( int i = 0 ; i < nJet; ++i ) {
    /// step 4
    JetCorrector->setJetEta(jetEta[i]);
    JetCorrector->setJetPt(jetRawPt[i]);
    JetCorrector->setJetA(jetArea[i]);
    JetCorrector->setRho(rho2012); 
    /// step 5
    double correction = JetCorrector->getCorrection();
    // vector<double> factors = JetCorrector->getSubCorrections();

    // cout << "   " << jetPt[i] << "   " << jetRawPt[i] << "   " << jetPt[i]/jetRawPt[i] << "   " << correction << endl;
    // cout << "   " << jetEn[i] << "   " << jetRawEn[i] << "   " << jetEn[i]/jetRawEn[i] << "   " << correction << endl;

    /// step 6
    //     correction *= 1 + jetJECUnc[i];
    jetPt[i] = jetRawPt[i] * correction;
    jetEn[i] = jetRawEn[i] * correction;
    /// --------------------------------------------------------------------------------------//

    //Do the regression if it is in place.
    jetPt[i] = correctJetPt(i);
    jetEn[i] *= jetPt[i]/(jetRawPt[i]*correction);

    //  if(i>20)continue;
    hAllJetPt->Fill(  jetPt[i],wei);
    hAllJetEta->Fill(jetEta[i],wei);

    float dR_jetg1 = deltaR(phoEta[iPho1], phoPhi[iPho1],jetEta[i],jetPhi[i]);
    float dR_jetg2 = deltaR(phoEta[iPho2], phoPhi[iPho2],jetEta[i],jetPhi[i]);
    //cout <<   " dR_jetg1  "<<  dR_jetg1  <<" dR_jetg2  "<<  dR_jetg2   <<endl;

    hDRAllJetg1->Fill(dR_jetg1,wei);
    hDRAllJetg2->Fill(dR_jetg2,wei);

    if( jetPt[i] > 25 ){

      if( dR_jetg1 < 0.5 || dR_jetg2 < 0.5 ) continue;

      hJetNHF->Fill(jetNHF[i],wei);
      hJetNEF->Fill(jetNEF[i],wei);
      hJetNConst->Fill(jetNConstituents[i],wei);

      hJetCHF->Fill(jetCHF[i],wei);
      hJetNCH->Fill(jetNCH[i],wei);
      hJetCEF->Fill(jetCEF[i],wei);

      // 2011
      if( nvtx <= 0 ) {
        if ( fabs(jetEta[i]) > 4.7 ) continue;
      }
      // PU Jet ID 2012
      else if( nvtx > 0 ) {
        if ( fabs(jetEta[i]) > 4.7 ) continue;
        // if ( (jetWPLevels[i][2] & 4) == 0 ) continue; // [0] simple MVA, [1] full MVA, [2] cutbased, [3] PhilV1; 1(= 1<<0): tight, 2(= 1<<1): medium, 4(= 2<<2): loose /// !! for ntuples after May30
        if (      fabs(jetEta[i]) < 2.5  ) { if (jetBetaStarClassicExt[i][selVtx]/log(nvtx-0.64) > 0.2 || jetDR2Mean[i] > 0.06 ) continue; }
        else if ( fabs(jetEta[i]) < 2.75 ) { if (jetBetaStarClassicExt[i][selVtx]/log(nvtx-0.64) > 0.3 || jetDR2Mean[i] > 0.05 ) continue; }
        else if ( fabs(jetEta[i]) < 3.0  ) { if (                                                         jetDR2Mean[i] > 0.05 ) continue; }
        else                               { if (                                                         jetDR2Mean[i] > 0.055) continue; }
      }

      jetIndex.push_back( i );

    }
  }

  delete ResJetPar;
  delete L3JetPar;
  delete L2JetPar;
  delete L1JetPar;
  delete JetCorrector;
 
  int leadJetId    = -99;
  int subleadJetId = -99;
  Float_t leadJetPt    = -99.;
  Float_t subleadJetPt = -99.;
  for ( unsigned ijet =0; ijet < jetIndex.size(); ijet++) {
    if      (jetPt[jetIndex[ijet]] > leadJetPt)    { subleadJetPt = leadJetPt;             subleadJetId=leadJetId; 
                                                     leadJetPt    = jetPt[jetIndex[ijet]]; leadJetId=ijet; }
    else if (jetPt[jetIndex[ijet]] > subleadJetPt) { subleadJetPt = jetPt[jetIndex[ijet]]; subleadJetId=ijet;}
  }

  if (leadJetId    > -1) jetSortedIndex.push_back(jetIndex[leadJetId]);
  if (subleadJetId > -1) jetSortedIndex.push_back(jetIndex[subleadJetId]);

  return jetSortedIndex;
}

vector<int> xAna::selectJetsRadion(Int_t nvtx, Int_t selVtx) {
  // The code is borrowed from xAna::selectJetsJEC() from above.
  // NOTE: jetPt and jetEn were already recalculated in xAna::selectJetsJEC().

  vector<int> jetIndex;
  vector<int> jetSortedIndex;

  // jet loop
  for(int i = 0; i < nJet; i++) {
    if (jetPt[i] < 25) continue;

    // 2011
    if( nvtx <= 0 ) {
      if ( fabs(jetEta[i]) > 4.7 ) continue;
    }
    // PU Jet ID 2012
    else if ( nvtx > 0 ) {
      if ( fabs(jetEta[i]) > 4.7 ) continue;
      // if ( (jetWPLevels[i][2] & 4) == 0 ) continue; // [0] simple MVA, [1] full MVA, [2] cutbased, [3] PhilV1; 1(= 1<<0): tight, 2(= 1<<1): medium, 4(= 2<<2): loose /// !! for ntuples after May30
      if (      fabs(jetEta[i]) < 2.5  ) { if (jetBetaStarClassicExt[i][selVtx]/log(nvtx-0.64) > 0.2 || jetDR2Mean[i] > 0.06 ) continue; }
      else if ( fabs(jetEta[i]) < 2.75 ) { if (jetBetaStarClassicExt[i][selVtx]/log(nvtx-0.64) > 0.3 || jetDR2Mean[i] > 0.05 ) continue; }
      else if ( fabs(jetEta[i]) < 3.0  ) { if (                                                         jetDR2Mean[i] > 0.05 ) continue; }
      else                               { if (                                                         jetDR2Mean[i] > 0.055) continue; }
    }

    jetIndex.push_back(i);
  }

  // sort by jetPt: slow but simple algorithm
  while (jetIndex.size() > 0) {
    int leadJetId = -99;
    Float_t leadJetPt = -99.;

    // find jet with highest pt
    for (unsigned ijet = 0; ijet < jetIndex.size(); ijet++)
      if (jetPt[jetIndex[ijet]] > leadJetPt) {
        leadJetPt = jetPt[jetIndex[ijet]];
        leadJetId = ijet;
      }

    // protection against negative jetPt and thus infinite loop
    if (leadJetId < 0) {
      std::cerr << "selectJetsRadion[ERROR]: jet with negative jetPt" << std::endl;
      break;
    }

    jetSortedIndex.push_back(jetIndex[leadJetId]);
    jetIndex.erase(jetIndex.begin() + leadJetId); // inefficient for vectors, use forward_list?
  } // jet sorting loop

  return jetSortedIndex;
}

void xAna::dijetSelection( const TLorentzVector & corPho1, const TLorentzVector & corPho2,
					 const vector<int> & goodJetsIndex,
					 const Double_t &wei, int ivtx,
					 int &vbftag, int &hstratag, int &catjet ) {
  vbftag   = 0;
  hstratag = 0;
  catjet   =-1;
 
  if( goodJetsIndex.size() < 2 ) return;
  int iJet1 = goodJetsIndex[0];
  int iJet2 = goodJetsIndex[1];

  TLorentzVector l_corr_gg;
  l_corr_gg = corPho1 + corPho2 ;

  TLorentzVector jet1;
  TLorentzVector jet2;
  TLorentzVector ljetjet;

  jet1.SetPtEtaPhiE(jetPt[iJet1], jetEta[iJet1], jetPhi[iJet1], jetEn[iJet1]);
  jet2.SetPtEtaPhiE(jetPt[iJet2], jetEta[iJet2], jetPhi[iJet2], jetEn[iJet2]);
  ljetjet = jet1 + jet2;

  float  deltaEta_jj = jetEta[iJet1]- jetEta[iJet2];
  float  massJetJet = ljetjet.M();
  // hggMass_twoJetTag->Fill(selectedMass,wei);
  // hggPt_twoJetTag->Fill(pt,wei);

  /* float dR_jet1g1 = deltaR(phoSCEta[iPho1], phoSCPhi[iPho1],jet1.Eta(),jet1.Phi()); */
  /* float dR_jet1g2 = deltaR(phoSCEta[iPho2], phoSCPhi[iPho2],jet1.Eta(),jet1.Phi()); */
  /* float dR_jet2g1 = deltaR(phoSCEta[iPho1], phoSCPhi[iPho1],jet2.Eta(),jet2.Phi()); */
  /* float dR_jet2g2 = deltaR(phoSCEta[iPho2], phoSCPhi[iPho2],jet2.Eta(),jet2.Phi()); */
  float dR_jet1g1 = deltaR(corPho1.Eta(), corPho1.Phi(),jet1.Eta(),jet1.Phi());
  float dR_jet1g2 = deltaR(corPho1.Eta(), corPho1.Phi(),jet1.Eta(),jet1.Phi());
  float dR_jet2g1 = deltaR(corPho1.Eta(), corPho1.Phi(),jet2.Eta(),jet2.Phi());
  float dR_jet2g2 = deltaR(corPho1.Eta(), corPho1.Phi(),jet2.Eta(),jet2.Phi());

  /// outside CMSSW do not have ROOT::Math::VectorUtil
  float dPhi_jet1jet2 = jet1.DeltaPhi(jet2);

  hDRjet1g1->Fill(dR_jet1g1);
  hDRjet1g2->Fill(dR_jet1g2);
  hDRjet2g1->Fill(dR_jet2g1);
  hDRjet2g2->Fill(dR_jet2g2);
  hDPhijet1jet2->Fill(dPhi_jet1jet2);

  hJet1Pt->Fill(jetPt[iJet1],wei);
  hJet2Pt->Fill(jetPt[iJet2],wei);
  hJet1Eta->Fill(jetEta[iJet1],wei);
  hJet2Eta->Fill(jetEta[iJet2],wei);

  hDeltaEtaJets->Fill(fabs(deltaEta_jj),wei);

  double zeppenfeld = l_corr_gg.Eta() - (jetEta[iJet1]+ jetEta[iJet2])/2.;

  hZeppenfeld->Fill(zeppenfeld,wei);

  double dPhiDiJetDiPho = ljetjet.DeltaPhi(l_corr_gg);
  hDeltaPhiDiJetDiPho->Fill(dPhiDiJetDiPho,wei);

  hMassJetJet->Fill(massJetJet,wei);

  hJetEtaProd_bc->Fill(jetEta[iJet1]*jetEta[iJet2]);
  //ProdEtaJets=jetEta[iJet1]*jetEta[iJet2];
  //dEtaJets=fabs(deltaEta_jj);
  //ZapJets=zeppenfeld;
  //MJets=massJetJet;

  // Di-Jet VBF Selection
  // if (corPho1.Et()>l_corr_gg.M()*55./120. && corPho2.Et()>25.) { // dijet VBF cut1
  if (jetPt[iJet1] > 30.) { // dijet VBF cut2
    hMgg_twoJetTag_jet1Pt->Fill(l_corr_gg.M(),wei);
    hDeltaEtaJets_bc->Fill(fabs(deltaEta_jj),wei);

    if( fabs(deltaEta_jj) > 3.5){ // dijet VBF cut3
      hMgg_twoJetTag_deltaEta->Fill(l_corr_gg.M(),wei);
      hZeppenfeld_bc->Fill(zeppenfeld,wei);
      
      if(fabs(zeppenfeld) < 2.5){ // dijet VBF cut4
	hMgg_twoJetTag_zep->Fill(l_corr_gg.M(),wei);
	hMassJetJet_bc->Fill(massJetJet,wei);

	if(massJetJet > 350){ // dijet VBF cut5
	  hMgg_twoJetTag_Mjets->Fill(l_corr_gg.M(),wei);
	  hDeltaPhiDiJetDiPho_bc->Fill(dPhiDiJetDiPho,wei);
	  
	  if (fabs(dPhiDiJetDiPho) > 2.6) { // dijet VBF cut6
	    hMgg_twoJetTag_dPhi->Fill(l_corr_gg.M(),wei);
	    
	  } // End of dijet VBF cut6 (dPhiDiJetDiPho > 2.6)
	} // End of dijet VBF cut5 (M_jj>350)
      } // End of dijet VBF cut4 (|Zeppenfeld|<2.5)
    } // End of dijet VBF cut3 (|deltaEta|>3.5)
  } // End of dijet VBF cut2 (lead Jet Pt > 30.)
  // } // End of dijet VBF cut1 (lead Pho Et > 55. * mass/120. && trail Pho Et > 25.)
  
  // Di-Jet VH Selection
  // if (corPho1.Et()>l_corr_gg.M()*60./120. && corPho2.Et()>25.) { // dijet VH cut1
  if (jetPt[iJet1] > 30.) { // dijet VH cut2
    hMgg_HSTRA_twoJetTag_jet1Pt->Fill(l_corr_gg.M() ,wei);
    hDeltaEtaJets_HSTRA_bc->Fill(fabs(deltaEta_jj),wei);
    
    if( fabs(deltaEta_jj) < 2.5){ // dijet VH cut3
      hMgg_HSTRA_twoJetTag_deltaEta->Fill(l_corr_gg.M(),wei);
      hZeppenfeld_HSTRA_bc->Fill(zeppenfeld,wei);
      
      if(fabs(zeppenfeld)<1.5){ // dijet VH cut4
	hMgg_HSTRA_twoJetTag_zep->Fill(l_corr_gg.M(),wei);
	hMassJetJet_HSTRA_bc->Fill(massJetJet,wei);

	if(fabs(massJetJet-85.)<30.){ // dijet VH cut5
          hMgg_HSTRA_twoJetTag_Mjets->Fill(l_corr_gg.M(),wei);

          if (corPho1.Pt()/l_corr_gg.M() > 0.5 && corPho2.Pt() > 25.) hstratag = 1;
	  
	} // End of dijet VH cut5 (|M_jj-85.|<30.)
      } // End of dijet VH cut4 (|Zeppenfeld|<1.5)
    } // End of dijet VH cut3 (|deltaEta|<2.5)
  } // End of dijet VH cut2 (lead Jet Pt > 30.)
  // } // End of dijet VH cut1 (lead Pho Et > 60. * mass/120. && trail Pho Et > 25.)

  // Daniele VBF jet category
  vbftag = 0;
  if( corPho1.Pt()/l_corr_gg.M() >60./120. && corPho2.Pt() > 25.) {

    if     ((corPho1.Pt()/l_corr_gg.M()>0.5)&&(corPho2.Pt()>25.)&&(jetPt[iJet1]>30.)&&(jetPt[iJet2]>30.)
	    &&(fabs(deltaEta_jj)>3.0)&&(fabs(zeppenfeld)< 2.5)&&(massJetJet>500)&&(fabs(dPhiDiJetDiPho)>2.6)) catjet = 0;
    else if((corPho1.Pt()/l_corr_gg.M()>0.5)&&(corPho2.Pt()>25.)&&(jetPt[iJet1]>30.)&&(jetPt[iJet2]>20.)
	    &&(fabs(deltaEta_jj)>3.0)&&(fabs(zeppenfeld)< 2.5)&&(massJetJet>250)&&(fabs(dPhiDiJetDiPho)>2.6)) catjet = 1;
  }


  _minitree->dijet_dEtaJJ    = deltaEta_jj;
  _minitree->dijet_zeppenfeld = zeppenfeld;
  _minitree->dijet_mJJ        = massJetJet;
  _minitree->dijet_dphiGGJJ   = dPhiDiJetDiPho;

  for( int ij = 0 ; ij < 2; ij++ ) {
    _minitree->dijet_ptj[ij] = jetPt[goodJetsIndex[ij]];
    _minitree->dijet_jec[ij] = jetPt[goodJetsIndex[ij]]/jetRawPt[goodJetsIndex[ij]];
    _minitree->dijet_etj[ij] = jetEt[goodJetsIndex[ij]];
    _minitree->dijet_eta[ij] = jetEta[goodJetsIndex[ij]];
    _minitree->dijet_phi[ij] = jetPhi[goodJetsIndex[ij]];
    _minitree->dijet_en[ij]  = jetEn[goodJetsIndex[ij]];
    _minitree->dijet_rms[ij] = jetDR2Mean[goodJetsIndex[ij]];
    _minitree->dijet_betaStar[ij]  = jetBetaStarClassicExt[goodJetsIndex[ij]][ivtx];
  }
  _minitree->dijet_ptj1 = _minitree->dijet_ptj[0];
  _minitree->dijet_ptj2 = _minitree->dijet_ptj[1];


  /// compute VBF mva discriminant
  myVBFLeadJPt= jetPt[goodJetsIndex[0]];
  myVBFSubJPt = jetPt[goodJetsIndex[1]];
  myVBF_Mjj   = massJetJet;
  myVBFdEta   = fabs(deltaEta_jj);
  myVBFZep    = fabs(zeppenfeld);
  myVBFdPhi   = fabs(dPhiDiJetDiPho);
  myVBF_Mgg   = l_corr_gg.M();
  myVBFDiPhoPtOverM   = l_corr_gg.Pt() / myVBF_Mgg;
  myVBFLeadPhoPtOverM = corPho1.Pt()   / myVBF_Mgg;
  myVBFSubPhoPtOverM  = corPho2.Pt()   / myVBF_Mgg;
  
  /// reset catjet for mva VBF selection
  if( DiscriVBF_useMvaSel ) catjet = -1;
  if( myVBFLeadPhoPtOverM > _config->diJetPtg1M() && 
      myVBFSubPhoPtOverM  > _config->diJetPtg2M() &&
      myVBFLeadJPt>30. && myVBFSubJPt>20. &&
      myVBF_Mjj           > _config->diJetMjj() ) {
    _minitree->mtree_fillDijetTree = true;
    _minitree->dijet_mvaVbf = DiscriVBF->EvaluateMVA(DiscriVBF_Method);
    _minitree->dijet_mvaVbf_sig    = 1-DiscriVBF->EvaluateMulticlass(DiscriVBF_Method)[0];
    _minitree->dijet_mvaVbf_dipho  = 1-DiscriVBF->EvaluateMulticlass(DiscriVBF_Method)[1];
    _minitree->dijet_mvaVbf_gluglu = 1-DiscriVBF->EvaluateMulticlass(DiscriVBF_Method)[2];
    
    if( DiscriVBF_useMvaSel )  {      
      if( myVBFLeadPhoPtOverM>1./3. && myVBFSubPhoPtOverM>1./4. && myVBFLeadJPt>30. && myVBFSubJPt>20. && myVBF_Mjj > 250. ) 
      for( unsigned icat = 0; icat < DiscriVBF_cat.size(); icat++ ) 
	if(  _minitree->dijet_mvaVbf >= DiscriVBF_cat[icat] ) { catjet = icat; break;  }
    }
  }

  if( catjet >= 0 ) vbftag = 1;
  if( vbftag == 1 ) _minitree->mtree_fillDijetTree = true;

  if( _minitree->mtree_fillDijetTree && catjet < 0 ) catjet = 999;
  

}

double xAna::correctJetPt(int i) {
  if(doJetRegression==0 || jetPt[i]<20 || fabs(jetEta[i])>2.5) return jetPt[i];
  else return JetRegression(i);
}

#endif
