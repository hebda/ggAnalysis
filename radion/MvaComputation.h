#include "xAna_allAna.h"
#include "VertexSelection.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <string>

float xAna::sumTrackPtInCone(int phoID, int vtxID, float minPt, 
					   float outerConeRadius, float innerConeRadius,
					   float etaStripHalfWidth, float dzMax, float d0Max) {
  if(vtxID < 0)return -999;
  TVector3 vtxpos(vtx[vtxID][0],vtx[vtxID][1],vtx[vtxID][2]);
  float sum = 0;
  for(int i = 0; i < nTrk; i++){
    TVector3 trackp(trkP[i][0],trkP[i][1],trkP[i][2]);
    if(trackp.Pt() < minPt)continue;
    TVector3 trkVtxVec(trkVtx[i][0],trkVtx[i][1],trkVtx[i][2]);
    double deltaz = fabs((vtxpos.Z() - trkVtxVec.Z()) - 
			 ( (trkVtxVec.X()-vtxpos.X())*trackp.Px() + (trkVtxVec.Y()-vtxpos.Y())*trackp.Py() )/trackp.Pt() * trackp.Pz()/trackp.Pt() );
    if(deltaz > dzMax)continue;
    double dxy = (-(trkVtxVec.X() - vtxpos.X())*trackp.Py() + (trkVtxVec.Y() - vtxpos.Y())*trackp.Px())/trackp.Pt();
    if(fabs(dxy) > d0Max)continue;
    double deta = fabs(phoEtaVtx[phoID][vtxID] - trackp.Eta());
    double dphi = fabs(phoPhiVtx[phoID][vtxID] - trackp.Phi());
    if(dphi > TMath::Pi())dphi = TMath::TwoPi() - dphi;
    double dR = sqrt(deta*deta + dphi*dphi);
    if(dR < outerConeRadius && dR >= innerConeRadius && deta >= etaStripHalfWidth)sum += trackp.Pt();
  }
  return sum;
}

float xAna::worstSumTrackPtInCone(int phoID, int &vtxID, float minPt, 
						float outerConeRadius, float innerConeRadius, 
						float etaStripHalfWidth, float dzMax, float d0Max) {
  int worstvtxID = -1;
  float maxisosum = -1;
  for(int i = 0; i < nVtx; i++){
    float isosum = sumTrackPtInCone(phoID,i,minPt,outerConeRadius,innerConeRadius,etaStripHalfWidth,dzMax,d0Max);
    if(isosum > maxisosum){
      maxisosum = isosum;
      worstvtxID = i;
    }
  }
  vtxID = worstvtxID;
  return maxisosum;
}


//////////////////////
//Setup PhotonID MVA//
//////////////////////
void xAna::Setup_MVA( ) {

  /// ------ setup 2011 mva phoID ------ ///
  string mvamethod = "BDT";
  cout << " Booking 2011 MVA Barrel photon ID with method: " << mvamethod << endl;
  phoID_2011[0]->AddVariable("HoE",&phoID_HoE);
  phoID_2011[0]->AddVariable("covIEtaIEta",&phoID_covIEtaIEta);
  phoID_2011[0]->AddVariable("tIso1abs",&phoID_tIso1abs);
  phoID_2011[0]->AddVariable("tIso3abs",&phoID_tIso3abs);
  phoID_2011[0]->AddVariable("tIso2abs",&phoID_tIso2abs);
  phoID_2011[0]->AddVariable("R9",&phoID_R9);
  phoID_2011[0]->AddVariable("absIsoEcal",&phoID_absIsoEcal);
  phoID_2011[0]->AddVariable("absIsoHcal",&phoID_absIsoHcal);
  phoID_2011[0]->AddVariable("NVertexes",&phoID_NVertexes);
  phoID_2011[0]->AddVariable("ScEta",&phoID_ScEta);
  phoID_2011[0]->AddVariable("EtaWidth",&phoID_EtaWidth);
  phoID_2011[0]->AddVariable("PhiWidth",&phoID_PhiWidth);
  phoID_2011[0]->BookMVA(mvamethod.c_str(),"mvaDiscriInputs/PhotonID_Barrel_Weights.xml");
  cout << " Booking 2011 MVA Endcap photon ID with method: " << mvamethod << endl;
  phoID_2011[1]->AddVariable("HoE",&phoID_HoE);
  phoID_2011[1]->AddVariable("covIEtaIEta",&phoID_covIEtaIEta);
  phoID_2011[1]->AddVariable("tIso1abs",&phoID_tIso1abs);
  phoID_2011[1]->AddVariable("tIso3abs",&phoID_tIso3abs);
  phoID_2011[1]->AddVariable("tIso2abs",&phoID_tIso2abs);
  phoID_2011[1]->AddVariable("R9",&phoID_R9);
  phoID_2011[1]->AddVariable("absIsoEcal",&phoID_absIsoEcal);
  phoID_2011[1]->AddVariable("absIsoHcal",&phoID_absIsoHcal);
  phoID_2011[1]->AddVariable("NVertexes",&phoID_NVertexes);
  phoID_2011[1]->AddVariable("ScEta",&phoID_ScEta);
  phoID_2011[1]->AddVariable("EtaWidth",&phoID_EtaWidth);
  phoID_2011[1]->AddVariable("PhiWidth",&phoID_PhiWidth);
  phoID_2011[1]->BookMVA(mvamethod.c_str(),"mvaDiscriInputs/PhotonID_Endcap_Weights.xml");
  /// ---------------------------- ///
  /// ------ setup 2012 mva phoID ------ ///
  cout << " Booking 2012 MVA Barrel photon ID with method: " << mvamethod << endl;
  phoID_2012[0] = new TMVA::Reader("!Color:Silent");
  phoID_2012[0]->AddVariable("ph.r9",                &phoID_R9 );
  phoID_2012[0]->AddVariable("ph.sigietaieta",       &phoID_covIEtaIEta );
  phoID_2012[0]->AddVariable("ph.scetawidth",        &phoID_EtaWidth);
  phoID_2012[0]->AddVariable("ph.scphiwidth",        &phoID_PhiWidth);
  phoID_2012[0]->AddVariable("ph.idmva_CoviEtaiPhi", &phoID_covIEtaIPhi);
  phoID_2012[0]->AddVariable("ph.idmva_s4ratio",     &phoID_S4);
  phoID_2012[0]->AddVariable("ph.idmva_GammaIso",    &phoID_pfPhoIso03 );
  phoID_2012[0]->AddVariable("ph.idmva_ChargedIso_selvtx",  &phoID_pfChIso03 );
  phoID_2012[0]->AddVariable("ph.idmva_ChargedIso_worstvtx",&phoID_pfChIso03worst );
  phoID_2012[0]->AddVariable("ph.sceta",             &phoID_ScEta );
  phoID_2012[0]->AddVariable("rho",                  &phoID_rho);
  phoID_2012[0]->BookMVA(mvamethod.c_str(),"mvaDiscriInputs/2012ICHEP_PhotonID_Barrel_BDT.weights.xml");
  cout << " Booking 2012 MVA Endcap photon ID with method: " << mvamethod << endl;
  phoID_2012[1]->AddVariable("ph.r9",                &phoID_R9 );
  phoID_2012[1]->AddVariable("ph.sigietaieta",       &phoID_covIEtaIEta );
  phoID_2012[1]->AddVariable("ph.scetawidth",        &phoID_EtaWidth);
  phoID_2012[1]->AddVariable("ph.scphiwidth",        &phoID_PhiWidth);
  phoID_2012[1]->AddVariable("ph.idmva_CoviEtaiPhi", &phoID_covIEtaIPhi);
  phoID_2012[1]->AddVariable("ph.idmva_s4ratio",     &phoID_S4);
  phoID_2012[1]->AddVariable("ph.idmva_GammaIso",    &phoID_pfPhoIso03 );
  phoID_2012[1]->AddVariable("ph.idmva_ChargedIso_selvtx",   &phoID_pfChIso03 );
  phoID_2012[1]->AddVariable("ph.idmva_ChargedIso_worstvtx", &phoID_pfChIso03worst );
  phoID_2012[1]->AddVariable("ph.sceta",             &phoID_ScEta );
  phoID_2012[1]->AddVariable("rho",                  &phoID_rho);
  phoID_2012[1]->AddVariable("ph.idmva_PsEffWidthSigmaRR",   &phoID_ESEffSigmaRR );
  phoID_2012[1]->BookMVA(mvamethod.c_str(),"mvaDiscriInputs/2012ICHEP_PhotonID_Endcap_BDT.weights_PSCorr.xml");
  /// ---------------------------- ///
  /// --- setup 2011 MVA dipho ---///
  cout << " Booking 2011 MVA diphoton dicriminant" << endl;
  DiscriDiPho_2011->AddVariable("masserrsmeared/mass",&Ddipho_masserr);
  DiscriDiPho_2011->AddVariable("masserrsmearedwrongvtx/mass",&Ddipho_masserrwrongvtx);
  DiscriDiPho_2011->AddVariable("vtxprob",&Ddipho_vtxprob);
  DiscriDiPho_2011->AddVariable("ph1.pt/mass",&Ddipho_piT1);
  DiscriDiPho_2011->AddVariable("ph2.pt/mass",&Ddipho_piT2);
  DiscriDiPho_2011->AddVariable("ph1.eta",&Ddipho_eta1);
  DiscriDiPho_2011->AddVariable("ph2.eta",&Ddipho_eta2);
  DiscriDiPho_2011->AddVariable("TMath::Cos(ph1.phi-ph2.phi)",&Ddipho_cosdPhi);
  DiscriDiPho_2011->AddVariable("ph1.idmva",&Ddipho_id1);
  DiscriDiPho_2011->AddVariable("ph2.idmva",&Ddipho_id2);
  //  DiscriDiPho_2011->BookMVA("BDTG","mvaDiscriInputs/MITDiPhotonWeights.xml");
  DiscriDiPho_2011->BookMVA("BDTG","mvaDiscriInputs/HggBambu_SMDipho_Jan16_BDTG.weights.xml");

  cout << " Booking 2012 MVA diphoton dicriminant" << endl;
  DiscriDiPho_2012->AddVariable("masserrsmeared/mass",&Ddipho_masserr);
  DiscriDiPho_2012->AddVariable("masserrsmearedwrongvtx/mass",&Ddipho_masserrwrongvtx);
  DiscriDiPho_2012->AddVariable("vtxprob",&Ddipho_vtxprob);
  DiscriDiPho_2012->AddVariable("ph1.pt/mass",&Ddipho_piT1);
  DiscriDiPho_2012->AddVariable("ph2.pt/mass",&Ddipho_piT2);
  DiscriDiPho_2012->AddVariable("ph1.eta",&Ddipho_eta1);
  DiscriDiPho_2012->AddVariable("ph2.eta",&Ddipho_eta2);
  DiscriDiPho_2012->AddVariable("TMath::Cos(ph1.phi-ph2.phi)",&Ddipho_cosdPhi);
  DiscriDiPho_2012->AddVariable("ph1.idmva",&Ddipho_id1);
  DiscriDiPho_2012->AddVariable("ph2.idmva",&Ddipho_id2);
  //  DiscriDiPho_2012->BookMVA("BDTG","mvaDiscriInputs/HggBambu_SMDipho_Jun19_BDTG.weights.xml");
  DiscriDiPho_2012->BookMVA("BDTG","mvaDiscriInputs/HggBambu_SMDipho_Jul07_BDTG.weights.xml");


  /// ------------------------------///
  /// ----- setup 2012 MVA VBF -----///
  DiscriVBF = new TMVA::Reader("!Color:Silent");
  string tmvaVbfMvaWeights = "mvaDiscriInputs/TMVA_vbf_6var_mjj100_diphopt_phopt_BDTG.weights.xml";
  DiscriVBF_Method = "BDTG";

  cout << " Booking 2012 MVA vbf dicriminant" << endl;
  DiscriVBF = new TMVA::Reader( "!Color:!Silent" );
  DiscriVBF->AddVariable("jet1pt", &myVBFLeadJPt);
  DiscriVBF->AddVariable("jet2pt", &myVBFSubJPt);
  DiscriVBF->AddVariable("abs(jet1eta-jet2eta)", &myVBFdEta);
  DiscriVBF->AddVariable("mj1j2", &myVBF_Mjj);
  DiscriVBF->AddVariable("zepp" , &myVBFZep);
  DiscriVBF->AddVariable("dphi" , &myVBFdPhi);
  if( DiscriVBF_UseDiPhoPt )  DiscriVBF->AddVariable("diphopt/diphoM", &myVBFDiPhoPtOverM); 
  if( DiscriVBF_UsePhoPt   ) {
      DiscriVBF->AddVariable("pho1pt/diphoM", &myVBFLeadPhoPtOverM);
      DiscriVBF->AddVariable("pho2pt/diphoM", &myVBFSubPhoPtOverM);
    }
  DiscriVBF->BookMVA( DiscriVBF_Method, tmvaVbfMvaWeights );
  // end of MVA VBF 


  //setup of MVA for b-jet energy regression
  if(doJetRegression==1){
    cout<<" Booking b-jet energy regression" << endl;
    jetRegres->AddVariable( "jet_eta", &jetRegEta);
    jetRegres->AddVariable( "jet_emfrac", &jetRegEMFrac);
    jetRegres->AddVariable( "jet_hadfrac", &jetRegHadFrac);
    jetRegres->AddVariable( "jet_nconstituents", &jetRegNConstituents);
    jetRegres->AddVariable( "jet_vtx3dL", &jetRegVtx3dL);
    jetRegres->AddVariable( "MET", &jetRegMET);
    jetRegres->AddVariable( "jet_dPhiMETJet", &jetRegDPhiMETJet);
    jetRegres->BookMVA("BDTG0","jetRegressionWeights/TMVARegression_Cat0_BDTG.weights.xml");  
    jetRegres->BookMVA("BDTG1","jetRegressionWeights/TMVARegression_Cat1_BDTG.weights.xml");  
  }
  //end of setup for b-jet energy regression


  /// -------------------------------///
  /// ----- setup vertex finder -----///
  vector<string> tmvaPerVtxVariables;
  tmvaPerVtxVariables.push_back("ptbal");
  tmvaPerVtxVariables.push_back("ptasym");
  tmvaPerVtxVariables.push_back("logsumpt2");  
  tmvaPerVtxVariables.push_back("limPullToConv");
  tmvaPerVtxVariables.push_back("nConv");
  string perVtxMvaWeights = "mvaDiscriInputs/TMVAClassification_BDTCat_conversions_tmva_407.weights.xml";
  _perVtxMvaMethod  = "BDTCat_conversions";
  string perEvtMvaWeights = "mvaDiscriInputs/TMVAClassification_evtBDTG_conversions_tmva_407.weights.xml";
  _perEvtMvaMethod  = "evtBDTG";
  string vertexsetup = "configFilesDir/vertex_selection_7TeV.dat";
  if( _config->setup().find("2012") != string::npos ) {
    perEvtMvaWeights = "mvaDiscriInputs/TMVAClassification_BDTvtxprob2012.weights.xml";
    vertexsetup      = "configFilesDir/vertex_selection_hcp2012.dat";
    _perEvtMvaMethod  = "BDTvtxprob2012";
  }
  _vtxAlgoParams = VertexAlgoParametersReader(vertexsetup);
  _vtxAna  = new HggVertexAnalyzer(_vtxAlgoParams);

  _perVtxReader = 0;
  _perEvtReader = 0;
  if(  perVtxMvaWeights != "" ) { 
    _perVtxReader = new TMVA::Reader( "!Color:!Silent" );
    HggVertexAnalyzer::bookVariables( *_perVtxReader, tmvaPerVtxVariables );
    _perVtxReader->BookMVA( _perVtxMvaMethod, perVtxMvaWeights );
  }
  if( perEvtMvaWeights != "" ) {
    _perEvtReader = new TMVA::Reader( "!Color:!Silent" );
    HggVertexAnalyzer::bookPerEventVariables( *_perEvtReader );
    _perEvtReader->BookMVA( _perEvtMvaMethod, perEvtMvaWeights );
  }
  cout << " Booking Vertex discriminants" << endl;
  cout << "       - use conversions     : " << _vtxAlgoParams.useAllConversions << endl;
  cout << "       - single leg sigma1Pix: " << _vtxAlgoParams.singlelegsigma1Pix << endl;
  cout << "       - check that if use single leg conversion (option=3) then sigma1Pix has meaningful value "
       << endl;

  // end of Vertex MVA 


}


float xAna::PhoID_MVA( int i, int selVtx ) {
  ////////////////////
  //do photon id mva//
  ////////////////////
  phoID_HoE         = phoHoverE[i];
  phoID_covIEtaIEta = phoSigmaIEtaIEta[i];
  phoID_R9          = phoR9[i];
  phoID_absIsoEcal  = phoEcalIsoDR03[i] - 0.17*rho25;
  phoID_absIsoHcal  = phoHcalIsoDR04[i] - 0.17*rho25;
  phoID_NVertexes   = nVtx;
  phoID_ScEta       = phoSCEta[i];
  phoID_EtaWidth    = phoSCEtaWidth[i];
  phoID_PhiWidth    = phoSCPhiWidth[i];


  if( _config->setup() == "Prompt2012_ichep" ) {
    phoID_covIEtaIPhi =  phoSigmaIEtaIPhi[i];
    phoID_S4 = phoS4ratio[i];
    phoID_ESEffSigmaRR = phoESEffSigmaRR[i][0];
    phoID_rho = rho2012;
    phoID_pfPhoIso03 = phoCiCPF4phopfIso03[i];
    phoID_pfChIso03  = phoCiCPF4chgpfIso03[i][selVtx];
    phoID_pfChIso03worst = -99;
    for( int iv=0; iv<nVtxBS ; iv++) 
      if( phoCiCPF4chgpfIso03[i][iv] >  phoID_pfChIso03worst )
	phoID_pfChIso03worst = phoCiCPF4chgpfIso03[i][iv];
  } else {
    int badvtx;
    float pho_tkiso_badvtx = worstSumTrackPtInCone(i,badvtx,0,0.4,0.02,0.0,1.0,0.1);
    float pho_tkiso_goodvtx = sumTrackPtInCone(i,selVtx,0,0.3,0.02,0.0,1.0,0.1);
    phoID_tIso1abs = (pho_tkiso_goodvtx + phoEcalIsoDR03[i] + phoHcalIsoDR04[i] - 0.17*rho25);
    phoID_tIso2abs = (pho_tkiso_badvtx + phoEcalIsoDR04[i] + phoHcalIsoDR04[i] - 0.52*rho25);
    phoID_tIso3abs = pho_tkiso_goodvtx;
  }

  return -1;
}


float xAna::DiPhoID_MVA( const TLorentzVector &glead, const TLorentzVector &gtrail, const TLorentzVector &hcand, 
				       const MassResolution & massResoCalc, float vtxProb,
				       const float &idlead, const float &idtrail) {
  //  Ddipho_masserr         = massResoCalc.relativeMassResolutionCorrVertex();
  Ddipho_masserr         = massResoCalc.relativeMassResolutionFab_energy();
  Ddipho_masserrwrongvtx = massResoCalc.relativeMassResolutionWrongVertex();
  Ddipho_vtxprob = vtxProb;
  Ddipho_piT1 = glead.Pt()/hcand.M();
  Ddipho_piT2 = gtrail.Pt()/hcand.M();
  Ddipho_eta1 = glead.Eta();
  Ddipho_eta2 = gtrail.Eta();
  Ddipho_cosdPhi = TMath::Cos(glead.Phi()-gtrail.Phi());
  Ddipho_id1 = idlead;
  Ddipho_id2 = idtrail;
  
  float diphomva = -1;
  if(  Ddipho_id1 > -0.3 &&  Ddipho_id2 > -0.3  ) 
    diphomva = Ddipho_mva->EvaluateMVA("BDTG");
  
  
  /// fill minitree with dipho inputs
  _minitree->mtree_massResoRightVtx = Ddipho_masserr;
  _minitree->mtree_massResoRandVtx  = Ddipho_masserrwrongvtx;
  _minitree->mtree_vtxProb          = Ddipho_vtxprob;
  
  _minitree->mtree_diphoMva = diphomva;

  return diphomva;
}

float xAna::JetRegression(int iJet){

  jetRegPt = jetPt[iJet];
  jetRegEta = jetEta[iJet];
  jetRegEMFrac = jetCEF[iJet]+jetNEF[iJet]+jetMEF[iJet];
  jetRegHadFrac = jetCHF[iJet]+jetNHF[iJet];
  jetRegNConstituents = jetNConstituents[iJet];
  jetRegVtx3dL = jetVtx3dL[iJet];
  jetRegMET = recoPfMET;
  float tmpPi = 3.1415927, tmpDPhi=fabs(jetPhi[iJet]-recoPfMETPhi);
  if(tmpDPhi>tmpPi) tmpDPhi=2*tmpPi-tmpDPhi;
  jetRegDPhiMETJet = tmpDPhi;

  if(doJetRegression==1){
    if(jetPt[iJet]<90) return jetPt[iJet]*jetRegres->EvaluateRegression("BDTG0")[0];
    if(jetPt[iJet]>90) return jetPt[iJet]*jetRegres->EvaluateRegression("BDTG1")[0];
  }
  else return jetPt[iJet];
}
