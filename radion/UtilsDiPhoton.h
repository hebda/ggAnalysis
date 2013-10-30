#include "xAna_allAna.h"

#include <TMath.h>
#include <stdexcept>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

bool xAna::selectTwoPhotons( int &i, int &j, int selVtx,
					   TLorentzVector &g1, TLorentzVector &g2, bool vetoElec[2] ) {

  if( selVtx < 0 || selVtx > nVtx ) 
    cout << " selected vertex no defined; ivtx : " << selVtx << " / " << nVtx << endl;

  /// do photonid preselection
  float etai = phoEtaVtx[i][selVtx];
  float phii = phoPhiVtx[i][selVtx];
  float etaj = phoEtaVtx[j][selVtx];
  float phij = phoPhiVtx[j][selVtx];
  g1.SetPtEtaPhiM( phoE[i]/cosh(etai),etai,phii, 0 ); 
  g2.SetPtEtaPhiM( phoE[j]/cosh(etaj),etaj,phij, 0 ); 
  int ilead = i, itrail = j;	  
  if( g1.Pt() < g2.Pt() ) {
    ilead = j; itrail = i;
    g1.SetPtEtaPhiM( phoE[j]/cosh(etaj),etaj,phij, 0 ); 
    g2.SetPtEtaPhiM( phoE[i]/cosh(etai),etai,phii, 0 ); 
  }
  
  int gID1 = photonID(ilead , selVtx, g1.Pt(), !vetoElec[0] );
  int gID2 = photonID(itrail, selVtx, g2.Pt(), !vetoElec[1] );
  if ( doControlSample == 0 && ( gID1 == 0 || gID2 == 0 ) ) return false; 
  else if ( doControlSample == 1 ){//redo one photon ID for the control sample
    int gID1_inv=0, gID2_inv=0;
    if( !mvaPhotonID(ilead,selVtx,g1.Pt(),!vetoElec[0]) ) gID1_inv=0;
    else gID1_inv = 1 - cic4PFPhotonID(ilead,selVtx,g1.Pt(),!vetoElec[0]);//inverted
    if( !mvaPhotonID(itrail,selVtx,g2.Pt(),!vetoElec[1]) ) gID1_inv=0;
    else gID1_inv = 1 - cic4PFPhotonID(itrail,selVtx,g2.Pt(),!vetoElec[1]);//inverted

    //Now ensure that there is one good photon and one good inverted photon
    if( !((gID1==1 && gID2_inv==1) || (gID1_inv==1 && gID2==1)) ) return false;
  }


  /// apply loose mva photon id cut inside the loop for mva selection
  if ( _config->analysisType() == "MVA" ) 
  for( int ii = 0 ; ii < 2; ii++ ) {
      int index = -1;
      if( ii == 0 ) index = ilead ;
      if( ii == 1 ) index = itrail;
      
      PhoID_MVA(index,selVtx);
      float phoid = -1;
      if( fabs(phoSCEta[index]) < 1.45 ) phoid = phoID_mva[0]->EvaluateMVA("BDT");
      else                               phoid = phoID_mva[1]->EvaluateMVA("BDT");
      
      if( phoid < _config->photonIdMvaCut() ) return false;
  }


  i = ilead;
  j = itrail;  

  return true;
}


void xAna::ReweightMC_phoIdVar( int i ) {    
  
  if       ( _config->setup() == "ReReco2011" && _config->analysisType() == "baseline" ) {
    if(fabs(phoSCEta[i])<1.479)  phoR9[i] = phoR9[i]*1.0048;
    else                         phoR9[i] = phoR9[i]*1.00492;
  } else if( _config->setup() == "ReReco2011" && _config->analysisType() == "MVA"      ) {
    phoR9[i]=phoR9[i]*1.0035;
    if( fabs(phoSCEta[i])<1.479  ) phoSigmaIEtaIEta[i] = phoSigmaIEtaIEta[i]*0.87+0.0011;
    else                           phoSigmaIEtaIEta[i] = phoSigmaIEtaIEta[i]*0.99;
    phoSCEtaWidth[i] *= 0.98;
    phoSCPhiWidth[i] *= 0.99;
  } else if( _config->setup() == "Prompt2012_ichep" && _config->analysisType() == "baselineCiC4PF" ) {
    /// need to be clarify: gglobe ue 
    if(fabs(phoSCEta[i])<1.479)   {
      phoR9[i]=phoR9[i]*1.0053;
      phoSigmaIEtaIEta[i]=phoSigmaIEtaIEta[i]*0.998;
    } else{
      phoR9[i]=phoR9[i]*1.0075;
      phoSigmaIEtaIEta[i]=phoSigmaIEtaIEta[i]*1.00;
    }
  } else if( _config->setup() == "Prompt2012_ichep" && _config->analysisType() == "MVA" ) {
    if(fabs(phoSCEta[i])<1.479) {
      phoR9[i]            = 1.0045*phoR9[i] + 0.0010;
      phoSigmaIEtaIEta[i] = 0.891832*phoSigmaIEtaIEta[i] + 0.0009133;
      phoSCPhiWidth[i]    = 1.00002 *phoSCPhiWidth[i] - 0.000371;
      phoSCEtaWidth[i]    = 1.04302 *phoSCEtaWidth[i] - 0.000618;
      phoS4ratio[i]       = 1.01894 *phoS4ratio[i]    - 0.01034;
    } else {
      phoR9[i]            = 1.0086*phoR9[i]- 0.0007 ;
      phoSigmaIEtaIEta[i] = 0.99470 *phoSigmaIEtaIEta[i]+ 0.00003;
      phoSCPhiWidth[i]    = 0.99992 *phoSCPhiWidth[i] - 0.00000048;
      phoSCEtaWidth[i]    = 0.903254*phoSCEtaWidth[i] + 0.001346;
      phoS4ratio[i]       = 1.04969 *phoS4ratio[i]    - 0.03642;
    }
  } else {
    throw std::runtime_error(string("Supported configs (setup,analysisType): \n") +
			     " (ReReco2011,baseline) ; (ReReco2011,MVA) ; (Prompt2012_ichep,baselineCiC4PF) ; (Prompt2012_ichep,MVA) ") ;
  }
}

void xAna::ReweightMC_phoIdVarBaseline( int i ) {
  if(fabs(phoSCEta[i])<1.479)
    phoR9[i]=phoR9[i]*1.0048;
  else
    phoR9[i]=phoR9[i]*1.00492;

}







bool xAna::isAGoodConv( int iconv  ) {
  /// take this part of the code from CMSSW ConvertionTools

  float probMin = 1e-6;
  float lxymin  = 2.0;

  /// if it is matched to a conversion, check its quality
  /// vertex validity is there by default
  //  if( !phoConvValidVtx[ipho]          ) return false;

  if( convChi2Probability[iconv] < probMin ) return false; 

  //compute transverse decay length
  //  XYZVector mom(convRefittedMomentum[iconv][0],convRefittedMomentum[iconv][1],convRefittedMomentum[iconv][2]);
  double xx = convRefittedMomentum[iconv][0];
  double yy = convRefittedMomentum[iconv][1];
  double dbsx = convVtx[iconv][0] - bspotPos[0];
  double dbsy = convVtx[iconv][1] - bspotPos[1];
  double lxy = (xx*dbsx + yy*dbsy)/sqrt(xx*xx+yy*yy);
  
  if( lxy < lxymin ) return false;

  if( convNHitsBeforeVtx[iconv][0] > 0 ) return false;
  if( convNHitsBeforeVtx[iconv][1] > 0 ) return false;
    
  return true;
}

bool xAna::isMatchedToAPromptElec( int ipho, float &dRclosestElec ) {
  /// take this part of the code from CMSSW ConvertionTools
  TLorentzVector pele;  
  bool match = false;
  dRclosestElec = 999;
  for( int iele = 0 ; iele < nEle; iele++ ) {
    if( eleMissHits[iele] > 0 ) continue;

    float deta = fabs( eleSCEta[iele] - phoSCEta[ipho] );
    float dphi = acos( cos( eleSCPhi[iele] - phoSCPhi[ipho] ) );
    float dR   = sqrt(deta*deta + dphi*dphi );
    if( dR < dRclosestElec ) dRclosestElec = dR;

    /// matched the photon and the electron from SC position
    if( deta > 0.01 ) continue;
    if( dphi > 0.01 ) continue;
    
    /// if this is matched to an electron then check if it is conversion
    /// this is slightly different from the original code
    /// but I have only these information here (should bre really close though)
    /* if( !phoIsConv[ipho] ) { */
    /*   /// if the photon has no conversion track, then match */
    /*   match = true; */
    /*   continue; */
    /* } else {       */
    /*   /// else try to find a conversion that would match the electron and check it */
    /*   pele.SetPtEtaPhiM( elePt[iele], eleEta[iele], elePhi[iele], 0 ); */
    /*   bool isConv = false; */
    /*   for( int iconv = 0; iconv < nConv; iconv++ ) { */
    /* 	TLorentzVector pconv; */
    /* 	pconv.SetPxPyPzE( convP4[iconv][0],convP4[iconv][1],convP4[iconv][2],convP4[iconv][3]); */
    /* 	if( pconv.Pt() < 1 ) continue; */
    /* 	if( pele.DeltaR( pconv ) < 0.1 && isAGoodConv( iconv ) ) {isConv = true; break;}  */
    /*   } */
    /*   if( isConv ) continue; */
    /*   match = true; */
    /* } */
    match = true;
  }

  if( match ) {
    /// now check this is not a conversion
    /// try something different: see conversion tools
    TVector3 calpos(phoCaloPos[ipho][0], phoCaloPos[ipho][1],phoCaloPos[ipho][2]);
    bool isConv = false;
    for( int iconv = 0; iconv < nConv; iconv++ ) { 
      TVector3 convP3( convRefittedMomentum[iconv][0], convRefittedMomentum[iconv][1], convRefittedMomentum[iconv][2] );
      TVector3 convVTX(convVtx[iconv][0], convVtx[iconv][1],convVtx[iconv][2]);
      TVector3 phoFromConvVtx = calpos - convVTX;
      float drConv = phoFromConvVtx.DeltaR( convP3 ); 
      if( drConv < 0.1 && isAGoodConv( iconv ) ) { isConv = true; break;}
    }
    if( isConv ) match = false;
  }

  return match ;
}
 
bool xAna::isInGAP_EB( int ipho ) {
  if( abs(phoSCEta[ipho]) > 1.45 ) return false;
  //  if( phoPos[ipho] == 2 || phoPos[ipho] == 4 ) return true;
  int detid1 = abs(phoSeedDetId1[ipho]);
  int detid2 = abs(phoSeedDetId2[ipho]);
  if( fabs(detid1-  0.5) < 5 || fabs(detid1- 25.5) < 5 || fabs(detid1- 45.5) < 5 ||
      fabs(detid1- 65.5) < 5 || fabs(detid1- 85.5) < 5 ) return true;

  if( fabs( detid2%20-10.5 ) > 5 ) return true;
  return false;
}

#include <string>
#include <sstream>
std::string itostr(int i ) {
  std::stringstream out;
  out << i;
  return out.str();
}

