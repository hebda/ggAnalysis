#ifndef MCtruth_functions_hh__
#define MCtruth_functions_hh__
#include "xAna_allAna.h"



bool xAna::isZgamma(void) {

  for ( int i=0; i<nMC; ++i) { 
    if( (mcPID[i] == 22 && mcMomPID[i] == 22 &&  mcGMomPID[i] == 22) || 
	(mcPID[i] == 22 && mcMomPID[i] == 22 && ( mcGMomPID[i] == 21 ||  fabs(mcGMomPID[i])<6) ) ||
	(mcPID[i] == 22 && mcMomPID[i] == 22 &&  ( mcGMomPID[i] == 23 || fabs(mcGMomPID[i])==24)) ||
	(mcPID[i] == 22 &&( mcMomPID[i] == 23 || fabs(mcMomPID[i])==24)) ||
	(mcPID[i] == 22 && (fabs(mcMomPID[i]) == 11 || fabs(mcMomPID[i]) == 13 ) && (fabs(mcGMomPID[i]) == 24 || fabs(mcGMomPID[i]) == 23 ) ) 
	){
      cout << "----- ientry = " << -1 <<" removing  Zgamma  -------------" <<endl;
      cout <<i<<"  mcPID[i] "<<   mcPID[i]  <<   "  mcMomPID[i] "<<   mcMomPID[i]   <<"  mcGMomPID[i]  "<< mcGMomPID[i] << endl;
      return true;
    } // if DC
  } // loop over mc
  // mode Zjets
  return false;
}


void xAna::fillMCtruthInfo_HiggsSignal( void ) {
  int hind(-1);
  vector<int> gind;
  for ( int i=0; i<nMC; ++i) { 
    if(  mcPID[i]     == 25 && mcMomPID[i] == 25 ) hind = i;         //higgs 
    if(  mcGMomPID[i] == 25 && mcPID[i]    == 22 ) gind.push_back(i); //gammas
    /// for JHU, the Higgs is not there
    else if ( _weight_manager->getCrossSection()->getHiggsMass() < 10 && mcPID[i]    == 22 ) gind.push_back(i); //gammas
  }
 
  
  if( gind.size() != unsigned(2) ) {
    cerr << " fillMCtruthInfo_HiggsSignal[ERROR]: # of photons != 2, nphotons = "
	 << gind.size() << endl;
    if( gind.size() != 2 ) return;
  }

  TLorentzVector g[2],h;
  for( unsigned ip = 0 ; ip <  gind.size(); ip++ ) {
    g[ip].SetPtEtaPhiE( mcPt[ gind[ip]], mcEta[gind[ip]], 
			mcPhi[gind[ip]], mcE[  gind[ip]] );
  }
  if( hind >= 0 )  h.SetPtEtaPhiE( mcPt[hind], mcEta[hind], mcPhi[hind], mcE[hind] );
  else             h = g[0]+g[1];

  if( hind >=0 && fabs((g[0]+g[1]).M() - h.M() ) > 0.01 )  {    
    cerr << " fillMCtruthInfo_HiggsSignal[ERROR]: the diphoton mass is not the same mass as the Higgs mass " << endl
	 << "  - higgs mass   : " << h.M() << endl
	 << "  - diphoton mass: " << (g[0]+g[1]).M() << endl;
    return;
  }
  _minitree->mc_cThetaStar_CS = 
    2*(g[0].E()*g[1].Pz() - g[1].E()*g[0].Pz())/(h.M()*sqrt(h.M2()+h.Pt()*h.Pt()));
  _minitree->mc_mH = h.M();
  
}



#endif
