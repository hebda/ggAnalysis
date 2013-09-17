#ifndef __energycorrection_hh__
#define __energycorrection_hh__

#include <TRandom3.h> 
 



class photonOverSmearing {
 public:
  enum category { EBlowEtaGold_CM = 0, EBlowEtaGold_GAP = 1, 
		  EBlowEtaBad_CM  = 2, EBlowEtaBad_GAP  = 3,
		  EBhighEtaGold   = 4, EBhighEtaBad     = 5,
		  EElowEtaGold    = 6, EElowEtaBad      = 7,
		  EEhighEtaGold   = 8, EEhighEtaBad     = 9,
		  UNKNOWN = 10 };


  photonOverSmearing(string setupType);
  photonOverSmearing::category photonOverSmearingCategory( float scEta, float r9, bool isEBGap );

  float meanOverSmearing( float scEta, float r9, bool isEBGap, int syst = 0 );
  float randOverSmearing( float scEta, float r9, bool isEBGap, int syst = 0 );

 private:
  vector<float> _oversmearing;
  vector<float> _oversmearing_err;

  TRandom3 overSmearingGen;
};

photonOverSmearing::photonOverSmearing(string setupType) {
  _oversmearing.resize(10,-1);
  _oversmearing_err.resize(10,-1);

  /// given in percents

  /// Jan16 smearing numbers
  /*
  _oversmearing[photonOverSmearing::EBlowEtaGold_CM ] = 0.67;
  _oversmearing[photonOverSmearing::EBlowEtaGold_GAP] = 0.77;
  _oversmearing[photonOverSmearing::EBlowEtaBad_CM  ] = 1.00;
  _oversmearing[photonOverSmearing::EBlowEtaBad_GAP ] = 0.89;
  _oversmearing[photonOverSmearing::EBhighEtaGold   ] = 1.37;
  _oversmearing[photonOverSmearing::EBhighEtaBad    ] = 1.88;
  _oversmearing[photonOverSmearing::EElowEtaGold    ] = 2.68;
  _oversmearing[photonOverSmearing::EElowEtaBad     ] = 2.79;
  _oversmearing[photonOverSmearing::EEhighEtaGold   ] = 2.93;
  _oversmearing[photonOverSmearing::EEhighEtaBad    ] = 3.01;
*/
  _oversmearing_err[photonOverSmearing::EBlowEtaGold_CM ] = 0.36;
  _oversmearing_err[photonOverSmearing::EBlowEtaGold_GAP] = 0.26;
  _oversmearing_err[photonOverSmearing::EBlowEtaBad_CM  ] = 0.25;
  _oversmearing_err[photonOverSmearing::EBlowEtaBad_GAP ] = 0.25;
  _oversmearing_err[photonOverSmearing::EBhighEtaGold   ] = 0.73;
  _oversmearing_err[photonOverSmearing::EBhighEtaBad    ] = 0.61;
  _oversmearing_err[photonOverSmearing::EElowEtaGold    ] = 0.95;
  _oversmearing_err[photonOverSmearing::EElowEtaBad     ] = 0.34;
  _oversmearing_err[photonOverSmearing::EEhighEtaGold   ] = 0.37;
  _oversmearing_err[photonOverSmearing::EEhighEtaBad    ] = 0.55;

  if( setupType == "Prompt2012_ichep" ) {
    /// numbers from sync http: https://twiki.cern.ch/twiki/bin/viewauth/CMS/UnblindingChecklist#Photon_energy_error 
    /* _oversmearing[photonOverSmearing::EBlowEtaGold_GAP] = 1.03; */
    /* _oversmearing[photonOverSmearing::EBlowEtaBad_GAP ] = 1.13; */
    /* _oversmearing[photonOverSmearing::EBlowEtaGold_CM ] = 0.81; */
    /* _oversmearing[photonOverSmearing::EBlowEtaBad_CM  ] = 1.26; */
    /* _oversmearing[photonOverSmearing::EBhighEtaGold   ] = 1.95; */
    /* _oversmearing[photonOverSmearing::EBhighEtaBad    ] = 2.22; */

    /* _oversmearing[photonOverSmearing::EElowEtaGold    ] = 3.66; */
    /* _oversmearing[photonOverSmearing::EElowEtaBad     ] = 3.34; */
    /* _oversmearing[photonOverSmearing::EEhighEtaGold   ] = 5.28; */
    /* _oversmearing[photonOverSmearing::EEhighEtaBad    ] = 5.58; */

    /// new numbers from HN https://hypernews.cern.ch/HyperNews/CMS/get/higgs2g/829.html
    /* _oversmearing[photonOverSmearing::EBlowEtaGold_GAP] = 1.08; */
    /* _oversmearing[photonOverSmearing::EBlowEtaGold_CM ] = 1.06; */
    /* _oversmearing[photonOverSmearing::EBlowEtaBad_GAP ] = 1.19; */
    /* _oversmearing[photonOverSmearing::EBlowEtaBad_CM  ] = 1.19; */

    /* _oversmearing[photonOverSmearing::EBhighEtaGold   ] = 1.49; */
    /* _oversmearing[photonOverSmearing::EBhighEtaBad    ] = 2.40; */

    /* _oversmearing[photonOverSmearing::EElowEtaGold    ] = 3.75; */
    /* _oversmearing[photonOverSmearing::EEhighEtaGold   ] = 6.07; */
    /* _oversmearing[photonOverSmearing::EElowEtaBad     ] = 3.30; */
    /* _oversmearing[photonOverSmearing::EEhighEtaBad    ] = 6.02; */

    /// new numbers for 53X: see https://twiki.cern.ch/twiki/pub/CMS/EcalEnergyScaleWithZeeForHgg/Hgg-scale-RUN2012ABC-53-1.pdf
    /* _oversmearing[photonOverSmearing::EBlowEtaGold_GAP] = 0.98; */
    /* _oversmearing[photonOverSmearing::EBlowEtaGold_CM ] = 0.99; */
    /* _oversmearing[photonOverSmearing::EBlowEtaBad_GAP ] = 0.94; */
    /* _oversmearing[photonOverSmearing::EBlowEtaBad_CM  ] = 1.13; */

    /* _oversmearing[photonOverSmearing::EBhighEtaGold   ] = 1.95; */
    /* _oversmearing[photonOverSmearing::EBhighEtaBad    ] = 1.98; */

    /* _oversmearing[photonOverSmearing::EElowEtaGold    ] = 2.91; */
    /* _oversmearing[photonOverSmearing::EElowEtaBad     ] = 2.75; */
    /* _oversmearing[photonOverSmearing::EEhighEtaGold   ] = 3.25;     */
    /* _oversmearing[photonOverSmearing::EEhighEtaBad    ] = 3.56; */

    //FOR MORIOND RUN ABCD SMEARING 
    _oversmearing[photonOverSmearing::EBlowEtaGold_GAP] = 1.11;
    _oversmearing[photonOverSmearing::EBlowEtaGold_CM ] = 1.11;
    _oversmearing[photonOverSmearing::EBlowEtaBad_GAP ] = 1.07;
    _oversmearing[photonOverSmearing::EBlowEtaBad_CM  ] = 1.07;
    _oversmearing[photonOverSmearing::EBhighEtaGold   ] = 1.55;
    _oversmearing[photonOverSmearing::EBhighEtaBad    ] = 1.94;
    
    _oversmearing[photonOverSmearing::EElowEtaGold    ] = 2.95;
    _oversmearing[photonOverSmearing::EElowEtaBad     ] = 2.76;
    _oversmearing[photonOverSmearing::EEhighEtaGold   ] = 3.70;    
    _oversmearing[photonOverSmearing::EEhighEtaBad    ] = 3.71;
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_GAP] = sqrt(0.07*0.07+0.22*0.22);
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_CM ] = sqrt(0.07*0.07+0.22*0.22);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_GAP ] = sqrt(0.06*0.06+0.24*0.24);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_CM  ] = sqrt(0.06*0.06+0.24*0.24);
    _oversmearing_err[photonOverSmearing::EBhighEtaGold   ] = sqrt(0.40*0.40+0.60*0.60);
    _oversmearing_err[photonOverSmearing::EBhighEtaBad    ] = sqrt(0.11*0.11+0.59*0.59);
    _oversmearing_err[photonOverSmearing::EElowEtaGold    ] = sqrt(0.25*0.25+0.90*0.90);
    _oversmearing_err[photonOverSmearing::EElowEtaBad     ] = sqrt(0.13*0.13+0.30*0.30);
    _oversmearing_err[photonOverSmearing::EEhighEtaGold   ] = sqrt(0.11*0.11+0.34*0.34);
    _oversmearing_err[photonOverSmearing::EEhighEtaBad    ] = sqrt(0.16*0.16+0.52*0.52);

  }

  if( setupType == "oversmear_hcp2012" ) {
    /// HCP2012 numbers
    _oversmearing[photonOverSmearing::EBlowEtaGold_GAP] = 1.10;
    _oversmearing[photonOverSmearing::EBlowEtaGold_CM ] = 1.00;
    _oversmearing[photonOverSmearing::EBlowEtaBad_GAP ] = 1.06;
    _oversmearing[photonOverSmearing::EBlowEtaBad_CM  ] = 1.06;
    _oversmearing[photonOverSmearing::EBhighEtaGold   ] = 1.86;
    _oversmearing[photonOverSmearing::EBhighEtaBad    ] = 1.96;

    _oversmearing[photonOverSmearing::EElowEtaGold    ] = 2.83;
    _oversmearing[photonOverSmearing::EElowEtaBad     ] = 2.67;
    _oversmearing[photonOverSmearing::EEhighEtaGold   ] = 3.43;    
    _oversmearing[photonOverSmearing::EEhighEtaBad    ] = 3.45;
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_GAP] = sqrt(0.14*0.14+0.22*0.22);
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_CM ] = sqrt(0.28*0.28+0.22*0.22);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_GAP ] = sqrt(0.08*0.08+0.24*0.24);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_CM  ] = sqrt(0.08*0.08+0.24*0.24);
    _oversmearing_err[photonOverSmearing::EBhighEtaGold   ] = sqrt(0.42*0.42+0.60*0.60);
    _oversmearing_err[photonOverSmearing::EBhighEtaBad    ] = sqrt(0.14*0.14+0.59*0.59);
    _oversmearing_err[photonOverSmearing::EElowEtaGold    ] = sqrt(0.31*0.31+0.90*0.90);
    _oversmearing_err[photonOverSmearing::EElowEtaBad     ] = sqrt(0.17*0.17+0.30*0.30);
    _oversmearing_err[photonOverSmearing::EEhighEtaGold   ] = sqrt(0.15*0.15+0.34*0.34);
    _oversmearing_err[photonOverSmearing::EEhighEtaBad    ] = sqrt(0.18*0.18+0.52*0.52);

  }



  if( setupType == "oversmear_moriond2013" ) {
    _oversmearing[photonOverSmearing::EBlowEtaGold_GAP] = 1.11;
    _oversmearing[photonOverSmearing::EBlowEtaGold_CM ] = 1.11;
    _oversmearing[photonOverSmearing::EBlowEtaBad_GAP ] = 1.07;
    _oversmearing[photonOverSmearing::EBlowEtaBad_CM  ] = 1.07;
    _oversmearing[photonOverSmearing::EBhighEtaGold   ] = 1.55;
    _oversmearing[photonOverSmearing::EBhighEtaBad    ] = 1.94;
    
    _oversmearing[photonOverSmearing::EElowEtaGold    ] = 2.95;
    _oversmearing[photonOverSmearing::EElowEtaBad     ] = 2.76;
    _oversmearing[photonOverSmearing::EEhighEtaGold   ] = 3.70;    
    _oversmearing[photonOverSmearing::EEhighEtaBad    ] = 3.71;
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_GAP] = sqrt(0.07*0.07+0.22*0.22);
    _oversmearing_err[photonOverSmearing::EBlowEtaGold_CM ] = sqrt(0.07*0.07+0.22*0.22);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_GAP ] = sqrt(0.06*0.06+0.24*0.24);
    _oversmearing_err[photonOverSmearing::EBlowEtaBad_CM  ] = sqrt(0.06*0.06+0.24*0.24);
    _oversmearing_err[photonOverSmearing::EBhighEtaGold   ] = sqrt(0.40*0.40+0.60*0.60);
    _oversmearing_err[photonOverSmearing::EBhighEtaBad    ] = sqrt(0.11*0.11+0.59*0.59);
    _oversmearing_err[photonOverSmearing::EElowEtaGold    ] = sqrt(0.25*0.25+0.90*0.90);
    _oversmearing_err[photonOverSmearing::EElowEtaBad     ] = sqrt(0.13*0.13+0.30*0.30);
    _oversmearing_err[photonOverSmearing::EEhighEtaGold   ] = sqrt(0.11*0.11+0.34*0.34);
    _oversmearing_err[photonOverSmearing::EEhighEtaBad    ] = sqrt(0.16*0.16+0.52*0.52);

  }




  

  /// absolute:
  for( unsigned ios = 0 ; ios < _oversmearing.size()    ; ios++ ) _oversmearing[ios]     /= 100.;
  for( unsigned ios = 0 ; ios < _oversmearing_err.size(); ios++ ) _oversmearing_err[ios] /= 100.;

}




photonOverSmearing::category photonOverSmearing::photonOverSmearingCategory( float scEta, float r9, bool isEBGap  ) {
  scEta = fabs(scEta);
  if     ( scEta < 1.0 ) {
    if( r9 > 0.94 ) {
      if( isEBGap ) return photonOverSmearing::EBlowEtaGold_GAP;
      else          return photonOverSmearing::EBlowEtaGold_CM;
    } else {
      if( isEBGap ) return photonOverSmearing::EBlowEtaBad_GAP;
      else          return photonOverSmearing::EBlowEtaBad_CM;
    }
  } else if( scEta < 1.45 )  {
    if( r9 > 0.94 ) return photonOverSmearing::EBhighEtaGold;
    else            return photonOverSmearing::EBhighEtaBad;
  } else if( scEta < 2.0 ) {
    if( r9 > 0.94 ) return photonOverSmearing::EElowEtaGold;
    else            return photonOverSmearing::EElowEtaBad;
  } else if( scEta < 2.5 ) {
    if( r9 > 0.94 ) return photonOverSmearing::EEhighEtaGold;
    else            return photonOverSmearing::EEhighEtaBad;
  }

  return  photonOverSmearing::UNKNOWN;
}


float photonOverSmearing::meanOverSmearing( float scEta, float r9, bool isEBGap, int syst) {
  photonOverSmearing::category cat = photonOverSmearingCategory(scEta, r9, isEBGap);
  if( cat == photonOverSmearing::UNKNOWN ) return 0;
  return _oversmearing[cat] + syst * _oversmearing_err[cat];
}


float photonOverSmearing::randOverSmearing( float scEta, float r9,  bool isEBGap, int syst) {
  float oversmearing = meanOverSmearing( scEta, r9, isEBGap, syst );
  return overSmearingGen.Gaus( 0, oversmearing);
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// FC: DiPhotonMassResolution estimator //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MassResolution {
 public:
  MassResolution(void);
  void addSmearingTerm( bool addSmearing = true );
  void setP4CalPosVtxResoSmear(const TLorentzVector & p1, const TLorentzVector &p2,
			       const TVector3   &calpos1, const TVector3  &calpos2,
			       const TVector3 &vtx,
			       float relReso[2], float relSmear[2] );


  /// vtxType:  
  /// - -1: no angular term
  /// - 0: use correct vertex dZ RMS
  /// - 1: use random vertex dZ RMS
  double angularTerm(int vtxType = 0) const ;  
  double relativeMassResolution( int vtxType ) const ;
  double relativeMassResolutionCorrVertex( ) const ;
  double relativeMassResolutionWrongVertex( ) const ;

  /// FC calc.
  double relativeMassResolutionFab_total( double pCorrectVtx )  const ;
  double relativeMassResolutionFab_angular( void ) const;
  double relativeMassResolutionFab_energy( void ) const;

 private:
  TLorentzVector _p4[2];
  TVector3 _cal[2];
  TVector3 _vtx;
  double _relResE[2];
  double _relSmearE[2];
  double _vtxSigZ[2];
  
  bool  _addSmearingTerm;
  
  double SecH(double x) const ;
  double TanH(double x) const ;

};


MassResolution::MassResolution(void) {
  _addSmearingTerm = false;
  _vtxSigZ[0] = 0.10;        /// correct vertex
  _vtxSigZ[1] = sqrt(2)*5.8; /// random  vertex
  _vtxSigZ[1] = sqrt(2)*5.0; /// random  vertex (update 15Oct2012)
  
}

void  MassResolution::addSmearingTerm( bool addSmearing  ) {
  _addSmearingTerm = addSmearing;
}

void MassResolution::setP4CalPosVtxResoSmear(const TLorentzVector & p1, const TLorentzVector &p2,
					     const TVector3   &calpos1,	const TVector3  &calpos2,
					     const TVector3 &vtx,
					     float relReso[2], float relSmear[2] ) {

  _vtx = vtx;
  for( int ie = 0 ; ie < 2 ;ie++ ) {
    _relResE[ie] = relReso[ie];
    _relSmearE[ie] = relSmear[ie];
    if( ie == 0 ) {
      _cal[0] = calpos1;
      _p4[0]  = p1;
    } else {
      _cal[1] = calpos2;
      _p4[1]  = p2;
    }       
  }
}

double MassResolution::angularTerm(int vtxType)  const  {
  if( vtxType < 0 || vtxType > 1 ) return 0;
  double sigma_dz = _vtxSigZ[vtxType];

  TVector3 epos[2];
  double r[2];
  for( int ie = 0 ; ie < 2; ie++ ) {
    epos[ie] = _cal[ie] - _vtx;
    r[ie]    = epos[ie].Mag(); 
  } 
  
  /// take this from gglobe
  /// need to cross check the actual formula
  double cos_term = TMath::Cos(epos[0].Phi()-epos[1].Phi());
  double sech1 = SecH(epos[0].Eta());
  double sech2 = SecH(epos[1].Eta());
  double tanh1 = TanH(epos[0].Eta());
  double tanh2 = TanH(epos[1].Eta());

  double numerator1 = sech1*(sech1*tanh2-tanh1*sech2*cos_term);
  double numerator2 = sech2*(sech2*tanh1-tanh2*sech1*cos_term);
  double denominator = 1. - tanh1*tanh2 - sech1*sech2*cos_term;

  double ResTerm = (-1.*sigma_dz/denominator)*(numerator1/r[0] + numerator2/r[1]);
  return ResTerm;
}


// utility functions
double MassResolution::SecH(double x)  const {
  return 1.0/TMath::CosH(x);
}

double MassResolution::TanH(double x)  const {
  return TMath::TanH(x);
}
 

double MassResolution::relativeMassResolution( int vtxType )  const {

  double massResolution2 = 0;

  //// energy term
  for( int ie = 0 ; ie < 2; ie++ ) {
    massResolution2 += _relResE[ie]*_relResE[ie];
    if( _addSmearingTerm ) 
          massResolution2 += _relSmearE[ie]*_relSmearE[ie];
  }
  /// angular term
  double alpha = _p4[0].Angle( _p4[1].Vect() );
  double angular_derivative = TMath::Sin(alpha)/(1.-TMath::Cos(alpha));
  /// FC: I don't understand why the derivative is set to 1 for wrong vertex but follow gglobe...
  if( vtxType == 1 ) angular_derivative = 1;
  double sigma_dalpha = angularTerm( vtxType ) * angular_derivative;

  massResolution2 += (sigma_dalpha*sigma_dalpha);

  return 0.5*sqrt( massResolution2 );
  
}

double MassResolution::relativeMassResolutionCorrVertex(void)  const {
  return relativeMassResolution(0);
}

double MassResolution::relativeMassResolutionWrongVertex(void)  const {
  return relativeMassResolution(1);
}


double MassResolution::relativeMassResolutionFab_total( double pCorrectVtx ) const {
  

  double sigma2_PVz = pCorrectVtx*_vtxSigZ[0]*_vtxSigZ[0] + (1-pCorrectVtx)*_vtxSigZ[1]*_vtxSigZ[1];
  
  /// M = sqrt ( 2 * E1/chEta1 * E2/chEta2 * (cosh Deta - cos Dphi ) )
  /// (sigmaM / M)^2 = 1/4 * { (sigmaE1 / E1)^2 + (sigmaE2 / E2)^2
  ///                          + ( sinhDeta / (cosh Deta - cos Dphi) * (dEta1/dPVz - dEta2/dPVz )
  ///                              -tanhEta1 *dEta1/dPVz -  tanhEta2 *dEta2/dPVz )^2 * sigma_PVz^2  }
  ///  dEta / dPVz = -1/ sqrt( Rcalo^2 + (zSC- PVz)^2 )

  //// energy term
  double energyTerm2 = 0;
  for( int ie = 0 ; ie < 2; ie++ ) {
    energyTerm2 += _relResE[ie]*_relResE[ie];
    if( _addSmearingTerm ) 
          energyTerm2 += _relSmearE[ie]*_relSmearE[ie];
  }

  /// add angular term
  TVector3 epos[2];
  for( int ie = 0 ; ie < 2; ie++ ) 
    epos[ie] = _cal[ie] - _vtx;
  
  float Deta = epos[0].Eta() - epos[1].Eta();
  float cos_dphi = cos( epos[0].Phi() - epos[1].Phi() ); 
  float deriv = 
    + sinh( Deta ) /(cosh(Deta)-cos_dphi) * (-1./epos[0].Mag() + 1./epos[1].Mag()) 
    - tanh( epos[0].Eta() )* (-1./epos[0].Mag() )
    - tanh( epos[1].Eta() )* (-1./epos[1].Mag() );
  
  double angularTerm2 = deriv*deriv*sigma2_PVz;
  
  return 0.5 *sqrt( energyTerm2 +  angularTerm2);

}

double MassResolution::relativeMassResolutionFab_angular( void ) const {
  
  double sigma_PVz = _vtxSigZ[1]; /// only wrong vertex

  /// add angular term
  TVector3 epos[2];
  for( int ie = 0 ; ie < 2; ie++ ) 
    epos[ie] = _cal[ie] - _vtx;
  
  float Deta = epos[0].Eta() - epos[1].Eta();
  float cos_dphi = cos( epos[0].Phi() - epos[1].Phi() ); 
  float deriv = 
    + sinh( Deta ) /(cosh(Deta)-cos_dphi) * (-1./epos[0].Mag() + 1./epos[1].Mag()) 
    - tanh( epos[0].Eta() )* (-1./epos[0].Mag() )
    - tanh( epos[1].Eta() )* (-1./epos[1].Mag() );
  
  
  return 0.5 *fabs(deriv)*sigma_PVz;

}

double MassResolution::relativeMassResolutionFab_energy( void ) const {
  //// energy term
  double energyTerm2 = 0;
  for( int ie = 0 ; ie < 2; ie++ ) {
    energyTerm2 += _relResE[ie]*_relResE[ie];
    if( _addSmearingTerm ) 
          energyTerm2 += _relSmearE[ie]*_relSmearE[ie];
  }
  return 0.5 *sqrt( energyTerm2 );
}



#endif
