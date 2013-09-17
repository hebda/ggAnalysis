#ifndef energyscalereader_HH_
#define energyscalereader_HH_

#include <TSystem.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

class EnergyScaleReader {
 public:
  inline EnergyScaleReader( void );
  ~EnergyScaleReader(void);

  inline int EcalPart( float r9, float scEta );
  inline int EcalPart( string ecal );
  inline string EcalPartString( float r9, float scEta );

  inline float energyScale(    float r9, float scEta, int run, int systShift = 0 );

  inline int runRange( int run, int iEcal );

  inline bool setup( string file );

 private:

  bool _isSetup;
  vector<int>    _runStart[8];
  vector<int>    _runStop[8];
  vector<string> _ecalPart;
  vector<float>  _eScale[8];
  vector<float>  _err_eScale[8];

  int _nWarnings;
};
inline int EnergyScaleReader::EcalPart( string ecal ) {
  
  /// internal code convention
  if( ecal == "EBlowEtaBad"   ) return 0;
  if( ecal == "EBlowEtaGold"  ) return 1;
  if( ecal == "EBhighEtaBad"  ) return 2;
  if( ecal == "EBhighEtaGold" ) return 3;
  if( ecal == "EElowEtaBad"   ) return 4;
  if( ecal == "EElowEtaGold"  ) return 5;
  if( ecal == "EEhighEtaBad"  ) return 6;
  if( ecal == "EEhighEtaGold" ) return 7;

  /// read other shervin's table format
  if( ecal == "EB-absEta_0_1-bad"      ) return 0;
  if( ecal == "EB-absEta_0_1-gold"     ) return 1;
  if( ecal == "EB-absEta_1_1.4442-bad" ) return 2;
  if( ecal == "EB-absEta_1_1.4442-gold") return 3;
  if( ecal == "EE-absEta_1.566_2-bad"  ) return 4;
  if( ecal == "EE-absEta_1.566_2-gold" ) return 5;
  if( ecal == "EE-absEta_2_2.5-bad"    ) return 6;
  if( ecal == "EE-absEta_2_2.5-gold"   ) return 7;
  
  return -1;
}

inline string EnergyScaleReader::EcalPartString( float r9, float scEta ) {
  string ecal = "Unknown";
  if     ( fabs(scEta) <= 1.00 ) ecal = "EBlowEta";
  else if( fabs(scEta) <= 1.45 ) ecal = "EBhighEta";
  else if( fabs(scEta) <= 2.00 ) ecal = "EElowEta";
  else if( fabs(scEta) <= 2.50 ) ecal = "EEhighEta";
  
  if( r9 > 0.94 ) ecal += "Gold";
  else            ecal += "Bad";
  
  return ecal;
}

inline int EnergyScaleReader::EcalPart(float r9, float scEta ) {
  return EcalPart( EcalPartString( r9, scEta ) );
}


bool EnergyScaleReader::setup( string energyscaleName ) {
   ifstream input(energyscaleName.c_str());

   //// default energyscale values


   cout << " Opening energy scale file: " << energyscaleName<< endl;
   while ( input.good() && !input.eof() ) {
     string line;
    getline(input,line,'\n');
    
    /// comment
    if( line.find("#") != string::npos ) continue;

    istringstream isstream1(line);
    istringstream isstream2(line);
    string ecal;
    int runMin,runMax;
    float eScale, err_eScale;
    isstream1 >> ecal >> runMin >> runMax >> eScale >> err_eScale ;

    if( runMin <= 0 ) {
      string dummy;
      /// most likely we are reading the new Format
      isstream2 >> ecal >> dummy >> runMin >> runMax >> eScale >> err_eScale ;
    }

    int iEcal = EcalPart(ecal);
    if( iEcal < 0 ) continue;
    _runStart[  iEcal].push_back( runMin );
    _runStop[   iEcal].push_back( runMax );
    _eScale[    iEcal].push_back( eScale );
    _err_eScale[iEcal].push_back( err_eScale );
   }
   _isSetup = true;
   return true;
}
EnergyScaleReader::EnergyScaleReader(void) {_isSetup = false; _nWarnings = 0;}
EnergyScaleReader::~EnergyScaleReader(void) {}


float EnergyScaleReader::energyScale(    float r9, float scEta, int run, int systShift ) {
  
  if( !_isSetup ) {
    cout << " EnergyScaleReader is not setup! use function setup with the ad hoc setup file" << endl;
    return 1;
  }

  int iEcal = EcalPart( r9, scEta );
  if( iEcal < 0 ) {
    cout << " EnergyScaleReader unknown Ecal part for: r9 = " << r9 << " eta = " << scEta << endl;
    return 1;
  }
    
  int irun = runRange( run, iEcal );

  if( irun < 0 ) {
    cout << " EnergyScaleReader: can not find run: " << run << " for ECAL: " << EcalPartString( r9, scEta ) << endl;
    return 1;
  }
  
  return _eScale[iEcal][irun] + systShift * _err_eScale[iEcal][irun];
}


inline int EnergyScaleReader::runRange( int run, int iEcal ) {
  for( unsigned irun = 0; irun < _runStart[iEcal].size(); irun++ )
    if( run >= _runStart[iEcal][irun] && run <= _runStop[iEcal][irun] ) return int(irun);


  if( run > _runStop[iEcal][_runStart[iEcal].size()-1] ) {
    if( _nWarnings < 20 ) {
      cout << "run " << run << " > last run known: " << _runStop[iEcal][_runStart[iEcal].size()-1] << " for iEcal " << iEcal << endl;
      _nWarnings++;
    }
    return _runStart[iEcal].size()-1;
  }

  return -1;
}


void testEnergyScale(int run,  string ecalCalibInput = "ecalCalibFiles/EnergyScale2012_prompt.txt") {

  
  EnergyScaleReader myEnergyScale;
  myEnergyScale.setup( ecalCalibInput );

  float r9[]  = { 0.8, 0.8, 0.83, 0.85, 0.98, 0.96, 0.99, 1.05, 0.7 };
  float eta[] = { 0.5, 1.4, 1.8, 2.4, 0.6, 1.2, 1.9, 2.23, 3.0 };
  for( int itest = 0 ; itest < 9; itest++ )
    cout << " eScale( r9 = " << r9[itest] << ", eta = "<< eta[itest] << ") = " << myEnergyScale.energyScale(r9[itest],eta[itest],run) << endl;
  
}
  

#endif
