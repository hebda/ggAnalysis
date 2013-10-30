#ifndef configreader_HH_
#define configreader_HH_

#include <TSystem.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

class MiniTree;

using namespace std;

void removeSpaces(string& str);

class ConfigReader {
 public:
  inline ConfigReader( string configName );
  inline ~ConfigReader(void);

  /// accessors

  inline bool   isLepTag(void) const {return _runLepTag;}
  inline bool   isOnlyLepTag(void) const {return _runOnlyLepTag;}

  inline bool   isEleGamma(void) const {return _runEleGamma;}
  inline bool   isGammaEle(void) const {return _runGammaEle;}


  inline bool   isTest(void) const {return _doTest;}
  inline string outDir(void) const {return _outDir;}
  inline string analysisType(void) const {return _analysisType;}
  inline string configName(void)   const {return _configName;}
  inline string energyScaleFile(void) const {return _energyScaleCorrectionFile; }
  inline string setup(void) const {return _setup;}
  inline bool   invertElectronVeto(void) const {return _invertElectronVeto;}
  inline bool   noElectronVeto(void) const {return _noElectronVeto;}
  inline bool   doSkimming(void) const {return _doSkimming;}
  inline bool   doVBFmvaCat(void) const {return _vbfMvaCat; }
  inline string storageElementDir(void) const {return _se;}
  inline float  pt1Cut(void) const {return _pt1Cut;}
  inline float  pt2Cut(void) const {return _pt2Cut;}
  inline float  loosePt1MCut() const { return _loosePt1MCut; }
  inline float  loosePt2MCut() const { return _loosePt2MCut; }

  inline float  mggCut(void) const {return _mggCut;}
  inline float  photonIdMvaCut(void) const { return _phoIdMvaCut; }

  inline float  luminosity(void) const {return _luminosity;}

  inline int    nEvtsPerJob(void) const {return _nEvtsPerJob;}
  inline int    getEnergyOverSmearingSyst(void) const {return _systEnergyResolution;}
  inline int    doPUJetID(void) const {return _doPUJetID;}

  //  inline string setupType(void) const {return _setupType;}
  /// dijetTree cuts
  inline float  diJetPtg1M(void) const { return _diJetPtg1M; }
  inline float  diJetPtg2M(void) const { return _diJetPtg2M; }
  inline float  diJetMjj(void)   const { return _diJetMjj;   }

  inline int  vtxSelType(void)   const { return _vtxSelType;   }
  inline int  getDoJetRegression(void)   const { return _doJetRegression;   }
  inline int  getDoControlSample(void)   const { return _doControlSample;   }
  

  inline ConfigReader& operator=( const ConfigReader &config );




 private:
  MiniTree *_minitree;
  string _configName;


  bool _runLepTag;
  bool _runOnlyLepTag;
  bool _runEleGamma;
  bool _runGammaEle;

  /// test: - switch the outputdir to test_out/
  bool _doTest;

  /// VBF category, use Mva?
  bool _vbfMvaCat;

  /// analysis type: baseline, MVA
  /// this mostly changes the photonID selection
  string _analysisType;
  string _setup;

  /// EnergyScale file
  string _energyScaleCorrectionFile;
  
  /// output dir
  string _outDir;
  string _se; /// for skimming store on eos

  /// electron veto: Z analysis: invert the electron veto
  bool _invertElectronVeto;

  /// do not cut on electron veto, usefull if on needs to do skimming
  bool _noElectronVeto;
  bool _doSkimming;

  /// in case one want to cut on pt.
  float _pt1Cut;
  float _pt2Cut;
  float _mggCut;
  float _phoIdMvaCut;
  float _loosePt1MCut; 
  float _loosePt2MCut; 
  

  /// # of evts per job (allows to split one file in several jobs)
  int _nEvtsPerJob;

  /// luminosity
  float _luminosity;

  /// energy resolution over smearing systematic
  int _systEnergyResolution;

  /// Jet ID for 2012 high PU
  bool  _doPUJetID;

  /// dijet cuts to fill dijet tree
  float _diJetPtg1M,_diJetPtg2M,_diJetMjj;

  // method to select vtx. 0 for MVA, 1 for sumpt2
  int _vtxSelType;

  /// regression and control sample
  int _doJetRegression,_doControlSample;

};


ConfigReader::ConfigReader( string myConfigName ) {
  _configName = myConfigName;
   ifstream input(_configName.c_str());

   //// default config values
   /// analysis: baseline / MVA / ...
   _analysisType = "baseline";
   /// setup: ReReco2011 / prompt2012_ichep /...
   _setup = "ReReco2011";
   
   _diJetPtg1M = 1./3.;
   _diJetPtg2M = 1./4.;
   _diJetMjj   = 250;

   _outDir = "./output_test_baseline/";
   _se = "unknown";
   _runEleGamma = 0;
   _runGammaEle = 0;
   _runLepTag = 0;
   _runOnlyLepTag = 0;
   _doTest = 1;
   _invertElectronVeto = 0;
   _pt1Cut = 25;
   _pt2Cut = 25;
   _loosePt1MCut = 1./3.;
   _loosePt2MCut = 1./4.;
   _mggCut = 90;
   _nEvtsPerJob = 500000;
   _systEnergyResolution = 0;
   _doPUJetID = 0;
   _noElectronVeto = 0;

   _phoIdMvaCut = -0.2; /// 2012 setup (-0.3 in 2011).
   _doSkimming = false;
   _luminosity = 5089;
   _vbfMvaCat = false;
   _energyScaleCorrectionFile = "unknown";
   _vtxSelType = 0;
   _doJetRegression = 0;
   _doControlSample = 0;
   cout << " Opening config file: " << _configName<< endl;
   while ( input.good() && !input.eof() ) {
     string line;
    getline(input,line,'\n');
    
    /// comment
    if( line.find("#") != string::npos ) continue;
    
    int posSeparator = line.find(":");
    
    if( posSeparator < 0 ) continue;
    
    /// split the line in key : val
    string key = line.substr(0,posSeparator);
    string val = line.substr(posSeparator+1,line.size());
    removeSpaces(key);
    removeSpaces(val);
    
    if( key.find( "analysisType" ) != string::npos ) _analysisType = val; 

    if( key.find("pt1_cut" ) != string::npos ) _pt1Cut = atof(val.c_str());
    if( key.find("pt2_cut" ) != string::npos ) _pt2Cut = atof(val.c_str());
    if( key.find("loosePt1M_cut") != string::npos ) _loosePt1MCut = atof(val.c_str());
    if( key.find("loosePt2M_cut") != string::npos ) _loosePt2MCut = atof(val.c_str());

    if( key.find("mgg_cut" ) != string::npos ) _mggCut = atof(val.c_str());
    if( key.find("phoIdMva_cut" ) != string::npos ) _phoIdMvaCut = atof(val.c_str());

    if( key.find("luminosity" ) != string::npos ) _luminosity = atof(val.c_str());

    if( key.find("doTest") != string::npos ) _doTest = atoi(val.c_str());
    if( key.find("runLepTag") != string::npos ) _runLepTag = atoi(val.c_str());
    if( key.find("runOnlyLepTag") != string::npos ) _runOnlyLepTag = atoi(val.c_str());

    if( key.find("outputDir") != string::npos ||
	key.find("outDir") != string::npos ||
	key.find("outdir") != string::npos ) _outDir = val;

   if( key.find("invertElectronVeto" ) != string::npos ||
       key.find("invertEleVeto")       != string::npos ||
       key.find("invEleVeto")          != string::npos ) _invertElectronVeto = atoi(val.c_str());

   if( key.find("noElectronVeto" ) != string::npos ) _noElectronVeto = atoi(val.c_str());
   if( key.find("doSkimming"     ) != string::npos ) _doSkimming = atoi(val.c_str());
   if( key.find("storageElement" ) != string::npos ) _se = val;

   if( key.find("VBFmvaCategory" ) != string::npos ) _vbfMvaCat = atoi(val.c_str());

   if( key.find("nEvtsPerJob") != string::npos ) _nEvtsPerJob = atoi(val.c_str());
   if( key.find("systEnergyResolution") != string::npos ) _systEnergyResolution = atoi(val.c_str());
   if( key.find("doPUJetID") != string::npos ) _doPUJetID = atoi(val.c_str());

   if( key.find("runEleGamma") != string::npos ) _runEleGamma = atoi(val.c_str());
   if( key.find("runGammaEle") != string::npos ) _runGammaEle = atoi(val.c_str());

   if( key.find("energyScaleCorrection") != string::npos ) _energyScaleCorrectionFile = val;

   if( key.find("diJetPtg1M_cut") != string::npos ) _diJetPtg1M = atof(val.c_str());
   if( key.find("diJetPtg2M_cut") != string::npos ) _diJetPtg2M = atof(val.c_str());
   if( key.find("diJetMjj_cut"  ) != string::npos ) _diJetMjj   = atof(val.c_str());

   if( key.find("vtxSelType") != string::npos ) _vtxSelType = atoi(val.c_str());
   if( key.find("doJetRegression") != string::npos ) _doJetRegression = atoi(val.c_str());
   if( key.find("doControlSample") != string::npos ) _doControlSample = atoi(val.c_str());

   if( key.find("setup") != string::npos || key.find("Setup") != string::npos ) {
     /// first verify the setup is allowed
     if( ! ( val == "ReReco2011" || val == "Prompt2012_ichep")  ) throw std::runtime_error("setup can only be: \n" +
											   string("    - ReReco2011 \n"      ) +
											   string("    - Prompt2012_ichep \n") );     
     _setup = val;
   }  
   }
   
   if( _diJetPtg1M < _loosePt1MCut ) _loosePt1MCut = _diJetPtg1M;
   if( _diJetPtg2M < _loosePt2MCut ) _loosePt2MCut = _diJetPtg2M;
	
   
   if( _doTest ) {
     _outDir = "./output_test/";
     _outDir = "./output_test_baseline/";
   }
   string mkdir = "mkdir -p " + _outDir;
   gSystem->Exec( mkdir.c_str() );
   
   ///// summary of the config file content...
   cout << " ============= Config file Content ================" << endl;
   cout << " -------- config file   : " << _configName << endl;
   cout << "    - is lepton tag?    : " << _runLepTag << endl; 
   cout << "    - is only lepton tag?    : " << _runOnlyLepTag << endl; 
   cout << "    - analysis type     : " << _analysisType << endl; 
   cout << "    - setup             : " << _setup << endl; 
   cout << "    - output directory  : " << _outDir << endl; 
   cout << "    - invert elec veto  : " << _invertElectronVeto << endl; 
   cout << "    - no elec veto      : " << _noElectronVeto << " (do not use for real analysis)" <<endl; 
   cout << "    - perform test      : " << _doTest << endl; 
   cout << "    - pt1 cut           : " << _pt1Cut << endl;
   cout << "    - pt2 cut           : " << _pt2Cut << endl;
   cout << "    - loose pt1/m cut   : " << _loosePt1MCut << endl;
   cout << "    - loose pt2/m cut   : " << _loosePt2MCut << endl;
   cout << "    - mgg cut           : " << _mggCut << endl;
   cout << "    - lumiosity         : " << _luminosity << endl;
   cout << "    - #evts per job     : " << _nEvtsPerJob << endl;
   cout << "    - oversmearing syst : " << _systEnergyResolution << endl;
   cout << "    - Jet PU ID cut     : " << _doPUJetID  << endl;
   cout << "    - DiJet Tree pt1/M  : " << _diJetPtg1M << endl;
   cout << "    - DiJet Tree pt2/M  : " << _diJetPtg2M << endl;
   cout << "    - DiJet Tree Mjj    : " << _diJetMjj   << endl;
   cout << "    - vtxSelType        : " << _vtxSelType << endl;
   cout << "    - doJetRegression   : " << _doJetRegression << endl;
   cout << "    - doControlSample   : " << _doControlSample << endl;
   cout << "    - skimming          : " << _doSkimming << endl;
   cout << "    - storage dir (eos) : " << _se << "(used for skimming only)" << endl;

   
   cout << " WARNING: if doTest ==> output dir: ./output_test/"  << endl;
   cout << " ==================================================" << endl;
   

   input.close();
   
}


ConfigReader& ConfigReader::operator=( const ConfigReader &c ) {
  _analysisType = c._analysisType;
  _outDir       = c._outDir;
  _doTest       = c._doTest;
  _pt1Cut       = c._pt1Cut;
  _pt2Cut       = c._pt2Cut;
  _loosePt1MCut = c._loosePt1MCut;
  _loosePt2MCut = c._loosePt2MCut;

  _nEvtsPerJob  = c._nEvtsPerJob;
  _invertElectronVeto = c._invertElectronVeto;
  _doPUJetID    = c._doPUJetID;

  _runEleGamma = c._runEleGamma;
  _runGammaEle = c._runGammaEle;
  _runOnlyLepTag = c._runOnlyLepTag; 

  _diJetPtg1M = c._diJetPtg1M;
  _diJetPtg2M = c._diJetPtg2M;
  _diJetMjj  = c._diJetMjj;
  _vtxSelType = c._vtxSelType;
  _doJetRegression = c._doJetRegression;
  _doControlSample = c._doControlSample;
  return *this;
}
  

ConfigReader::~ConfigReader(void) {}


void removeSpaces(string &str) {
  unsigned long ipos = -1;
  while( (ipos = str.find(" ")) != string::npos ) str.replace(ipos,1,"");
}

void testConfigFileReader(string configFile ) {
  ConfigReader config( configFile );

}

////////////////////// Higgs cross section file reader //////////////////////////
#include <map>
map<float,float> ReadHiggsCrossSecionFile(string file);
map<float,float> ReadHiggsBranchingRatioFile(string file);
class SMHiggsCrossSection {
 public:
  inline SMHiggsCrossSection( void );
  inline ~SMHiggsCrossSection(void);

  inline void is7TeV( void ) { _xsec_at_nTeV =  7; }
  inline void is8TeV( void ) { _xsec_at_nTeV =  8; }
  inline void is14TeV(void ) { _xsec_at_nTeV = 14; }

  inline float HiggsSMxsec_ggh( float mass );
  inline float HiggsSMxsec_vbf( float mass );
  inline float HiggsSMxsec_wh(  float mass );
  inline float HiggsSMxsec_zh(  float mass );
  inline float HiggsSMxsec_tth( float mass );
  inline float HiggsBR(  float mass );

 private:
  int     _xsec_at_nTeV;
  string _file_7TeV_ggh;
  string _file_8TeV_ggh;
  string _file_14TeV_ggh;
  string _file_7TeV_vbf;
  string _file_8TeV_vbf;
  string _file_14TeV_vbf;
  string _file_7TeV_wh;
  string _file_8TeV_wh;
  string _file_14TeV_wh;
  string _file_7TeV_zh;
  string _file_8TeV_zh;
  string _file_14TeV_zh;
  string _file_7TeV_tth;
  string _file_8TeV_tth;
  string _file_14TeV_tth;

  string _file_br;

  map<float,float> _xsec_ggh;
  map<float,float> _xsec_vbf;
  map<float,float> _xsec_wh;
  map<float,float> _xsec_zh;
  map<float,float> _xsec_tth;
  map<float,float> _br;


};


SMHiggsCrossSection::SMHiggsCrossSection(void) {
  _file_7TeV_ggh  = "sigCrossSections/higgsCrossSection_7TeV_ggH.txt";
  _file_8TeV_ggh  = "sigCrossSections/higgsCrossSection_8TeV_ggH.txt";
  _file_14TeV_ggh = "sigCrossSections/higgsCrossSection_14TeV_ggH.txt";
  _file_7TeV_vbf  = "sigCrossSections/higgsCrossSection_7TeV_VBF.txt";
  _file_8TeV_vbf  = "sigCrossSections/higgsCrossSection_8TeV_VBF.txt";
  _file_14TeV_vbf = "sigCrossSections/higgsCrossSection_14TeV_VBF.txt";
  _file_7TeV_wh   = "sigCrossSections/higgsCrossSection_7TeV_WH.txt";
  _file_8TeV_wh   = "sigCrossSections/higgsCrossSection_8TeV_WH.txt";
  _file_14TeV_wh  = "sigCrossSections/higgsCrossSection_14TeV_WH.txt";
  _file_7TeV_zh   = "sigCrossSections/higgsCrossSection_7TeV_ZH.txt";
  _file_8TeV_zh   = "sigCrossSections/higgsCrossSection_8TeV_ZH.txt";
  _file_14TeV_zh  = "sigCrossSections/higgsCrossSection_14TeV_ZH.txt";
  _file_7TeV_tth  = "sigCrossSections/higgsCrossSection_7TeV_TTH.txt";
  _file_8TeV_tth  = "sigCrossSections/higgsCrossSection_8TeV_TTH.txt";
  _file_14TeV_tth = "sigCrossSections/higgsCrossSection_14TeV_TTH.txt";
  
  _file_br        = "sigCrossSections/higgsBranchingRatio.txt";

  _xsec_at_nTeV = 7;
}
SMHiggsCrossSection::~SMHiggsCrossSection(void) {};

float SMHiggsCrossSection::HiggsBR( float mass ) {
  if( _br.size() == 0 ) {
    string filein = _file_br;
    _br = ReadHiggsBranchingRatioFile( filein );
  }

  /// now the map should be correctly filled
  if( _br.size() == 0 ) {
    cout << " ========== Higgs SM branching ratios are not filled properly... =========" << endl;
    return -1;
  }

  map<float,float>::iterator massit = _br.find( mass );
  if( massit == _br.end() ) {
    cout << "  --- did not find mass: " << mass << " for branching ratio (return +1)" << endl;
    return +1;
  }

  return massit->second;
}



float SMHiggsCrossSection::HiggsSMxsec_ggh( float mass ) {
  if( _xsec_ggh.size() == 0 ) {
    string filein = "unknown";
    if( _xsec_at_nTeV ==  7 ) filein = _file_7TeV_ggh;
    if( _xsec_at_nTeV ==  8 ) filein = _file_8TeV_ggh;
    if( _xsec_at_nTeV == 14 ) filein = _file_14TeV_ggh;
    _xsec_ggh = ReadHiggsCrossSecionFile( filein );
  }

  /// now the map should be correctly filled
  if( _xsec_ggh.size() == 0 ) {
    cout << " ========== Higgs SM xsec are not filled properly... =========" << endl;
    return -1;
  }

  map<float,float>::iterator massit = _xsec_ggh.find( mass );
  if( massit == _xsec_ggh.end() ) {
    cout << "  --- did not find mass: " << mass << " for process ggh (return -1)" << endl;
    return -1;
  }

  return massit->second;
}



float SMHiggsCrossSection::HiggsSMxsec_vbf( float mass ) {
  if( _xsec_vbf.size() == 0 ) {
    string filein = "unknown";
    if( _xsec_at_nTeV ==  7 ) filein = _file_7TeV_vbf;
    if( _xsec_at_nTeV ==  8 ) filein = _file_8TeV_vbf;
    if( _xsec_at_nTeV == 14 ) filein = _file_14TeV_vbf;
    _xsec_vbf = ReadHiggsCrossSecionFile( filein );
  }

  /// now the map should be correctly filled
  if( _xsec_vbf.size() == 0 ) {
    cout << " ========== Higgs SM xsec are not filled properly... =========" << endl;
    return -1;
  }

  map<float,float>::iterator massit = _xsec_vbf.find( mass );
  if( massit == _xsec_vbf.end() ) {
    cout << "  --- did not find mass: " << mass << " for process vbf (return -1)" << endl;
    return -1;
  }

  return massit->second;
}


float SMHiggsCrossSection::HiggsSMxsec_wh( float mass ) {
  if( _xsec_wh.size() == 0 ) {
    string filein = "unknown";
    if( _xsec_at_nTeV ==  7 ) filein = _file_7TeV_wh;
    if( _xsec_at_nTeV ==  8 ) filein = _file_8TeV_wh;
    if( _xsec_at_nTeV == 14 ) filein = _file_14TeV_wh;
    _xsec_wh = ReadHiggsCrossSecionFile( filein );
  }

  /// now the map should be correctly filled
  if( _xsec_wh.size() == 0 ) {
    cout << " ========== Higgs SM xsec are not filled properly... =========" << endl;
    return -1;
  }

  map<float,float>::iterator massit = _xsec_wh.find( mass );
  if( massit == _xsec_wh.end() ) {
    cout << "  --- did not find mass: " << mass << " for process wh (return -1)" << endl;
    return -1;
  }

  return massit->second;
}


float SMHiggsCrossSection::HiggsSMxsec_zh( float mass ) {
  if( _xsec_zh.size() == 0 ) {
    string filein = "unknown";
    if( _xsec_at_nTeV ==  7 ) filein = _file_7TeV_zh;
    if( _xsec_at_nTeV ==  8 ) filein = _file_8TeV_zh;
    if( _xsec_at_nTeV == 14 ) filein = _file_14TeV_zh;
    _xsec_zh = ReadHiggsCrossSecionFile( filein );
  }

  /// now the map should be correctly filled
  if( _xsec_zh.size() == 0 ) {
    cout << " ========== Higgs SM xsec are not filled properly... =========" << endl;
    return -1;
  }

  map<float,float>::iterator massit = _xsec_zh.find( mass );
  if( massit == _xsec_zh.end() ) {
    cout << "  --- did not find mass: " << mass << " for process zh (return -1)" << endl;
    return -1;
  }

  return massit->second;
}



float SMHiggsCrossSection::HiggsSMxsec_tth( float mass ) {
  if( _xsec_tth.size() == 0 ) {
    string filein = "unknown";
    if( _xsec_at_nTeV ==  7 ) filein = _file_7TeV_tth;
    if( _xsec_at_nTeV ==  8 ) filein = _file_8TeV_tth;
    if( _xsec_at_nTeV == 14 ) filein = _file_14TeV_tth;
    _xsec_tth = ReadHiggsCrossSecionFile( filein );
  }

  /// now the map should be correctly filled
  if( _xsec_tth.size() == 0 ) {
    cout << " ========== Higgs SM xsec are not filled properly... =========" << endl;
    return -1;
  }

  map<float,float>::iterator massit = _xsec_tth.find( mass );
  if( massit == _xsec_tth.end() ) {
    cout << "  --- did not find mass: " << mass << " for process tth (return -1)" << endl;
    return -1;
  }

  return massit->second;
}


map<float,float> ReadHiggsCrossSecionFile(string file) {
  ifstream input(file.c_str());
  map<float,float> higgsCrossSection;

   cout << " --- read higgs cross section file: " << file << endl;
   if( !input.good()) {
     cout << "    ===> can not open file: " << file << endl;
     return higgsCrossSection;
   }

   while ( input.good() && !input.eof() ) {
     string line;
     getline(input,line,'\n');
    
     /// comment
     if( line.find("#") != string::npos ) continue;
     
     /// read line
     float mass, xsec;
     istringstream isstream(line);
     isstream >> mass >> xsec;
     higgsCrossSection.insert( make_pair( mass, xsec ) );
   }

   return higgsCrossSection;
}




map<float,float> ReadHiggsBranchingRatioFile(string file) {
  ifstream input(file.c_str());
  map<float,float> higgsBranchingRatio;

   cout << " --- read higgs cross section file: " << file << endl;
   if( !input.good()) {
     cout << "    ===> can not open file: " << file << endl;
     return higgsBranchingRatio;
   }

   while ( input.good() && !input.eof() ) {
     string line;
     getline(input,line,'\n');
    
     /// comment
     if( line.find("#") != string::npos ) continue;
     
     /// read line
     float mass, br_gluglu, br_gamgam;
     istringstream isstream(line);
     isstream >> mass >> br_gluglu >> br_gamgam;
     higgsBranchingRatio.insert( make_pair( mass, br_gamgam ) );
   }

   return higgsBranchingRatio;
}


void testXsecReader( int eLHC = 8) {
  SMHiggsCrossSection signalXsecReader;
  if( eLHC ==  7 )  signalXsecReader.is7TeV();
  if( eLHC ==  8 )  signalXsecReader.is8TeV();
  if( eLHC == 14 )  signalXsecReader.is14TeV();

  float mass[5] = {100 ,114.5, 120, 125, 140 };
  for( int im = 0; im < 5; im++ ) 
   cout << " - mass = " << mass[im] 
	<< "; ggh : " << signalXsecReader.HiggsSMxsec_ggh( mass[im] )
	<< "; vbf : " << signalXsecReader.HiggsSMxsec_vbf( mass[im] )
	<< "; wh  : " << signalXsecReader.HiggsSMxsec_wh( mass[im] )
	<< "; zh  : " << signalXsecReader.HiggsSMxsec_zh( mass[im] )
	<< "; tth : " << signalXsecReader.HiggsSMxsec_tth( mass[im] )
	<< "  ---> BR: " << signalXsecReader.HiggsBR( mass[im] )
	<< endl;

}



#endif
