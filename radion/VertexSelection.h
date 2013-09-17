#ifndef vertexSelection__hh__
#define vertexSelection__hh__

#include "xAna_allAna.h"
#include "SetupReader.h"

VertexAlgoParameters VertexAlgoParametersReader( string file );


int  xAna::matchPhotonToConversion(int lpho) {
  
  if(DoDebugEvent)  _xcheckTextFile << " ------- xAna::matchPhotonToConversion (pho=" << lpho << " ) --------" << endl; 
  
  TVector3 pho_xyz; 
  pho_xyz.SetXYZ(phoCaloPos[lpho][0],
		 phoCaloPos[lpho][1],
		 phoCaloPos[lpho][2]);
  
  float dRMin = 999.;
  int iMatch=-1;     
  
  //// defining a good conversion
  float probMin = 1e-6;
  for(int iconv=0; iconv<nConv; iconv++) {
    TVector3 qconv; 
    if     ( convNTracks[iconv] == 2 ) qconv.SetXYZ( convRefittedMomentum[iconv][0],
						     convRefittedMomentum[iconv][1],
						     convRefittedMomentum[iconv][2] );
    else if( convNTracks[iconv] == 1 ) qconv.SetXYZ( convPairMomentum[iconv][0],
						     convPairMomentum[iconv][1],
						     convPairMomentum[iconv][2] ); 

    TVector3 vtx_xyz( convVtx[iconv][0], convVtx[iconv][1], convVtx[iconv][2] );
    TVector3 qpho = pho_xyz - vtx_xyz;        
    float dR = qpho.DeltaR( qconv ); 
    if(DoDebugEvent)    {
      _xcheckTextFile  << " check conv[" << iconv << "]: pt = " << qconv.Pt() << endl
		       << "              : nTrk = " <<  convNTracks[iconv]<< endl
		       << "              : vtxValid? = " << convValidVtx[iconv] << endl
		       << "              : chi2 = " << convChi2Probability[iconv] << endl
		       << "              : z   = " << vtx_xyz.Z() << " vs PV: " << vtxbs[0][2] << endl
		       << "              : phoDR = " << dR << endl;	
    }

    if( qconv.Pt() < 10 ) continue;
    if( _vtxAlgoParams.useAllConversions == 1 && convNTracks[iconv] != 1 ) continue;
    if( _vtxAlgoParams.useAllConversions == 2 && convNTracks[iconv] != 2 ) continue;
    if( _vtxAlgoParams.useAllConversions == 3 && convNTracks[iconv] != 1 && convNTracks[iconv] != 2 ) continue;
    if( convNTracks[iconv] == 2 &&
	! ( convValidVtx[iconv] && convChi2Probability[iconv] > probMin ) ) continue;
    if(DoDebugEvent)  _xcheckTextFile  << "              : pass quality cut. " << endl;

    if ( dR < dRMin )  {
      dRMin=dR;
      iMatch=iconv;
    }
  }
  

  if ( dRMin < 0.1 )  {
    if(DoDebugEvent)    _xcheckTextFile  << "       =====>  matched conversion index " << iMatch << endl;
    return iMatch;
  } 
 
  return -1;
}


PhotonInfo xAna::fillPhotonInfo( int iPho, bool useAllConvs) {

    bool pho_isEB_1 = false;
    if(fabs(phoSCEta[iPho]) < 1.4442) pho_isEB_1 = true;

    int iConv1 = -1;
    if( _vtxAlgoParams.useAllConversions > 0 ) 
      iConv1 = matchPhotonToConversion(iPho);


    if(DoDebugEvent)  _xcheckTextFile << "iPho = " << iPho << "  iConv1 = " << iConv1 << endl;
	
    TVector3 qpho; qpho.SetXYZ( phoCaloPos[iPho][0], phoCaloPos[iPho][1], phoCaloPos[iPho][2] );
    TVector3 bs  ; bs  .SetXYZ( bspotPos[0], bspotPos[1], bspotPos[2] );


    if(iConv1>=0) {
      TVector3 vtxconv; vtxconv.SetXYZ(convVtx[iConv1][0],convVtx[iConv1][1],convVtx[iConv1][2]);
   
      if(DoDebugEvent)
	_xcheckTextFile << " fillPhotonInfo: conversion match, convZofPrimVtxFromTrks["<<iConv1<<"] = " 
			<<  convZofPrimVtxFromTrks[iConv1] << " : defVtx = " << vtxbs[0][2] << endl;
      
      
      PhotonInfo pho(iPho, qpho, bs, vtxconv,
		     convNTracks[iConv1] == 1 ? 
		     TVector3(convPairMomentum[iConv1][0], convPairMomentum[iConv1][1], convPairMomentum[iConv1][2]):
		     TVector3(convRefittedMomentum[iConv1][0], convRefittedMomentum[iConv1][1], convRefittedMomentum[iConv1][2]),
		     phoE[iPho],  pho_isEB_1,
		     convNTracks[iConv1], convValidVtx[iConv1], convChi2Probability[iConv1], convEoverP[iConv1] ); 
      
      return pho;      
    } 
    float phoConvV_X = 0;//phoConvVtx[iPho][0];
    float phoConvV_Y = 0;//phoConvVtx[iPho][1];
    float phoConvV_Z = 0;//phoConvVtx[iPho][2];
    
    /* PhotonInfo pho(iPho, qpho, bs, TVector3(phoConvV_X,phoConvV_Y,phoConvV_Z), */
    /* 		   TVector3(phoConvRefittedMomentum[iPho][0], phoConvRefittedMomentum[iPho][1], phoConvRefittedMomentum[iPho][2]), */
    /* 		   phoE[iPho], pho_isEB_1, */
    /* 		   phoConvNTrks[iPho]   , phoConvValidVtx[iPho], */
    /* 		   phoConvChi2Prob[iPho], phoConvEoverP[iPho] );  */
    
    PhotonInfo pho(iPho, qpho, bs, TVector3(phoConvV_X,phoConvV_Y,phoConvV_Z),
		   TVector3(phoConvRefittedMomentum[iPho][0], phoConvRefittedMomentum[iPho][1], phoConvRefittedMomentum[iPho][2]),
		   phoE[iPho], pho_isEB_1,
		   phoConvNTrks[iPho]   , 0,
		   -1, phoConvEoverP[iPho] ); 

    return pho;
}


vector<int>  xAna::getSelectedVertex(int iPho1, int iPho2, bool useAllConvs, bool useMva) {
  PhotonInfo pho1 = fillPhotonInfo(iPho1,useAllConvs);
  PhotonInfo pho2 = fillPhotonInfo(iPho2,useAllConvs);


  ///// FC: try to fill it differently since it does not give the same results as gglobe
  vector<float> vtx_x,vtx_y,vtx_z;
  vector<int>   vtx_ntks;
  vector<vector<unsigned short> > vtx_tkind;
  vector<vector<float> > vtx_tkweight;
  vector<float> tk_px, tk_py, tk_pz, tk_pterr;
  vector<int>   tk_quality;
  
  /// fill tracks
  for( int t = 0; t < nTrk; ++t) {
    tk_px.push_back( trkP[t][0] );
    tk_py.push_back( trkP[t][1] );
    tk_pz.push_back( trkP[t][2] );
    tk_quality.push_back( trkQuality[t] );
    tk_pterr  .push_back( trkPtErr[t] );
  }
  
  /// fix in case the track branch does not exist
  bool trkBranchExist = fChain->GetBranchStatus("nTrk");

  /// fill vtx
  vtx_tkind.resize(nVtxBS);
  vtx_tkweight.resize(nVtxBS);

  for( int v = 0; v < nVtxBS; ++v) {
    vtx_x.push_back( vtxbs[v][0] );
    vtx_y.push_back( vtxbs[v][1] );
    vtx_z.push_back( vtxbs[v][2] );

    unsigned tksize = 0;
    if( trkBranchExist ) tksize = (*vtxbsTkIndex)[v].size();


    vtx_ntks.push_back( tksize );

    /// thenindex should be already ok...    
    for( unsigned it = 0 ; it <  tksize; it++ ) {
      vtx_tkind[v].push_back((unsigned short)(*vtxbsTkIndex)[v][it]);
      vtx_tkweight[v].push_back( (*vtxbsTkWeight)[v][it]);
    }
  }
  
  ggVertexInfo vinfo(vtx_x, vtx_y,vtx_z,
		     vtx_ntks,
		     tk_px,tk_py,tk_pz,
		     tk_pterr, tk_quality,
		     vtx_tkind, vtx_tkweight
		     ) ;


  _vtxAna->clear();
  _vtxAna->analyze(vinfo, pho1, pho2);
  //int result = -1;
  vector<int> result;
  result.push_back(-1);
  if (useMva) { 
    vector<int> rankresults = _vtxAna->rank(*_perVtxReader, _perVtxMvaMethod);
    vertMVA  = _vtxAna->perEventMva(*_perEvtReader ,_perEvtMvaMethod,rankresults);
    vertProb = _vtxAna->vertexProbability(vertMVA,nVtxBS);
    if( ! trkBranchExist ) {
      vertProb = 0.7; 
      if( rankresults.size() > 0 ) rankresults[0];
      else  rankresults.push_back(0);
    }

    if( _minitree->addSyncVariables ) {
      for( unsigned iv = 0; iv < 3; iv++ ) {
	_minitree->vtxId[ iv] = iv < rankresults.size() ? rankresults[iv] : -1;
	_minitree->vtxMva[iv] = iv < rankresults.size() ?_vtxAna->mva(rankresults[iv]) : -2;
	_minitree->vtxDz[iv]  = iv < rankresults.size() ?_vtxAna->vertexz(rankresults[iv]) - _vtxAna->vertexz(rankresults[0]) : -999;
      }
      _minitree->vtxPtBal      = _vtxAna->ptbal(rankresults[0]);
      _minitree->vtxAsym       = _vtxAna->ptasym(rankresults[0]);
      _minitree->vtxLogSumPt2  = _vtxAna->logsumpt2(rankresults[0]);
      _minitree->vtxPullToConv = _vtxAna->pulltoconv(rankresults[0]);
      _minitree->vtxNConv      = _vtxAna->nconv(rankresults[0]);
 
      if( DoDebugEvent ) {
	_xcheckTextFile << "  vertex analysis   nConversion: " << _minitree->vtxNConv << endl;
	_xcheckTextFile << "  vertex analysis:  pullToConv: " << _minitree->vtxPullToConv << endl;
      }
   }

    return rankresults;
  }


  return result; 
}



//////////// vertex adaptator /////////////////
ggVertexInfo::ggVertexInfo( const vector<float> &vtx_x, const vector<float> &vtx_y, const vector<float> &vtx_z,
			    const vector<int> &vtx_ntks,
			    const vector<float> &tk_px, const vector<float> &tk_py, const vector<float> &tk_pz,
			    const vector<float> &tk_pterr, const vector<int> &tk_quality,
			    vector<vector<unsigned short> > & vtx_tkind,
			    vector<vector<float> >    & vtx_tkweight
			    ) 
  
 {
   _vtx_x = vtx_x;
   _vtx_y = vtx_y;
   _vtx_z = vtx_z;
   _vtx_ntks = vtx_ntks;
   
   _tk_px = tk_px;
   _tk_py = tk_py;
   _tk_pz = tk_pz;
   _tk_pterr   = tk_pterr;
   _tk_quality = tk_quality;
   
   _vtx_tkind    = vtx_tkind;
   _vtx_tkweight = vtx_tkweight;

   
  _tk_n  = int(tk_px.size());
  _vtx_n = int(vtx_x.size());
}


//////////// vertex selection algorithm parameters /////////////////
VertexAlgoParameters VertexAlgoParametersReader( string file ){
  VertexAlgoParameters params;
  SetupReader reader( file );
  bool btmp(false);
  int   itmp;
  float ftmp;
  string stmp;
  
  itmp = reader.getInt( "fixTkIndex"        , btmp); if( btmp ) params.fixTkIndex         = bool(itmp);
  itmp = reader.getInt( "rescaleTkPtByError", btmp); if( btmp ) params.rescaleTkPtByError = bool(itmp);
  itmp = reader.getInt( "highPurityOnly"    , btmp); if( btmp ) params.highPurityOnly     = bool(itmp);
  itmp = reader.getInt( "removeTracksInCone", btmp); if( btmp ) params.removeTracksInCone = bool(itmp);
  itmp = reader.getInt( "useAllConversions" , btmp); if( btmp ) params.useAllConversions  = itmp;
  ftmp = reader.getFloat( "trackCountThr", btmp); if( btmp ) params.trackCountThr = ftmp;
  ftmp = reader.getFloat( "maxD0Signif"  , btmp); if( btmp ) params.maxD0Signif = ftmp;
  ftmp = reader.getFloat( "maxDzSignif"  , btmp); if( btmp ) params.maxDzSignif = ftmp;
  ftmp = reader.getFloat( "coneSize"     , btmp); if( btmp ) params.coneSize = ftmp;
  ftmp = reader.getFloat( "sigma1Pix"    , btmp); if( btmp ) params.sigma1Pix = ftmp;
  ftmp = reader.getFloat( "sigma1Tib"    , btmp); if( btmp ) params.sigma1Tib = ftmp;
  ftmp = reader.getFloat( "sigma1Tob"    , btmp); if( btmp ) params.sigma1Tob = ftmp;
  ftmp = reader.getFloat( "sigma1PixFwd" , btmp); if( btmp ) params.sigma1PixFwd = ftmp;
  ftmp = reader.getFloat( "sigma1Tid"    , btmp); if( btmp ) params.sigma1Tid = ftmp;
  ftmp = reader.getFloat( "sigma1Tec"    , btmp); if( btmp ) params.sigma1Tec = ftmp;
  ftmp = reader.getFloat( "sigma2Pix"    , btmp); if( btmp ) params.sigma2Pix = ftmp;
  ftmp = reader.getFloat( "sigma2Tib"    , btmp); if( btmp ) params.sigma2Tib = ftmp;
  ftmp = reader.getFloat( "sigma2Tob"    , btmp); if( btmp ) params.sigma2Tob = ftmp;
  ftmp = reader.getFloat( "sigma2PixFwd" , btmp); if( btmp ) params.sigma2PixFwd = ftmp;
  ftmp = reader.getFloat( "sigma2Tid"    , btmp); if( btmp ) params.sigma2Tid = ftmp;
  ftmp = reader.getFloat( "sigma2Tec"    , btmp); if( btmp ) params.sigma2Tec = ftmp;
  ftmp = reader.getFloat( "singlelegsigma1Pix"   , btmp); if( btmp ) params.singlelegsigma1Pix = ftmp;
  ftmp = reader.getFloat( "singlelegsigma1Tib"   , btmp); if( btmp ) params.singlelegsigma1Tib = ftmp;
  ftmp = reader.getFloat( "singlelegsigma1Tob"   , btmp); if( btmp ) params.singlelegsigma1Tob = ftmp;
  ftmp = reader.getFloat( "singlelegsigma1PixFwd", btmp); if( btmp ) params.singlelegsigma1PixFwd = ftmp;
  ftmp = reader.getFloat( "singlelegsigma1Tid"   , btmp); if( btmp ) params.singlelegsigma1Tid = ftmp;
  ftmp = reader.getFloat( "singlelegsigma1Tec"   , btmp); if( btmp ) params.singlelegsigma1Tec = ftmp;
  ftmp = reader.getFloat( "singlelegsigma2Pix"   , btmp); if( btmp ) params.singlelegsigma2Pix = ftmp;
  ftmp = reader.getFloat( "singlelegsigma2Tib"   , btmp); if( btmp ) params.singlelegsigma2Tib = ftmp;
  ftmp = reader.getFloat( "singlelegsigma2Tob"   , btmp); if( btmp ) params.singlelegsigma2Tob = ftmp;
  ftmp = reader.getFloat( "singlelegsigma2PixFwd", btmp); if( btmp ) params.singlelegsigma2PixFwd = ftmp;
  ftmp = reader.getFloat( "singlelegsigma2Tid"   , btmp); if( btmp ) params.singlelegsigma2Tid = ftmp;
  ftmp = reader.getFloat( "singlelegsigma2Tec"   , btmp); if( btmp ) params.singlelegsigma2Tec = ftmp;

  stmp = reader.getString("vtxProbFormula", btmp); if( btmp ) params.vtxProbFormula = stmp;

  return params;
}



#endif
