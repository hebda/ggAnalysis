#include "xAna_allAna.h"

bool xAna::PassTriggerSelection( void ) {
  int runNum = run;
  Int_t passHLT = 1;
  
  if (mode_ == -1) { // 2012 data
    // if(runNum>=160404 && runNum <= 161176){      
    if       ( runNum >= 190456 && runNum <= 194269 ) {
      // HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v*  --> [0]
      // HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v* -- > 8 
      // HLTIndex[0]
      if(HLTIndex[0]==-1 || HLTIndex[8]==-1 )
	{cout <<" run "<<  runNum <<"  event "<< event   <<" HLTIndex[0]  "<<  HLTIndex[0] <<"     HLTIndex[8]  "<<  HLTIndex[8] <<endl;}
      else{	
	if ( HLT[HLTIndex[0]] == 0 &&  HLT[HLTIndex[8]] == 0)passHLT=0;
      }
    } else if( runNum >= 194270 && runNum < 203894 ) {
      //  HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v*  --> 0 
      //  HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v*  --> 7
      //  HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v* --> 8
      if(HLTIndex[0]==-1 || HLTIndex[7]==-1 ||  HLTIndex[8]==-1)
	{cout <<" run "<<  runNum <<"  event "<< event   <<" HLTIndex[0]  "<<  HLTIndex[0] <<
	    " HLTIndex[7]  "<<  HLTIndex[7] <<"  HLTIndex[8]  "<<  HLTIndex[8] <<endl;}
      else {	
	if ( HLT[HLTIndex[0]] == 0 && HLT[HLTIndex[7]] == 0 &&  HLT[HLTIndex[8]] == 0 ) passHLT=0;  }
     
    } else if( runNum >= 203894 ) {
      //  HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v*  --> dropped for run2012D
      //  HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v*  --> 7
      //  HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v* --> 8
      if( HLTIndex[7]==-1 ||  HLTIndex[8]== -1 ){ 
	cout <<" run "<<  runNum <<"  event "<< event   <<" HLTIndex[0]  "<<  HLTIndex[0] <<
	  " HLTIndex[7]  "<<  HLTIndex[7] <<"  HLTIndex[8]  "<<  HLTIndex[8] <<endl;
      } else {	
	if ( HLT[HLTIndex[7]] == 0 &&  HLT[HLTIndex[8]] == 0 ) passHLT=0;  
      }
    }
  
  } else {   // is MC
    passHLT=1;   //==================== NO HLT applied to MC ===============================
  }
  
  return passHLT;
  
}
