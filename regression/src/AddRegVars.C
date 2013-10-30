#include "TFile.h"
#include "TTree.h"
#include "TMVA/Reader.h"
#include "TLorentzVector.h"
using namespace TMVA;

void AddRegVars(int, char*, char*, int);

int main(int argc, char *argv[]){
  if(argc==1){
    printf("AddRegVars [optimIndex] [inputFileName] [outputFileName] [mass]\n");
    return 1;
  }

  int optimIndex=445;
  int mass=300;
  bool splitFile=1;
  char *newFileBase="";
  if(argc>1) optimIndex=atoi(argv[1]);
  if(argc>2) mass=atoi(argv[2]);
  if(argc>3) splitFile=atoi(argv[3]);
  if(argc>4) newFileBase=argv[4];

  char inputFile[160], inputFileEven[160], inputFileOdd[160];
  sprintf(inputFile,"/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_5_2_5/src/ggAnalysis/regression/TrainingFiles/tree_M%i_RD.root",mass);
  sprintf(inputFileEven,"/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_5_2_5/src/ggAnalysis/regression/TrainingFiles/tree_M%i_RD_even.root",mass);
  sprintf(inputFileOdd,"/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_5_2_5/src/ggAnalysis/regression/TrainingFiles/tree_M%i_RD_odd.root",mass);

  char outputFile[160], outputFileEven[160], outputFileOdd[160];
  sprintf(outputFile,"RegOutputFiles/tree_M%i_%i.root",mass,optimIndex);
  sprintf(outputFileEven,"RegOutputFiles/tree_M%i_%i_even.root",mass,optimIndex);
  sprintf(outputFileOdd,"RegOutputFiles/tree_M%i_%i_odd.root",mass,optimIndex);

  char newInputFile[160], newOutputFile[160];
  sprintf(newInputFile,"/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_5_2_5/src/ggAnalysis/regression/TrainingFiles/%s_regInput.root",newFileBase);
  sprintf(newOutputFile,"RegOutputFiles/%s_regOutput_%i.root",newFileBase,optimIndex);

  if(newFileBase[0]!='\0')
    AddRegVars(optimIndex, newInputFile, newOutputFile, mass);
  else if(splitFile){
    AddRegVars(optimIndex, inputFileEven, outputFileEven, mass);
    AddRegVars(optimIndex, inputFileOdd, outputFileOdd, mass);
  }
  else
    AddRegVars(optimIndex, inputFile, outputFile, mass);

  return 0;
}

void AddRegVars(int optimIndex, char* inputFileName, char* outputFileName, int mass){

  TFile *inFile = new TFile(inputFileName,"READ");
  TTree *inTree = (TTree*)inFile->Get("Events");
  TFile *outFile = new TFile(outputFileName,"RECREATE");
  TTree *outTree = new TTree("Events","Regression data with results");

  Long64_t evtNum;
  float ijet_genJetPt,ijet_pt,ijet_eta,ijet_phi,ijet_En,ijet_cef,ijet_nef,ijet_mef,ijet_emfrac,ijet_hadfrac,ijet_nconstituents,ijet_chf,ijet_JECUnc,ijet_ptLeadTrack,ijet_vtxPt,ijet_vtx3dL,ijet_SoftLeptPtCut,ijet_dPhiMETJet,ijet_nch,ijet_vtx3deL,ijet_vtxMass,ijet_ptRaw,ijet_EnRaw,ijet_SoftLeptptRelCut,ijet_SoftLeptdRCut,ijet_partonID,ijet_dRJetGenJet;
  float MET, rho25;

  inTree->SetBranchAddress("evtNum", &evtNum);
  inTree->SetBranchAddress("jet_genJetPt",     &ijet_genJetPt          );
  inTree->SetBranchAddress("jet_pt",           &ijet_pt                );
  inTree->SetBranchAddress("jet_eta",		 &ijet_eta		);     
  inTree->SetBranchAddress("jet_phi",		 &ijet_phi		);     
  inTree->SetBranchAddress("jet_En",		 &ijet_En		);     
  inTree->SetBranchAddress("jet_cef",		 &ijet_cef		);     
  inTree->SetBranchAddress("jet_nef",		 &ijet_nef		);     
  inTree->SetBranchAddress("jet_mef",		 &ijet_mef		);     
  inTree->SetBranchAddress("jet_emfrac",	 &ijet_emfrac		);     
  inTree->SetBranchAddress("jet_hadfrac",	 &ijet_hadfrac		);     
  inTree->SetBranchAddress("jet_nconstituents",&ijet_nconstituents	);     
  inTree->SetBranchAddress("jet_chf",		 &ijet_chf		);     
  inTree->SetBranchAddress("jet_JECUnc",	 &ijet_JECUnc		);     
  inTree->SetBranchAddress("jet_ptLeadTrack",	 &ijet_ptLeadTrack	);     
  inTree->SetBranchAddress("jet_vtxPt",	 &ijet_vtxPt		);     
  inTree->SetBranchAddress("jet_vtx3dL",	 &ijet_vtx3dL		);     
  inTree->SetBranchAddress("jet_SoftLeptPtCut",&ijet_SoftLeptPtCut	);     
  inTree->SetBranchAddress("MET",		 &MET		   	);
  inTree->SetBranchAddress("jet_dPhiMETJet",	 &ijet_dPhiMETJet	);     
  inTree->SetBranchAddress("jet_nch",		 &ijet_nch		);     
  inTree->SetBranchAddress("jet_vtx3deL",	 &ijet_vtx3deL	   	);
  inTree->SetBranchAddress("jet_vtxMass",	 &ijet_vtxMass	   	);
  inTree->SetBranchAddress("jet_ptRaw",	 &ijet_ptRaw		);     
  inTree->SetBranchAddress("jet_EnRaw",	 &ijet_EnRaw		);     
  inTree->SetBranchAddress("jet_SoftLeptptRelCut",&ijet_SoftLeptptRelCut  );
  inTree->SetBranchAddress("jet_SoftLeptdRCut",&ijet_SoftLeptdRCut	);     
  inTree->SetBranchAddress("rho25",		 &rho25		     	);
  inTree->SetBranchAddress("jet_partonID",	 &ijet_partonID	   	);
  inTree->SetBranchAddress("jet_dRJetGenJet",  &ijet_dRJetGenJet      );

  const int nhJet = 2;
  int nhJet_=nhJet, hJet_hasSoftLept[nhJet];
  float hJet_pt[nhJet], hJet_eta[nhJet], hJet_phi[nhJet], hJet_En[nhJet], hJet_genJetPt[nhJet], hJet_dRJetGenJet[nhJet];
  float Hmass, HmassCorr, HmassGen,  hJet_ptCorr[nhJet];

  outTree->Branch("nhJet", &nhJet_, "nhJet/I");
  outTree->Branch("hJet_pt", hJet_pt, "hJet_pt[nhJet]/F");
  outTree->Branch("hJet_eta", hJet_eta, "hJet_eta[nhJet]/F");
  outTree->Branch("hJet_genJetPt", hJet_genJetPt, "hJet_genJetPt[nhJet]/F");
  outTree->Branch("hJet_ptCorr", hJet_ptCorr, "hJet_ptCorr[nhJet]/F");
  outTree->Branch("hJet_hasSoftLept", hJet_hasSoftLept, "hJet_hasSoftLept[nhJet]/I");
  outTree->Branch("hJet_dRJetGenJet", hJet_dRJetGenJet, "hJet_dRJetGenJet[nhJet]/F");
  outTree->Branch("Hmass", &Hmass, "Hmass/F");
  outTree->Branch("HmassCorr", &HmassCorr, "HmassCorr/F");

  float var1,var2,var3,var4,var5,var6,var7;
  Reader *readerRegres = new Reader( "!Color:!Silent" );

  readerRegres->AddVariable( "jet_eta", &ijet_eta);
  readerRegres->AddVariable( "jet_emfrac", &ijet_emfrac);
  readerRegres->AddVariable( "jet_hadfrac", &ijet_hadfrac);
  readerRegres->AddVariable( "jet_nconstituents", &ijet_nconstituents);
  readerRegres->AddVariable( "jet_vtx3dL", &ijet_vtx3dL);
  readerRegres->AddVariable( "MET", &MET);
  readerRegres->AddVariable( "jet_dPhiMETJet", &ijet_dPhiMETJet);
  readerRegres->BookMVA("BDTG0",TString::Format("weights/TMVARegression_%i_Cat0_BDTG.weights.xml",optimIndex));
  readerRegres->BookMVA("BDTG1",TString::Format("weights/TMVARegression_%i_Cat1_BDTG.weights.xml",optimIndex));

  unsigned nEntries = inTree->GetEntries();
  inTree->GetEntry(0);
  Long_t oldEventNum=evtNum-1;//different event numbers to begin.
  int iHiggsJet=0;
  for(unsigned i=0; i<nEntries; ++i){
    inTree->GetEntry(i);
    if(evtNum!=oldEventNum){//new event
      oldEventNum=evtNum;
      iHiggsJet=0;
      hJet_pt[0]=ijet_pt;
      hJet_eta[0]=ijet_eta;
      hJet_phi[0]=ijet_phi;
      hJet_En[0]=ijet_En;
      hJet_genJetPt[0]=ijet_genJetPt;
      hJet_hasSoftLept[0]=0;//for now
      hJet_dRJetGenJet[0]=ijet_dRJetGenJet;
      float corr_factor = hJet_pt[0] < 90 ? (float)(readerRegres->EvaluateRegression("BDTG0")[0]) : (float)(readerRegres->EvaluateRegression("BDTG1")[0]);
      hJet_ptCorr[0] = hJet_pt[0]*corr_factor;
      ++iHiggsJet;
      continue;
    }
    else if(iHiggsJet==2)//3 or more jets
      continue;
    else if(iHiggsJet==1){//2nd of 2 jets. Fill tree.
      hJet_pt[1]=ijet_pt;
      hJet_eta[1]=ijet_eta;
      hJet_phi[1]=ijet_phi;
      hJet_En[1]=ijet_En;
      hJet_genJetPt[1]=ijet_genJetPt;
      hJet_hasSoftLept[1]=0;//for now
      hJet_dRJetGenJet[1]=ijet_dRJetGenJet;
      float corr_factor = hJet_pt[1] < 90 ? (float)(readerRegres->EvaluateRegression("BDTG0")[0]) : (float)(readerRegres->EvaluateRegression("BDTG1")[0]);
      hJet_ptCorr[1] = hJet_pt[1]*corr_factor;
      ++iHiggsJet;
    }

    TLorentzVector hJet41,hJet42;
    hJet41.SetPtEtaPhiE(hJet_pt[0],hJet_eta[0],hJet_phi[0],hJet_En[0]);
    hJet42.SetPtEtaPhiE(hJet_pt[1],hJet_eta[1],hJet_phi[1],hJet_En[1]);
    Hmass = (hJet41+hJet42).M();

    hJet41.SetPtEtaPhiE(hJet_ptCorr[0],hJet_eta[0],hJet_phi[0],hJet_En[0]*hJet_ptCorr[0]/hJet_pt[0]);
    hJet42.SetPtEtaPhiE(hJet_ptCorr[1],hJet_eta[1],hJet_phi[1],hJet_En[1]*hJet_ptCorr[1]/hJet_pt[1]);
    HmassCorr = (hJet41+hJet42).M();

    outTree->Fill();
  }
  outTree->Write();
  outFile->Close();
  inFile->Close();
}
