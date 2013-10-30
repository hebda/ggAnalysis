
/*______________________________________________________________________________
 * Subtree producer.
 *
 * List of cuts is given below in the code.
 *
 * Usage: check MASS variable, check cuts, then run
 *    root -l mktree.cc+
 */

#include <TFile.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "../KinematicFit/DiJetKinFitter.cc"
#include "../KinematicFit/TKinFitter.cc"
#include "../KinematicFit/TFitParticleEtEtaPhi.cc"
#include "../KinematicFit/TAbsFitParticle.cc"
#include "../KinematicFit/TFitConstraintM.cc"
#include "../KinematicFit/TAbsFitConstraint.cc"

// radion mass to process (defines input and output root file names)
Int_t MASS = 1500;

// photon/jet pt cuts (in GeV)
Double_t minPtGamma = 20;
Double_t minPtJet = 25;

// jet eta cut
Double_t maxEtaJet = 2.5;

// jet b-tag cuts
const char* bJetTagsBranchName = "bJetTags";  // name of tree branch with bJetTag values
Double_t minBJetTag = 0.679;                  // b-jet tagging threshold value

// dR(gamma, gamma) and dR(gammas, jets) cuts
Double_t maxDeltaRGG = -1;  // -1 means unset
Double_t minDeltaRGJ = 0.4;  // -1 means unset

// allowed gamma+gamma mass region (in GeV)
// Double_t minMassGG = 125 - 1.5*3;
// Double_t maxMassGG = 125 + 1.5*3;
Double_t minMassGG = 0;
Double_t maxMassGG = 1e+10;

// allowed jet+jet mass region (in GeV)
// Double_t minMassJJ = 132 - 14.5*3;
// Double_t maxMassJJ = 132 + 14.5*3;
Double_t minMassJJ = 0;
Double_t maxMassJJ = 1e+10;

Bool_t includeWeightVars = kTRUE;
Bool_t doKinFit =kFALSE;
Bool_t includeTrainingVars = kFALSE;

//method for TArray class
void RemoveElement(TArrayF *arr, int index){
  for(int j=index+1; j<arr->GetSize(); ++j){
    arr->SetAt(arr->GetAt(j),j-1);
  }
  arr->Set(arr->GetSize()-1);
}

//______________________________________________________________________________
void mktree_one(const char* infile, const char* treename, int eventCut)
{
   // current working directory
   TDirectory* wd = gDirectory;

   // open input tree
   TFile* fi = TFile::Open(infile);
   TTree* tree = (TTree*) fi->Get("RToHH");

   // cd back into previous current working directory
   if (wd) wd->cd();

   // variables for input tree
   Long64_t evtNum;
   float wei;
   TLorentzVector *gamma1, *gamma2;
   TClonesArray* jets;
   TArrayF* bJetTags;
   TArrayF *jet_genJetPt;
   TArrayF *jet_eta;		   
   TArrayF *jet_cef;		   
   TArrayF *jet_nef;		   
   TArrayF *jet_mef;		   
   TArrayF *jet_nconstituents;	   
   TArrayF *jet_chf;		   
   TArrayF *jet_JECUnc;		   
   TArrayF *jet_ptLeadTrack;	   
   TArrayF *jet_vtxPt;		   
   TArrayF *jet_vtx3dL;		   
   TArrayF *jet_SoftLeptPtCut;	   
   float    MET;		   
   TArrayF *jet_dPhiMETJet;	   
   TArrayF *jet_nch;		   
   TArrayF *jet_vtx3deL;	   
   TArrayF *jet_vtxMass;	   
   TArrayF *jet_ptRaw;		   
   TArrayF *jet_EnRaw;		   
   TArrayF *jet_SoftLeptptRelCut;  
   TArrayF *jet_SoftLeptdRCut;	   
   float    rho25;		   
   TArrayF *jet_partonID;	   
   TArrayF *jet_dRJetGenJet;           

   float ijet_genJetPt,ijet_pt,ijet_eta,ijet_phi,ijet_En,ijet_cef,ijet_nef,ijet_mef,ijet_emfrac,ijet_hadfrac,ijet_nconstituents,ijet_chf,ijet_JECUnc,ijet_ptLeadTrack,ijet_vtxPt,ijet_vtx3dL,ijet_SoftLeptPtCut,ijet_dPhiMETJet,ijet_nch,ijet_vtx3deL,ijet_vtxMass,ijet_ptRaw,ijet_EnRaw,ijet_SoftLeptptRelCut,ijet_SoftLeptdRCut,ijet_partonID,ijet_dRJetGenJet;

   gamma1 = gamma2 = NULL;
   jets = NULL;
   bJetTags = NULL;
   jet_genJetPt = NULL;
   jet_eta = NULL;
   jet_cef = NULL;
   jet_nef = NULL;
   jet_mef = NULL;
   jet_nconstituents = NULL;
   jet_chf = NULL;
   jet_JECUnc = NULL;
   jet_ptLeadTrack = NULL;
   jet_vtxPt = NULL;
   jet_vtx3dL = NULL;
   jet_SoftLeptPtCut = NULL;
   jet_dPhiMETJet = NULL;
   jet_nch = NULL;
   jet_vtx3deL = NULL;
   jet_vtxMass = NULL;
   jet_ptRaw = NULL;
   jet_EnRaw = NULL;
   jet_SoftLeptptRelCut = NULL;
   jet_SoftLeptdRCut = NULL;
   jet_partonID = NULL;
   jet_dRJetGenJet = NULL;

   // associate input tree branches with variables
   tree->SetBranchAddress("evtNum", &evtNum);
   tree->SetBranchAddress("wei", &wei);
   tree->SetBranchAddress("gamma1", &gamma1);
   tree->SetBranchAddress("gamma2", &gamma2);
   tree->SetBranchAddress("jets", &jets);
   tree->SetBranchAddress(bJetTagsBranchName, &bJetTags);
   tree->SetBranchAddress("jet_genJetPt",        &jet_genJetPt          );
   tree->SetBranchAddress("jet_eta",		 &jet_eta		);     
   tree->SetBranchAddress("jet_cef",		 &jet_cef		);     
   tree->SetBranchAddress("jet_nef",		 &jet_nef		);     
   tree->SetBranchAddress("jet_mef",		 &jet_mef		);     
   tree->SetBranchAddress("jet_nconstituents",	 &jet_nconstituents	);     
   tree->SetBranchAddress("jet_chf",		 &jet_chf		);     
   tree->SetBranchAddress("jet_JECUnc",		 &jet_JECUnc		);     
   tree->SetBranchAddress("jet_ptLeadTrack",	 &jet_ptLeadTrack	);     
   tree->SetBranchAddress("jet_vtxPt",		 &jet_vtxPt		);     
   tree->SetBranchAddress("jet_vtx3dL",		 &jet_vtx3dL		);     
   tree->SetBranchAddress("jet_SoftLeptPtCut",	 &jet_SoftLeptPtCut	);     
   tree->SetBranchAddress("MET",		 &MET		   	);
   tree->SetBranchAddress("jet_dPhiMETJet",	 &jet_dPhiMETJet	);     
   tree->SetBranchAddress("jet_nch",		 &jet_nch		);     
   tree->SetBranchAddress("jet_vtx3deL",	 &jet_vtx3deL	   	);
   tree->SetBranchAddress("jet_vtxMass",	 &jet_vtxMass	   	);
   tree->SetBranchAddress("jet_ptRaw",		 &jet_ptRaw		);     
   tree->SetBranchAddress("jet_EnRaw",		 &jet_EnRaw		);     
   tree->SetBranchAddress("jet_SoftLeptptRelCut",&jet_SoftLeptptRelCut  );
   tree->SetBranchAddress("jet_SoftLeptdRCut",	 &jet_SoftLeptdRCut	);     
   tree->SetBranchAddress("rho25",		 &rho25		     	);
   tree->SetBranchAddress("jet_partonID",	 &jet_partonID	   	);
   tree->SetBranchAddress("jet_dRJetGenJet",         &jet_dRJetGenJet           );

   // create output tree
   TTree* treeOut = new TTree(treename, "Radion tree");

   // variables for output tree
   Double_t wei2;
   Double_t mGG, mJJ, mRad, ptRad, ptG1, ptG2;
   Double_t deltaRGJ, deltaRGG;
   Int_t bJetTagCategory;

   // associate output tree branches with variables
   treeOut->Branch("evtNum", &evtNum);
   treeOut->Branch("wei", &wei2);
   if(includeTrainingVars){
     treeOut->Branch("jet_genJetPt",     &ijet_genJetPt          );
     treeOut->Branch("jet_pt",           &ijet_pt                );
     treeOut->Branch("jet_eta",		 &ijet_eta		);     
     treeOut->Branch("jet_phi",		 &ijet_phi		);     
     treeOut->Branch("jet_En",		 &ijet_En		);     
     treeOut->Branch("jet_cef",		 &ijet_cef		);     
     treeOut->Branch("jet_nef",		 &ijet_nef		);     
     treeOut->Branch("jet_mef",		 &ijet_mef		);     
     treeOut->Branch("jet_emfrac",	 &ijet_emfrac		);     
     treeOut->Branch("jet_hadfrac",	 &ijet_hadfrac		);     
     treeOut->Branch("jet_nconstituents",&ijet_nconstituents	);     
     treeOut->Branch("jet_chf",		 &ijet_chf		);     
     treeOut->Branch("jet_JECUnc",	 &ijet_JECUnc		);     
     treeOut->Branch("jet_ptLeadTrack",	 &ijet_ptLeadTrack	);     
     treeOut->Branch("jet_vtxPt",	 &ijet_vtxPt		);     
     treeOut->Branch("jet_vtx3dL",	 &ijet_vtx3dL		);     
     treeOut->Branch("jet_SoftLeptPtCut",&ijet_SoftLeptPtCut	);     
     treeOut->Branch("MET",		 &MET		   	);
     treeOut->Branch("jet_dPhiMETJet",	 &ijet_dPhiMETJet	);     
     treeOut->Branch("jet_nch",		 &ijet_nch		);     
     treeOut->Branch("jet_vtx3deL",	 &ijet_vtx3deL	   	);
     treeOut->Branch("jet_vtxMass",	 &ijet_vtxMass	   	);
     treeOut->Branch("jet_ptRaw",	 &ijet_ptRaw		);     
     treeOut->Branch("jet_EnRaw",	 &ijet_EnRaw		);     
     treeOut->Branch("jet_SoftLeptptRelCut",&ijet_SoftLeptptRelCut  );
     treeOut->Branch("jet_SoftLeptdRCut",&ijet_SoftLeptdRCut	);     
     treeOut->Branch("rho25",		 &rho25		     	);
     treeOut->Branch("jet_partonID",	 &ijet_partonID	   	);
     treeOut->Branch("jet_dRJetGenJet",  &ijet_dRJetGenJet      );
   }
   else{
     treeOut->Branch("mGG", &mGG);
     treeOut->Branch("mJJ", &mJJ);
     treeOut->Branch("mRad", &mRad);
     treeOut->Branch("ptRad", &ptRad);
     if(includeWeightVars){
       treeOut->Branch("ptG1", &ptG1);
       treeOut->Branch("ptG2", &ptG2);
     }
     treeOut->Branch("deltaRGJ", &deltaRGJ);
     treeOut->Branch("deltaRGG", &deltaRGG);
     treeOut->Branch("bJetTagCategory", &bJetTagCategory);
   }

   // event loop
   for (Int_t i = 0; i < tree->GetEntriesFast(); i++) {
      tree->GetEntry(i);
      if(eventCut==0 && evtNum%2==1) continue; //keep even events
      if(eventCut==1 && evtNum%2==0) continue; //keep odd events

      wei2 = wei;

      //if training vars are included, do as few quality cuts as possible.
      if(includeTrainingVars){

	for (Int_t j = 0; j<jets->GetEntriesFast(); ++j) {
	  TLorentzVector* jet = (TLorentzVector*) jets->At(j);

	  // select jets coming from b quarks
	  if (abs(jet_partonID->At(j)) == 5){
	  
	    ijet_genJetPt         = jet_genJetPt         ->At(j);
	    ijet_pt               = jet->Pt();
	    ijet_eta		   = jet_eta		  ->At(j);
	    ijet_phi		   = jet->Phi();
	    ijet_En		   = jet->E();
	    ijet_cef		   = jet_cef		  ->At(j);
	    ijet_nef		   = jet_nef		  ->At(j);
	    ijet_mef		   = jet_mef		  ->At(j);
	    ijet_nconstituents	   = jet_nconstituents	  ->At(j);
	    ijet_chf		   = jet_chf		  ->At(j);
	    ijet_JECUnc	   = jet_JECUnc		  ->At(j);
	    ijet_ptLeadTrack	   = jet_ptLeadTrack	  ->At(j);
	    ijet_vtxPt	   = jet_vtxPt		  ->At(j);
	    ijet_vtx3dL	   = jet_vtx3dL		  ->At(j);
	    ijet_SoftLeptPtCut	   = jet_SoftLeptPtCut	  ->At(j);
	    ijet_dPhiMETJet	   = jet_dPhiMETJet	  ->At(j);
	    ijet_nch		   = jet_nch		  ->At(j);
	    ijet_vtx3deL	   = jet_vtx3deL	  ->At(j);
	    ijet_vtxMass	   = jet_vtxMass	  ->At(j);
	    ijet_ptRaw	   = jet_ptRaw		  ->At(j);
	    ijet_EnRaw	   = jet_EnRaw		  ->At(j);
	    ijet_SoftLeptptRelCut= jet_SoftLeptptRelCut ->At(j);
	    ijet_SoftLeptdRCut	   = jet_SoftLeptdRCut	  ->At(j);
	    ijet_partonID	   = jet_partonID	  ->At(j);
	    ijet_dRJetGenJet      = jet_dRJetGenJet      ->At(j);
	    ijet_emfrac = ijet_cef+ijet_nef+ijet_mef;
	    ijet_hadfrac = ijet_chf+ijet_nch;
	    
	    treeOut->Fill(); //MET and rho25 do not change when looping over jets.

	  }
	}//jet loop

	continue;
      }


      // if less than two photons or jets
      if (gamma1->E() <= 0 || gamma2->E() <= 0 || jets->GetEntriesFast() < 2)
         continue;

      // cut on photons' pt
      if (gamma1->Pt() < minPtGamma || gamma2->Pt() < minPtGamma)
         continue;

      // arrays of accepted + b-tagged/not b-tagged jets
      TObjArray bJets, otherJets;

      // fill the two jet arrays
      for (Int_t j = 0; j < jets->GetEntriesFast(); j++) {
         TLorentzVector* jet = (TLorentzVector*) jets->At(j);

         // cuts on jet's pt and eta
         if (jet->Pt() < minPtJet) continue;
         if (TMath::Abs(jet->Eta()) > maxEtaJet) continue;

         // skip jets pointing to directions of the photons
         Double_t v1 = gamma1->DeltaR(*jet);
         Double_t v2 = gamma2->DeltaR(*jet);
         Double_t deltaR = TMath::Min(v1, v2);
         if (deltaR < 0.3) continue;

         if ((*bJetTags)[j] >= minBJetTag)
            bJets.Add(jet);
         else
            otherJets.Add(jet);
      }

      // require to have at least two jets
      if (bJets.GetEntriesFast() + otherJets.GetEntriesFast() < 2)
         continue;

      // set event's b-tag category
      if (bJets.GetEntriesFast() >= 2)
         bJetTagCategory = 2;
      else
         bJetTagCategory = bJets.GetEntriesFast();

      // take two jets which give highest two-body pt
      TLorentzVector *jet1, *jet2;
      TLorentzVector *tmpJet1, *tmpJet2;
      float tmpHpt=-999;
      int tmpIndex1=-1, tmpIndex2=-1;

      if (bJets.GetEntriesFast() >= 2) {
	for(int itmp=0; itmp<bJets.GetEntriesFast(); ++itmp){
	  for(int jtmp=itmp+1; jtmp<bJets.GetEntriesFast(); ++jtmp){
	    tmpJet1 = (TLorentzVector*)bJets.At(itmp);
	    tmpJet2 = (TLorentzVector*)bJets.At(jtmp);
	    if(tmpHpt < ((TLorentzVector)(*tmpJet1+*tmpJet2)).M() ){
	      tmpHpt = ((TLorentzVector)(*tmpJet1+*tmpJet2)).M();
	      tmpIndex1=itmp;
	      tmpIndex2=jtmp;
	    }
	  }
	}
	jet1 = (TLorentzVector*) bJets.At(tmpIndex1);
	jet2 = (TLorentzVector*) bJets.At(tmpIndex2);
      }
      else if (bJets.GetEntriesFast() == 1) {
	tmpIndex1=0;
	tmpJet1 = (TLorentzVector*)bJets.At(tmpIndex1);
	for(int itmp=0; itmp<otherJets.GetEntriesFast(); ++itmp){
	  tmpJet2 = (TLorentzVector*)otherJets.At(itmp);
	    if(tmpHpt < ((TLorentzVector)(*tmpJet1+*tmpJet2)).M() ){
	      tmpHpt = ((TLorentzVector)(*tmpJet1+*tmpJet2)).M();
	      tmpIndex2=itmp;
	    }	  
	}
	jet1 = (TLorentzVector*) bJets.At(tmpIndex1);
	jet2 = (TLorentzVector*) otherJets.At(tmpIndex2);
      }
      else {
	for(int itmp=0; itmp<otherJets.GetEntriesFast(); ++itmp){
	  for(int jtmp=itmp+1; jtmp<otherJets.GetEntriesFast(); ++jtmp){
	    tmpJet1 = (TLorentzVector*)otherJets.At(itmp);
	    tmpJet2 = (TLorentzVector*)otherJets.At(jtmp);
	    if(tmpHpt < ((TLorentzVector)(*tmpJet1+*tmpJet2)).M() ){
	      tmpHpt = ((TLorentzVector)(*tmpJet1+*tmpJet2)).M();
	      tmpIndex1=itmp;
	      tmpIndex2=jtmp;
	    }
	  }
	}
	jet1 = (TLorentzVector*) otherJets.At(tmpIndex1);
	jet2 = (TLorentzVector*) otherJets.At(tmpIndex2);
      }

      // apply kinematic fitter
      if(doKinFit){
	DiJetKinFitter* fitter_jetsH = new DiJetKinFitter( "fitter_jetsH", "fitter_jets", 125.0 );
	TLorentzVector _jet1, _jet2;
	_jet1.SetPtEtaPhiE(jet1->Pt(),jet1->Eta(),jet1->Phi(),jet1->E());
	_jet2.SetPtEtaPhiE(jet2->Pt(),jet2->Eta(),jet2->Phi(),jet2->E());
	std::pair<TLorentzVector,TLorentzVector> jets_kinfitH = fitter_jetsH->fit(_jet1, _jet2);
	jet1 = &jets_kinfitH.first;
	jet2 = &jets_kinfitH.second;
      }

      // calculate deltaR
      Double_t v1 = gamma1->DeltaR(*jet1);
      Double_t v2 = gamma1->DeltaR(*jet2);
      Double_t v3 = gamma2->DeltaR(*jet1);
      Double_t v4 = gamma2->DeltaR(*jet2);
      deltaRGJ = TMath::Min(TMath::Min(v1, v2), TMath::Min(v3, v4));
      deltaRGG = gamma1->DeltaR(*gamma2);

      // deltaR cuts
      if (minDeltaRGJ >= 0 && deltaRGJ < minDeltaRGJ) continue;
      if (maxDeltaRGG >= 0 && deltaRGG > maxDeltaRGG) continue;

      // GG and JJ part
      TLorentzVector higgsGG = *gamma1 + *gamma2;
      TLorentzVector higgsJJ = *jet1 + *jet2;

      mGG = higgsGG.M();
      mJJ = higgsJJ.M();

      // cuts on reconstructed Higgs mass position
      if (mGG <= 0 || mJJ <= 0) continue;
      if (mGG < minMassGG || mGG > maxMassGG || mJJ < minMassJJ || mJJ > maxMassJJ)
         continue;

      // radion part
      TLorentzVector radion = higgsGG + higgsJJ;

      mRad = radion.M();
      if (mRad < 0) continue;
      ptRad = radion.Pt();

      //gamma pt for control sample
      ptG1=gamma1->Pt();
      ptG2=gamma2->Pt();

      treeOut->Fill();
   } // event loop

   treeOut->Write("", TObject::kOverwrite);

   // cleanup
   fi->Close();
   delete fi;
}

//______________________________________________________________________________
void mktree(char *input, char *output, int eventCut=-1)
{

   TFile* fo = TFile::Open(output, "RECREATE");

   mktree_one(input, "Events", eventCut);

   fo->Close();
   delete fo;
}

