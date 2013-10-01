
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
Double_t minDeltaRGJ = -1;  // -1 means unset

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
Bool_t doKinFit = kTRUE;

// correct TLorentzVectors of gamma+gamma and jet+jet by putting the Higgs mass
// to its place
Bool_t correctHiggsMass = kFALSE;

//______________________________________________________________________________
void mktree_one(const char* infile, const char* treename)
{
   // current working directory
   TDirectory* wd = gDirectory;

   // open input tree
   TFile* fi = TFile::Open(infile);
   TTree* tree = (TTree*) fi->Get("RToHH");

   // cd back into previous current working directory
   if (wd) wd->cd();

   // variables for input tree
   int evtNum;
   float wei;
   TLorentzVector *gamma1, *gamma2;
   TClonesArray* jets;
   TArrayF* bJetTags;

   gamma1 = gamma2 = NULL;
   jets = NULL;
   bJetTags = NULL;

   // associate input tree branches with variables
   tree->SetBranchAddress("evtNum", &evtNum);
   tree->SetBranchAddress("wei", &wei);
   tree->SetBranchAddress("gamma1", &gamma1);
   tree->SetBranchAddress("gamma2", &gamma2);
   tree->SetBranchAddress("jets", &jets);
   tree->SetBranchAddress(bJetTagsBranchName, &bJetTags);

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

   // event loop
   for (Int_t i = 0; i < tree->GetEntriesFast(); i++) {
      tree->GetEntry(i);

      wei2 = wei;

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

      // take two jets with highest pt
      TLorentzVector *jet1, *jet2;

      if (bJets.GetEntriesFast() >= 2) {
         jet1 = (TLorentzVector*) bJets.At(0);
         jet2 = (TLorentzVector*) bJets.At(1);
      }
      else if (bJets.GetEntriesFast() == 1) {
         jet1 = (TLorentzVector*) bJets.At(0);
         jet2 = (TLorentzVector*) otherJets.At(0);
      }
      else {
         jet1 = (TLorentzVector*) otherJets.At(0);
         jet2 = (TLorentzVector*) otherJets.At(1);
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

      // put Higgs to its mass; in principle, other algorithms are possible
      if (correctHiggsMass) {
         higgsGG *= 125.5/mGG;
         higgsJJ *= 125.5/mJJ;
      }

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
void mktree(char *input, char *output)
{

   TFile* fo = TFile::Open(output, "RECREATE");

   mktree_one(input, "Events");

   fo->Close();
   delete fo;
}

