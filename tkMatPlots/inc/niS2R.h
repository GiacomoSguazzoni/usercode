//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun  5 17:13:05 2010 by ROOT version 5.22/00b
// from TTree ntupleS2R/sim2reco
// found on file: ../7TeV/ntuple_nuclint_CMSSW358p3_goodcoll7TeV_2010-06-04.root
//////////////////////////////////////////////////////////

#ifndef niS2R_h
#define niS2R_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class niS2R {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           isAssoc;
   Char_t          isNuclSim;
   Char_t          isKSim;
   Float_t         pt;
   Float_t         phi;
   Float_t         theta;
   Float_t         ptOut;
   Float_t         phiOut;
   Float_t         thetaOut;
   Float_t         mOut;
   Int_t           nOut;
   Float_t         x;
   Float_t         y;
   Float_t         z;
   Float_t         deltapt;
   Float_t         deltaphi;
   Float_t         deltatheta;
   Float_t         deltax;
   Float_t         deltay;
   Float_t         deltaz;
   Char_t          isTherePrimaryTrack;
   Char_t          isThereMergedTrack;
   Char_t          isNucl;
   Char_t          isNuclLoose;
   Char_t          isNuclKink;
   Char_t          isFake;
   Char_t          isK0;
   Char_t          isLambda;
   Char_t          isLambdaBar;
   Char_t          isKminusLoose;
   Char_t          isKplusLoose;
   Char_t          isLooper;
   Char_t          isConvLoose;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_isAssoc;   //!
   TBranch        *b_isNuclSim;   //!
   TBranch        *b_isKSim;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_ptOut;   //!
   TBranch        *b_phiOut;   //!
   TBranch        *b_thetaOut;   //!
   TBranch        *b_mOut;   //!
   TBranch        *b_nOut;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_deltapt;   //!
   TBranch        *b_deltaphi;   //!
   TBranch        *b_deltatheta;   //!
   TBranch        *b_deltax;   //!
   TBranch        *b_deltay;   //!
   TBranch        *b_deltaz;   //!
   TBranch        *b_isTherePrimaryTrack;   //!
   TBranch        *b_isThereMergedTrack;   //!
   TBranch        *b_isNucl;   //!
   TBranch        *b_isNuclLoose;   //!
   TBranch        *b_isNuclKink;   //!
   TBranch        *b_isFake;   //!
   TBranch        *b_isK0;   //!
   TBranch        *b_isLambda;   //!
   TBranch        *b_isLambdaBar;   //!
   TBranch        *b_isKminusLoose;   //!
   TBranch        *b_isKplusLoose;   //!
   TBranch        *b_isLooper;   //!
   TBranch        *b_isConvLoose;   //!

   niS2R(TTree *tree=0);
   virtual ~niS2R();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef niS2R_cxx
niS2R::niS2R(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../7TeV/ntuple_nuclint_CMSSW358p3_goodcoll7TeV_2010-06-04.root");
      if (!f) {
         f = new TFile("../7TeV/ntuple_nuclint_CMSSW358p3_goodcoll7TeV_2010-06-04.root");
      }
      tree = (TTree*)gDirectory->Get("ntupleS2R");

   }
   Init(tree);
}

niS2R::~niS2R()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t niS2R::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t niS2R::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void niS2R::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("isAssoc", &isAssoc, &b_isAssoc);
   fChain->SetBranchAddress("isNuclSim", &isNuclSim, &b_isNuclSim);
   fChain->SetBranchAddress("isKSim", &isKSim, &b_isKSim);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("ptOut", &ptOut, &b_ptOut);
   fChain->SetBranchAddress("phiOut", &phiOut, &b_phiOut);
   fChain->SetBranchAddress("thetaOut", &thetaOut, &b_thetaOut);
   fChain->SetBranchAddress("mOut", &mOut, &b_mOut);
   fChain->SetBranchAddress("nOut", &nOut, &b_nOut);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("deltapt", &deltapt, &b_deltapt);
   fChain->SetBranchAddress("deltaphi", &deltaphi, &b_deltaphi);
   fChain->SetBranchAddress("deltatheta", &deltatheta, &b_deltatheta);
   fChain->SetBranchAddress("deltax", &deltax, &b_deltax);
   fChain->SetBranchAddress("deltay", &deltay, &b_deltay);
   fChain->SetBranchAddress("deltaz", &deltaz, &b_deltaz);
   fChain->SetBranchAddress("isTherePrimaryTrack", &isTherePrimaryTrack, &b_isTherePrimaryTrack);
   fChain->SetBranchAddress("isThereMergedTrack", &isThereMergedTrack, &b_isThereMergedTrack);
   fChain->SetBranchAddress("isNucl", &isNucl, &b_isNucl);
   fChain->SetBranchAddress("isNuclLoose", &isNuclLoose, &b_isNuclLoose);
   fChain->SetBranchAddress("isNuclKink", &isNuclKink, &b_isNuclKink);
   fChain->SetBranchAddress("isFake", &isFake, &b_isFake);
   fChain->SetBranchAddress("isK0", &isK0, &b_isK0);
   fChain->SetBranchAddress("isLambda", &isLambda, &b_isLambda);
   fChain->SetBranchAddress("isLambdaBar", &isLambdaBar, &b_isLambdaBar);
   fChain->SetBranchAddress("isKminusLoose", &isKminusLoose, &b_isKminusLoose);
   fChain->SetBranchAddress("isKplusLoose", &isKplusLoose, &b_isKplusLoose);
   fChain->SetBranchAddress("isLooper", &isLooper, &b_isLooper);
   fChain->SetBranchAddress("isConvLoose", &isConvLoose, &b_isConvLoose);
   Notify();
}

Bool_t niS2R::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void niS2R::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t niS2R::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef niS2R_cxx
