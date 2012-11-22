//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct  8 15:00:59 2012 by ROOT version 5.32/00
// from TTree sT/sT
// found on file: /afs/cern.ch/user/s/sguazzon/myWorkarea/split/minBias_nH4_v2.root
//////////////////////////////////////////////////////////

#ifndef sT_h
#define sT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class sT {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nRun;
   Int_t           nEvent;
   Int_t           iEvent;
   Int_t           NTrack;
   Int_t           iTrack;
   Float_t         p;
   Float_t         pErr;
   Float_t         pAve;
   Float_t         pAveErr;
   Float_t         pt;
   Float_t         ptErr;
   Float_t         eta;
   Float_t         the;
   Float_t         phi;
   Int_t           nHitVal;
   Int_t           nHitMis;
   Int_t           nHitIna;
   Int_t           nHitBad;
   Int_t           nEffHit;
   Int_t           q;
   Float_t         dxy;
   Float_t         dz;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         chi2;
   Float_t         maxchi2;
   Float_t         T;
   Float_t         rIn;
   Float_t         rOut;
   Float_t         xIn;
   Float_t         yIn;
   Float_t         zIn;
   Float_t         xOut;
   Float_t         yOut;
   Float_t         zOut;
   Int_t           iTrackSim;
   Float_t         pSim;
   Float_t         ptSim;
   Float_t         etaSim;
   Float_t         theSim;
   Float_t         phiSim;
   Float_t         pLossSim;
   Float_t         dpdxSim;
   Int_t           nHitSim;
   Float_t         hitFrac;
   Float_t         TSim;
   Float_t         rInSim;
   Float_t         rOutSim;
   Float_t         xInSim;
   Float_t         yInSim;
   Float_t         zInSim;
   Float_t         xOutSim;
   Float_t         yOutSim;
   Float_t         zOutSim;
   Int_t           NSplit;
   Float_t         TInSplit;
   Float_t         TOutSplit;
   Float_t         dpdxSplit;
   Float_t         dpdxErrSplit;
   Float_t         chi2Split;
   Float_t         dpdxTSplit;
   Float_t         dpdxTErrSplit;
   Float_t         chi2TSplit;
   Float_t         dpdxZSplit;
   Float_t         dpdxZErrSplit;
   Float_t         chi2ZSplit;

   // List of branches
   TBranch        *b_nRun;   //!
   TBranch        *b_nEvent;   //!
   TBranch        *b_iEvent;   //!
   TBranch        *b_NTrack;   //!
   TBranch        *b_iTrack;   //!
   TBranch        *b_p;   //!
   TBranch        *b_pErr;   //!
   TBranch        *b_pAve;   //!
   TBranch        *b_pAveErr;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_ptErr;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_the;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_nHitVal;   //!
   TBranch        *b_nHitMis;   //!
   TBranch        *b_nHitIna;   //!
   TBranch        *b_nHitBad;   //!
   TBranch        *b_nEffHit;   //!
   TBranch        *b_q;   //!
   TBranch        *b_dxy;   //!
   TBranch        *b_dz;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_maxchi2;   //!
   TBranch        *b_T;   //!
   TBranch        *b_rIn;   //!
   TBranch        *b_rOut;   //!
   TBranch        *b_xIn;   //!
   TBranch        *b_yIn;   //!
   TBranch        *b_zIn;   //!
   TBranch        *b_xOut;   //!
   TBranch        *b_yOut;   //!
   TBranch        *b_zOut;   //!
   TBranch        *b_iTrackSim;   //!
   TBranch        *b_pSim;   //!
   TBranch        *b_ptSim;   //!
   TBranch        *b_etaSim;   //!
   TBranch        *b_theSim;   //!
   TBranch        *b_phiSim;   //!
   TBranch        *b_pLossSim;   //!
   TBranch        *b_dpdxSim;   //!
   TBranch        *b_nHitSim;   //!
   TBranch        *b_hitFrac;   //!
   TBranch        *b_TSim;   //!
   TBranch        *b_rInSim;   //!
   TBranch        *b_rOutSim;   //!
   TBranch        *b_xInSim;   //!
   TBranch        *b_yInSim;   //!
   TBranch        *b_zInSim;   //!
   TBranch        *b_xOutSim;   //!
   TBranch        *b_yOutSim;   //!
   TBranch        *b_zOutSim;   //!
   TBranch        *b_NSplit;   //!
   TBranch        *b_TInSplit;   //!
   TBranch        *b_TOutSplit;   //!
   TBranch        *b_dpdxSplit;   //!
   TBranch        *b_dpdxSplitErr;   //!
   TBranch        *b_chi2Split;   //!
   TBranch        *b_dpdxTSplit;   //!
   TBranch        *b_dpdxTSplitErr;   //!
   TBranch        *b_chi2TSplit;   //!
   TBranch        *b_dpdxZSplit;   //!
   TBranch        *b_dpdxZSplitErr;   //!
   TBranch        *b_chi2ZSplit;   //!

   sT(TTree *tree=0);
   virtual ~sT();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef sT_cxx
sT::sT(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/user/s/sguazzon/myWorkarea/split/minBias_nH4_v2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/afs/cern.ch/user/s/sguazzon/myWorkarea/split/minBias_nH4_v2.root");
      }
      f->GetObject("sT",tree);

   }
   Init(tree);
}

sT::~sT()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t sT::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t sT::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void sT::Init(TTree *tree)
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

   fChain->SetBranchAddress("nRun", &nRun, &b_nRun);
   fChain->SetBranchAddress("nEvent", &nEvent, &b_nEvent);
   fChain->SetBranchAddress("iEvent", &iEvent, &b_iEvent);
   fChain->SetBranchAddress("NTrack", &NTrack, &b_NTrack);
   fChain->SetBranchAddress("iTrack", &iTrack, &b_iTrack);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("pErr", &pErr, &b_pErr);
   fChain->SetBranchAddress("pAve", &pAve, &b_pAve);
   fChain->SetBranchAddress("pAveErr", &pAveErr, &b_pAveErr);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("ptErr", &ptErr, &b_ptErr);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("the", &the, &b_the);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("nHitVal", &nHitVal, &b_nHitVal);
   fChain->SetBranchAddress("nHitMis", &nHitMis, &b_nHitMis);
   fChain->SetBranchAddress("nHitIna", &nHitIna, &b_nHitIna);
   fChain->SetBranchAddress("nHitBad", &nHitBad, &b_nHitBad);
   fChain->SetBranchAddress("nEffHit", &nEffHit, &b_nEffHit);
   fChain->SetBranchAddress("q", &q, &b_q);
   fChain->SetBranchAddress("dxy", &dxy, &b_dxy);
   fChain->SetBranchAddress("dz", &dz, &b_dz);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("maxchi2", &maxchi2, &b_maxchi2);
   fChain->SetBranchAddress("T", &T, &b_T);
   fChain->SetBranchAddress("rIn", &rIn, &b_rIn);
   fChain->SetBranchAddress("rOut", &rOut, &b_rOut);
   fChain->SetBranchAddress("xIn", &xIn, &b_xIn);
   fChain->SetBranchAddress("yIn", &yIn, &b_yIn);
   fChain->SetBranchAddress("zIn", &zIn, &b_zIn);
   fChain->SetBranchAddress("xOut", &xOut, &b_xOut);
   fChain->SetBranchAddress("yOut", &yOut, &b_yOut);
   fChain->SetBranchAddress("zOut", &zOut, &b_zOut);
   fChain->SetBranchAddress("iTrackSim", &iTrackSim, &b_iTrackSim);
   fChain->SetBranchAddress("pSim", &pSim, &b_pSim);
   fChain->SetBranchAddress("ptSim", &ptSim, &b_ptSim);
   fChain->SetBranchAddress("etaSim", &etaSim, &b_etaSim);
   fChain->SetBranchAddress("theSim", &theSim, &b_theSim);
   fChain->SetBranchAddress("phiSim", &phiSim, &b_phiSim);
   fChain->SetBranchAddress("pLossSim", &pLossSim, &b_pLossSim);
   fChain->SetBranchAddress("dpdxSim", &dpdxSim, &b_dpdxSim);
   fChain->SetBranchAddress("nHitSim", &nHitSim, &b_nHitSim);
   fChain->SetBranchAddress("hitFrac", &hitFrac, &b_hitFrac);
   fChain->SetBranchAddress("TSim", &TSim, &b_TSim);
   fChain->SetBranchAddress("rInSim", &rInSim, &b_rInSim);
   fChain->SetBranchAddress("rOutSim", &rOutSim, &b_rOutSim);
   fChain->SetBranchAddress("xInSim", &xInSim, &b_xInSim);
   fChain->SetBranchAddress("yInSim", &yInSim, &b_yInSim);
   fChain->SetBranchAddress("zInSim", &zInSim, &b_zInSim);
   fChain->SetBranchAddress("xOutSim", &xOutSim, &b_xOutSim);
   fChain->SetBranchAddress("yOutSim", &yOutSim, &b_yOutSim);
   fChain->SetBranchAddress("zOutSim", &zOutSim, &b_zOutSim);
   fChain->SetBranchAddress("NSplit", &NSplit, &b_NSplit);
   fChain->SetBranchAddress("TInSplit", &TInSplit, &b_TInSplit);
   fChain->SetBranchAddress("TOutSplit", &TOutSplit, &b_TOutSplit);
   fChain->SetBranchAddress("dpdxSplit", &dpdxSplit, &b_dpdxSplit);
   fChain->SetBranchAddress("dpdxErrSplit", &dpdxErrSplit, &b_dpdxSplitErr);
   fChain->SetBranchAddress("chi2Split", &chi2Split, &b_chi2Split);
   fChain->SetBranchAddress("dpdxTSplit", &dpdxTSplit, &b_dpdxTSplit);
   fChain->SetBranchAddress("dpdxTErrSplit", &dpdxTErrSplit, &b_dpdxTSplitErr);
   fChain->SetBranchAddress("chi2TSplit", &chi2TSplit, &b_chi2TSplit);
   fChain->SetBranchAddress("dpdxZSplit", &dpdxZSplit, &b_dpdxZSplit);
   fChain->SetBranchAddress("dpdxZErrSplit", &dpdxZErrSplit, &b_dpdxZSplitErr);
   fChain->SetBranchAddress("chi2ZSplit", &chi2ZSplit, &b_chi2ZSplit);
   Notify();
}

Bool_t sT::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void sT::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t sT::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef sT_cxx
