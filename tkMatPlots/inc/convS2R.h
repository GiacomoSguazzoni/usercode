//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb  4 15:15:52 2013 by ROOT version 5.32/00
// from TTree ntupleS2R/sim2reco
// found on file: /raid/sguazz/maps/minbias_mc_reco_conv.root
//////////////////////////////////////////////////////////

#ifndef convS2R_h
#define convS2R_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class convS2R {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           isAssoc;
   Int_t           q1;
   Float_t         pt1;
   Float_t         phi1;
   Float_t         theta1;
   Float_t         qseed1;
   Float_t         ptseed1;
   Float_t         phiseed1;
   Float_t         thetaseed1;
   Float_t         rseed1;
   Float_t         zseed1;
   Int_t           r_algo1;
   Float_t         r_pt1;
   Float_t         r_d01;
   Float_t         r_theta1;
   Float_t         r_phi1;
   Float_t         r_chi21;
   Int_t           r_hits1;
   Int_t           r_before1;
   Int_t           r_missHits1;
   Int_t           q2;
   Float_t         pt2;
   Float_t         phi2;
   Float_t         theta2;
   Float_t         qseed2;
   Float_t         ptseed2;
   Float_t         phiseed2;
   Float_t         thetaseed2;
   Float_t         rseed2;
   Float_t         zseed2;
   Int_t           r_algo2;
   Float_t         r_pt2;
   Float_t         r_d02;
   Float_t         r_theta2;
   Float_t         r_phi2;
   Float_t         r_chi22;
   Int_t           r_hits2;
   Int_t           r_before2;
   Int_t           r_missHits2;
   Float_t         pt;
   Float_t         phi;
   Float_t         theta;
   Float_t         deltapt;
   Float_t         deltaphi;
   Float_t         deltatheta;
   Float_t         x;
   Float_t         y;
   Float_t         z;
   Float_t         deltax;
   Float_t         deltay;
   Float_t         deltaz;
   Float_t         r_minapp;
   Float_t         r_chi2prob;
   Int_t           q_gto;
   Int_t           q_aes;
   Int_t           q_am;
   Int_t           q_ameg;
   Int_t           q_hp;
   Int_t           q_he;
   Int_t           q_em1t;
   Int_t           q_em2t;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_isAssoc;   //!
   TBranch        *b_q1;   //!
   TBranch        *b_pt1;   //!
   TBranch        *b_phi1;   //!
   TBranch        *b_theta1;   //!
   TBranch        *b_qseed1;   //!
   TBranch        *b_ptseed1;   //!
   TBranch        *b_phiseed1;   //!
   TBranch        *b_thetaseed1;   //!
   TBranch        *b_rseed1;   //!
   TBranch        *b_zseed1;   //!
   TBranch        *b_r_algo1;   //!
   TBranch        *b_r_pt1;   //!
   TBranch        *b_r_d01;   //!
   TBranch        *b_r_theta1;   //!
   TBranch        *b_r_phi1;   //!
   TBranch        *b_r_chi21;   //!
   TBranch        *b_r_hits1;   //!
   TBranch        *b_r_before1;   //!
   TBranch        *b_r_missHits1;   //!
   TBranch        *b_q2;   //!
   TBranch        *b_pt2;   //!
   TBranch        *b_phi2;   //!
   TBranch        *b_theta2;   //!
   TBranch        *b_qseed2;   //!
   TBranch        *b_ptseed2;   //!
   TBranch        *b_phiseed2;   //!
   TBranch        *b_thetaseed2;   //!
   TBranch        *b_rseed2;   //!
   TBranch        *b_zseed2;   //!
   TBranch        *b_r_algo2;   //!
   TBranch        *b_r_pt2;   //!
   TBranch        *b_r_d02;   //!
   TBranch        *b_r_theta2;   //!
   TBranch        *b_r_phi2;   //!
   TBranch        *b_r_chi22;   //!
   TBranch        *b_r_hits2;   //!
   TBranch        *b_r_before2;   //!
   TBranch        *b_r_missHits2;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_deltapt;   //!
   TBranch        *b_deltaphi;   //!
   TBranch        *b_deltatheta;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_deltax;   //!
   TBranch        *b_deltay;   //!
   TBranch        *b_deltaz;   //!
   TBranch        *b_r_minapp;   //!
   TBranch        *b_r_chi2prob;   //!
   TBranch        *b_q_gto;   //!
   TBranch        *b_q_aes;   //!
   TBranch        *b_q_am;   //!
   TBranch        *b_q_ameg;   //!
   TBranch        *b_q_hp;   //!
   TBranch        *b_q_he;   //!
   TBranch        *b_q_em1t;   //!
   TBranch        *b_q_em2t;   //!

   convS2R(TTree *tree=0);
   virtual ~convS2R();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef convS2R_cxx
convS2R::convS2R(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/raid/sguazz/maps/minbias_mc_reco_conv.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/raid/sguazz/maps/minbias_mc_reco_conv.root");
      }
      f->GetObject("ntupleS2R",tree);

   }
   Init(tree);
}

convS2R::~convS2R()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t convS2R::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t convS2R::LoadTree(Long64_t entry)
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

void convS2R::Init(TTree *tree)
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
   fChain->SetBranchAddress("q1", &q1, &b_q1);
   fChain->SetBranchAddress("pt1", &pt1, &b_pt1);
   fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
   fChain->SetBranchAddress("theta1", &theta1, &b_theta1);
   fChain->SetBranchAddress("qseed1", &qseed1, &b_qseed1);
   fChain->SetBranchAddress("ptseed1", &ptseed1, &b_ptseed1);
   fChain->SetBranchAddress("phiseed1", &phiseed1, &b_phiseed1);
   fChain->SetBranchAddress("thetaseed1", &thetaseed1, &b_thetaseed1);
   fChain->SetBranchAddress("rseed1", &rseed1, &b_rseed1);
   fChain->SetBranchAddress("zseed1", &zseed1, &b_zseed1);
   fChain->SetBranchAddress("r_algo1", &r_algo1, &b_r_algo1);
   fChain->SetBranchAddress("r_pt1", &r_pt1, &b_r_pt1);
   fChain->SetBranchAddress("r_d01", &r_d01, &b_r_d01);
   fChain->SetBranchAddress("r_theta1", &r_theta1, &b_r_theta1);
   fChain->SetBranchAddress("r_phi1", &r_phi1, &b_r_phi1);
   fChain->SetBranchAddress("r_chi21", &r_chi21, &b_r_chi21);
   fChain->SetBranchAddress("r_hits1", &r_hits1, &b_r_hits1);
   fChain->SetBranchAddress("r_before1", &r_before1, &b_r_before1);
   fChain->SetBranchAddress("r_missHits1", &r_missHits1, &b_r_missHits1);
   fChain->SetBranchAddress("q2", &q2, &b_q2);
   fChain->SetBranchAddress("pt2", &pt2, &b_pt2);
   fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fChain->SetBranchAddress("theta2", &theta2, &b_theta2);
   fChain->SetBranchAddress("qseed2", &qseed2, &b_qseed2);
   fChain->SetBranchAddress("ptseed2", &ptseed2, &b_ptseed2);
   fChain->SetBranchAddress("phiseed2", &phiseed2, &b_phiseed2);
   fChain->SetBranchAddress("thetaseed2", &thetaseed2, &b_thetaseed2);
   fChain->SetBranchAddress("rseed2", &rseed2, &b_rseed2);
   fChain->SetBranchAddress("zseed2", &zseed2, &b_zseed2);
   fChain->SetBranchAddress("r_algo2", &r_algo2, &b_r_algo2);
   fChain->SetBranchAddress("r_pt2", &r_pt2, &b_r_pt2);
   fChain->SetBranchAddress("r_d02", &r_d02, &b_r_d02);
   fChain->SetBranchAddress("r_theta2", &r_theta2, &b_r_theta2);
   fChain->SetBranchAddress("r_phi2", &r_phi2, &b_r_phi2);
   fChain->SetBranchAddress("r_chi22", &r_chi22, &b_r_chi22);
   fChain->SetBranchAddress("r_hits2", &r_hits2, &b_r_hits2);
   fChain->SetBranchAddress("r_before2", &r_before2, &b_r_before2);
   fChain->SetBranchAddress("r_missHits2", &r_missHits2, &b_r_missHits2);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("deltapt", &deltapt, &b_deltapt);
   fChain->SetBranchAddress("deltaphi", &deltaphi, &b_deltaphi);
   fChain->SetBranchAddress("deltatheta", &deltatheta, &b_deltatheta);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("deltax", &deltax, &b_deltax);
   fChain->SetBranchAddress("deltay", &deltay, &b_deltay);
   fChain->SetBranchAddress("deltaz", &deltaz, &b_deltaz);
   fChain->SetBranchAddress("r_minapp", &r_minapp, &b_r_minapp);
   fChain->SetBranchAddress("r_chi2prob", &r_chi2prob, &b_r_chi2prob);
   fChain->SetBranchAddress("q_gto", &q_gto, &b_q_gto);
   fChain->SetBranchAddress("q_aes", &q_aes, &b_q_aes);
   fChain->SetBranchAddress("q_am", &q_am, &b_q_am);
   fChain->SetBranchAddress("q_ameg", &q_ameg, &b_q_ameg);
   fChain->SetBranchAddress("q_hp", &q_hp, &b_q_hp);
   fChain->SetBranchAddress("q_he", &q_he, &b_q_he);
   fChain->SetBranchAddress("q_em1t", &q_em1t, &b_q_em1t);
   fChain->SetBranchAddress("q_em2t", &q_em2t, &b_q_em2t);
   Notify();

}

Bool_t convS2R::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void convS2R::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t convS2R::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef convS2R_cxx
