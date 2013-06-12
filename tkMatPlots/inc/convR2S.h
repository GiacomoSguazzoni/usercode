//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb  4 09:35:51 2013 by ROOT version 5.32/00
// from TTree ntupleR2S/reco2sim
// found on file: /raid/sguazz/maps/minbias_run2012D_aod_conv.root
//////////////////////////////////////////////////////////

#ifndef convR2S_h
#define convR2S_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class convR2S {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           isAssoc;
   Int_t           isDouble;
   Int_t           q1;
   Int_t           algo1;
   Float_t         d01;
   Float_t         dz1;
   Float_t         pt1;
   Float_t         phi1;
   Float_t         iphi1;
   Float_t         theta1;
   Float_t         x1;
   Float_t         y1;
   Float_t         z1;
   Float_t         chi21;
   Int_t           nHits1;
   Int_t           beforeHits1;
   Int_t           npHits1;
   Int_t           nsHits1;
   Int_t           missHits1;
   Float_t         lambdaError1;
   Float_t         r_pt1;
   Float_t         r_d01;
   Float_t         r_phi1;
   Float_t         r_theta1;
   Float_t         r_px1;
   Float_t         r_py1;
   Float_t         r_pz1;
   Float_t         ipx1;
   Float_t         ipy1;
   Float_t         ipz1;
   Float_t         ix1;
   Float_t         iy1;
   Float_t         iz1;
   Float_t         ox1;
   Float_t         oy1;
   Float_t         oz1;
   Int_t           q2;
   Int_t           algo2;
   Float_t         d02;
   Float_t         dz2;
   Float_t         pt2;
   Float_t         phi2;
   Float_t         iphi2;
   Float_t         theta2;
   Float_t         x2;
   Float_t         y2;
   Float_t         z2;
   Float_t         chi22;
   Int_t           nHits2;
   Int_t           beforeHits2;
   Int_t           npHits2;
   Int_t           nsHits2;
   Int_t           missHits2;
   Float_t         lambdaError2;
   Float_t         r_pt2;
   Float_t         r_d02;
   Float_t         r_phi2;
   Float_t         r_theta2;
   Float_t         r_px2;
   Float_t         r_py2;
   Float_t         r_pz2;
   Float_t         ipx2;
   Float_t         ipy2;
   Float_t         ipz2;
   Float_t         ix2;
   Float_t         iy2;
   Float_t         iz2;
   Float_t         ox2;
   Float_t         oy2;
   Float_t         oz2;
   Float_t         pt;
   Float_t         phi;
   Float_t         theta;
   Float_t         mass;
   Float_t         r_pt;
   Float_t         r_phi;
   Float_t         r_theta;
   Float_t         r_mass;
   Float_t         chi2;
   Float_t         deltapt;
   Float_t         deltaphi;
   Float_t         deltatheta;
   Float_t         x;
   Float_t         y;
   Float_t         z;
   Float_t         deltax;
   Float_t         deltay;
   Float_t         deltaz;
   Float_t         minapp;
   Int_t           refit;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         tx;
   Float_t         ty;
   Float_t         tz;
   Float_t         bsx;
   Float_t         bsy;
   Float_t         bsz;
   Float_t         chi2prob;
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
   TBranch        *b_isDouble;   //!
   TBranch        *b_q1;   //!
   TBranch        *b_algo1;   //!
   TBranch        *b_d01;   //!
   TBranch        *b_dz1;   //!
   TBranch        *b_pt1;   //!
   TBranch        *b_phi1;   //!
   TBranch        *b_iphi1;   //!
   TBranch        *b_theta1;   //!
   TBranch        *b_x1;   //!
   TBranch        *b_y1;   //!
   TBranch        *b_z1;   //!
   TBranch        *b_chi21;   //!
   TBranch        *b_nHits1;   //!
   TBranch        *b_beforeHits1;   //!
   TBranch        *b_npHits1;   //!
   TBranch        *b_nsHits1;   //!
   TBranch        *b_missHits1;   //!
   TBranch        *b_lambdaError1;   //!
   TBranch        *b_r_pt1;   //!
   TBranch        *b_r_d01;   //!
   TBranch        *b_r_phi1;   //!
   TBranch        *b_r_theta1;   //!
   TBranch        *b_r_px1;   //!
   TBranch        *b_r_py1;   //!
   TBranch        *b_r_pz1;   //!
   TBranch        *b_ipx1;   //!
   TBranch        *b_ipy1;   //!
   TBranch        *b_ipz1;   //!
   TBranch        *b_ix1;   //!
   TBranch        *b_iy1;   //!
   TBranch        *b_iz1;   //!
   TBranch        *b_ox1;   //!
   TBranch        *b_oy1;   //!
   TBranch        *b_oz1;   //!
   TBranch        *b_q2;   //!
   TBranch        *b_algo2;   //!
   TBranch        *b_d02;   //!
   TBranch        *b_dz2;   //!
   TBranch        *b_pt2;   //!
   TBranch        *b_phi2;   //!
   TBranch        *b_iphi2;   //!
   TBranch        *b_theta2;   //!
   TBranch        *b_x2;   //!
   TBranch        *b_y2;   //!
   TBranch        *b_z2;   //!
   TBranch        *b_chi22;   //!
   TBranch        *b_nHits2;   //!
   TBranch        *b_beforeHits2;   //!
   TBranch        *b_npHits2;   //!
   TBranch        *b_nsHits2;   //!
   TBranch        *b_missHits2;   //!
   TBranch        *b_lambdaError2;   //!
   TBranch        *b_r_pt2;   //!
   TBranch        *b_r_d02;   //!
   TBranch        *b_r_phi2;   //!
   TBranch        *b_r_theta2;   //!
   TBranch        *b_r_px2;   //!
   TBranch        *b_r_py2;   //!
   TBranch        *b_r_pz2;   //!
   TBranch        *b_ipx2;   //!
   TBranch        *b_ipy2;   //!
   TBranch        *b_ipz2;   //!
   TBranch        *b_ix2;   //!
   TBranch        *b_iy2;   //!
   TBranch        *b_iz2;   //!
   TBranch        *b_ox2;   //!
   TBranch        *b_oy2;   //!
   TBranch        *b_oz2;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_r_pt;   //!
   TBranch        *b_r_phi;   //!
   TBranch        *b_r_theta;   //!
   TBranch        *b_r_mass;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_deltapt;   //!
   TBranch        *b_deltaphi;   //!
   TBranch        *b_deltatheta;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_deltax;   //!
   TBranch        *b_deltay;   //!
   TBranch        *b_deltaz;   //!
   TBranch        *b_minapp;   //!
   TBranch        *b_refit;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_tx;   //!
   TBranch        *b_ty;   //!
   TBranch        *b_tz;   //!
   TBranch        *b_bsx;   //!
   TBranch        *b_bsy;   //!
   TBranch        *b_bsz;   //!
   TBranch        *b_chi2prob;   //!
   TBranch        *b_q_gto;   //!
   TBranch        *b_q_aes;   //!
   TBranch        *b_q_am;   //!
   TBranch        *b_q_ameg;   //!
   TBranch        *b_q_hp;   //!
   TBranch        *b_q_he;   //!
   TBranch        *b_q_em1t;   //!
   TBranch        *b_q_em2t;   //!

   convR2S(TTree *tree=0);
   virtual ~convR2S();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef convR2S_cxx
convR2S::convR2S(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/raid/sguazz/maps/minbias_run2012D_aod_conv.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/raid/sguazz/maps/minbias_run2012D_aod_conv.root");
      }
      f->GetObject("ntupleR2S",tree);

   }
   Init(tree);
}

convR2S::~convR2S()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t convR2S::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t convR2S::LoadTree(Long64_t entry)
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

void convR2S::Init(TTree *tree)
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
   fChain->SetBranchAddress("isDouble", &isDouble, &b_isDouble);
   fChain->SetBranchAddress("q1", &q1, &b_q1);
   fChain->SetBranchAddress("algo1", &algo1, &b_algo1);
   fChain->SetBranchAddress("d01", &d01, &b_d01);
   fChain->SetBranchAddress("dz1", &dz1, &b_dz1);
   fChain->SetBranchAddress("pt1", &pt1, &b_pt1);
   fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
   fChain->SetBranchAddress("iphi1", &iphi1, &b_iphi1);
   fChain->SetBranchAddress("theta1", &theta1, &b_theta1);
   fChain->SetBranchAddress("x1", &x1, &b_x1);
   fChain->SetBranchAddress("y1", &y1, &b_y1);
   fChain->SetBranchAddress("z1", &z1, &b_z1);
   fChain->SetBranchAddress("chi21", &chi21, &b_chi21);
   fChain->SetBranchAddress("nHits1", &nHits1, &b_nHits1);
   fChain->SetBranchAddress("beforeHits1", &beforeHits1, &b_beforeHits1);
   fChain->SetBranchAddress("npHits1", &npHits1, &b_npHits1);
   fChain->SetBranchAddress("nsHits1", &nsHits1, &b_nsHits1);
   fChain->SetBranchAddress("missHits1", &missHits1, &b_missHits1);
   fChain->SetBranchAddress("lambdaError1", &lambdaError1, &b_lambdaError1);
   fChain->SetBranchAddress("r_pt1", &r_pt1, &b_r_pt1);
   fChain->SetBranchAddress("r_d01", &r_d01, &b_r_d01);
   fChain->SetBranchAddress("r_phi1", &r_phi1, &b_r_phi1);
   fChain->SetBranchAddress("r_theta1", &r_theta1, &b_r_theta1);
   fChain->SetBranchAddress("r_px1", &r_px1, &b_r_px1);
   fChain->SetBranchAddress("r_py1", &r_py1, &b_r_py1);
   fChain->SetBranchAddress("r_pz1", &r_pz1, &b_r_pz1);
   fChain->SetBranchAddress("ipx1", &ipx1, &b_ipx1);
   fChain->SetBranchAddress("ipy1", &ipy1, &b_ipy1);
   fChain->SetBranchAddress("ipz1", &ipz1, &b_ipz1);
   fChain->SetBranchAddress("ix1", &ix1, &b_ix1);
   fChain->SetBranchAddress("iy1", &iy1, &b_iy1);
   fChain->SetBranchAddress("iz1", &iz1, &b_iz1);
   fChain->SetBranchAddress("ox1", &ox1, &b_ox1);
   fChain->SetBranchAddress("oy1", &oy1, &b_oy1);
   fChain->SetBranchAddress("oz1", &oz1, &b_oz1);
   fChain->SetBranchAddress("q2", &q2, &b_q2);
   fChain->SetBranchAddress("algo2", &algo2, &b_algo2);
   fChain->SetBranchAddress("d02", &d02, &b_d02);
   fChain->SetBranchAddress("dz2", &dz2, &b_dz2);
   fChain->SetBranchAddress("pt2", &pt2, &b_pt2);
   fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fChain->SetBranchAddress("iphi2", &iphi2, &b_iphi2);
   fChain->SetBranchAddress("theta2", &theta2, &b_theta2);
   fChain->SetBranchAddress("x2", &x2, &b_x2);
   fChain->SetBranchAddress("y2", &y2, &b_y2);
   fChain->SetBranchAddress("z2", &z2, &b_z2);
   fChain->SetBranchAddress("chi22", &chi22, &b_chi22);
   fChain->SetBranchAddress("nHits2", &nHits2, &b_nHits2);
   fChain->SetBranchAddress("beforeHits2", &beforeHits2, &b_beforeHits2);
   fChain->SetBranchAddress("npHits2", &npHits2, &b_npHits2);
   fChain->SetBranchAddress("nsHits2", &nsHits2, &b_nsHits2);
   fChain->SetBranchAddress("missHits2", &missHits2, &b_missHits2);
   fChain->SetBranchAddress("lambdaError2", &lambdaError2, &b_lambdaError2);
   fChain->SetBranchAddress("r_pt2", &r_pt2, &b_r_pt2);
   fChain->SetBranchAddress("r_d02", &r_d02, &b_r_d02);
   fChain->SetBranchAddress("r_phi2", &r_phi2, &b_r_phi2);
   fChain->SetBranchAddress("r_theta2", &r_theta2, &b_r_theta2);
   fChain->SetBranchAddress("r_px2", &r_px2, &b_r_px2);
   fChain->SetBranchAddress("r_py2", &r_py2, &b_r_py2);
   fChain->SetBranchAddress("r_pz2", &r_pz2, &b_r_pz2);
   fChain->SetBranchAddress("ipx2", &ipx2, &b_ipx2);
   fChain->SetBranchAddress("ipy2", &ipy2, &b_ipy2);
   fChain->SetBranchAddress("ipz2", &ipz2, &b_ipz2);
   fChain->SetBranchAddress("ix2", &ix2, &b_ix2);
   fChain->SetBranchAddress("iy2", &iy2, &b_iy2);
   fChain->SetBranchAddress("iz2", &iz2, &b_iz2);
   fChain->SetBranchAddress("ox2", &ox2, &b_ox2);
   fChain->SetBranchAddress("oy2", &oy2, &b_oy2);
   fChain->SetBranchAddress("oz2", &oz2, &b_oz2);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("mass", &mass, &b_mass);
   fChain->SetBranchAddress("r_pt", &r_pt, &b_r_pt);
   fChain->SetBranchAddress("r_phi", &r_phi, &b_r_phi);
   fChain->SetBranchAddress("r_theta", &r_theta, &b_r_theta);
   fChain->SetBranchAddress("r_mass", &r_mass, &b_r_mass);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("deltapt", &deltapt, &b_deltapt);
   fChain->SetBranchAddress("deltaphi", &deltaphi, &b_deltaphi);
   fChain->SetBranchAddress("deltatheta", &deltatheta, &b_deltatheta);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("deltax", &deltax, &b_deltax);
   fChain->SetBranchAddress("deltay", &deltay, &b_deltay);
   fChain->SetBranchAddress("deltaz", &deltaz, &b_deltaz);
   fChain->SetBranchAddress("minapp", &minapp, &b_minapp);
   fChain->SetBranchAddress("refit", &refit, &b_refit);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("tx", &tx, &b_tx);
   fChain->SetBranchAddress("ty", &ty, &b_ty);
   fChain->SetBranchAddress("tz", &tz, &b_tz);
   fChain->SetBranchAddress("bsx", &bsx, &b_bsx);
   fChain->SetBranchAddress("bsy", &bsy, &b_bsy);
   fChain->SetBranchAddress("bsz", &bsz, &b_bsz);
   fChain->SetBranchAddress("chi2prob", &chi2prob, &b_chi2prob);
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

Bool_t convR2S::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void convR2S::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t convR2S::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef convR2S_cxx
