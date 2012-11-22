#include <cmath>
#include <iostream>
#include "TROOT.h"
#include "TString.h"
#include "TH1F.h"
#include "TChain.h"
#include "TCanvas.h"

int studyCut(TChain* chain, TString cut, int nBefore, int nmax){

  TH1F * h = new TH1F("h","h",1,0.,2.);
  
  chain->Draw("1>>h",cut,"goff",nmax);
  int n = h->GetEntries();
  std::cout << ">>" << n << " events (" << 1.*n/(1.*nBefore) << " efficiency) for the cut " << cut << std::endl;

  delete h;

  return n;

}

void cutFlow(){

  gROOT->SetStyle("Plain");
  
  TChain *uno = new TChain("sT");
  uno->Add("/afs/cern.ch/user/s/sguazzon/myWorkarea/split/minBias_nH4_v2.root"); 
  //  uno->Add("/afs/cern.ch/user/s/sguazzon/myWorkarea/split/minbias_run2012A_nH4_v2.root"); 
  //  uno->Add("/afs/cern.ch/user/s/sguazzon/myWorkarea/split/qcd1530_nH4_v2.root"); 

  std::cout << " Sample:"  << uno->GetEntries() << std::endl;
    
  // Cuts
  std::vector<TString> cuts;
  cuts.push_back(TString("pt>0.9&&pt<1.1"));
  cuts.push_back(TString("rIn<8.&&((rOut>100.)||(zOut>260.)||(zOut<-260.))"));
  cuts.push_back(TString("(chi2Split/(NSplit-2.))>0."));
  cuts.push_back(TString("dz<5.&&dz>-5."));
  cuts.push_back(TString("chi2<2"));
  cuts.push_back(TString("maxchi2<4."));
  cuts.push_back(TString("((nHitIna+nHitMis+nHitBad)==0)"));
  cuts.push_back(TString("(chi2Split/(NSplit-2.))<4."));

  int maxentries = 1000000;
  int nbefore = maxentries;
  TString cut("");

  std::vector<TString>::iterator it;

  for ( it=cuts.begin() ; it < cuts.end(); it++ ){

    cut = cut + *it;
    nbefore = studyCut(uno, cut, nbefore, maxentries);
    cut = cut + "&&";
    

  }

}
