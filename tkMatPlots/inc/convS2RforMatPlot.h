#ifndef convS2RforMatPlot_h
#define convS2RforMatPlot_h

#include <iostream>
#include <vector>
#include <math.h>

#include <TH2.h>
#include <TString.h>

#ifdef UNFOLD
#include <RooUnfoldResponse.h>
#endif //#ifdef UNFOLD

#include "EffVsUV.h"
#include "GeoCut.h"

#include "convS2R.h"

class convS2RforMatPlot : public convS2R {
public :

  std::vector<GeoCut>* geoCuts;

  EffVsUV* effUV;

  Int_t evRangeMin;
  Int_t evRangeMax;

  Double_t x0;
  Double_t y0;
  Double_t z0;

  //
  
  Int_t uIndex;
  Int_t vIndex;
  Int_t uCutIndex;
  Int_t vCutIndex;
  
  void SetUIndex(Int_t uI, Int_t uCutI){ uIndex = uI; uCutIndex = uCutI; };
  void SetVIndex(Int_t vI, Int_t vCutI){ vIndex = vI; vCutIndex = vCutI; };

  //

  void SetEvRangeMin(Int_t val){ evRangeMin = val; };
  void SetEvRangeMax(Int_t val){ evRangeMax = val; };
  void SetGeoCuts(std::vector<GeoCut>* gcut, EffVsUV* eff){ geoCuts = gcut; effUV = eff;};
  void LoopForFill(TH1*);
  void LoopForFill(TH2*);
#ifdef UNFOLD
  void LoopForTrain(RooUnfoldResponse*);
#endif //#ifdef UNFOLD
  Int_t QualityCut() { return 1; };

  convS2RforMatPlot(const char* filename);
  convS2RforMatPlot(TTree* tree);
  ~convS2RforMatPlot();

};

#endif

#ifdef convS2RforMatPlot_cxx
convS2RforMatPlot::convS2RforMatPlot(const char* filename)
{
  
  std::cout << " Istanzio convS2RforMatPlot..." << std::endl;
  std::cout << " Cerco di aprire questo file: " << filename << std::endl;

  x0=0.;
  y0=0.;
  z0=0.;

  TFile * f = new TFile(filename);
  
  TTree * tree = (TTree*)f->Get("ntupleS2R");

  Init(tree);

}

convS2RforMatPlot::convS2RforMatPlot(TTree* tree)
{
  
  std::cout << " Istanzio convS2RforMatPlot tramite un tree..." << std::endl;

  x0=0.;
  y0=0.;
  z0=0.;

  Init(tree);

  fChain->SetBranchStatus("event",1);
  fChain->SetBranchStatus("isAssoc",1);
  fChain->SetBranchStatus("x",1);
  fChain->SetBranchStatus("y",1);
  fChain->SetBranchStatus("z",1);

}

convS2RforMatPlot::~convS2RforMatPlot()
{

  std::cout << " Distruggo convS2RforMatPlot..." << std::endl;

}

#endif // #ifdef convS2RforMatPlot_cxx
