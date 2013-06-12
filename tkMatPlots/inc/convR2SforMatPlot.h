#ifndef convR2SforMatPlot_h
#define convR2SforMatPlot_h

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

#include "convR2S.h"

class convR2SforMatPlot : public convR2S {
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
  void LoopForFill(TH1*, TH1*, TH1*);
  void LoopForFill(TH2*, TH2*, TH2*);
#ifdef UNFOLD
  void LoopForTrain(RooUnfoldResponse*);
#endif //#ifdef UNFOLD
  Int_t QualityCut();
  void SetCenterCoord(Double_t XX, Double_t YY, Double_t ZZ){ x0 = XX; y0 = YY; z0 = ZZ;};

  convR2SforMatPlot(const char* filename);
  convR2SforMatPlot(TTree *tree);
  ~convR2SforMatPlot();

};

#endif

#ifdef convR2SforMatPlot_cxx
convR2SforMatPlot::convR2SforMatPlot(const char* filename)
{
  
  std::cout << " Istanzio convR2SforMatPlot..." << std::endl;
  std::cout << " Cerco di aprire questo file: " << filename << std::endl;

  x0=0.;
  y0=0.;
  z0=0.;

  TFile * f = new TFile(filename);
  
  TTree * tree = (TTree*)f->Get("ntupleR2S");

  Init(tree);

}

convR2SforMatPlot::convR2SforMatPlot(TTree *tree)
{
  
  std::cout << " Istanzio convR2SforMatPlot tramite un tree..." << std::endl;

  x0=0.;
  y0=0.;
  z0=0.;

  Init(tree);

  //Select only used branches

  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("event",1);
  fChain->SetBranchStatus("isAssoc",1);
  fChain->SetBranchStatus("pt1",1);
  fChain->SetBranchStatus("pt2",1);
  fChain->SetBranchStatus("nHits1",1);
  fChain->SetBranchStatus("nHits2",1);
  fChain->SetBranchStatus("x",1);
  fChain->SetBranchStatus("y",1);
  fChain->SetBranchStatus("z",1);
  //  fChain->SetBranchStatus("deltax",1);
  //  fChain->SetBranchStatus("deltay",1);
  //  fChain->SetBranchStatus("deltaz",1);
  fChain->SetBranchStatus("minapp",1);
  fChain->SetBranchStatus("chi2prob",1);
  fChain->SetBranchStatus("q_hp",1);
  fChain->SetBranchStatus("d01",1);
  fChain->SetBranchStatus("d02",1);
  fChain->SetBranchStatus("theta1",1);
  fChain->SetBranchStatus("theta2",1);

}

convR2SforMatPlot::~convR2SforMatPlot()
{

  std::cout << " Distruggo convR2SforMatPlot..." << std::endl;

}

Int_t convR2SforMatPlot::QualityCut()
{

  Int_t iCut = 0;

#include "ConvQualityCut.cxx"

  return iCut;

}

#endif // #ifdef convR2SforMatPlot_cxx
