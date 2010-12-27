#ifndef convS2RforMatPlot_h
#define convS2RforMatPlot_h

#include <iostream>
#include <vector>
#include <math.h>

#include <TH2.h>
#include <TString.h>

#include <RooUnfoldResponse.h>

#include "EffVsRadius.h"
#include "GeoCut.h"

#include "convS2R.h"

class convS2RforMatPlot : public convS2R {
public :

  std::vector<GeoCut>* geoCuts;

  EffVsRadius* effR;

  Int_t evRangeMin;
  Int_t evRangeMax;

  Double_t x0;
  Double_t y0;
  Double_t z0;

  //
  
  Int_t uIndex;
  Int_t vIndex;
  
  void SetUIndex(Int_t uI){ uIndex = uI; };
  void SetVIndex(Int_t vI){ vIndex = vI; };

  //

  void SetEvRangeMin(Int_t val){ evRangeMin = val; };
  void SetEvRangeMax(Int_t val){ evRangeMax = val; };
  void SetGeoCuts(std::vector<GeoCut>* gcut){ geoCuts = gcut; };
  void SetEffRadius(EffVsRadius* eR){ effR = eR; };
  void LoopForFill(TH1*, TH1*);
  void LoopForFill(TH2*, TH2*);
  void LoopForTrain(RooUnfoldResponse*);
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

}

convS2RforMatPlot::~convS2RforMatPlot()
{

  std::cout << " Distruggo convS2RforMatPlot..." << std::endl;

}

#endif // #ifdef convS2RforMatPlot_cxx
