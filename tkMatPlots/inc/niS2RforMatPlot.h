#ifndef niS2RforMatPlot_h
#define niS2RforMatPlot_h

#include <iostream>
#include <vector>
#include <math.h>

#include <TH2.h>
#include <TString.h>

#include <RooUnfoldResponse.h>

#include "EffVsRadius.h"
#include "GeoCut.h"

#include "niS2R.h"

class niS2RforMatPlot : public niS2R {
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
  Int_t QualityCut();

  niS2RforMatPlot(const char* filename);
  niS2RforMatPlot(TTree* tree);
  ~niS2RforMatPlot();

};

#endif

#ifdef niS2RforMatPlot_cxx
niS2RforMatPlot::niS2RforMatPlot(const char* filename)
{
  
  std::cout << " Istanzio niS2RforMatPlot..." << std::endl;
  std::cout << " Cerco di aprire questo file: " << filename << std::endl;

  x0=0.;
  y0=0.;
  z0=0.;

  TFile * f = new TFile(filename);
  
  TTree * tree = (TTree*)f->Get("ntupleS2R");

  Init(tree);

}

niS2RforMatPlot::niS2RforMatPlot(TTree* tree)
{
  
  std::cout << " Istanzio niS2RforMatPlot tramite un tree..." << std::endl;

  x0=0.;
  y0=0.;
  z0=0.;

  Init(tree);

}

niS2RforMatPlot::~niS2RforMatPlot()
{

  std::cout << " Distruggo niS2RforMatPlot..." << std::endl;

}

Int_t niS2RforMatPlot::QualityCut()
{

  Int_t iCut = 0;

#include "NiSimQualityCut.cxx"

  return iCut;

}

#endif // #ifdef niS2RforMatPlot_cxx
