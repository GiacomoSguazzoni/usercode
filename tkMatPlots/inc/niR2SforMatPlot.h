#ifndef niR2SforMatPlot_h
#define niR2SforMatPlot_h

#include <iostream>
#include <vector>
#include <math.h>

#include <TH2.h>
#include <TString.h>

#include <RooUnfoldResponse.h>

#include "EffVsRadius.h"
#include "GeoCut.h"

#include "niR2S.h"

class niR2SforMatPlot : public niR2S {
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
  void SetCenterCoord(Double_t XX, Double_t YY, Double_t ZZ){ x0 = XX; y0 = YY; z0 = ZZ;};

  niR2SforMatPlot(const char* filename);
  niR2SforMatPlot(TTree* tree);
  ~niR2SforMatPlot();

};

#endif

#ifdef niR2SforMatPlot_cxx
niR2SforMatPlot::niR2SforMatPlot(const char* filename)
{
  
  std::cout << " Istanzio niR2SforMatPlot..." << std::endl;
  std::cout << " Cerco di aprire questo file: " << filename << std::endl;

  x0=0.;
  y0=0.;
  z0=0.;

  TFile * f = new TFile(filename);
  
  TTree * tree = (TTree*)f->Get("ntupleR2S");

  Init(tree);

}

niR2SforMatPlot::niR2SforMatPlot(TTree *tree)
{
  
  std::cout << " Istanzio niR2SforMatPlot tramite un tree..." << std::endl;

  x0=0.;
  y0=0.;
  z0=0.;

  Init(tree);

}

niR2SforMatPlot::~niR2SforMatPlot()
{

  std::cout << " Distruggo niR2SforMatPlot..." << std::endl;

}

Int_t niR2SforMatPlot::QualityCut()
{

  Int_t iCut = 0;

#include "NiQualityCut.cxx"

  return iCut;

}

#endif // #ifdef niR2SforMatPlot_cxx
