#ifndef niS2RforMatPlot_h
#define niS2RforMatPlot_h

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

#include "NtupleReaderNuclearInteractions.h"

class niS2RforMatPlot : public NtupleReaderNuclearInteractions {
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
  //  void LoopForFill(TH2*);
#ifdef UNFOLD
  void LoopForTrain(RooUnfoldResponse*);
#endif //#ifdef UNFOLD
  Int_t QualityCut(int);

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

  beginJob(tree);

  //  Init(tree);

}

niS2RforMatPlot::niS2RforMatPlot(TTree* tree)
{
  
  std::cout << " Istanzio niS2RforMatPlot tramite un tree..." << std::endl;

  x0=0.;
  y0=0.;
  z0=0.;

  beginJob(tree);

  //  Init(tree);

}

niS2RforMatPlot::~niS2RforMatPlot()
{

  std::cout << " Distruggo niS2RforMatPlot..." << std::endl;

}

Int_t niS2RforMatPlot::QualityCut(int i)
{

  Int_t iCut = 0;

#include "NiSimQualityCut.cxx"

  return iCut;

}

#endif // #ifdef niS2RforMatPlot_cxx
