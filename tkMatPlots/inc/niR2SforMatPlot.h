#ifndef niR2SforMatPlot_h
#define niR2SforMatPlot_h

#include <iostream>
#include <vector>
#include <math.h>

#include <TH2.h>
#include <TString.h>
#include <TVector3.h>
#include <Math/VectorUtil.h>

#ifdef UNFOLD
#include <RooUnfoldResponse.h>
#endif //#ifdef UNFOLD

#include "EffVsUV.h"
#include "GeoCut.h"

#include "NtupleReaderNuclearInteractions.h"

class niR2SforMatPlot : public NtupleReaderNuclearInteractions {
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
  //  void LoopForFill(TH2*, TH2*, TH2*);
#ifdef UNFOLD
  void LoopForTrain(RooUnfoldResponse*);
#endif //#ifdef UNFOLD
  Int_t QualityCut(int);
  void SetCenterCoord(Double_t XX, Double_t YY, Double_t ZZ){ x0 = XX; y0 = YY; z0 = ZZ;};

  //  niR2SforMatPlot(const char* filename);
  niR2SforMatPlot(TTree* tree);
  ~niR2SforMatPlot();

};

#endif

#ifdef niR2SforMatPlot_cxx
/*
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
*/

niR2SforMatPlot::niR2SforMatPlot(TTree *tree)
{
  
  std::cout << " Istanzio niR2SforMatPlot tramite un tree..." << std::endl;

  x0=0.;
  y0=0.;
  z0=0.;
  
  beginJob(tree);

  //  Init(tree);

}

niR2SforMatPlot::~niR2SforMatPlot()
{

  std::cout << " Distruggo niR2SforMatPlot..." << std::endl;

}

Int_t niR2SforMatPlot::QualityCut(int i)
{

  Int_t iCut = 0;

#include "NiQualityCut.cxx"

  return iCut;

}

#endif // #ifdef niR2SforMatPlot_cxx
