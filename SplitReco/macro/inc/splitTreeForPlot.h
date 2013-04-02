#ifndef splitTreeForPlot_h
#define splitTreeForPlot_h

#include <iostream>
#include <vector>
#include <math.h>

#include <TH1F.h>
#include <TH2.h>
#include <TString.h>

#include "sT.h"
#include "elPlot.h"

class splitTreeForPlot : public sT {
public :

  int evRangeMin;
  int evRangeMax;
  int evRangeMaxNorm;

  void fillNormHisto();
  float getWeiFromNormHisto(float&);
  void SetEvRangeMax(int val){ evRangeMax = val; };
  void SetEvRangeMaxNorm(int val){ evRangeMaxNorm = val; };
  void LoopForFill(elPlot*);
  int QualityCut();
  void setDPdxQuantities();
  int AnalysisCut();
  double normalizedPhi(double phi);
  double sguazzNormalizedPhi(double phi);

  splitTreeForPlot(TTree *tree = 0, int evMax=-1, int evMaxNorm=-1);
  //  splitTreeForPlot(const char* filename);
  ~splitTreeForPlot();

  bool iMC;

  TH1F * normHisto;
  TH1F * histoPtWei;
  TH1F * histoPt2;
  
  float ptmin, ptmax;
  double dPdxTrue, dPdx, dPdxErr, dPdxChi2;

};

#endif // #ifdef splitTreeForPlot
