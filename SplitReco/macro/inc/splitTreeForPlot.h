#ifndef splitTreeForPlot_h
#define splitTreeForPlot_h

#include <iostream>
#include <vector>
#include <math.h>

#include <TH2.h>
#include <TString.h>

#include "sT.h"
#include "elPlot.h"

class splitTreeForPlot : public sT {
public :

  int evRangeMin;
  int evRangeMax;

  void SetEvRangeMax(int val){ evRangeMax = val; };
  void LoopForFill(elPlot*);
  int QualityCut();
  double normalizedPhi(double phi);
  double sguazzNormalizedPhi(double phi);

  splitTreeForPlot(const char* filename);
  splitTreeForPlot(TTree *tree);
  ~splitTreeForPlot();

};

#endif // #ifdef splitTreeForPlot
