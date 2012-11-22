#ifndef elPlot_h
#define elPlot_h

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

#include "binHistos.h"
#include "dpdxGraph.h"

class elPlot {
 public :
  
  TString TSName;
  
  int uNBin, vNBin;
  double uMin, uMax, vMin, vMax, uBinW, vBinW;
  
  //Vector of bin histograms
  std::vector<binHistos*> bHVector;
  //

  elPlot(const char*, double, double, int, double, double, int);
  elPlot(const char*);
  ~elPlot();
  void PlotAndWriteAll();
  void FitAll();
  void Fill(double, double, double, double, double);

 private:

  void otherInits();

  TFile * theFile;
  int FindBinHisto(double, double);

  /*
  std::vector<double> dPdxMeasInBinVect;
  std::vector<double> dPdxGausInBinVect;
  std::vector<double> dPdxSimInBinVect;
  */

  dpdxGraph * measGraph;
  dpdxGraph * gausGraph;
  dpdxGraph * meanGraph;
  dpdxGraph * medianGraph;
  dpdxGraph * meanTrunGraph;

  dpdxGraph * gausVsEtaGraph;
  dpdxGraph * meanVsEtaGraph;
  dpdxGraph * medianVsEtaGraph;
  dpdxGraph * simVsEtaGraph;

  dpdxGraph * gausVsPhiGraph;
  dpdxGraph * simVsPhiGraph;

  int nBins;

};

#endif
