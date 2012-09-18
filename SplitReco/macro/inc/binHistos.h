#ifndef binHistos_H
#define binHistos_H

#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TFile.h>
#include <RooRealVar.h>

class binHistos {
 public:

  binHistos(const char*, const char*, double, double);
  binHistos(const char*, TFile*);
  ~binHistos();
  void Fill(double, double, double);
  void PlotAll();
  void FitAll(double&, double&, double&, double&, double&, double&);
  void WriteAll(TFile*);
  TString getName(){ return TSId; };

  TH1D * getHistoFromFile(TString, TFile*);
  double GetBinUCenter(){return uBinCenter;};
  double GetBinVCenter(){return vBinCenter;};

 private:
  void Plot(TH1*);
  void Plot(TH2*);
  double FitMeas(TH1*);
  double FitGaus(TH1*);
  double FitRoot(TH1*);
  double FitSim(TH1*);
  void Write(TH1*,TFile*);

  double setMyRange(TH1*, RooRealVar*);
  void estMyValues(TH1*, double&, double&);
  double Median(const TH1D*);
  double meanInRange(const TH1D*, double&, double&);

  TString TSId;

  //Histograms
  TH1D* dPdxSim;
  TH1D* dPdxMeas;
  TH1D* dPdxPull;
  TH2D* dPdxMeasVSdPdxSim;

  double uBinCenter;
  double vBinCenter;

};

#endif 
