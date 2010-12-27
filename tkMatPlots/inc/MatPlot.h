#ifndef MatPlot_h
#define MatPlot_h

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>

#include <TROOT.h>
#include <TNamed.h>
#include <TH1D.h>
#include <TH2D.h>

#ifdef UNFOLD
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>
#endif

#include "convR2SforMatPlot.h"
#include "convS2RforMatPlot.h"
#include "niR2SforMatPlot.h"
#include "niS2RforMatPlot.h"

#include "EffVsRadius.h"
#include "GeoCut.h"

class MatPlot {
public :

  TString TSName;

  //2D

  TH2D * rawSim2DH;
  TH2D * rawMC2DH;
  TH2D * rawMCFs2DH;
  TH2D * rawFake2DH;
  TH2D * rawData2DH;
  TH2D * rawDataFs2DH;

  TH2D * Sim2DH;
  TH2D * MC2DH;
  TH2D * MCFs2DH;
  TH2D * Fake2DH;
  TH2D * Data2DH;
  TH2D * DataFs2DH;

#ifdef UNFOLD
  TH2D * trainTrue2DH;
  TH2D * train2DH;

  TH2D * unfoldBay1MCFs2DH;
  TH2D * unfoldBay2MCFs2DH;
  TH2D * unfoldBay3MCFs2DH;
  TH2D * unfoldBay4MCFs2DH;

  TH2D * unfoldSVD1MCFs2DH;
  TH2D * unfoldSVD2MCFs2DH;
  TH2D * unfoldSVD3MCFs2DH;
  TH2D * unfoldSVD4MCFs2DH;

  TH2D * unfoldBay1DataFs2DH;
  TH2D * unfoldBay2DataFs2DH;
  TH2D * unfoldBay3DataFs2DH;
  TH2D * unfoldBay4DataFs2DH;

  TH2D * unfoldSVD1DataFs2DH;
  TH2D * unfoldSVD2DataFs2DH;
  TH2D * unfoldSVD3DataFs2DH;
  TH2D * unfoldSVD4DataFs2DH;
#endif

  //1D

  TH1D * rawSim1DH;
  TH1D * rawMC1DH;
  TH1D * rawMCFs1DH;
  TH1D * rawFake1DH;
  TH1D * rawData1DH;
  TH1D * rawDataFs1DH;

  TH1D * Sim1DH;
  TH1D * MC1DH;
  TH1D * MCFs1DH;
  TH1D * Fake1DH;
  TH1D * Data1DH;
  TH1D * DataFs1DH;

#ifdef UNFOLD
  TH1D * unfoldBay1MCFs1DH;
  TH1D * unfoldBay2MCFs1DH;
  TH1D * unfoldBay3MCFs1DH;
  TH1D * unfoldBay4MCFs1DH;

  TH1D * unfoldSVD1MCFs1DH;
  TH1D * unfoldSVD2MCFs1DH;
  TH1D * unfoldSVD3MCFs1DH;
  TH1D * unfoldSVD4MCFs1DH;

  TH1D * unfoldBay1DataFs1DH;
  TH1D * unfoldBay2DataFs1DH;
  TH1D * unfoldBay3DataFs1DH;
  TH1D * unfoldBay4DataFs1DH;

  TH1D * unfoldSVD1DataFs1DH;
  TH1D * unfoldSVD2DataFs1DH;
  TH1D * unfoldSVD3DataFs1DH;
  TH1D * unfoldSVD4DataFs1DH;
#endif

  //

  Float_t MCScaleFact;
  Float_t SimScaleFact;
  Float_t DataScaleFact;

  //

#ifdef UNFOLD
  RooUnfoldResponse * response;
  RooUnfoldBayes * unfoldBay;
  RooUnfoldBayes * unfoldSVD;
#endif

  EffVsRadius* effRadius;
  void SetEffRadii(std::vector<Double_t>* radii){ effRadius = new EffVsRadius(radii);};

  std::vector<GeoCut>* geoCuts;
  void SetGeoCuts(std::vector<GeoCut>* gcut){ geoCuts = gcut; };

  //
  
  Int_t uIndex;
  Int_t vIndex;
  
  void SetUIndex(Int_t uI){ uIndex = uI; };
  void SetVIndex(Int_t vI){ vIndex = vI; };

  //

  void test(Double_t val = 0.);

  MatPlot(const char*, Double_t, Double_t, Double_t, Double_t);
  MatPlot(const char*, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
  ~MatPlot();
  void Plot(TH2*);
  void Plot(TH1*);
  //
  void WBScale1D(TH1*);
  void EffScale1D(TH1*);
  //
  void WBScale2D(TH2*);
  void EffScale2D(TH2*);
  //
  void WBScaleRZ(TH2*);
  void EffScaleRZ(TH2*);
  //
  void WBScaleXY(TH2*);
  void EffScaleXY(TH2*);
  //
  void WBScaleR(TH1*);
  void EffScaleR(TH1*);
  //
  void WBScalePhi(TH1*);
  void EffScalePhi(TH1*);
  //
  Double_t WBWeiCompDxDyDz(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
  Double_t WBWeiCompDrDzDphi(Double_t, Double_t, Double_t, Double_t, Double_t phi1=-3.141592653589793, Double_t phi2=3.141592653589793);
  void PlotAll();
  //
  void Fill(convR2SforMatPlot*, TH2*, TH2* fake = 0);
  void Fill(convS2RforMatPlot*, TH2*);
  void Fill(convR2SforMatPlot*, TH1*, TH1* fake = 0);
  void Fill(convS2RforMatPlot*, TH1*);
  void FillData(convR2SforMatPlot*, Int_t minE = 1, Int_t maxE = 999999999);
  void FillMC(convR2SforMatPlot*, Int_t minE = 1, Int_t maxE = 999999999);
  void FillSim(convS2RforMatPlot*, Int_t minE = 1, Int_t maxE = 999999999);
  //
  void Fill(niR2SforMatPlot*, TH2*, TH2* fake = 0);
  void Fill(niS2RforMatPlot*, TH2*);
  void Fill(niR2SforMatPlot*, TH1*, TH1* fake = 0);
  void Fill(niS2RforMatPlot*, TH1*);
  void FillData(niR2SforMatPlot*, Int_t minE = 1, Int_t maxE = 999999999);
  void FillMC(niR2SforMatPlot*, Int_t minE = 1, Int_t maxE = 999999999);
  void FillSim(niS2RforMatPlot*, Int_t minE = 1, Int_t maxE = 999999999);
  //
  void SetMCScaleFact(Double_t sf = 1.){ MCScaleFact = sf; };
  void SetSimScaleFact(Double_t sf = 1.){ SimScaleFact = sf; };
  void SetDataScaleFact(Double_t sf = 1.){ DataScaleFact = sf; };
  Double_t FindRMax(Double_t x1, Double_t x2, Double_t y1, Double_t y2){
    //    return 1.;
    return sqrt(std::max(fabs(x1),fabs(x2))*std::max(fabs(x1),fabs(x2))+std::max(fabs(y1),fabs(y2))*std::max(fabs(y1),fabs(y2)));
  };

  
#ifdef UNFOLD
  void responseTrain(convS2RforMatPlot*, convR2SforMatPlot*, Int_t minE = 1, Int_t maxE = 999999999);
  void responseTrain(niS2RforMatPlot*, niR2SforMatPlot*, Int_t minE = 1, Int_t maxE = 999999999);
  void doUnfold();
  TH1D * doUnfoldHist(const char *, Int_t, TH1D *);
  TH2D * doUnfoldHist(const char *, Int_t, TH2D *);
  // S2R dx = SIM - RECO ==> RECO = SIM - dx
  //
#endif

};

#endif

#ifdef MatPlot_cxx
MatPlot::MatPlot(const char* name, Double_t uMin, Double_t uMax, Double_t vMin, Double_t vMax, Double_t binW, Double_t binWSim)
{

  std::cout << " Istanzio 2D... " << name << std::endl;

  //Default
  uIndex = 1;
  vIndex = 2;
  
  MCScaleFact = 1.; 
  SimScaleFact = 1.;
  DataScaleFact = 1.; 

  Int_t uNBinSim = (uMax-uMin)/binWSim;
  Int_t vNBinSim = (vMax-vMin)/binWSim;

  Int_t uNBin = (uMax-uMin)/binW;
  Int_t vNBin = (vMax-vMin)/binW;

  TSName.Append(name);

  rawSim1DH    = 0;
  rawMC1DH     = 0;
  rawMCFs1DH   = 0;
  rawFake1DH   = 0;
  rawData1DH   = 0;
  rawDataFs1DH = 0;

  Sim1DH    = 0;
  MC1DH     = 0;
  MCFs1DH   = 0;
  Fake1DH   = 0;
  Data1DH   = 0;
  DataFs1DH = 0;

  rawSim2DH    = new TH2D("rSim_"+TSName,"rSimH" , uNBinSim, uMin, uMax, vNBinSim, vMin, vMax);
  rawMC2DH     = new TH2D("rMC_"+TSName,"rMCH"  , uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  rawMCFs2DH   = new TH2D("rMCFs_"+TSName,"rMCFsH"  , uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  rawFake2DH   = new TH2D("rFake_"+TSName,"rFakeH", uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  rawData2DH   = new TH2D("rData_"+TSName,"rDataH", uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  rawDataFs2DH = new TH2D("rDataFs_"+TSName,"rDataFsH", uNBin,    uMin, uMax, vNBin,    vMin, vMax);

  /*
  Sim2DH    = new TH2D("Sim_"+TSName,"SimH" , uNBinSim, uMin, uMax, vNBinSim, vMin, vMax);
  MC2DH     = new TH2D("MC_"+TSName,"MCH"  , uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  MCFs2DH   = new TH2D("MCFs_"+TSName,"MCFsH"  , uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  Fake2DH   = new TH2D("Fake_"+TSName,"FakeH", uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  Data2DH   = new TH2D("Data_"+TSName,"DataH", uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  DataFs2DH = new TH2D("DataFs_"+TSName,"DataFsH", uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  */

#ifdef UNFOLD
  trainTrue2DH = new TH2D("trainTrue_"+TSName,"trainTrue", uNBin, uMin, uMax, vNBin, vMin, vMax);
  train2DH = new TH2D("train_"+TSName,"train", uNBin, uMin, uMax, vNBin, vMin, vMax);

  unfoldBay1MCFs2DH = 0;
  unfoldBay2MCFs2DH = 0;
  unfoldBay3MCFs2DH = 0;
  unfoldBay4MCFs2DH = 0;

  unfoldSVD1MCFs2DH = 0;
  unfoldSVD2MCFs2DH = 0;
  unfoldSVD3MCFs2DH = 0;
  unfoldSVD4MCFs2DH = 0;

  unfoldBay1DataFs2DH = 0;
  unfoldBay2DataFs2DH = 0;
  unfoldBay3DataFs2DH = 0;
  unfoldBay4DataFs2DH = 0;

  unfoldSVD1DataFs2DH = 0;
  unfoldSVD2DataFs2DH = 0;
  unfoldSVD3DataFs2DH = 0;
  unfoldSVD4DataFs2DH = 0;
  
  response= new RooUnfoldResponse("response2D","response2D");
  //  response->Setup(train2DH, trainTrue2DH);
#endif

}

MatPlot::MatPlot(const char* name, Double_t uMin, Double_t uMax, Double_t binW, Double_t binWSim)
{

  std::cout << " Istanzio 1D... " << name << std::endl;

  //Default
  uIndex = 4;
  vIndex = 0;

  MCScaleFact = 1.; 
  SimScaleFact = 1.;
  DataScaleFact = 1.; 

  Int_t uNBinSim = (uMax-uMin)/binWSim;

  Int_t uNBin = (uMax-uMin)/binW;

  TSName.Append(name);

  rawSim2DH    = 0;
  rawMC2DH     = 0;
  rawMCFs2DH   = 0;
  rawFake2DH   = 0;
  rawData2DH   = 0;
  rawDataFs2DH = 0;

  Sim2DH    = 0;
  MC2DH     = 0;
  MCFs2DH   = 0;
  Fake2DH   = 0;
  Data2DH   = 0;
  DataFs2DH = 0;

#ifdef UNFOLD
  trainTrue2DH = 0;
  train2DH = 0;
#endif

  rawSim1DH  = new TH1D( "rSim_"+TSName,"rSimH" , uNBinSim, uMin, uMax);
  rawMC1DH   = new TH1D( "rMC_"+TSName,"rMCH"  , uNBin,    uMin, uMax);
  rawMCFs1DH = new TH1D(  "rMCFs_"+TSName,"rMCFsH"  , uNBin,    uMin, uMax);
  rawFake1DH = new TH1D("rFake_"+TSName,"rFakeH", uNBin,    uMin, uMax);
  rawData1DH = new TH1D("rData_"+TSName,"rDataH", uNBin,    uMin, uMax);
  rawDataFs1DH = new TH1D("rDataFs_"+TSName,"rDataFsH", uNBin,    uMin, uMax);

  /*
  Sim1DH  = new TH1D("Sim_"+TSName,"SimH" , uNBinSim, uMin, uMax);
  MC1DH   = new TH1D("MC_"+TSName,"MCH"  , uNBin,    uMin, uMax);
  MCFs1DH = new TH1D("MCFs_"+TSName,"MCFsH"  , uNBin,    uMin, uMax);
  Fake1DH = new TH1D("Fake_"+TSName,"FakeH", uNBin,    uMin, uMax);
  Data1DH = new TH1D("Data_"+TSName,"DataH", uNBin,    uMin, uMax);
  DataFs1DH = new TH1D("DataFs_"+TSName,"DataFsH", uNBin,    uMin, uMax);
  */

#ifdef UNFOLD
  unfoldBay1MCFs1DH = 0;
  unfoldBay2MCFs1DH = 0;
  unfoldBay3MCFs1DH = 0;
  unfoldBay4MCFs1DH = 0;

  unfoldSVD1MCFs1DH = 0;
  unfoldSVD2MCFs1DH = 0;
  unfoldSVD3MCFs1DH = 0;
  unfoldSVD4MCFs1DH = 0;

  unfoldBay1DataFs1DH = 0;
  unfoldBay2DataFs1DH = 0;
  unfoldBay3DataFs1DH = 0;
  unfoldBay4DataFs1DH = 0;

  unfoldSVD1DataFs1DH = 0;
  unfoldSVD2DataFs1DH = 0;
  unfoldSVD3DataFs1DH = 0;
  unfoldSVD4DataFs1DH = 0;
  
  response = new RooUnfoldResponse(uNBin, uMin, uMax);
#endif

}

MatPlot::~MatPlot()
{

  std::cout << " Distruggo..." << std::endl;

}

#endif // #ifdef MatPlot_cxx
