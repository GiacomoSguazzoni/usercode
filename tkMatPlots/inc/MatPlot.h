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

#include "EffVsUV.h"
#include "GeoCut.h"

class MatPlot {
public :

  TString TSName;

  //2D

  TH2D * rawSim2DH;              // Simulated Monte Carlo True vertices
  TH2D * rawMC2DH;               // Reco Monte Carlo vertices 
  TH2D * rawMCFakeSub2DH;        // (Reco Monte Carlo - Fake Monte Carlo)
  TH2D * rawMCFakeSubWrtSimPosition2DH; // (Reco Monte Carlo - Fake Monte Carlo) wrt to sim position
  TH2D * rawMCFake2DH;           // Fake Monte Carlo
  TH2D * rawData2DH;             // All Reco Data vertices
  TH2D * rawDataFakeSub2DH;      // All Reco Fake Subtracted vertices 

  // Categories
  // ---------------------------------
  // Succeed(RecoPos) => Within acceptance and quality (cuts), within the same bin (sim and reco position are in the same bin), and Associated
  // Xtalk(RecoPos) => Within acceptance and quality (cuts), outside the bin (sim and reco position are not in the same bin),  and Associated
  //                 // Xtalk has been implemented as XtalkU (cross talk in the U coordinate) and XtalkV (cross talk in the V coordinate)
  // Fake(RecoPos) => Within acceptance and quality (cuts) but not Associated
  // Fail(SimPos) => Outside acceptance and quality (cuts) and Associated || Sim Not Associated
  // 
  // Cand(RecoPos) => Candidates in data

  // For the new efficiency infrastructure
  // Fail(SimPos) => Outside acceptance and quality (cuts) and Associated
  // NotReco(SimPos) => Sim Not Associated //This category has been added for the new efficiency infrastructure

  // Efficiency is defined: (Succeed+Xtalk) / (Succeed+XTalk+Fail+NotReco) 

  TH2D * Sim2DH;
  TH2D * MC2DH;
  TH2D * MCFakeSub2DH;
  TH2D * MCFake2DH;
  TH2D * Data2DH;
  TH2D * DataFakeSub2DH;

#ifdef UNFOLD
  TH2D * trainTrue2DH;
  TH2D * train2DH;

  TH2D * unfoldBay1MCFakeSub2DH;
  TH2D * unfoldBay2MCFakeSub2DH;
  TH2D * unfoldBay3MCFakeSub2DH;
  TH2D * unfoldBay4MCFakeSub2DH;

  TH2D * unfoldSVD1MCFakeSub2DH;
  TH2D * unfoldSVD2MCFakeSub2DH;
  TH2D * unfoldSVD3MCFakeSub2DH;
  TH2D * unfoldSVD4MCFakeSub2DH;

  TH2D * unfoldBay1DataFakeSub2DH;
  TH2D * unfoldBay2DataFakeSub2DH;
  TH2D * unfoldBay3DataFakeSub2DH;
  TH2D * unfoldBay4DataFakeSub2DH;

  TH2D * unfoldSVD1DataFakeSub2DH;
  TH2D * unfoldSVD2DataFakeSub2DH;
  TH2D * unfoldSVD3DataFakeSub2DH;
  TH2D * unfoldSVD4DataFakeSub2DH;
#endif

  //1D

  TH1D * rawSim1DH;
  TH1D * rawMC1DH;
  TH1D * rawMCFakeSub1DH;
  TH1D * rawMCFakeSubWrtSimPosition1DH;
  TH1D * rawMCFake1DH;
  TH1D * rawData1DH;
  TH1D * rawDataFakeSub1DH;

  TH1D * Sim1DH;
  TH1D * MC1DH;
  TH1D * MCFakeSub1DH;
  TH1D * MCFake1DH;
  TH1D * Data1DH;
  TH1D * DataFakeSub1DH;

#ifdef UNFOLD
  TH1D * unfoldBay1MCFakeSub1DH;
  TH1D * unfoldBay2MCFakeSub1DH;
  TH1D * unfoldBay3MCFakeSub1DH;
  TH1D * unfoldBay4MCFakeSub1DH;

  TH1D * unfoldSVD1MCFakeSub1DH;
  TH1D * unfoldSVD2MCFakeSub1DH;
  TH1D * unfoldSVD3MCFakeSub1DH;
  TH1D * unfoldSVD4MCFakeSub1DH;

  TH1D * unfoldBay1DataFakeSub1DH;
  TH1D * unfoldBay2DataFakeSub1DH;
  TH1D * unfoldBay3DataFakeSub1DH;
  TH1D * unfoldBay4DataFakeSub1DH;

  TH1D * unfoldSVD1DataFakeSub1DH;
  TH1D * unfoldSVD2DataFakeSub1DH;
  TH1D * unfoldSVD3DataFakeSub1DH;
  TH1D * unfoldSVD4DataFakeSub1DH;
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

  EffVsUV* effUV;

  std::vector<GeoCut>* geoCuts;
  void SetGeoCuts(std::vector<GeoCut>* gcut){ geoCuts = gcut; };

  //
  
  Int_t uIndex;
  Int_t vIndex;
  Int_t uCutIndex;
  Int_t vCutIndex;
  
  void SetUIndex(Int_t uI, Int_t uCutI){ uIndex = uI; uCutIndex = uCutI; };
  void SetVIndex(Int_t vI, Int_t vCutI){ vIndex = vI; vCutIndex = vCutI; };

  //

  void test(Double_t val = 0.);

  MatPlot(const char*, Double_t, Double_t, Double_t, Double_t, Double_t);
  MatPlot(const char*, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
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
  void WBScaleXY(TH2*);
  void EffScaleUV(TH2*);
  void EffScaleUV(TH1*);
  //
  void WBScaleR(TH1*);
  void EffScaleR(TH1*);
  //
  void WBScalePhi(TH1*);
  //
  void WBScaleTheta(TH1*);
  //
  Double_t WBWeiCompDxDyDz(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
  Double_t WBWeiCompDrDzDphi(Double_t, Double_t, Double_t, Double_t, Double_t phi1=-3.141592653589793, Double_t phi2=3.141592653589793);
  Double_t WBWeiCompDrDtethaDphi(Double_t, Double_t, Double_t, Double_t, Double_t phi1=-3.141592653589793, Double_t phi2=3.141592653589793);
  void PlotAll();
  //
  void Fill(convR2SforMatPlot*, TH2*, TH2* fake = 0, TH2* histTrue = 0);
  void Fill(convS2RforMatPlot*, TH2*);
  void Fill(convR2SforMatPlot*, TH1*, TH1* fake = 0, TH1* histTrue = 0);
  void Fill(convS2RforMatPlot*, TH1*);
  void FillData(convR2SforMatPlot*, Int_t minE = 1, Int_t maxE = 0);
  void FillMC(convR2SforMatPlot*, Int_t minE = 1, Int_t maxE = 0);
  void FillSim(convS2RforMatPlot*, Int_t minE = 1, Int_t maxE = 0);
  //
  void Fill(niR2SforMatPlot*, TH2*, TH2* fake = 0, TH2* histTrue = 0);
  void Fill(niS2RforMatPlot*, TH2*);
  void Fill(niR2SforMatPlot*, TH1*, TH1* fake = 0, TH1* histTrue = 0);
  void Fill(niS2RforMatPlot*, TH1*);
  void FillData(niR2SforMatPlot*, Int_t minE = 1, Int_t maxE = 0);
  void FillMC(niR2SforMatPlot*, Int_t minE = 1, Int_t maxE = 0);
  void FillSim(niS2RforMatPlot*, Int_t minE = 1, Int_t maxE = 0);
  //
  void SetMCScaleFact(Double_t sf = 1.){ MCScaleFact = sf; };
  void SetSimScaleFact(Double_t sf = 1.){ SimScaleFact = sf; };
  void SetDataScaleFact(Double_t sf = 1.){ DataScaleFact = sf; };
  Double_t FindRMax(Double_t x1, Double_t x2, Double_t y1, Double_t y2){
    //    return 1.;
    return sqrt(std::max(fabs(x1),fabs(x2))*std::max(fabs(x1),fabs(x2))+std::max(fabs(y1),fabs(y2))*std::max(fabs(y1),fabs(y2)));
  };

  
#ifdef UNFOLD
  void responseTrain(convS2RforMatPlot*, convR2SforMatPlot*, Int_t minE = 1, Int_t maxE = 0);
  void responseTrain(niS2RforMatPlot*, niR2SforMatPlot*, Int_t minE = 1, Int_t maxE = 0);
  void doUnfold();
  TH1D * doUnfoldHist(const char *, Int_t, TH1D *);
  TH2D * doUnfoldHist(const char *, Int_t, TH2D *);
  // S2R dx = SIM - RECO ==> RECO = SIM - dx
  //
#endif

};

#endif

#ifdef MatPlot_cxx
MatPlot::MatPlot(const char* name, Double_t uMin, Double_t uMax, Double_t vMin, Double_t vMax, Double_t binW, Double_t binWSim, Double_t binWEff)
{

  std::cout << " Istanzio 2D... " << name << std::endl;

  //Default
  uIndex = 1;
  vIndex = 2;
  
  MCScaleFact = 1.; 
  SimScaleFact = 1.;
  DataScaleFact = 1.; 

  Int_t uNBin = round((uMax-uMin)/binW);
  Int_t vNBin = round((vMax-vMin)/binW);

  Int_t uNBinSim = round((uMax-uMin)/binWSim);
  Int_t vNBinSim = round((vMax-vMin)/binWSim);

  Int_t uNBinEff = round((uMax-uMin)/binWEff);
  Int_t vNBinEff = round((vMax-vMin)/binWEff);

  if ( binWEff < 0. ){
    uNBinEff = 1;
    vNBinEff = 1;
  }

  std::cout << " NBins: Reco U " << uNBin << " Reco V " << vNBin << std::endl;
  std::cout << " NBins: Sim U " << uNBinSim << " Sim V " << vNBinSim << std::endl;
  std::cout << " NBins: Eff U " << uNBinEff << " Eff V " << vNBinEff << std::endl;

  TSName.Append(name);

  rawSim1DH    = 0;
  rawMC1DH     = 0;
  rawMCFakeSub1DH   = 0;
  rawMCFakeSubWrtSimPosition1DH   = 0;
  rawMCFake1DH   = 0;
  rawData1DH   = 0;
  rawDataFakeSub1DH = 0;

  Sim1DH    = 0;
  MC1DH     = 0;
  MCFakeSub1DH   = 0;
  MCFake1DH   = 0;
  Data1DH   = 0;
  DataFakeSub1DH = 0;

  rawSim2DH    = new TH2D("rSim_"+TSName,"rSimH" , uNBinSim, uMin, uMax, vNBinSim, vMin, vMax);
  rawMC2DH     = new TH2D("rMC_"+TSName,"rMCH"  , uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  rawMCFakeSub2DH   = new TH2D("rMCFakeSub_"+TSName,"rMCFakeSubH"  , uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  rawMCFakeSubWrtSimPosition2DH   = new TH2D("rMCFakeSubWrtSimPosition_"+TSName,"rMCFakeSubWrtSimPositionH"  , uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  rawMCFake2DH   = new TH2D("rFake_"+TSName,"rFakeH", uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  rawData2DH   = new TH2D("rData_"+TSName,"rDataH", uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  rawDataFakeSub2DH = new TH2D("rDataFakeSub_"+TSName,"rDataFakeSubH", uNBin,    uMin, uMax, vNBin,    vMin, vMax);

  /*
  Sim2DH    = new TH2D("Sim_"+TSName,"SimH" , uNBinSim, uMin, uMax, vNBinSim, vMin, vMax);
  MC2DH     = new TH2D("MC_"+TSName,"MCH"  , uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  MCFakeSub2DH   = new TH2D("MCFakeSub_"+TSName,"MCFakeSubH"  , uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  MCFake2DH   = new TH2D("MCFake_"+TSName,"MCFakeH", uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  Data2DH   = new TH2D("Data_"+TSName,"DataH", uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  DataFakeSub2DH = new TH2D("DataFakeSub_"+TSName,"DataFakeSubH", uNBin,    uMin, uMax, vNBin,    vMin, vMax);
  */

  //Efficiency
  effUV = new EffVsUV(TSName, uNBinEff, uMin, uMax, vNBinEff, vMin, vMax);

#ifdef UNFOLD
  trainTrue2DH = new TH2D("trainTrue_"+TSName,"trainTrue", uNBin, uMin, uMax, vNBin, vMin, vMax);
  train2DH = new TH2D("train_"+TSName,"train", uNBin, uMin, uMax, vNBin, vMin, vMax);

  unfoldBay1MCFakeSub2DH = 0;
  unfoldBay2MCFakeSub2DH = 0;
  unfoldBay3MCFakeSub2DH = 0;
  unfoldBay4MCFakeSub2DH = 0;

  unfoldSVD1MCFakeSub2DH = 0;
  unfoldSVD2MCFakeSub2DH = 0;
  unfoldSVD3MCFakeSub2DH = 0;
  unfoldSVD4MCFakeSub2DH = 0;

  unfoldBay1DataFakeSub2DH = 0;
  unfoldBay2DataFakeSub2DH = 0;
  unfoldBay3DataFakeSub2DH = 0;
  unfoldBay4DataFakeSub2DH = 0;

  unfoldSVD1DataFakeSub2DH = 0;
  unfoldSVD2DataFakeSub2DH = 0;
  unfoldSVD3DataFakeSub2DH = 0;
  unfoldSVD4DataFakeSub2DH = 0;
  
  response= new RooUnfoldResponse("response2D","response2D");
  //  response->Setup(train2DH, trainTrue2DH);
#endif

}

MatPlot::MatPlot(const char* name, Double_t uMin, Double_t uMax, Double_t binW, Double_t binWSim, Double_t binWEff)
{

  std::cout << " Istanzio 1D... " << name << std::endl;

  //Default
  uIndex = 4;
  vIndex = 0;

  MCScaleFact = 1.; 
  SimScaleFact = 1.;
  DataScaleFact = 1.; 

  Int_t uNBin = round((uMax-uMin)/binW);
  Int_t uNBinSim = round((uMax-uMin)/binWSim);
  Int_t uNBinEff = round((uMax-uMin)/binWEff);

  if ( binWEff < 0. ){
    uNBinEff = 1;
  }

  TSName.Append(name);

  rawSim2DH    = 0;
  rawMC2DH     = 0;
  rawMCFakeSub2DH   = 0;
  rawMCFakeSubWrtSimPosition2DH   = 0;
  rawMCFake2DH   = 0;
  rawData2DH   = 0;
  rawDataFakeSub2DH = 0;

  Sim2DH    = 0;
  MC2DH     = 0;
  MCFakeSub2DH   = 0;
  MCFake2DH   = 0;
  Data2DH   = 0;
  DataFakeSub2DH = 0;

#ifdef UNFOLD
  trainTrue2DH = 0;
  train2DH = 0;
#endif

  rawSim1DH  = new TH1D( "rSim_"+TSName,"rSimH" , uNBinSim, uMin, uMax);
  rawMC1DH   = new TH1D( "rMC_"+TSName,"rMCH"  , uNBin,    uMin, uMax);
  rawMCFakeSub1DH = new TH1D(  "rMCFakeSub_"+TSName,"rMCFakeSubH"  , uNBin,    uMin, uMax);
  rawMCFakeSubWrtSimPosition1DH = new TH1D(  "rMCFakeSubWrtSimPosition_"+TSName,"rMCFakeSubWrtSimPositionH"  , uNBin,    uMin, uMax);
  rawMCFake1DH = new TH1D("rFake_"+TSName,"rFakeH", uNBin,    uMin, uMax);
  rawData1DH = new TH1D("rData_"+TSName,"rDataH", uNBin,    uMin, uMax);
  rawDataFakeSub1DH = new TH1D("rDataFakeSub_"+TSName,"rDataFakeSubH", uNBin,    uMin, uMax);

  /*
  Sim1DH  = new TH1D("Sim_"+TSName,"SimH" , uNBinSim, uMin, uMax);
  MC1DH   = new TH1D("MC_"+TSName,"MCH"  , uNBin,    uMin, uMax);
  MCFakeSub1DH = new TH1D("MCFakeSub_"+TSName,"MCFakeSubH"  , uNBin,    uMin, uMax);
  MCFake1DH = new TH1D("MCFake_"+TSName,"MCFakeH", uNBin,    uMin, uMax);
  Data1DH = new TH1D("Data_"+TSName,"DataH", uNBin,    uMin, uMax);
  DataFakeSub1DH = new TH1D("DataFakeSub_"+TSName,"DataFakeSubH", uNBin,    uMin, uMax);
  */

  //Efficiency
  effUV = new EffVsUV(TSName, uNBinEff, uMin, uMax, 1, 0., 2.);

#ifdef UNFOLD
  unfoldBay1MCFakeSub1DH = 0;
  unfoldBay2MCFakeSub1DH = 0;
  unfoldBay3MCFakeSub1DH = 0;
  unfoldBay4MCFakeSub1DH = 0;

  unfoldSVD1MCFakeSub1DH = 0;
  unfoldSVD2MCFakeSub1DH = 0;
  unfoldSVD3MCFakeSub1DH = 0;
  unfoldSVD4MCFakeSub1DH = 0;

  unfoldBay1DataFakeSub1DH = 0;
  unfoldBay2DataFakeSub1DH = 0;
  unfoldBay3DataFakeSub1DH = 0;
  unfoldBay4DataFakeSub1DH = 0;

  unfoldSVD1DataFakeSub1DH = 0;
  unfoldSVD2DataFakeSub1DH = 0;
  unfoldSVD3DataFakeSub1DH = 0;
  unfoldSVD4DataFakeSub1DH = 0;
  
  response = new RooUnfoldResponse(uNBin, uMin, uMax);
#endif

}

MatPlot::~MatPlot()
{

  std::cout << " Distruggo..." << std::endl;

}

#endif // #ifdef MatPlot_cxx
