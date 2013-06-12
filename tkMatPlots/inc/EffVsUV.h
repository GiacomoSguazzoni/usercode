#ifndef EffVsUV_h
#define EffVsUV_h

#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TH2D.h>
#include <TString.h>

class EffVsUV {
public :

  Int_t uNBin, vNBin;
  Double_t uMin, uMax, vMin, vMax;

  TH2D * Success2DH;
  TH2D * XTalkU2DH;
  TH2D * XTalkV2DH;
  TH2D * Fail2DH;
  TH2D * Cand2DH;
  TH2D * Fake2DH;
  TH2D * NotReco2DH;
  
  EffVsUV(TString, Int_t, Double_t, Double_t, Int_t, Double_t, Double_t);
  ~EffVsUV();
  
  void Success(Double_t, Double_t, Double_t, Double_t);
  void Fail(Double_t, Double_t);
  void Cand(Double_t, Double_t);
  void Fake(Double_t, Double_t);
  void NotReco(Double_t, Double_t);

  Double_t EffForUV(Double_t, Double_t);
  Double_t EffErrForUV(Double_t, Double_t);
  Double_t EffForBin(Int_t, Int_t);
  Double_t EffErrForBin(Int_t, Int_t);

  void Plot(TH2*);

};

#endif

#ifdef EffVsUV_cxx
EffVsUV::EffVsUV(TString TSName, Int_t un, Double_t umi, Double_t uma, Int_t vn, Double_t vmi, Double_t vma)
{

  uNBin = un;
  uMin  = umi;
  uMax  = uma;
  vNBin = vn;
  vMin  = vmi;
  vMax  = vma;

  std::cout << ">>>>Building efficiency interface with U: " << uNBin << " " << uMin << " " << uMax <<
    " " << vNBin << " " << vMin << " " << vMax << std::endl;

  Success2DH = new TH2D("effSuccess_"+TSName,"effSuccessH", uNBin, uMin, uMax, vNBin, vMin, vMax);
  XTalkU2DH  = new TH2D("effXTalkU_" +TSName,"effXTalkUH" , uNBin, uMin, uMax, vNBin, vMin, vMax);
  XTalkV2DH  = new TH2D("effXTalkV_" +TSName,"effXTalkVH" , uNBin, uMin, uMax, vNBin, vMin, vMax);
  Fail2DH    = new TH2D("effFail_"   +TSName,"effFailH"   , uNBin, uMin, uMax, vNBin, vMin, vMax);
  Cand2DH    = new TH2D("effCand_"   +TSName,"effCandH"   , uNBin, uMin, uMax, vNBin, vMin, vMax);
  Fake2DH    = new TH2D("effFake_"   +TSName,"effFakeH"   , uNBin, uMin, uMax, vNBin, vMin, vMax);  
  NotReco2DH = new TH2D("effNotReco_"+TSName,"effNotRecoH", uNBin, uMin, uMax, vNBin, vMin, vMax);

}

EffVsUV::~EffVsUV()
{

  Plot(Success2DH);
  Plot(XTalkU2DH);
  Plot(XTalkV2DH); 
  Plot(Fail2DH);   
  Plot(Cand2DH);   
  Plot(Fake2DH);   
  Plot(NotReco2DH);

}

#endif // #ifdef EffVsUV_cxx
