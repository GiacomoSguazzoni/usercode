#define EffVsUV_cxx
#include "EffVsUV.h"

void EffVsUV::Success(Double_t uTrue, Double_t vTrue, Double_t uMeas, Double_t vMeas){

  Int_t uTrueBin = Success2DH->GetXaxis()->FindBin(uTrue);
  Int_t vTrueBin = Success2DH->GetYaxis()->FindBin(vTrue);
  Int_t uMeasBin = Success2DH->GetXaxis()->FindBin(uMeas);
  Int_t vMeasBin = Success2DH->GetYaxis()->FindBin(vMeas);

  if ( uTrueBin == uMeasBin && vTrueBin == vMeasBin ) {
    Success2DH->Fill(uTrue, vTrue);
  } else {
    if ( ! (uTrueBin == uMeasBin) ) XTalkU2DH->Fill(uTrue, vTrue);
    if ( ! (vTrueBin == vMeasBin) ) XTalkV2DH->Fill(uTrue, vTrue); 
  }
}

void EffVsUV::Fail(Double_t uTrue, Double_t vTrue){

  Fail2DH->Fill(uTrue, vTrue);
  //  fail.at(radiusBin(rTrue))++;

}

void EffVsUV::NotReco(Double_t uTrue, Double_t vTrue){

  NotReco2DH->Fill(uTrue, vTrue);
  //  fail.at(radiusBin(rTrue))++;

}

void EffVsUV::Cand(Double_t uMeas, Double_t vMeas){

  Cand2DH->Fill(uMeas, vMeas);
  //  cand.at(radiusBin(rMeas))++;

}

void EffVsUV::Fake(Double_t uMeas, Double_t vMeas){

  Fake2DH->Fill(uMeas, vMeas);
  //  fake.at(radiusBin(rMeas))++;

}

Double_t EffVsUV::EffForUV(Double_t u, Double_t v){

  Int_t uBin = Success2DH->GetXaxis()->FindBin(u);
  Int_t vBin = Success2DH->GetYaxis()->FindBin(v);
  return EffForBin(uBin, vBin);

}

Double_t EffVsUV::EffErrForUV(Double_t u, Double_t v){

  Int_t uBin = Success2DH->GetXaxis()->FindBin(u);
  Int_t vBin = Success2DH->GetYaxis()->FindBin(v);
  return EffErrForBin(uBin, vBin);

}

Double_t EffVsUV::EffForBin(Int_t uBin, Int_t vBin){

  Double_t suc = 1.*Success2DH->GetBinContent(uBin, vBin);
  Double_t fai = 1.*Fail2DH->GetBinContent(uBin, vBin);
  Double_t xta = 1.*(XTalkU2DH->GetBinContent(uBin, vBin)+XTalkV2DH->GetBinContent(uBin, vBin)) ;
  Double_t nor = 1.*NotReco2DH->GetBinContent(uBin, vBin);

  Double_t den = suc+xta+fai+nor;

  if ( den ) {
    return 1.*(suc+xta)/den;
  } else {
    return 0.;
  }

}

#include <math.h>

Double_t EffVsUV::EffErrForBin(Int_t uBin, Int_t vBin){

  Double_t eff = EffForBin(uBin, vBin);
  Double_t suc = 1.*Success2DH->GetBinContent(uBin, vBin);
  Double_t fai = 1.*Fail2DH->GetBinContent(uBin, vBin);
  Double_t xta = 1.*(XTalkU2DH->GetBinContent(uBin, vBin)+XTalkV2DH->GetBinContent(uBin, vBin)) ;
  Double_t nor = 1.*NotReco2DH->GetBinContent(uBin, vBin);

  Double_t den = suc+xta+fai+nor;

  if ( den ) {
    return sqrt(eff*(1.-eff)/den);
  } else {
    return 0.;
  }

}

#include <TCanvas.h>

void EffVsUV::Plot(TH2* hist)
{

  std::cout << " E ora plotto " << hist->GetName() << std::endl;

  TCanvas * can = new TCanvas(hist->GetName(),"",1000,1000);

  hist->Draw();

  can->SaveAs(TString(hist->GetName())+".png");
  hist->SaveAs(TString(hist->GetName())+".root");

}
