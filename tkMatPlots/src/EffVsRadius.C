#define EffVsRadius_cxx
#include "EffVsRadius.h"

Int_t EffVsRadius::radiusBin(Double_t radius){
  
  std::vector<Double_t>::iterator rIt;
  Int_t nBin = 0;

  for ( rIt=rBoundaries->begin() ; rIt < rBoundaries->end()-1; rIt++ ){
    //    std::cout << " " << nBin << " radius " << radius << " rmin rmax " << *(rIt) << " " << *(rIt+1) << std::endl;

    if ( radius < *(rIt+1) ) return nBin;
    nBin++;
    
  }
  
  return nBin-1;
    
}

void EffVsRadius::Succeed(Double_t rTrue, Double_t rMeas){

  Int_t nTrue = radiusBin(rTrue);
  Int_t nMeas = radiusBin(rMeas);

  if ( nTrue == nMeas ) {
    succeed.at(nTrue)++;
  } else {
    //    fail.at(nTrue)++;   // E` questo il problema?!?!?
    xtalk.at(nMeas)++;
  }
}

void EffVsRadius::Fail(Double_t rTrue){

  fail.at(radiusBin(rTrue))++;

}

void EffVsRadius::Cand(Double_t rMeas){

  cand.at(radiusBin(rMeas))++;

}

void EffVsRadius::Fake(Double_t rMeas){

  fake.at(radiusBin(rMeas))++;

}

void EffVsRadius::Xtalk(Double_t rMeas){

  xtalk.at(radiusBin(rMeas))++;

}

Double_t EffVsRadius::EffForRadius(Double_t radius){

  Int_t nBin = radiusBin(radius);
  return EffForBin(nBin);

}

Double_t EffVsRadius::EffErrForRadius(Double_t radius){

  Int_t nBin = radiusBin(radius);
  return EffErrForBin(nBin);

}

Double_t EffVsRadius::EffForBin(Int_t nBin){

  Double_t suc = 1.*succeed.at(nBin);
  Double_t fai = 1.*fail.at(nBin);
  Double_t xt = 1.*xtalk.at(nBin);
  if (suc) {
    return 1.*(suc+xt)/(suc+xt+fai);
  } else {
    return 0.;
  }

}

#include <math.h>

Double_t EffVsRadius::EffErrForBin(Int_t nBin){

  Double_t eff = EffForBin(nBin);
  Double_t suc = 1.*succeed.at(nBin);
  Double_t fai = 1.*fail.at(nBin);
  Double_t xt = 1.*xtalk.at(nBin);

  if (suc+xt+fai) {
    return sqrt(eff*(1.-eff)/(suc+xt+fai));
  } else {
    return 0.;
  }

}

void EffVsRadius::DumpEfficiency(){
  //
  std::vector<Double_t>::iterator rIt;
  Int_t nBin = 0;
  for ( rIt=rBoundaries->begin() ; rIt < rBoundaries->end()-1; rIt++ ){
    Double_t eff = EffForBin(nBin);
    Double_t eeff = EffErrForBin(nBin);
    std::cout << " " << nBin << " rMin " << *(rIt) << " rMax " << *(rIt+1) << " eff " << eff << " +- " << eeff << std::endl;
    nBin++;
  }
}

void EffVsRadius::DumpRawNumbers(){

  std::vector<Double_t>::iterator rIt;
  Int_t nBin = 0;
  Int_t mcCand = 0;
  Int_t Cand = 0;
  for ( rIt=rBoundaries->begin(); rIt < rBoundaries->end()-1; rIt++ ){
    mcCand += succeed.at(nBin) + xtalk.at(nBin) + fake.at(nBin);;
    Cand += cand.at(nBin);
    nBin++;
  }

  Float_t norm = Cand/(1.*mcCand);

  nBin = 0;
  for ( rIt=rBoundaries->begin(); rIt < rBoundaries->end()-1; rIt++ ){
    //
    //
    Float_t suc = 1.*succeed.at(nBin);
    Float_t fai = 1.*fail.at(nBin);
    Float_t can = 1.*cand.at(nBin);
    Float_t fak = 1.*fake.at(nBin);
    Float_t xt = 1.*xtalk.at(nBin);

    //
    //    Float_t ratio = (norm*suc)/(can-norm*(xt+fak));
    Float_t rawratio = can/(norm*(suc+fak+xt));
    Float_t erawratio = rawratio*sqrt(1/can + 1/(suc+fak+xt));
    Float_t ratio = (can-norm*(xt+fak))/(norm*suc);
    Float_t eratioOverRatio2 = 1./suc + ( can/norm/norm + xt + fak )/( can/norm - xt - fak )/( can/norm - xt - fak );
    Float_t eratio = ratio*sqrt(eratioOverRatio2);

    std::cout << " " << nBin << " rMin " << *(rIt) << " rMax " << *(rIt+1) << 
      " norm " << norm <<
      " suc " << suc <<
      " fai " << fai <<
      " can " << can <<
      " fak " << fak <<
      " xt " << xt <<
      " rawr " << rawratio <<
      " erawr " << erawratio <<
      " r " << ratio <<
      " er " << eratio <<
      std::endl;
    nBin++;
  }

}

