#ifndef EffVsRadius_h
#define EffVsRadius_h

#include <iostream>
#include <vector>

#include <TROOT.h>

class EffVsRadius {
public :

  std::vector<Double_t>* rBoundaries;
  
  std::vector<Int_t> succeed;
  std::vector<Int_t> xtalk;
  std::vector<Int_t> fail;
  std::vector<Int_t> cand;
  std::vector<Int_t> fake;

  EffVsRadius(std::vector<Double_t>*);
  ~EffVsRadius();
  
  Int_t radiusBin(Double_t);

  void Succeed(Double_t, Double_t);
  void Fail(Double_t);
  void Cand(Double_t);
  void Fake(Double_t);
  void Xtalk(Double_t);

  Double_t EffForRadius(Double_t);
  Double_t EffErrForRadius(Double_t);
  Double_t EffForBin(Int_t);
  Double_t EffErrForBin(Int_t);
  void DumpEfficiency();
  void DumpRawNumbers();

};

#endif

#ifdef EffVsRadius_cxx
EffVsRadius::EffVsRadius(std::vector<Double_t>* rb)
{

  rBoundaries = rb;

  succeed.resize(rBoundaries->size()-1,0);
  fail.resize(rBoundaries->size()-1,0);
  cand.resize(rBoundaries->size()-1,0);
  fake.resize(rBoundaries->size()-1,0);
  xtalk.resize(rBoundaries->size()-1,0);
  
}

EffVsRadius::~EffVsRadius()
{

}

#endif // #ifdef EffVsRadius_cxx
