#ifndef GeoCut_h
#define GeoCut_h

#include <iostream>

#include <TROOT.h>

class GeoCut {
public :

  Double_t uMin;
  Double_t uMax;
  Double_t vMin;
  Double_t vMax;

  GeoCut(Double_t, Double_t, Double_t, Double_t);
 ~GeoCut();

  Double_t GetUMin(){ return uMin; };
  Double_t GetUMax(){ return uMax; };
  Double_t GetVMin(){ return vMin; };
  Double_t GetVMax(){ return vMax; };

  Int_t GeoCutOk(Double_t, Double_t);

};

#endif

#ifdef GeoCut_cxx
GeoCut::GeoCut(Double_t umin, Double_t umax, Double_t vmin, Double_t vmax)
{

  uMin = umin;
  uMax = umax;
  vMin = vmin;
  vMax = vmax;

}

GeoCut::~GeoCut()
{

}

#endif // #ifdef GeoCut_cxx
