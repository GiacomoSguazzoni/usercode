#ifndef GeoCut_h
#define GeoCut_h

#include <iostream>

#include <TROOT.h>

class GeoCut {
public :

  Double_t zMin;
  Double_t zMax;
  Double_t rMin;
  Double_t rMax;

  GeoCut(Double_t, Double_t, Double_t, Double_t);
 ~GeoCut();

  Double_t GetZMin(){ return zMin; };
  Double_t GetZMax(){ return zMax; };
  Double_t GetRMin(){ return rMin; };
  Double_t GetRMax(){ return rMax; };

  Int_t GeoCutOk(Double_t, Double_t);

};

#endif

#ifdef GeoCut_cxx
GeoCut::GeoCut(Double_t zmin, Double_t zmax, Double_t rmin, Double_t rmax)
{

  zMin = zmin;
  zMax = zmax;
  rMin = rmin;
  rMax = rmax;

}

GeoCut::~GeoCut()
{

}

#endif // #ifdef GeoCut_cxx
