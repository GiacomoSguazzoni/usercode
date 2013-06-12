#define GeoCut_cxx
#include "GeoCut.h"

Int_t GeoCut::GeoCutOk(Double_t u, Double_t v){

  //  std::cout << " uMin:" << uMin << " u:" << u << " uMax:" << uMax << " || vMin:" << vMin << " v:" << v << " vMax:" << vMax << std::endl;  

  if ( uMin<u && uMax>u && vMin<v && vMax>v ) return 1;
  
  return 0;

}
