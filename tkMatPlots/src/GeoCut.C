#define GeoCut_cxx
#include "GeoCut.h"

Int_t GeoCut::GeoCutOk(Double_t r, Double_t z){

  if ( rMin<r && rMax>r && zMin<z && zMax>z ) return 1;
  
  return 0;

}
