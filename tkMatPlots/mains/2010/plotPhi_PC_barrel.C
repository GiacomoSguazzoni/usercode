#include <vector>

#include "convR2SforMatPlot.h"
#include "convS2RforMatPlot.h"
#include "MatPlot.h"
#include "GeoCut.h"

//
// Normalizzazione
//
// Fattori scala radiali
//
// Plot raggio

int main(){
  
  Double_t PI=3.141592653589793;

  std::vector<Double_t> radii;
  radii.push_back(2.0); // 1.
  radii.push_back(3.0); // 3.5
  radii.push_back(6.0);
  radii.push_back(9.0);
  radii.push_back(17.0);
  radii.push_back(22.0); // 20.0
  radii.push_back(25.5);
  radii.push_back(30.0);
  radii.push_back(34.0);
  radii.push_back(38.0);
  radii.push_back(45.5); // 46.0
  radii.push_back(54.5); // 54.0
  radii.push_back(57.0);
  radii.push_back(61.5); // 61.0
  radii.push_back(65.0);

  Double_t rMax = 65.;
  std::vector<GeoCut> geoCuts;
  //
  // Multiple geocuts work only for R plot, RZ, XY plots!!!
  //
  // geoCuts.push_back(GeoCut(-26., 26., 2., 3.5)); // BP
  geoCuts.push_back(GeoCut(-26., 26., 2., 20.)); // BP+PXB
  //  geoCuts.push_back(GeoCut(-66., 66., 20., rMax));

  //  MatPlot myPlot("R", 1., rMax, 0.5, 0.5);
  MatPlot myPlot("Phi", -PI, PI, 2.*PI/10., 2.*PI/10.);
  //  MatPlot myPlot("Phi", -PI, PI, 2.*PI/10., PI/180.);
  //  MatPlot myPlot("Phi", -PI, PI, PI/360., PI/360.);
  myPlot.SetEffRadii(&radii);
  myPlot.SetGeoCuts(&geoCuts);
  myPlot.SetUIndex(6);

  //
  myPlot.SetMCScaleFact(-1.);
  myPlot.SetSimScaleFact(-1.);
  myPlot.SetDataScaleFact(1.);

  //MC Reco
  //
  //#include "files_conv_mcreco_20100619.cxx"
#include "files_conv_mcreco_20100713.cxx"
  convR2SforMatPlot convMC(convmc);
//  convR2SforMatPlot convMC("../7TeV/ntuple_conversion_CMSSW358p3_minbias7TeV_2010-06-04_v2.root");
  myPlot.FillMC(&convMC, 1, 100000000); 

  //MC Sim
  //
  //#include "files_conv_mcsim_20100619.cxx"
#include "files_conv_mcsim_20100713.cxx"
  convS2RforMatPlot convSim(convsim);
  //  convS2RforMatPlot convSim("../7TeV/ntuple_conversion_CMSSW358p3_minbias7TeV_2010-06-04.root");
  myPlot.FillSim(&convSim, 1, 100000000);

  //Data
  //
  //#include "files_conv_data_20100619.cxx"
#include "files_conv_data_20100713.cxx"
  convData.SetCenterCoord(-0.1475, -0.3782, -0.4847);
  myPlot.FillData(&convData);

  myPlot.PlotAll();
  myPlot.test();
 
}
