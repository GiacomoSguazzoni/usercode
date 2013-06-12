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

  std::vector<GeoCut> geoCuts;
  geoCuts.push_back(GeoCut(-26., 26., 1., 20.));
  geoCuts.push_back(GeoCut(-66., 66., 20., 65.));

  Double_t boxSize = 60.;
  //MatPlot myPlot("XY", -boxSize, boxSize, -boxSize, boxSize, 1., 0.1);
  MatPlot myPlot("XY", -boxSize, boxSize, -boxSize, boxSize, 0.5, 0.1);
  myPlot.SetEffRadii(&radii);
  myPlot.SetGeoCuts(&geoCuts);

  //
  myPlot.SetMCScaleFact(-1.);
  myPlot.SetSimScaleFact(-1.);
  myPlot.SetDataScaleFact(1.);

  //MC Reco
  //
  //#include "files_conv_20100619.cxx"
  //#include "files_conv_mcreco_20100713.cxx"
#include "src/files_conv_mcreco_386p1_20101124.cxx"
  convR2SforMatPlot convMC(mc);
  myPlot.FillMC(&convMC); // All events

  //MC Sim
  //
  //#include "files_conv_20100619.cxx"
  //#include "files_conv_mcsim_20100713.cxx"
#include "src/files_conv_mcsim_386p1_20101124.cxx"
  convS2RforMatPlot convSim(sim);
  myPlot.FillSim(&convSim); // All events

  //Data
  //
  //#include "files_conv_20100619.cxx"
  //#include "files_conv_data_20100713.cxx"
#include "src/files_conv_data_386p1_20101124.cxx"
  convR2SforMatPlot convData(data);
  convData.SetCenterCoord(-0.1475, -0.3782, -0.4847);
  myPlot.FillData(&convData);

  myPlot.test();
  myPlot.PlotAll();

}
