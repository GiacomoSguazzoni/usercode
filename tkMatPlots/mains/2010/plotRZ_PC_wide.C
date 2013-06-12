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

  Double_t boxSize = 150.;
  std::vector<GeoCut> geoCuts;
  //Endcap
  geoCuts.push_back(GeoCut(-10., boxSize, 1., 60.));

  //No End Cap
  //  geoCuts.push_back(GeoCut(-32.5, 32.5, 1., 17.5));
  //  geoCuts.push_back(GeoCut(-boxSize, boxSize, 17., 54.));

  MatPlot myPlot("RZ", -10., boxSize, 0., 60., .5, 0.05);

  myPlot.SetEffRadii(&radii);
  myPlot.SetGeoCuts(&geoCuts);
  // vAr[0]=1.;
  // vAr[1]=X;
  // vAr[2]=Y;
  // vAr[3]=Z;
  // vAr[4]=sqrt(X*X+Y*Y);
  myPlot.SetUIndex(3);
  myPlot.SetVIndex(4);

  //
  myPlot.SetMCScaleFact(-1.);
  myPlot.SetSimScaleFact(-1.);
  myPlot.SetDataScaleFact(1.);

  //MC Reco
  //
  //#include "files_conv_20100619.cxx"
#include "files_conv_mcreco_20100713.cxx"
  convR2SforMatPlot convMC(convmc);
  myPlot.FillMC(&convMC); // All events

  //MC Sim
  //
  //#include "files_conv_20100619.cxx"
#include "files_conv_mcsim_20100713.cxx"
  convS2RforMatPlot convSim(convsim);
  myPlot.FillSim(&convSim); // All events

  //Data
  //
  //#include "files_conv_20100619.cxx"
#include "files_conv_data_20100713.cxx"
  convData.SetCenterCoord(-0.1475, -0.3782, -0.4847);
  myPlot.FillData(&convData);

  myPlot.test();
  myPlot.PlotAll();

}
