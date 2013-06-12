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

  Double_t pi = 3.141592653589793;

  std::vector<GeoCut> geoCuts;
  geoCuts.push_back(GeoCut(13., 29., 26., 50.));

  Double_t binW = 10.*pi/180.;

  MatPlot myPlot("Phi", -pi, pi, binW, binW, -1.);
  myPlot.SetGeoCuts(&geoCuts);
  myPlot.SetUIndex(6, 3); //phi, Z
  myPlot.SetVIndex(0, 4); //--, R

  //
  myPlot.SetMCScaleFact(-1.);
  myPlot.SetSimScaleFact(-1.);
  myPlot.SetDataScaleFact(1.);

  //MC Reco
  //
#include "src/files_conv_minbias_reco.cxx"
  convR2SforMatPlot convMC(mc);
  myPlot.FillMC(&convMC);

  //MC Sim
  //
#include "src/files_conv_minbias_sim.cxx"
  convS2RforMatPlot convSim(sim);
  myPlot.FillSim(&convSim, 1, 10000000);

  //Data
  //
#include "src/files_conv_run2012ABCD.cxx"
  convR2SforMatPlot convData(data);
  //  convData.SetCenterCoord(-0.1475, -0.3782, -0.4847);
  //
  // run2012D
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TkAlignmentPixelPosition 
  convData.SetCenterCoord(-0.10822, -0.376969, -0.410475);
  myPlot.FillData(&convData);

  myPlot.test();
  myPlot.PlotAll();


 
}
