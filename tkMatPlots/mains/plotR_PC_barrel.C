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

  Double_t rMax = 65.;
  std::vector<GeoCut> geoCuts;
  geoCuts.push_back(GeoCut(-26., 26., 2., 20.));
  geoCuts.push_back(GeoCut(-66., 66., 20., rMax));

  MatPlot myPlot("R", 1., rMax, 0.5, 0.5, -1.);
  myPlot.SetGeoCuts(&geoCuts);
  myPlot.SetUIndex(4, 3); //Radius for plot, z for cut
  myPlot.SetVIndex(0, 4); //Plot is 1d, radius for cut

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
