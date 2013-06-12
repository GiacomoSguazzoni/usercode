#include <vector>

#include "convR2SforMatPlot.h"
#include "convS2RforMatPlot.h"
#include "MatPlot.h"
#include "GeoCut.h"

int main(){

  std::vector<GeoCut> geoCuts;
  geoCuts.push_back(GeoCut(-180., 180., 1., 100.));

  //No End Cap
  //  geoCuts.push_back(GeoCut(-32.5, 32.5, 1., 17.5));
  //  geoCuts.push_back(GeoCut(-boxSize, boxSize, 17., 54.));

  MatPlot myPlot("RZ", -180., 180., 0., 100., 0.5, 0.5, -1.);
  myPlot.SetGeoCuts(&geoCuts);
  myPlot.SetUIndex(3);
  myPlot.SetVIndex(4);

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
