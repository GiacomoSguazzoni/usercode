#include <vector>

#include "convR2SforMatPlot.h"
#include "convS2RforMatPlot.h"
#include "MatPlot.h"
#include "GeoCut.h"

int main(){

  std::vector<GeoCut> geoCuts;
  geoCuts.push_back(GeoCut(10., 30., 10., 80.));

  Double_t boxSize = 60.;
  //MatPlot myPlot("XY", -boxSize, boxSize, -boxSize, boxSize, 1., 0.1);
  //binsize, binsizesim, binsizeeff
  MatPlot myPlot("XY", -45., -5., 10., 50., 0.2, 0.1, 120.);
  myPlot.SetGeoCuts(&geoCuts);

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
