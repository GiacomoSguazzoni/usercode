#include <vector>

#include "niR2SforMatPlot.h"
#include "niS2RforMatPlot.h"
#include "MatPlot.h"
#include "GeoCut.h"

int main(){

  std::vector<GeoCut> geoCuts;
  geoCuts.push_back(GeoCut(-24., 24., 1., 20.));
  geoCuts.push_back(GeoCut(-66., 66., 20., 65.));

  Double_t boxSize = 22.;
  //MatPlot myPlot("XY", -boxSize, boxSize, -boxSize, boxSize, 1., 0.1);
  //binsize, binsizesim, binsizeeff
  MatPlot myPlot("XY", -boxSize, boxSize, -boxSize, boxSize, 0.1, 0.1, -1.); // era 2 
  myPlot.SetGeoCuts(&geoCuts);
  myPlot.SetUIndex(1, 3); //X for plot, zeta for cut
  myPlot.SetVIndex(2, 4); //Y for plot, radius for cut

  //
  myPlot.SetMCScaleFact(-1.);
  myPlot.SetSimScaleFact(-1.);
  myPlot.SetDataScaleFact(1.);

  //MC Reco
  //
#include "src/files_ni_dummy_reco.cxx"
  niR2SforMatPlot niMC(mc);
  myPlot.FillMC(&niMC);

  //MC Sim
  //
#include "src/files_ni_dummy_sim.cxx"
  niS2RforMatPlot niSim(sim);
  myPlot.FillSim(&niSim, 1, 2000000);

  //Data
  //
#include "src/files_ni_run2012D.cxx"
  niR2SforMatPlot niData(data);
  //  convData.SetCenterCoord(-0.1475, -0.3782, -0.4847);
  //
  // run2012D
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TkAlignmentPixelPosition 
  niData.SetCenterCoord(-0.10822, -0.376969, -0.410475);
  myPlot.FillData(&niData);

  myPlot.test();
  myPlot.PlotAll();

}
