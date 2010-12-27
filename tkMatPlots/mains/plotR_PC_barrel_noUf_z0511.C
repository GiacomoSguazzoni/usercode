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

  Double_t rMax = 65.;
  std::vector<GeoCut> geoCuts;
  geoCuts.push_back(GeoCut(5., 11., 2., rMax));

  MatPlot myPlotR("R", 1., rMax, 0.5, 0.5);
  myPlotR.SetEffRadii(&radii);
  myPlotR.SetGeoCuts(&geoCuts);

  //
  myPlotR.SetMCScaleFact(-1.);
  myPlotR.SetSimScaleFact(-1.);
  myPlotR.SetDataScaleFact(1.);

  //MC Reco
  //
  //#include "files_conv_mcreco_20100619.cxx"
#include "files_conv_mcreco_20100713.cxx"
  convR2SforMatPlot convMC(convmc);
//  convR2SforMatPlot convMC("../7TeV/ntuple_conversion_CMSSW358p3_minbias7TeV_2010-06-04_v2.root");
  myPlotR.FillMC(&convMC, 1, 100000000); 

  //MC Sim
  //
  //#include "files_conv_mcsim_20100619.cxx"
#include "files_conv_mcsim_20100713.cxx"
  convS2RforMatPlot convSim(convsim);
  //  convS2RforMatPlot convSim("../7TeV/ntuple_conversion_CMSSW358p3_minbias7TeV_2010-06-04.root");
  myPlotR.FillSim(&convSim, 1, 100000000);

  //Data
  //
  //#include "files_conv_data_20100619.cxx"
#include "files_conv_data_20100713.cxx"
  convData.SetCenterCoord(-0.1475, -0.3782, -0.4847);
  myPlotR.FillData(&convData);

  myPlotR.PlotAll();
  myPlotR.test();
 
}
