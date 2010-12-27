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

  /*
  std::vector<Double_t> radii;
  radii.push_back(1.0);
  radii.push_back(6.0);
  radii.push_back(9.0);
  radii.push_back(17.0);
  radii.push_back(20.0);
  radii.push_back(25.5);
  radii.push_back(30.0);
  radii.push_back(34.0);
  radii.push_back(38.0);
  radii.push_back(46.0);
  radii.push_back(54.0);
  radii.push_back(57.0);
  radii.push_back(61.0);
  radii.push_back(65.0);
  */

  std::vector<Double_t> radii;
  radii.push_back(2.0); // 1.
  radii.push_back(3.0); // 3.5
  radii.push_back(6.0);
  radii.push_back(9.0);
  radii.push_back(17.0);
  radii.push_back(20.0);
  radii.push_back(25.5);
  radii.push_back(30.0);
  radii.push_back(34.0);
  radii.push_back(38.0);
  radii.push_back(46.0);
  radii.push_back(54.5); // 54.0
  radii.push_back(57.0);
  radii.push_back(61.5); // 61.0
  radii.push_back(65.0);

  std::vector<GeoCut> geoCuts;
  geoCuts.push_back(GeoCut(-26., 26., 1., 20.));
  //  geoCuts.push_back(GeoCut(-66., 66., 20., 65.));

  Double_t boxSize = 12.;
  MatPlot myPlot("XY", -boxSize, boxSize, -boxSize, boxSize, 0.25, 0.05);
  myPlot.SetEffRadii(&radii);
  myPlot.SetGeoCuts(&geoCuts);

  /*
  // Normalization factors 
  Float_t nDataEv = 10809460;
  Float_t nSimEv = 6106513;
  // Float_t nSimEv = 6106513*5./11.;
  Float_t scaleFact = nDataEv/nSimEv;
  myPlot.SetMCScaleFact(scaleFact);
  myPlot.SetSimScaleFact(scaleFact);
  myPlot.SetDataScaleFact(1.);
  */

  myPlot.SetMCScaleFact(-1.);
  myPlot.SetSimScaleFact(-1.);
  myPlot.SetDataScaleFact(1.);

  convR2SforMatPlot convMC("../7TeV/ntuple_conversion_CMSSW358p3_minbias7TeV_2010-06-04.root");
  myPlot.FillMC(&convMC); // All events

  convS2RforMatPlot convSim("../7TeV/ntuple_conversion_CMSSW358p3_minbias7TeV_2010-06-04.root");
  //  myPlot.responseTrain(&convSim, &convMC, 5000000, 20000000); //Remaining events
  myPlot.FillSim(&convSim); // All events

  convR2SforMatPlot convData("../7TeV/ntuple_conversion_CMSSW358p3_goodcoll7TeV_2010-06-04.root");
  convData.SetCenterCoord(-0.1475, -0.3782, -0.4847);
  myPlot.FillData(&convData);
  //  myPlot.doUnfold();

  myPlot.test();
  myPlot.PlotAll();

}
