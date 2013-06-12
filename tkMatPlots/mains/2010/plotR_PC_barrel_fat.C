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
  geoCuts.push_back(GeoCut(-26., 26., 2., 20.));
  geoCuts.push_back(GeoCut(-66., 66., 20., rMax));

  MatPlot myPlotR("R", 1., rMax, 0.5, 0.5);
  myPlotR.SetEffRadii(&radii);
  myPlotR.SetGeoCuts(&geoCuts);
  Float_t nDataEv = 10809460;
  // Float_t nSimEv = 6106513;
  Float_t nSimEv = 6106513*5./11.;
  Float_t scaleFact = nDataEv/nSimEv;
  //
  // Assumo che il numero di eventi iniziale fosse 11000000
  //
  // I primi 5M si usano come pseudo dati
  // i secondi 6M per il training...
  //

  //
  myPlotR.SetMCScaleFact(-1.);
  myPlotR.SetSimScaleFact(-1.);
  myPlotR.SetDataScaleFact(1.);

  convR2SforMatPlot convMC("../7TeV/ntuple_conversion_CMSSW358p3_minbias7TeV_2010-06-04.root");
  myPlotR.FillMC(&convMC, 5000000, 100000000); // First 5M events

  convS2RforMatPlot convSim("../7TeV/ntuple_conversion_CMSSW358p3_minbias7TeV_2010-06-04.root");
  myPlotR.responseTrain(&convSim, &convMC, 1, 5000000); //Remaining events
  myPlotR.FillSim(&convSim, 5000000, 100000000); // First 5000000 events

  convR2SforMatPlot convData("../7TeV/TIB1Int_Plus15_2.root");
  //  convData.SetCenterCoord(-0.1475, -0.3782, -0.4847);
  myPlotR.FillData(&convData);
  myPlotR.doUnfold();

  /*
  myPlotR.SetMCScaleFact(scaleFact);
  myPlotR.SetSimScaleFact(scaleFact);
  myPlotR.SetDataScaleFact(1.);

  convR2SforMatPlot convMC("../7TeV/ntuple_conversion_CMSSW356_minbias7TeV_2010-05-10.root");
  myPlotR.FillMC(&convMC, 1, 5000000); // First 5M events

  convS2RforMatPlot convSim("../7TeV/ntuple_conversion_CMSSW356_minbias7TeV_2010-05-10.root");
  myPlotR.responseTrain(&convSim, &convMC, 5000000, 20000000); //Remaining events
  myPlotR.FillSim(&convSim, 1, 5000000); // First 5000000 events

  convR2SforMatPlot convData("../7TeV/ntuple_conversion_goodcoll7TeV_2010-05-10.root");
  convData.SetCenterCoord(-0.1475, -0.3782, -0.4847);
  myPlotR.FillData(&convData);
  myPlotR.doUnfold();
  */

  myPlotR.PlotAll();
  myPlotR.test();
 
}
