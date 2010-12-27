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

  std::vector<GeoCut> geoCuts;
  geoCuts.push_back(GeoCut(-5., 5., 1., 65.));

  MatPlot myPlotR("R", 1., 60., 1., 1.);
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
  myPlotR.SetMCScaleFact(scaleFact);
  myPlotR.SetSimScaleFact(scaleFact);
  myPlotR.SetDataScaleFact(1.);
  //  MatPlot myPlotEta("Eta", -3., 3., 0.1, 0.01);
  //  myPlotEta.SetGeoCuts(&geoCuts);

  //  Double_t boxSize = 60.;
  //  MatPlot myPlotXY("XY", -boxSize, boxSize, -boxSize, boxSize, 1., 0.1);
  //  myPlotXY.SetGeoCuts(&geoCuts);

  convR2SforMatPlot convMC("../7TeV/ntuple_conversion_CMSSW356_minbias7TeV_2010-05-10.root");
  myPlotR.FillMC(&convMC, 1, 5000000); // First 5M events
  //  myPlotXY.FillMC(&convMC);

  convS2RforMatPlot convSim("../7TeV/ntuple_conversion_CMSSW356_minbias7TeV_2010-05-10.root");
  myPlotR.responseTrain(&convSim, &convMC, 5000000, 20000000); //Remaining events
  myPlotR.FillSim(&convSim, 1, 5000000); // First 5000000 events
  //  myPlotXY.FillSim(&convSim);

  convR2SforMatPlot convData("../7TeV/ntuple_conversion_goodcoll7TeV_2010-05-10.root");
  myPlotR.FillData(&convData);
  myPlotR.doUnfold();
  //  myPlotXY.FillData(&convData);

  myPlotR.PlotAll();
  //  myPlotXY.PlotAll();
 
}
