#include <vector>

#include "convR2SforMatPlot.h"
#include "convS2RforMatPlot.h"
#include "niR2SforMatPlot.h"
#include "niS2RforMatPlot.h"
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
  /*
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
  radii.push_back(2.0); // 1.
  radii.push_back(3.0); // 3.5
  radii.push_back(6.0);
  radii.push_back(9.0);
  radii.push_back(17.0);
  radii.push_back(20.0);
  radii.push_back(25.5);
  radii.push_back(28.0); //Added!!!
  radii.push_back(30.0);
  radii.push_back(34.0);
  radii.push_back(36.5); //Added!!!
  radii.push_back(38.0);
  radii.push_back(46.0);
  radii.push_back(54.5); // 54.0
  radii.push_back(57.0);
  radii.push_back(61.5); // 61.0
  radii.push_back(65.0);


  Double_t rMax = 65.;
  std::vector<GeoCut> geoCuts;
  geoCuts.push_back(GeoCut(-26., 26., 1., 20.));
  geoCuts.push_back(GeoCut(-66., 66., 20., 61.));

  MatPlot myPlotR("R", 1., rMax, .5, .5);
  myPlotR.SetEffRadii(&radii);
  myPlotR.SetGeoCuts(&geoCuts);

  //
  myPlotR.SetMCScaleFact(-1.);
  myPlotR.SetSimScaleFact(-1.);
  myPlotR.SetDataScaleFact(1.);

  /*  
      TChain *convmc = new TChain("ntupleR2S");
      convmc->Add("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04.root");
      convmc->Add("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04_v2.root");
  */
#include "src/files_ni_mcreco_386p1_20101124.cxx"
  niR2SforMatPlot convMC(mc);
  //  myPlotR.FillMC(&convMC, 10000000, 100000000); 
  myPlotR.FillMC(&convMC, 1, 1000000000); 

  /*
    TChain *convsim = new TChain("ntupleS2R");
    convsim->Add("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04.root");
    convsim->Add("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04_v2.root");
  */
#include "src/files_ni_mcsim_386p1_20101124.cxx"
  niS2RforMatPlot convSim(sim);
  //  myPlotR.responseTrain(&convSim, &convMC, 1, 10000000); //First 10M events for filling
  //  myPlotR.FillSim(&convSim, 10000000, 100000000); 
  myPlotR.FillSim(&convSim, 1, 1000000000); 

  /*
    niR2SforMatPlot convData("../7TeV/ntuple_nuclint_CMSSW358p3_goodcoll7TeV_2010-06-04.root");
  */
#include "src/files_ni_data_386p1_20101124.cxx"
  niR2SforMatPlot convData(data);
  convData.SetCenterCoord(-0.1475, -0.3782, -0.4847);
  myPlotR.FillData(&convData);
  //  myPlotR.doUnfold();

  myPlotR.PlotAll();

  myPlotR.test();
}
