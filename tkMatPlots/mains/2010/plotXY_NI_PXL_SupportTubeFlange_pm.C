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

  std::vector<GeoCut> geoCuts;
  geoCuts.push_back(GeoCut(-123., -109.5, 1., 30.));
  geoCuts.push_back(GeoCut( 109.5,  123., 1., 30.)); //Plot is ok but I'm not sure correction factors are correctly computed

  Double_t boxSize = 30.;
  MatPlot myPlot("XY", -boxSize, boxSize, -boxSize, boxSize, 0.2, 0.05);
  myPlot.SetEffRadii(&radii);
  myPlot.SetGeoCuts(&geoCuts);

  //
  myPlot.SetMCScaleFact(-1.);
  myPlot.SetSimScaleFact(-1.);
  myPlot.SetDataScaleFact(1.);

  TChain *convmc = new TChain("ntupleR2S");
  convmc->Add("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04.root");
  convmc->Add("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04_v2.root");
  niR2SforMatPlot convMC(convmc);
  //  niR2SforMatPlot convMC("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04.root");
  myPlot.FillMC(&convMC); // All events

  TChain *convsim = new TChain("ntupleS2R");
  convsim->Add("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04.root");
  convsim->Add("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04_v2.root");
  niS2RforMatPlot convSim(convsim);
  //  niS2RforMatPlot convSim("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04.root");
  //  myPlot.responseTrain(&convSim, &convMC, 5000000, 20000000); //Remaining events
  myPlot.FillSim(&convSim); // All events

  niR2SforMatPlot convData("../7TeV/ntuple_nuclint_CMSSW358p3_goodcoll7TeV_2010-06-04.root");
  convData.SetCenterCoord(-0.1475, -0.3782, -0.4847);
  myPlot.FillData(&convData);
  //  myPlot.doUnfold();

  myPlot.test();
  myPlot.PlotAll();

}
