#include <vector>

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

  Double_t boxSize = 150.;
  std::vector<GeoCut> geoCuts;
  //Endcap
  geoCuts.push_back(GeoCut(-10., boxSize, 1., 60.));

  //No End Cap
  //  geoCuts.push_back(GeoCut(-32.5, 32.5, 1., 17.5));
  //  geoCuts.push_back(GeoCut(-boxSize, boxSize, 17., 54.));

  MatPlot myPlot("RZ", -10., boxSize, 0., 60., .5, 0.05);
  myPlot.SetEffRadii(&radii);
  myPlot.SetGeoCuts(&geoCuts);
  // vAr[0]=1.;
  // vAr[1]=X;
  // vAr[2]=Y;
  // vAr[3]=Z;
  // vAr[4]=sqrt(X*X+Y*Y);
  myPlot.SetUIndex(3);
  myPlot.SetVIndex(4);

  //
  myPlot.SetMCScaleFact(-1.);
  myPlot.SetSimScaleFact(-1.);
  myPlot.SetDataScaleFact(1.);

  TChain *nimc = new TChain("ntupleR2S");
  nimc->Add("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04.root");
  nimc->Add("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04_v2.root");
  niR2SforMatPlot niMC(nimc);
  myPlot.FillMC(&niMC); // All events

  TChain *nisim = new TChain("ntupleS2R");
  nisim->Add("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04.root");
  nisim->Add("../7TeV/ntuple_nuclint_CMSSW358p3_minbias7TeV_2010-06-04_v2.root");
  niS2RforMatPlot niSim(nisim);
  myPlot.FillSim(&niSim); // All events

  niR2SforMatPlot niData("../7TeV/ntuple_nuclint_CMSSW358p3_goodcoll7TeV_2010-06-04.root");
  niData.SetCenterCoord(-0.1475, -0.3782, -0.4847);
  myPlot.FillData(&niData);

  myPlot.test();
  myPlot.PlotAll();

}
