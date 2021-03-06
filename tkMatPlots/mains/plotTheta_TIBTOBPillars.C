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

  Double_t pi = 3.141592653589793;

  std::vector<GeoCut> geoCuts;

  Double_t thetaMax = 2.*atan(exp(1.));
  Double_t thetaMin = 2.*atan(exp(-1.));
			    
  std::cout << " thetaMin " << thetaMin << " thetaMax " << thetaMax << std::endl;

  // eta = -log(tan(theta/2.)) ==> theta= 2.*tan(exp(-eta))
  geoCuts.push_back(GeoCut(-pi, pi, 26., 67.));

  Double_t binW = (thetaMax-thetaMin)/40.;

  std::cout << " thetaMin " << thetaMin << " thetaMax " << thetaMax << " binW " << binW << std::endl;

  MatPlot myPlot("Theta", thetaMin, thetaMax, binW, binW, -1.);
  myPlot.SetGeoCuts(&geoCuts);
  myPlot.SetUIndex(5, 6); //Tetha for plot, phi for cut
  myPlot.SetVIndex(0, 4); //Plot is 1d, radius for cut

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
