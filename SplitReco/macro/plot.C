#include <vector>
#include <cmath>

#include "splitTreeForPlot.h"
#include "elPlot.h"

#include <TColor.h>
#include <TStyle.h>

void setPlotStyle() {

    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    gROOT->SetStyle("Plain");

}

int main(){

  elPlot * myPlot = 0;

  setPlotStyle();

  int iWrite = 1;
  
  if ( iWrite ) {
    
    //  elPlot myPlot("EtaPhi", -1.5, 1.5, 30, -PI, PI, 5);
    //  elPlot myPlot("EtaPhiOut_1", -1.5, 1.5, 15, 0., 2.*M_PI, 5);
    //  elPlot myPlot("EtaPhiOut_2", -1.5, 1.5, 15, 0., 2.*M_PI, 5);
    /////  elPlot myPlot("test", -1.5, 1.5, 5, 0., 2.*M_PI, 1);
    //  elPlot myPlot("test", -0.9, 0.9, 9, -0.125*M_PI, 2.*M_PI-0.125*M_PI, 8);
    
    myPlot = new elPlot("test", -2.1, 2.1, 21, -0.125*M_PI, 2.*M_PI-0.125*M_PI, 8);
    //    myPlot = new elPlot("test", -0.9, 0.9, 9, -0.125*M_PI, 2.*M_PI-0.125*M_PI, 8);
    //myPlot = new elPlot("test", -0.9, 0.9, 3, -0.125*M_PI, 2.*M_PI-0.125*M_PI, 4);
    
    //  elPlot myPlot("EtaPhi", 0.9, 1.0, 1, 1.88, 2.51, 10);
    //  elPlot myPlot("EtaPhi", 0.9, 1.0, 1, 1.88, 2.51, 10);
    
    //MC Reco
    //
    TChain *mc = new TChain("sT");
    //  mc->Add("./split_QCD_Pt-15to30_TuneZ2_7TeV_pythia6_GEN-SIM-RECODEBUG_PU_S6_START44_V5-v1.root");
    //    mc->Add("./split_qcd1530_v2.root");
    //    mc->Add("/raid/sguazz/split/minbias_run2012Cpart2c_nH4.root"); 
    //    mc->Add("/raid/sguazz/split/minbias_run2012Cpart1_nH4.root");
    // mc->Add("/raid/sguazz/split/minbias_run2012Cpart2a_nH4.root");
    // mc->Add("/raid/sguazz/split/minbias_run2012Cpart2b_nH4.root");
    // mc->Add("/raid/sguazz/split/minbias_run2012Cpart2c_nH4.root");
    mc->Add("/raid/sguazz/split/minbias_run2012A_nH4_v2.root"); 
    //mc->Add("/raid/sguazz/split/minBias_nH4_v2.root"); 
    //mc->Add("/raid/sguazz/split/MuGun_NomGeo.root"); 
    //mc->Add("/raid/sguazz/split/MuGun_NomGeo_part1.root"); 
    //mc->Add("/raid/sguazz/split/MuGun_NomGeo_part2.root"); 
    //mc->Add("/raid/sguazz/split/MuGun_ExtFlat10.root"); 
    //mc->Add("/raid/sguazz/split/MuGun_ExtFlat20.root"); 
    splitTreeForPlot split(mc, -1, -1); //entries for plots, entries for normalization
    //splitTreeForPlot split(mc, 10000000, 10000000);
    split.LoopForFill(myPlot);
    
    myPlot->PlotAndWriteAll();
    myPlot->FitAll();
    
  } else {

    myPlot = new elPlot("test");
    myPlot->FitAll();

  }

  delete myPlot;

}


