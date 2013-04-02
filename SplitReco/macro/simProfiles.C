#include <cmath>
#include <iostream>
#include "TROOT.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TStyle.h"

void palette(){


  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);


}



void simProfiles(){

  gROOT->SetStyle("Plain");
  
  TChain *uno = new TChain("sT");
  uno->Add("/raid/sguazz/split/MuGun_NomGeo.root"); 

  TH2F * dpdxSimVsEta = new TH2F("dpdxSimVsEta",";#eta;dP_{SIM}/dx [MeV/cm]",210,-2.1,2.1,120,0.,0.6);
  TH2F * plossSimVsEta = new TH2F("plossSimVsEta",";#eta;dP_{SIM} [MeV]",210,-2.1,2.1,100,0.,100.);
  TH2F * plossSimOverTVsEta = new TH2F("plossSimOverTVsEta",";#eta;dP_{SIM}/dx [MeV/cm]",210,-2.1,2.1,100,0.,0.6);
  TH2F * plossSimOverTSimVsEtaWCut = new TH2F("plossSimOverTSimVsEtaWCut",";#eta;dP_{SIM}/dx [MeV/cm]",210,-2.1,2.1,100,0.,0.6);
  TH2F * TSimVsEta = new TH2F("TSimVsEta",";#eta;T_{SIM} [cm]",210,-2.1,2.1,240,80.,320.);
  TH2F * TSimVsEtaWCut = new TH2F("TSimVsEtaWCut",";#eta;T_{SIM} [cm]",210,-2.1,2.1,240,80.,320.);
  TH2F * TVsEta = new TH2F("TVsEta",";#eta;T_{SIM} [cm]",210,-2.1,2.1,240,80.,320.);
  //  TH2F * RSimOutVsZSimOut = new TH2F("RSimOutVsZSimOut",";R_{OUT} Sim [cm];Z_{OUT} Sim [cm]",300,-600.,600.,150,0.,600.);
  
  // Cut
  TString cut("pt>0.9&&pt<1.1&&chi2<2.&&rIn<8.&&((rOut>100.)||(zOut>260.)||(zOut<-260.))&&pLossSim>0.&&dz<5.&&dz>-5.&&dpdxSim<0.001");
  TString cutReco("pt>0.9&&pt<1.1&&chi2<2.&&rIn<8.&&((rOut>100.)||(zOut>260.)||(zOut<-260.))&&(nHitIna+nHitMis+nHitBad)==0&&pLossSim>0.&&dz<5.&&dz>-5.&&dpdxSim<0.001");
  TString tsimcut("(((abs(eta)<0.3)&&(TSim<130.))||((abs(eta)>0.3)&&(TSim<154.545*abs(eta)+83.636)))");

  int maxentries = 10000000;

  uno->Draw("1000.*dpdxSim:eta>>dpdxSimVsEta",cut,"goff",maxentries);
  uno->Draw("1000.*pLossSim:eta>>plossSimVsEta",cut,"goff",maxentries);
  uno->Draw("1000.*pLossSim/T:eta>>plossSimOverTVsEta",cut,"goff",maxentries);
  uno->Draw("1000.*pLossSim/TSim:eta>>plossSimOverTSimVsEtaWCut",cut+"&&"+tsimcut,"goff",maxentries);
  uno->Draw("TSim:eta>>TSimVsEta",cut,"goff",maxentries);
  uno->Draw("TSim:eta>>TSimVsEtaWCut",cut+"&&"+tsimcut,"goff",maxentries);
  uno->Draw("T:eta>>TVsEta",cutReco,"goff",maxentries);
  //  uno->Draw("rOutSim:>>RSimOutVsZSimOut",cutReco,"goff",maxentries);

  palette();

  TCanvas *can = new TCanvas("can","",1100,800);
  can->SetRightMargin(0.15);
  //  can->SetLogz(1);

  dpdxSimVsEta->SetStats(kFALSE);
  dpdxSimVsEta->Draw("colz");
  can->SaveAs("dpdxSimVsEta.png");

  plossSimVsEta->SetStats(kFALSE);
  plossSimVsEta->Draw("colz");
  can->SaveAs("plossSimVsEta.png");

  plossSimOverTVsEta->SetStats(kFALSE);
  plossSimOverTVsEta->Draw("colz");
  can->SaveAs("plossSimOverTVsEta.png");

  plossSimOverTSimVsEtaWCut->SetStats(kFALSE);
  plossSimOverTSimVsEtaWCut->Draw("colz");
  can->SaveAs("plossSimOverTSimVsEtaWCut.png");

  TSimVsEta->SetStats(kFALSE);
  TSimVsEta->Draw("colz");
  can->SaveAs("TSimVsEta.png");

  TSimVsEtaWCut->SetStats(kFALSE);
  TSimVsEtaWCut->Draw("colz");
  can->SaveAs("TSimVsEtaWCut.png");

  TVsEta->SetStats(kFALSE);
  TVsEta->Draw("colz");
  can->SaveAs("TVsEta.png");

}
