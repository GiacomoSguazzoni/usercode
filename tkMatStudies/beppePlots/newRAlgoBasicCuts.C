#include <iostream>
#include <cstdlib>

#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TColor.h"
#include "TAxis.h"
#include "TRandom3.h"

TString generateRandomName(){
  //  TRandom3 ran;
  //  Int_t nomeNumero = (Int_t)(ran.Uniform(0.,99999.));
  Int_t nomeNumero = (Int_t)(std::rand());

  char nome[20];
  sprintf(nome,"h%d",nomeNumero);
  return TString(nome);
}


TCanvas* makeNiceCanvas(Int_t pixelPerBin, Int_t nbinx, Int_t nbiny, Double_t top, Double_t bottom, Double_t left, Double_t right) {

  Int_t rubaX = 4; //determinato sperimentalmente
  Int_t rubaY = 28; //determinato sperimentalmente

  TString name = generateRandomName();

  Int_t plotBaseDimX = pixelPerBin*nbinx;
  Int_t plotBaseDimY = pixelPerBin*nbiny;
  Int_t XX = (Int_t)(plotBaseDimX/(1.-left-right));
  Int_t YY = (Int_t)(plotBaseDimY/(1.-top-bottom));
  TCanvas* can = new TCanvas(name,name,XX+rubaX,YY+rubaY);
  can->SetTopMargin(top);
  can->SetBottomMargin(bottom);
  can->SetRightMargin(right);
  can->SetLeftMargin(left);
  std::cout << "Nice canvas " << XX << " * " << YY << std::endl;
  return can;

}

Double_t algoFun(Double_t algo){

  Double_t offset=0.5;
  if (algo==4) return 0+offset;
  if (algo==5) return 1+offset;
  if (algo==6) return 2+offset;
  if (algo==7) return 3+offset;
  if (algo==8) return 4+offset;
  if (algo==9) return 5+offset;
  if (algo==12) return 8+offset;
  if (algo==29) return 9+offset;

  return -1+offset;

}

Double_t algoCharge(Int_t q, Double_t algo1, Double_t algo2, Double_t q1, Double_t q2){

  if ( q==(Int_t)q1 ){
    return algoFun(algo1);
  } else {
    return algoFun(algo2);
  }

  return 0;

}

void SetFancyPalette(){
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void SetIterLabels(TAxis * axis, Int_t offset){

  axis->SetBinLabel(0+offset,"iter 0");
  axis->SetBinLabel(1+offset,"iter 1");
  axis->SetBinLabel(2+offset,"iter 2");
  axis->SetBinLabel(3+offset,"iter 3");
  axis->SetBinLabel(4+offset,"iter 4");
  axis->SetBinLabel(5+offset,"iter 5");
  //  axis->SetBinLabel(6,"SingleTk");
  //  axis->SetBinLabel(7,"SingleTk");
  axis->SetBinLabel(8+offset,"SingleTk");
  axis->SetBinLabel(9+offset,"GSF");

}

void SetSpecialLabels(TH2F * h){

  h->GetXaxis()->SetBinLabel(2, "Total");
  SetIterLabels(h->GetXaxis(), 4);
  SetIterLabels(h->GetYaxis(), 2);

}

void doAlgo2vs1Plot(const char * name, TChain *mc, TChain *mcEvt, TString cut){

  TString Name(name);

  Int_t nbinX = 13;
  Int_t nbinY = 11;

  TString algoName = generateRandomName();

  TH2F * algo2vs1 = new TH2F(algoName,algoName+";Tk 1;Tk 2;Conv/Evt/10^{-2}",nbinX,-3.,10.,nbinY,-1.,10.);
  mc->Draw("algoFun(algo2):algoFun(algo1)>>"+algoName,cut,"goff");
  Int_t nEvt = mcEvt->GetEntries();
  Double_t scaleFact = 1./nEvt;
  algo2vs1->Scale(scaleFact*100.);

  //
  //Fill cumulative values
  for (Int_t iY=2;iY<nbinY+1; iY++){
    Double_t sum = 0.;
    for (Int_t iX=4;iX<nbinX+1; iX++){
      sum+=algo2vs1->GetBinContent(iX, iY);
    }
    algo2vs1->SetBinContent(2, iY, sum);
  }


  TCanvas* can = makeNiceCanvas(65, nbinX, nbinY, 0.1, 0.15, 0.15, 0.2);

  can->SetGrid();
  can->cd();
  SetFancyPalette();
  gStyle->SetPaintTextFormat("2.2f");
  algo2vs1->SetStats(0);
  algo2vs1->SetTitle("");
  algo2vs1->SetTitleOffset(1.3,"XY");
  algo2vs1->SetMaximum(1.05);

  SetSpecialLabels(algo2vs1);

  algo2vs1->Draw("Zcol");
  algo2vs1->SetMarkerColor(kGray+1);
  algo2vs1->SetMarkerSize(1.4);
  algo2vs1->Draw("textsame");
  can->SaveAs("algo2vs1_"+Name+".png");

}

void doAlgoPosVsElePlot(const char * name, TChain *mc, TChain *mcEvt, TString cut){

  TString Name(name);

  Int_t nbinX = 11;
  Int_t nbinY = 11;

  TString algoName = generateRandomName();

  TH2F * algo2vs1 = new TH2F(algoName,algoName+";Ele;Pos;Conv/Evt/10^{-2}",nbinX,-1.,10.,nbinY,-1.,10.);
  mc->Draw("algoCharge(1,algo1,algo2,q1,q2):algoCharge(-1,algo1,algo2,q1,q2)>>"+algoName,cut,"goff");
  Int_t nEvt = mcEvt->GetEntries();
  Double_t scaleFact = 1./nEvt;
  algo2vs1->Scale(scaleFact*100.);

  TCanvas* can = makeNiceCanvas(65, nbinX, nbinY, 0.1, 0.15, 0.15, 0.2);

  can->SetGrid();
  can->cd();
  SetFancyPalette();
  gStyle->SetPaintTextFormat("2.2f");
  algo2vs1->SetStats(0);
  algo2vs1->SetTitle("");
  algo2vs1->SetTitleOffset(1.3,"XY");
  algo2vs1->SetMaximum(1.05);

  SetIterLabels(algo2vs1->GetXaxis(), 2);
  SetIterLabels(algo2vs1->GetYaxis(), 2);

  algo2vs1->Draw("Zcol");
  algo2vs1->SetMarkerColor(kGray+1);
  algo2vs1->SetMarkerSize(1.4);
  algo2vs1->Draw("textsame");
  can->SaveAs("algoPosVsEle_"+Name+".png");

}


TH1F * addToStack(TH1F * in, TChain *mc, TString var, TString cut, Int_t colorIndex){

  Int_t nbin = in->GetNbinsX();
  Double_t minbin = in->GetBinLowEdge(1);
  Double_t maxbin = in->GetBinLowEdge(nbin)+in->GetBinWidth(nbin);

  TString nome = generateRandomName();

  TH1F * newHtmp = new TH1F(nome+"tmp", nome+"tmp", nbin, minbin, maxbin);
  mc->Draw(var+">>"+nome+"tmp",cut,"goff");

  TH1F * newH = new TH1F(nome, nome, nbin, minbin, maxbin);

  for (Int_t j=1; j<nbin+1; j++){
      newH->SetBinContent(j, newHtmp->GetBinContent(j)+in->GetBinContent(j));
  }
  newH->SetFillColor(colorIndex);

  return newH;

}

void doVarPlot(const char * name, TChain *mc, TChain *mcEvt, TString cut, TString var, TString varName, Int_t nbin, Double_t minbin, Double_t maxbin){

  TString Name(name);

  TString canName = generateRandomName();
  TCanvas * c1 = new TCanvas(canName,canName,800+4,600+28);
  gStyle->SetOptStat(0);

  TString nome = generateRandomName();
  TH1F* dum = new TH1F(nome,nome,nbin, minbin, maxbin);
  TH1F* rA4 = addToStack(dum, mc, var, cut+"&&max(algo1,algo2)==4", kBlack);
  TH1F* rA5 = addToStack(rA4, mc, var, cut+"&&max(algo1,algo2)==5", kRed);
  TH1F* rA6 = addToStack(rA5, mc, var, cut+"&&max(algo1,algo2)==6", kGreen);
  TH1F* rA7 = addToStack(rA6, mc, var, cut+"&&max(algo1,algo2)==7", kMagenta);
  TH1F* rA8 = addToStack(rA7, mc, var, cut+"&&max(algo1,algo2)==8", kCyan);
  TH1F* rA9 = addToStack(rA8, mc, var, cut+"&&max(algo1,algo2)==9", kOrange);
  TH1F* rA10 = addToStack(rA9, mc, var, cut+"&&max(algo1,algo2)==10", kGray);
  TH1F* rA11 = addToStack(rA10, mc, var, cut+"&&max(algo1,algo2)==11", kBlue);
  TH1F* rA29 = addToStack(rA11, mc, var, cut+"&&max(algo1,algo2)==29", kSpring+3);
  TH1F* rA12 = addToStack(rA29, mc, var, cut+"&&max(algo1,algo2)==12", kViolet+1);

  std::cout << " Integralone " << (float)rA12->Integral() << std::endl;

  Double_t scaleFact = 1./mcEvt->GetEntries();

  rA4->Scale(scaleFact);
  rA5->Scale(scaleFact);
  rA6->Scale(scaleFact);
  rA7->Scale(scaleFact);
  rA8->Scale(scaleFact);
  rA9->Scale(scaleFact);
  rA10->Scale(scaleFact);
  rA11->Scale(scaleFact);
  rA12->Scale(scaleFact);
  rA29->Scale(scaleFact);

  //--------------------

  c1->cd();
  rA12->SetTitle("");
  rA12->GetXaxis()->SetTitle(varName);
  rA12->GetYaxis()->SetTitle("Fraction of #gamma_{conv.}/event");
  rA12->GetYaxis()->SetTitleOffset(1.25);
  rA12->SetMinimum(0.00000002);
  rA12->SetMaximum(0.02);
  rA12->Draw();
  rA29->Draw("same");
  rA11->Draw("same");
  rA10->Draw("same");
  rA9->Draw("same");
  rA8->Draw("same");
  rA7->Draw("same");
  rA6->Draw("same");
  rA5->Draw("same");
  rA4->Draw("same");

  TLegend * leg = new TLegend(0.1,0.9,0.9,0.98);
  leg->SetNColumns(5);
  leg->SetFillColor(kWhite);
  leg->AddEntry(rA4, "N_{algo}=0","f");
  leg->AddEntry(rA5, "N_{algo}=1","f");
  leg->AddEntry(rA6, "N_{algo}=2","f");
  leg->AddEntry(rA7, "N_{algo}=3","f");
  leg->AddEntry(rA8, "N_{algo}=4","f");
  leg->AddEntry(rA9, "N_{algo}=5","f");
  leg->AddEntry(rA10,"N_{algo}=6","f");
  leg->AddEntry(rA11,"N_{algo}=7","f");
  leg->AddEntry(rA12,"N_{algo}=8","f");
  leg->AddEntry(rA29,"GSF","f");
  leg->Draw();

  TLatex *t1a = new TLatex();
  t1a->SetText(0.6,0.76,"CMS Preliminary / "+Name);
  t1a->SetTextFont(42);
  t1a->SetTextSize(0.05);
  t1a->SetTextAlign(22);
  t1a->SetNDC(kTRUE);
  t1a->Draw();

  c1->Update();
  c1->SetLogy();
  c1->SetGridy();
  c1->Update();

  c1->SaveAs("algo_"+varName+"_"+Name+".png");

  //-------------------------------------------------------------------

}


void newRAlgoBasicCuts(){
  
  std::cout << " Let's rock..." << std::endl;
  
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  TString r="sqrt(x*x+y*y)";
  TString rSIM="sqrt((x+deltax)*(x+deltax)+(y+deltay)*(y+deltay))";
  TString pt="pt";
  TString ptSIM="pt+deltapt";

  TString cut="pt<10 && abs(-log(tan(theta/2)))<1.4 &&  nHits1>4 && nHits2>4 && chi2prob>0.0005 && minapp<1.";// && pt1>0.5 && pt2>0.5"
  //TString cut="pt<10 && abs(-log(tan(theta/2)))<1.4 &&"+r+">15.";// && pt1>0.5 && pt2>0.5"

  
  TChain *mc = new TChain("ntupleR2S");
  TChain *mcEvt = new TChain("ntupleEvt");
    //#include "fileInclude_PAS.icc"
    //#include "fileInclude_20101124.icc"
    //#include "fileInclude_MinBiasMC_386p1_altTrig_DG_10MEvt_1kSeedLimit.icc"
    //#include "fileInclude_MinBiasMC_386p1_altTrig_DG_10MEvt_10kSeedLimit.icc"
    //#include "fileInclude_MinBiasMC_386p1_altTrig_DG_10MEvt_10kSeedLimit_GiuFilter.icc"
    //#include "fileInclude_MinBiasMC_386p1_altTrig_DG_10MEvt_10kSeedLimit_2kmaxInputSeeds_GiuFilter.icc"
    //#include "fileInclude_MinBiasMC_386p1_altTrig_DG_10MEvt_10kSeedLimit_2kmaxInputSeeds.icc"
    //#include "fileInclude_MinBiasMC_386p1_altTrig_DG_10MEvt_10kSeedLimit_10kmaxInputSeeds_GiuFilter.icc"
    //#include "fileInclude_MinBiasMC_386p1_altTrig_DG_10MEvt_10kSeedLimit_10kmaxInputSeeds.icc"
    //#include "fileInclude_MinBiasMC_397_altTrig_DG_10MEvt_10kSeedLimit_10kmaxInputSeeds_GiuFilter.icc"
    //#include "fileInclude_MinBiasMC_397_HLiuConv_altTrig_DG_10MEvt_10kSeedLimit_10kmaxInputSeeds_GiuFilter.icc"
    //#include "fileInclude_MinBiasMC_397_altTrig_DG_10MEvt_seedSingleTrack.icc"
#include "fileInclude_MinBiasMC_397_altTrig_DG_10MEvt_seedSingleTrack_inCKFandGSF.icc"
    //#include "fileInclude_MinBiasMC_397_altTrig_DG_10MEvt_seedSingleTrack_onlyInGSF.icc"
    //#include "fileInclude_MinBiasMC_397_altTrig_DG_10MEvt_seedSingleTrack_STD.icc"

  TString Name("_inCKFandGSF"); 
  //  TString Name("_onlyInGSF"); 
  //  TString Name("_STD"); 
  TString cutRECO=cut;
  TString cutFAKE=cutRECO+"&& isAssoc==0";
  TString cutASSO=cutRECO+"&& isAssoc==1";
  doAlgo2vs1Plot("Reco"+Name, mc, mcEvt, cutRECO);
  doAlgo2vs1Plot("Fake"+Name, mc, mcEvt, cutFAKE);
  doAlgo2vs1Plot("Asso"+Name, mc, mcEvt, cutASSO);
  doAlgoPosVsElePlot("Reco"+Name, mc, mcEvt, cutRECO);
  doAlgoPosVsElePlot("Fake"+Name, mc, mcEvt, cutFAKE);
  doAlgoPosVsElePlot("Asso"+Name, mc, mcEvt, cutASSO);
  doVarPlot("Reco"+Name, mc, mcEvt, cutRECO, r, "r", 80, 0., 80.);
  doVarPlot("Fake"+Name, mc, mcEvt, cutFAKE, r, "r", 80, 0., 80.);
  doVarPlot("Asso"+Name, mc, mcEvt, cutASSO, r, "r", 80, 0., 80.);
  doVarPlot("Reco"+Name, mc, mcEvt, cutRECO, "pt", "pt", 80, 0., 10.);
  doVarPlot("Fake"+Name, mc, mcEvt, cutFAKE, "pt", "pt", 80, 0., 10.);
  doVarPlot("Asso"+Name, mc, mcEvt, cutASSO, "pt", "pt", 80, 0., 10.);


//------------------------------------------------------------------------------

}



