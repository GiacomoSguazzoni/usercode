#include <cmath>
#include <iostream>
#include "TROOT.h"
#include "TString.h"
#include "TH1F.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"

void draw(TH1F * uno, TH1F * due, TString &n1, TString &n2){

  TString name(uno->GetName());

  std::cout << " Plot " << name << " uno:"  << uno->GetEntries() << " due:" << due->GetEntries() << std::endl;
  std::cout << " Plot " << name << " uno:"  << uno->Integral() << " due:" << due->Integral() << std::endl;

  name = "Plot_"+name;

  TCanvas * can = new TCanvas(name,"",1000,750);

  uno->Sumw2();
  due->Sumw2();

  uno->Scale(1./(uno->Integral()));
  due->Scale(1./(due->Integral()));

  due->SetStats(kFALSE);
  uno->SetStats(kFALSE);

  due->SetMaximum(1.1*max(uno->GetMaximum(),due->GetMaximum()));
  due->SetMinimum(0.00001);

  TH1F *dueBis = (TH1F*)due->Clone("dueBis");

  dueBis->SetLineColor(kBlue);
  dueBis->SetLineWidth(2);
  dueBis->SetFillStyle(0);
  dueBis->Draw("HIST");

  due->SetMarkerStyle(20);
  due->SetMarkerSize(0.55);
  due->SetMarkerColor(kAzure+7);
  due->SetFillColor(kAzure+8);
  due->SetLabelSize(0.06,"xy");
  due->SetTitleSize(0.06,"xy");
  due->Draw("PE2same");
  uno->SetMinimum(0.01);
  uno->SetMarkerStyle(20);
  uno->SetMarkerSize(0.75);
  uno->Draw("PE0X0same");


  TLegend * leg = new TLegend(0.1,.92,0.9,0.98);
  leg->SetNColumns(2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(due,n2,"fp");
  leg->AddEntry(uno,n1,"p");
  leg->Draw();

  can->SaveAs(name+".png");

}

void sampleComparison(){

  gROOT->SetStyle("Plain");
  
  /*
  TChain *uno = new TChain("sT");
  uno->Add("/raid/sguazz/split/MuGun_NomGeo.root"); 
  TString n1("MuGun NomGeo");
  
  TChain *due = new TChain("sT");
  due->Add("/raid/sguazz/split/minBias_nH4_v2.root"); 
  TString n2("MinBias");
  */
  TChain *uno = new TChain("sT");
  uno->Add("/raid/sguazz/split/minBias_nH4_v2.root"); 
  TString n1("minBias MC");
  
  TChain *due = new TChain("sT");
  due->Add("/raid/sguazz/split/minbias_run2012Cpart1_nH4.root"); 
  TString n2("minBias run2012C");
    
  std::cout << " Sample uno:"  << uno->GetEntries() << " due:" << due->GetEntries() << std::endl;
  
  TH1F * chi2SplitUno = new TH1F("chi2Split",";#chi^{2}_{SPLIT}/d.o.f.;",40, 0.,20.);
  TH1F * chi2SplitDue = new TH1F("chi2SplitDue",";#chi^{2}_{SPLIT}/d.o.f.;",40, 0.,20.);
  
  TH1F * etaUno = new TH1F("eta",";#eta;",50,-2.5,2.5);
  TH1F * etaDue = new TH1F("etaDue",";#eta;",50,-2.5,2.5);
  
  TH1F * phiUno = new TH1F("phi",";#phi;",50,-M_PI,M_PI);
  TH1F * phiDue = new TH1F("phiDue",";#phi;",50,-M_PI,M_PI);
  
  TH1F * ptUno = new TH1F("pt",";P_{T} [GeV];",50,0.7,1.3);
  TH1F * ptDue = new TH1F("ptDue",";P_{T} [GeV];",50,0.7,1.3);
  
  TH1F * ptNormUno = new TH1F("ptNorm",";P_{T} [GeV];",200,0.9,1.1);
  TH1F * ptNormDue = new TH1F("ptNormDue",";P_{T} [GeV];",200,0.9,1.1);
  
  TH1F * nUno = new TH1F("n",";N_{SPLIT};",7,-0.5,6.5);
  TH1F * nDue = new TH1F("nDue",";N_{SPLIT};",7,-0.5,6.5);

  TH1F * dzUno = new TH1F("dz",";dz;",40,-20.,20.);
  TH1F * dzDue = new TH1F("dzDue",";dz;",40,-20.,20.);
  
  TH1F * routUno = new TH1F("rOut",";rOut;",30,100.,115.);
  TH1F * routDue = new TH1F("rOutDue",";rOut;",30,100.,115.);
  
  // MC
  //    TString cut("pt>0.9&&pt<1.1&&chi2<2.&&rIn<8.&&((rOut>100.)||(zOut>260.)||(zOut<-260.))&&pLossSim>0.&&(nHitIna+nHitMis+nHitBad)==0&&dz<5.&&dz>-5.&&maxchi2<4.&&chi2Split>0.");
  //    TString cutNoPt("chi2<2.&&rIn<8.&&((rOut>100.)||(zOut>260.)||(zOut<-260.))&&pLossSim>0.&&(nHitIna+nHitMis+nHitBad)==0&&dz<5.&&dz>-5.&&maxchi2<4.&&chi2Split>0.");
  // Data
  TString        cut("pt>0.9&&pt<1.1&&chi2<2.&&rIn<8.&&((rOut>100.)||(zOut>260.)||(zOut<-260.))&&(nHitIna+nHitMis+nHitBad)==0&&dz<5.&&dz>-5.&&maxchi2<4.&&(chi2Split/(NSplit-2.))>0.&&(chi2Split/(NSplit-2.))<4.");
  TString    cutROut("pt>0.9&&pt<1.1&&chi2<2.&&rIn<8.&&((rOut>100.)||(zOut>260.)||(zOut<-260.))&&(nHitIna+nHitMis+nHitBad)==0&&dz<5.&&dz>-5.&&maxchi2<4.&&(chi2Split/(NSplit-2.))>0.&&(chi2Split/(NSplit-2.))<4.&&eta<0.9&&eta>-0.9");
  TString    cutNoDz("pt>0.9&&pt<1.1&&chi2<2.&&rIn<8.&&((rOut>100.)||(zOut>260.)||(zOut<-260.))&&(nHitIna+nHitMis+nHitBad)==0&&maxchi2<4.&&(chi2Split/(NSplit-2.))>0.&&(chi2Split/(NSplit-2.))<4.");
  TString   cutNoChi("pt>0.9&&pt<1.1&&chi2<2.&&rIn<8.&&((rOut>100.)||(zOut>260.)||(zOut<-260.))&&(nHitIna+nHitMis+nHitBad)==0&&dz<5.&&dz>-5.&&maxchi2<4.&&(chi2Split/(NSplit-2.))>0.");
  TString    cutNoPt(                "chi2<2.&&rIn<8.&&((rOut>100.)||(zOut>260.)||(zOut<-260.))&&(nHitIna+nHitMis+nHitBad)==0&&dz<5.&&dz>-5.&&maxchi2<4.&&(chi2Split/(NSplit-2.))>0.&&(chi2Split/(NSplit-2.))<4.");
  TString cutNoSplit("pt>0.9&&pt<1.1&&chi2<2.&&rIn<8.&&((rOut>100.)||(zOut>260.)||(zOut<-260.))&&(nHitIna+nHitMis+nHitBad)==0&&dz<5.&&dz>-5.&&maxchi2<4.");
  
  int maxentries = 10000000;
  //int maxentries = 999999999;
  
  //TString assoc("&&pLossSim>0.");
  TString assoc("");
  
    uno->Draw("(chi2Split/(NSplit-2.))>>chi2Split",cutNoChi+assoc,"goff",maxentries);
    due->Draw("(chi2Split/(NSplit-2.))>>chi2SplitDue",cutNoChi+assoc,"goff",maxentries);
    draw(chi2SplitUno, chi2SplitDue, n1, n2);
    
    uno->Draw("eta>>eta",cut+assoc,"goff",maxentries);
    due->Draw("eta>>etaDue",cut+assoc,"goff",maxentries);
    draw(etaUno, etaDue, n1, n2);

    uno->Draw("pt>>pt",cutNoPt+assoc,"goff",maxentries);
    due->Draw("pt>>ptDue",cutNoPt+assoc,"goff",maxentries);
    draw(ptUno, ptDue, n1, n2);

    uno->Draw("pt>>pt",cut+assoc,"goff",maxentries);
    due->Draw("pt>>ptDue",cut+assoc,"goff",maxentries);
    draw(ptUno, ptDue, n1, n2);
  
  uno->Draw("dz>>dz",cutNoDz+assoc,"goff",maxentries);
  due->Draw("dz>>dzDue",cutNoDz+assoc,"goff",maxentries);
  draw(dzUno, dzDue, n1, n2);
  
  uno->Draw("rOut>>rOut",cutROut+assoc,"goff",maxentries);
  due->Draw("rOut>>rOutDue",cutROut+assoc,"goff",maxentries);
  draw(routUno, routDue, n1, n2);
  
    uno->Draw("phi>>phi",cut+assoc,"goff",maxentries);
    due->Draw("phi>>phiDue",cut+assoc,"goff",maxentries);
    draw(phiUno, phiDue, n1, n2);
    
    uno->Draw("NSplit>>n",cutNoSplit+assoc,"goff",maxentries);
    due->Draw("NSplit>>nDue",cutNoSplit+assoc,"goff",maxentries);
    draw(nUno, nDue, n1, n2);

}
