#include <iostream>
#include <cstdlib>

#include <TROOT.h>
#include <TChain.h>
#include <TProfile2D.h>
#include <TString.h>
#include <TCanvas.h>
#include <TArrow.h>

void SetFancyGrayscalePalette(){
  const Int_t NRGBs = 3;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.5, 1.00 };
  Double_t red[NRGBs]   = { .75, 0.1, 0.0 };
  Double_t green[NRGBs] = { .75, 0.1, 0.0 };
  Double_t blue[NRGBs]  = { .75, 0.1, 0.0 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}


TString generateRandomName(){
  //  TRandom3 ran;                                                                                                                                       
  //  Int_t nomeNumero = (Int_t)(ran.Uniform(0.,99999.));                                                                                                 
  Int_t nomeNumero = (Int_t)(std::rand());

  char nome[20];
  sprintf(nome,"h%d",nomeNumero);
  return TString(nome);
}


TCanvas* makeNiceCanvas(Int_t pixelPerBin, Int_t nbinx, Int_t nbiny, Int_t top, Int_t bottom, Int_t left, Int_t right) {
  
  Int_t rubaX = 4; //determinato sperimentalmente                                                                          
  Int_t rubaY = 28; //determinato sperimentalmente                                                                                                        
  TString name = generateRandomName();

  Int_t plotBaseDimX = pixelPerBin*nbinx;
  Int_t plotBaseDimY = pixelPerBin*nbiny;
  Int_t XX = (plotBaseDimX+left+right+rubaX);
  Int_t YY = (plotBaseDimY+top+bottom+rubaY);
  TCanvas* can = new TCanvas(name,name,XX,YY);
  can->SetTopMargin(1.*top/YY);
  can->SetBottomMargin(1.*bottom/YY);
  can->SetRightMargin(1.*right/XX);
  can->SetLeftMargin(1.*left/XX);
  can->SetBorderMode(0);
  std::cout << "Nice canvas " << XX << " * " << YY << std::endl;
  return can;

}

/*
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
  can->SetBorderMode(0);
  std::cout << "Nice canvas " << XX << " * " << YY << std::endl;
  return can;

}
*/

void drawArrow(Double_t u, Double_t v, Double_t du, Double_t due, Double_t dv, Double_t dve){

  if ( du+dv ){ 

    Double_t arrowLengthSqr = du*du+dv*dv;
    Double_t relErrorOnArrowLength = sqrt(du*du*due*due+dv*dv*dve*dve)/arrowLengthSqr;
    
    Int_t color = 0; //-10, -9, -7, -4, 0
    Int_t lwid = 1;
    if ( relErrorOnArrowLength>0.5 ) color=kBlue-10;
    if ( relErrorOnArrowLength<0.4 ) color=kBlue-9;
    if ( relErrorOnArrowLength<0.3 ) color=kBlue-7;
    if ( relErrorOnArrowLength<0.2 ) color=kBlue-4;
    if ( relErrorOnArrowLength<0.1 ) color=kBlue;
    if ( relErrorOnArrowLength<0.05 ) {
      color=kBlue;
      lwid = 2;
    }

    TArrow * ar = new TArrow(u-0.5*du,v-0.5*dv,u+0.5*du,v+0.5*dv,0.005,">");
    ar->SetLineWidth(lwid);
    ar->SetLineColor(color);
    ar->Draw();
  }

}

void drawArrowPlot(TChain *mc, TString varU, TString deltaU, TString varV, TString deltaV, Int_t nbinU, Double_t minU, Double_t maxU, Int_t nbinV, Double_t minV, Double_t maxV){

  TString uvname = generateRandomName();
  TString uname = generateRandomName();
  TString vname = generateRandomName();
  TH2D * UV = new TH2D(uvname,uvname,10*nbinU,minU,maxU,10*nbinV,minV,maxV);
  TProfile2D * dU = new TProfile2D(uname,uname,nbinU,minU,maxU,nbinV,minV,maxV);
  TProfile2D * dV = new TProfile2D(vname,vname,nbinU,minU,maxU,nbinV,minV,maxV);

  mc->Draw(varV+":"+varU+">>"+uvname,"(isAssoc==1)","goff");
  mc->Draw(deltaU+":"+varV+":"+varU+">>"+uname,"(isAssoc==1)","goffprof");
  mc->Draw(deltaV+":"+varV+":"+varU+">>"+vname,"(isAssoc==1)","goffprof");
  
  SetFancyGrayscalePalette();
  UV->Draw("colsame");


  for (Int_t iU=1; iU<nbinU+1; iU++){
    for (Int_t iV=1; iV<nbinV+1; iV++){
      
      Double_t uu = dU->GetXaxis()->GetBinCenter(iU);
      Double_t vv = dU->GetYaxis()->GetBinCenter(iV);
      Double_t du = dU->GetBinContent(iU,iV);
      Double_t dv = dV->GetBinContent(iU,iV);
      Double_t due = dU->GetBinError(iU,iV);
      Double_t dve = dV->GetBinError(iU,iV);
      
      drawArrow(uu,vv,du,due,dv,dve);

    }
  }
}

Double_t ComputeRelSize(Int_t pxl){

  std::cout << "vvvvvv" << gPad->GetWh() << " dsa " << gPad->GetAbsHNDC() << std::endl;

  return (1.*pxl)/(1.*gPad->GetWh()*gPad->GetAbsHNDC());

}

Double_t ComputeTitleOffset(Int_t fontpxl, Int_t pxl){
  //
  //Sembra che dipenda dall'aspect ratio del plot in pxl: distanza ~ size * larghezza/altezza 
  //
  // con plot quadrato e pxl = 25 ==> 45
  // con plot quadrato e pxl = 50 ==> 90
  
  Double_t y1 = 45.; 
  Double_t y2 = 90.;
  Double_t x1 = 25.; 
  Double_t x2 = 50.;

  Double_t nomDistance = (y2-y1)*(1.*fontpxl-x1)/(1.*(x2-x1))+y1;

  nomDistance=nomDistance*(1.*gPad->GetWw()*gPad->GetAbsWNDC())/(1.*gPad->GetWh()*gPad->GetAbsHNDC());

  
  std::cout << " nomDis " << nomDistance << " " << 1.*pxl/nomDistance << std::endl;
  return 1.*pxl/nomDistance;

  //return 1.;

}


void drawFrame(TString nameU, TString nameV, Double_t minU, Double_t maxU, Double_t minV, Double_t maxV, Int_t canH){

  TString name = generateRandomName();
  TH2F *frame = new TH2F(name,name+";"+nameU+";"+nameV,1,minU,maxU,1,minV,maxV); 
  frame->SetStats(0);
  frame->SetTitle("");
  // frame->GetYaxis()->SetDecimals(1);
  //  frame->GetXaxis()->SetNdivisions(520); //Default 510
  // frame->GetYaxis()->SetNdivisions(510); //Default 510
  // frame->GetXaxis()->SetMoreLogLabels();
  //  frame->GetXaxis()->SetNoExponent();
  frame->GetYaxis()->SetTitleOffset(ComputeTitleOffset(25, 55));
  //  frame->GetXaxis()->SetTitleOffset(1.);
  frame->SetLabelSize(ComputeRelSize(25),"xy");
  frame->SetTitleSize(ComputeRelSize(25),"xy");
  frame->Draw();

}

void arrowPlot(){

  gROOT->Reset();
  gROOT->SetStyle("Plain");

  TChain *mc = new TChain("ntupleR2S");
  //MC reco
  mc->Add("./ntuple_conversion_MinBiasMC_01.root");
  mc->Add("./ntuple_conversion_MinBiasMC_02.root");
  mc->Add("./ntuple_conversion_MinBiasMC_03.root");
  mc->Add("./ntuple_conversion_MinBiasMC_04.root");

 
  TString rSIM = "sqrt((x+deltax)*(x+deltax)+(y+deltay)*(y+deltay))";
  TString xSIM = "(x+deltax)";
  TString ySIM = "(y+deltay)";
  TString zSIM = "(z+deltaz)";
  TString deltar = "sqrt((deltax)*(deltax)+(deltay)*(deltay))";
  TString deltax = "(deltax)";
  TString deltay = "(deltay)";
  TString deltaz = "(deltaz)";

  //Z
  Int_t nbinU = 64;
  Double_t minU = -80.;
  Double_t maxU = 80.;

  //R
  Int_t nbinV = 26;
  Double_t minV = 0.;
  Double_t maxV = 65.;

  TString varU = zSIM;
  TString varV = rSIM;
  TString deltaU = deltaz;
  TString deltaV = deltar;

  //  TCanvas * canrz = makeNiceCanvas(15, nbinU, nbinV, 0.04, 0.11, 0.07, 0.02);
  TCanvas * canrz = makeNiceCanvas(15, nbinU, nbinV, 20, 60, 70, 20);
  canrz->cd();
  drawFrame("true Z (cm)", "true R (cm)", minU, maxU, minV, maxV, canrz->GetWh());
  drawArrowPlot(mc, varU, deltaU, varV, deltaV, nbinU, minU, maxU, nbinV, minV, maxV);
  canrz->SaveAs("RZ_Arrow.png");

  //X
  nbinU = 52;
  minU = -65.;
  maxU = 65.;
  
  //Y
  nbinV = 52;
  minV = -65.;
  maxV = 65.;
  
  varU = xSIM;
  varV = ySIM;
  deltaU = deltax;
  deltaV = deltay;

  TCanvas * canxy = makeNiceCanvas(15, nbinU, nbinV, 20, 60, 70, 20);
  canxy->cd();
  drawFrame("true X (cm)", "true Y (cm)", minU, maxU, minV, maxV, canxy->GetWh());
  drawArrowPlot(mc, varU, deltaU, varV, deltaV, nbinU, minU, maxU, nbinV, minV, maxV);
  canxy->SaveAs("XY_Arrow.png");

}
