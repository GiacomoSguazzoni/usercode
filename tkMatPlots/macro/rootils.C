#include <iostream>
#include <cmath>
#include <cstdlib>
#include <TROOT.h>
#include <TStyle.h>
#include <TChain.h>
#include <TProfile2D.h>
#include <TArrow.h>
#include <TColor.h>
#include <TLatex.h>

#include "rootils.h"

namespace rootils {

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

void SetFancyColorPalette(){
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };



  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void SetFancyLightGreenToLightRedPalette(){

    const Int_t NRGBs = 3;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.50, 1.00 };
    Double_t red[NRGBs]   = { 0.40, 1.00, 1.00 };
    Double_t green[NRGBs] = { 1.00, 1.00, 0.40 };
    Double_t blue[NRGBs]  = { 0.40, 1.00, 0.40 };

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


TCanvas* makeNiceCanvasByPixMargins(Int_t pixelPerBinX, Int_t pixelPerBinY, Int_t nbinx, Int_t nbiny, Int_t top, Int_t bottom, Int_t left, Int_t right) {
  
  Int_t rubaX = 4; //determinato sperimentalmente                                                                          
  Int_t rubaY = 28; //determinato sperimentalmente                                                                                                        

  TString name = generateRandomName();

  Int_t plotBaseDimX = pixelPerBinX*nbinx;
  Int_t plotBaseDimY = pixelPerBinY*nbiny;
  Int_t XX = (Int_t)(plotBaseDimX+left+right);
  Int_t YY = (Int_t)(plotBaseDimY+top+bottom);
  TCanvas* can = new TCanvas(name,name,XX+rubaX,YY+rubaY);
  can->SetTopMargin((1.*top)/(1.*YY));
  can->SetBottomMargin((1.*bottom)/(1.*YY));
  can->SetRightMargin(right/(1.*XX));
  can->SetLeftMargin(left/(1.*XX));
  can->SetBorderMode(0);
  std::cout << "Nice canvas " << XX << " * " << YY << " Margin: t " << can->GetTopMargin() << " b " << can->GetBottomMargin() << " l " << can->GetLeftMargin() << " r " << can->GetRightMargin() << std::endl;

  return can;

}

  TCanvas* makeNiceCanvasByFracMargins(Int_t pixelPerBin, Int_t nbinx, Int_t nbiny, Double_t top, Double_t bottom, Double_t left, Double_t right) {
  
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

void drawArrow(Double_t u, Double_t v, Double_t du, Double_t due, Double_t dv, Double_t dve){

  if ( du+dv ){ 

    Double_t arrowLengthSqr = du*du+dv*dv;
    Double_t relErrorOnArrowLength = sqrt(du*du*due*due+dv*dv*dve*dve)/arrowLengthSqr;
    
    Int_t color = 0; //-10, -9, -7, -4, 0
    Int_t lwid = 2;
    if ( relErrorOnArrowLength>0.5 ) color=kBlue-9;
    if ( relErrorOnArrowLength<0.4 ) color=kBlue-7;
    if ( relErrorOnArrowLength<0.3 ) color=kBlue-4;
    if ( relErrorOnArrowLength<0.2 ) color=kBlue;
    //    if ( relErrorOnArrowLength<0.1 ) color=kBlue;
    if ( relErrorOnArrowLength<0.1 ) {
      color=kBlue;
      lwid = 3;
    }
   if ( relErrorOnArrowLength<0.05 ) {
      color=kBlue;
      lwid = 4;
    }
   
   Double_t arrowSize = 0.006;

    TArrow * ar = new TArrow(u-0.5*du,v-0.5*dv,u+0.5*du,v+0.5*dv,arrowSize,"|>");
    ar->SetLineWidth(lwid);
    ar->SetFillColor(color);
    ar->SetLineColor(color);
    ar->Draw();
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
  
  
  TH2F* drawFrame(TString nameU, TString nameV, Double_t minU, Double_t maxU, Double_t minV, Double_t maxV){
    
    TString name = generateRandomName();
    TH2F *frame = new TH2F(name,name+";"+nameU+";"+nameV,1,minU,maxU,1,minV,maxV); 
    frame->SetStats(0);
    frame->SetTitle("");
    // frame->GetYaxis()->SetDecimals(1);
    //  frame->GetXaxis()->SetNdivisions(520); //Default 510
    // frame->GetYaxis()->SetNdivisions(510); //Default 510
    // frame->GetXaxis()->SetMoreLogLabels();
    //  frame->GetXaxis()->SetNoExponent();
    Int_t effCharSize = 50;
    frame->GetYaxis()->SetTitleOffset(ComputeTitleOffset(effCharSize, 2.1*effCharSize));
    //  frame->GetXaxis()->SetTitleOffset(1.);
    frame->SetLabelSize(ComputeRelSize(effCharSize),"xy");
    frame->SetTitleSize(ComputeRelSize(effCharSize),"xy");
    frame->Draw();
    
    return frame;
    
  }
  
  void drawEtaValues(Float_t etax, Float_t etay, Float_t etaxout, Float_t etayout){
    
    //Add eta labels
    Float_t etas[33] = {-3.4, -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3.0, 3.4};
    //  Float_t etax = 2940.;
    //  Float_t etay = 1240.;
    Float_t offT = (etayout - etay)*0.1;
    //  Float_t lineL = offT;
    
    for (Int_t ieta=0; ieta<33; ieta++){
      Float_t th = 2.0*atan(exp(-1.0*etas[ieta]));
      Int_t talign = 21;
      
      //IP
      if ( 0 ) {
	TLine *lineh = new TLine(-20.,0.,20.,0.); 
	lineh->Draw();  
	TLine *linev = new TLine(0.,-10.,0.,10.); 
	linev->Draw();  
      }
      
      Double_t tang=tan(th);
      Double_t x1 = etay/tang;
      Double_t y1 = etax*abs(tang);
      Double_t x2 = etayout/tang;
      Double_t y2 = etaxout*abs(tang);
      Double_t xt = (etayout+offT)/tang;
      Double_t yt = (etaxout+offT)*abs(tang);
      if ( x1 < -etax ) x1=-etax;  
      if ( x1 > etax ) x1=etax;  
      if ( y1 > etay ) y1=etay;  
      if ( x2 < -etaxout ) x2=-etaxout;  
      if ( x2 > etaxout ) x2=etaxout;  
      if ( y2 > etayout ) y2=etayout;  
      if ( xt < -(etaxout+offT) ) xt=-(etaxout+offT);  
      if ( xt > (etaxout+offT) ) xt=(etaxout+offT);  
      if ( yt > (etayout+offT) ) yt=(etayout+offT);  
      
      
      /*
	if ( etas[ieta]>-1.6 && etas[ieta]<1.6 ){
	Float_t x1 = etay/tan(th);
	Float_t y1 = etay;
	} else if ( etas[ieta]<=-1.6 ) {
	Float_t x1 = -etax;
	Float_t y1 = -etax*tan(th);
	talign = 11;
	} else if ( etas[ieta]>=1.6 ){
	Float_t x1 = etax;
	Float_t y1 = etax*tan(th);
	talign = 31;
	}
	Float_t x2 = x1+lineL*cos(th);
	Float_t y2 = y1+lineL*sin(th);
      */
      
      //    Float_t xt = x2;
      //    Float_t yt = y2+offT;
      //      cout << isign << " " << th*180./pi << " " << x1 << " " << y1 << "\n";
      TLine *line1 = new TLine(x1,y1,x2,y2); 
      line1->SetLineColor(kBlue);
      line1->Draw();  
      char text[20];
      int rc = sprintf(text, "%3.1f", etas[ieta]);
      TLatex *t1=0;
      if ( etas[ieta] == 0 ) {
	t1 = new TLatex(xt,yt,"#eta = 0"); 
      } else {
	t1 = new TLatex(xt,yt,text); 
      }
      t1->SetTextColor(kBlue);
      t1->SetTextSize(0.03);
      t1->SetTextAlign(talign);
      t1->Draw();     
      
    }
    
  }


  void drawEtaValuesOnThetaPlot(Float_t etay, Float_t etayout, Float_t etaLimit, Float_t textFact){
    
    //Add eta labels
    Float_t etas[33] = {-3.4, -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3.0, 3.4};
    //  Float_t etax = 2940.;
    //  Float_t etay = 1240.;
    Float_t offT = (etayout - etay);
    //  Float_t lineL = offT;
    
    for (Int_t ieta=0; ieta<33; ieta++){

      if ( abs(etas[ieta]) > etaLimit ) continue;

      Float_t th = 2.0*atan(exp(-1.0*etas[ieta]));
      Int_t talign = 21;
      
      Float_t xt = th;
      Float_t yt = etay+offT;
      Float_t ylu = etay+offT*0.8;
      Float_t yld = etay;

      
      TLine *line1 = new TLine(xt,yld,xt,ylu); 
      line1->SetLineColor(kBlue);
      line1->Draw();  
      char text[20];
      int rc = sprintf(text, "%3.1f", etas[ieta]);
      TLatex *t1=0;
      if ( etas[ieta] == 0 ) {
	t1 = new TLatex(xt,yt,"#eta = 0"); 
      } else {
	t1 = new TLatex(xt,yt,text); 
      }
      t1->SetTextColor(kBlue);
      t1->SetTextSize(0.03*textFact);
      t1->SetTextAlign(talign);
      t1->Draw();     
      
    }
    
  }

  Double_t max(Double_t arg1, Double_t arg2){

    if ( arg1 > arg2 ) return arg1;
    
    return arg2;

  }

  Double_t min(Double_t arg1, Double_t arg2){

    if ( arg1 < arg2 ) return arg1;
    
    return arg2;

  }

  //
}

