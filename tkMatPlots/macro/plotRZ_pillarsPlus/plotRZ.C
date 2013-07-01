#include <algorithm>

Double_t GetMaxWithinRadius(TH2D * H, Double_t radius){

  Int_t nbinx = H->GetNbinsX();
  Int_t nbiny = H->GetNbinsY();

  Double_t max = 0.;

  for(Int_t i=1; i<nbinx+1; i++){
    for(Int_t j=1; j<nbiny+1; j++){
   
      Double_t x = H->GetXaxis()->GetBinCenter(i);
      Double_t y = H->GetYaxis()->GetBinCenter(j);

      if ( (x*x+y*y) < (radius*radius) ){
	Double_t value = H->GetBinContent(i,j);
        if ( value > max ) max = value;
      }

    }
  }

  return max;

} 

Double_t GetIntWithinRZ(TH2D * H, Double_t radius, Double_t Z){

  Int_t nbinz = H->GetNbinsX(); //z
  Int_t nbinr = H->GetNbinsY(); //R

  Double_t Int = 0.;

  for(Int_t i=1; i<nbinz+1; i++){
    for(Int_t j=1; j<nbinr+1; j++){
   
      Double_t z = H->GetXaxis()->GetBinCenter(i);
      Double_t r = H->GetYaxis()->GetBinCenter(j);

      if ( r < radius && z > -Z && z < Z ) ){
	Int += H->GetBinContent(i,j);
      }

    }
  }

  return Int;

} 

void writePrel(TH2* frame, const char* suffix){
  
  TString Suff(suffix);
  
  Float_t x0 = frame->GetXaxis()->GetXmin();
  Float_t x1 = frame->GetXaxis()->GetXmax();
  Float_t y0 = frame->GetYaxis()->GetXmin();
  Float_t y1 = frame->GetYaxis()->GetXmax();
  
  Float_t par = -0.08;

  //  prel = new TLatex(x0+0.1*(x1-x0),(max-min)*par,"CMS Preliminary 2010"+Suff);
  prel = new TLatex(x0+par*(x1-x0),y0,"CMS Preliminary 2010");
  prel->SetTextSize(.055);
  prel->SetTextAngle(90);
  prel->Draw();

  par = -0.17;
  label = new TLatex(x0,y0+par*(y1-y0),Suff);
  label->SetTextSize(.065);
  label->Draw();

}

void Plot2D(TH2D * H, char* label, TString dir){
  
  //Canvas
  
  Float_t hSide =H->GetXaxis()->GetXmax()-H->GetXaxis()->GetXmin();
  Float_t vSide =H->GetYaxis()->GetXmax()-H->GetYaxis()->GetXmin();;
  Float_t asRatio = hSide/vSide;

  cout << hSide << " " << vSide << endl;

  Int_t sqX = 1200;
  Int_t sqY = 1000;

  TCanvas *xyCan = new TCanvas("xyCan","xyCan",((Int_t)(sqX*asRatio)),sqY);
  //  TCanvas *xyCan = new TCanvas("xyCan","xyCan",1200,((Int_t)(1000/asRatio)));
  xyCan->SetLeftMargin(0.1);
  xyCan->SetRightMargin(0.2);
    //xyCan->SetLeftMargin(0.);
  xyCan->SetTopMargin(0.04);
  xyCan->SetBottomMargin(0.16);

  H->SetTitle();
  H->GetXaxis()->SetTitle("z (cm)");
  H->GetYaxis()->SetTitle("R (cm)");
  H->GetZaxis()->SetTitle("#propto N_{#gamma} #times P/X_{0}/event");
  H->GetZaxis()->SetTitleOffset(0.9);
  H->SetTitleSize(0.08,"xy");
  H->SetTitleSize(0.07,"z");
  H->SetTitleOffset(0.85,"x");
  H->SetTitleOffset(0.5,"y");
  H->SetLabelSize(0.07,"xyz");
  H->SetMinimum(0.);
  H->SetContour(255);
  H->SetStats(kFALSE);
  H->Draw("zcol");

  writePrel(H,label);

  xyCan->SaveAs(dir+TString(H->GetName())+".png");

}

void makeColorTable(){


  /*
  const Int_t NRGBs = 3;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.30, 1.00 };
  Double_t red[NRGBs]   = { 0.96, 0.20, 0.00 };
  Double_t green[NRGBs] = { 0.96, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.96, 0.20, 0.00 };
  */

  const Int_t NRGBs = 2;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 1.00 };
  Double_t red[NRGBs]   = { 0.99, 0.00 };
  Double_t green[NRGBs] = { 0.99, 0.00 };
  Double_t blue[NRGBs]  = { 0.99, 0.00 };

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

}

void makeColorTableRB(){


  const Int_t NRGBs = 3;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.50, 1.00 };
  Double_t red[NRGBs]   = { 0.00, .90, 1.00 };
  Double_t green[NRGBs] = { 0.00, .90, 0.00 };
  Double_t blue[NRGBs]  = { 1.00, .90, 0.00 };

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

}

void plotRZ(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetHistMinimumZero(kFALSE);

  //  TString dir("allStat_6hit/");
  //  TString dir("allStat_6hit_noEvtBug_wEndCap/");
  //  TString dir("allStat_6hit_noWeiBug/");
  TString dir("./");

  //Agguanta gli istogrammi
  TFile * SimF = new TFile(dir+"Sim_RZ.root");
  TH2D * SimH = (TH2D*) SimF->Get("Sim_RZ");;
  //
  TFile * MCF = new TFile(dir+"MC_RZ.root");
  TH2D * MCH = (TH2D*) MCF->Get("MC_RZ");;
  //
  //  TFile * MCFsF = new TFile(dir+"rMCFs_RZ.root");
  //  TH2D * MCFsH = (TH2D*) MCFsF->Get("rMCFs_RZ");;
  //
  TFile * FakeF = new TFile(dir+"MCFake_RZ.root");
  TH2D * FakeH = (TH2D*) FakeF->Get("MCFake_RZ");;
  //
  TFile * DataF = new TFile(dir+"Data_RZ.root");
  TH2D * DataH = (TH2D*) DataF->Get("Data_RZ");;
  //
  //  TFile * DataFsF = new TFile(dir+"rDataFs_RZ.root");
  //  TH2D * DataFsH = (TH2D*) DataFsF->Get("rDataFs_RZ");;
  //
  
  cout << MCH->GetEntries() << " " << MCH->GetEffectiveEntries() << " " << MCH->Integral() << endl;
  cout << DataH->GetEntries() << " " << DataH->GetEffectiveEntries() << " " << DataH->Integral() << endl;

  Float_t scaleFact = 1./23876690; //2010-06-04
  //  Float_t scaleFact = 1./10809460; //2010-05-10
  //Float_t scaleFact = 1.;
  SimH->Scale(scaleFact); 
  MCH->Scale(scaleFact);
  //  MCFsH->Scale(scaleFact);
  FakeH->Scale(scaleFact);
  DataH->Scale(scaleFact);
  //  DataFsH->Scale(scaleFact);

  cout << MCH->GetEntries() << " " << MCH->GetEffectiveEntries() << " " << MCH->Integral() << endl;
  cout << DataH->GetEntries() << " " << DataH->GetEffectiveEntries() << " " << DataH->Integral() << endl;

  //Max
  Double_t simMax=SimH->GetMaximum();
  Double_t simPxlMax=GetMaxWithinRadius(SimH, 15.);
  Double_t simPxlInt=GetIntWithinRZ(SimH, 15., 26.);
  Double_t MCMax=MCH->GetMaximum();
  Double_t MCPxlMax=GetMaxWithinRadius(MCH, 15.);
  Double_t MCPxlInt=GetIntWithinRZ(MCH, 15., 26.);
  Double_t DataMax=DataH->GetMaximum();
  Double_t DataPxlMax=GetMaxWithinRadius(DataH, 15.);
  Double_t DataPxlInt=GetIntWithinRZ(DataH, 15., 26.);
  Double_t mAx=max(MCMax,DataMax);


  std::cout << " max     >>>> mc " << MCMax << " data " << DataMax << " Sim " << simMax << std::endl;
  std::cout << " pxl max >>>> mc " << MCPxlMax << " data " << DataPxlMax << " Sim " << simPxlMax << std::endl;
  std::cout << " pxl int >>>> mc " << MCPxlInt << " data " << DataPxlInt << " Sim " << simPxlInt << std::endl;

  // This would be need to normalize SimH to the pixel region material (taking into account that SimH and DataH/MCH have different binning!)
  //
  //  SimH->Scale(DataPxlInt*(DataH->GetXaxis()->GetBinWidth(1)*DataH->GetYaxis()->GetBinWidth(1))/simPxlInt/(SimH->GetXaxis()->GetBinWidth(1)*SimH->GetYaxis()->GetBinWidth(1)));
  // I prefer to have Sim plots and Data/MC plot to look similar by appropriately scaling the SimH range, not by scaling the histo

  //  mAx=0.03;
  //  mAx=0.011;
  mAx=0.7*0.5*(MCPxlMax+DataPxlMax); //should be 0.8*... but I put something different to make the plot similar to NI

  Double_t simRangeScale = 1./(DataPxlInt*(DataH->GetXaxis()->GetBinWidth(1)*DataH->GetYaxis()->GetBinWidth(1))/simPxlInt/(SimH->GetXaxis()->GetBinWidth(1)*SimH->GetYaxis()->GetBinWidth(1)));

  SimH->SetMaximum(mAx*simRangeScale); 
  MCH->SetMaximum(mAx);  
  DataH->SetMaximum(mAx);

  //

  makeColorTable();

  Plot2D(SimH, "#gamma conv., MC Truth", dir);
  Plot2D(MCH, "#gamma conv., MC Reco #sqrt{s}=7TeV", dir);
  Plot2D(DataH, "#gamma conv., Data #sqrt{s}=7TeV", dir);

  //

  /*
  TH2D *MCCoarse = MCH->Rebin2D(2,2,"MCCoarse");
  TH2D *DataCoarse = DataH->Rebin2D(2,2,"DataCoarse");
  */

  /*
  cout << MCCoarse->GetEntries() << " " << MCCoarse->GetEffectiveEntries() << endl;
  cout << DataCoarse->GetEntries() << " " << DataCoarse->GetEffectiveEntries() << endl;
  */

  /*
  TH2D *diff = new TH2D("diff","diff", MCCoarse->GetNbinsX(),-60.,60.,MCCoarse->GetNbinsY(),-60.,60.);
  for (Int_t iBin = 0; iBin < MCCoarse->GetNbinsX(); iBin++){
    for (Int_t jBin = 0; jBin < MCCoarse->GetNbinsY(); jBin++){
      Double_t mc = MCCoarse->GetBinContent(iBin+1, jBin+1); 
      Double_t data = DataCoarse->GetBinContent(iBin+1, jBin+1); 
      cout << mc << " " << data << endl;
      if ( mc ) diff->SetBinContent(iBin+1,jBin+1,0.5+(mc-data)/mc);
    }
  }
  */

  MCH->Rebin2D(2,2);
  DataH->Rebin2D(2,2);

  TH2D *diff = new TH2D("diff","diff", MCH->GetNbinsX(),-60.,60.,MCH->GetNbinsY(),-60.,60.);
  for (Int_t iBin = 0; iBin < MCH->GetNbinsX(); iBin++){
    for (Int_t jBin = 0; jBin < MCH->GetNbinsY(); jBin++){
      Double_t mc = MCH->GetBinContent(iBin+1, jBin+1); 
      Double_t data = DataH->GetBinContent(iBin+1, jBin+1); 
      if ( mc+data ) {
      if ( ! mc ) mc=0.0000001;
      Double_t val = data/mc;
      //      Double_t val = (mc-data)/mc;
      Double_t valup = 1.2;
      Double_t vallo = 0.8;
      Double_t offset = 0.;
      if ( val > valup ) diff->SetBinContent(iBin+1,jBin+1,offset+valup);
      if ( val < vallo ) diff->SetBinContent(iBin+1,jBin+1,offset+vallo);
      if ( val > vallo && val < valup ) diff->SetBinContent(iBin+1,jBin+1,offset+val);
      }
    }
  }


  diff->SetMinimum(0.);
  diff->SetMaximum(2.);

  makeColorTableRB();

  Plot2D(diff,"Ratio", dir);

}
