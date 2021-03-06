void writePrel(TH2* frame, Float_t par = 1.03){
  
  Float_t  x0 = frame->GetXaxis()->GetXmin();
  Float_t  x1 = frame->GetXaxis()->GetXmax();
  Float_t min = frame->GetYaxis()->GetXmin();
  Float_t max = frame->GetYaxis()->GetXmax();

  prel = new TLatex(x1-0.01*(x1-x0),(max-min)*par,"CMS Preliminary 2010, #sqrt{s}=7TeV");
  prel->SetTextAlign(31);
  prel->SetTextSize(0.08);
  prel->Draw();

}

void writeLabels(TH2* frame){
  
  Float_t  x0 = frame->GetXaxis()->GetXmin();
  Float_t  x1 = frame->GetXaxis()->GetXmax();
  Float_t min = frame->GetYaxis()->GetXmin();
  Float_t max = frame->GetYaxis()->GetXmax();

  cout << " Lab " << x0 << " " << x1 << " " << min << " " << max << endl;

  Float_t parlo = 0.5;
  Float_t parup = 0.6;

  Int_t n=11;
  Float_t off = 0.25;
  Float_t  binR[11]={29.,43.,71.,102.,200.,256.,339.,419.,500.,550.,609.};
  Float_t  paroff[11]={off,off,off,0.2*off,0.,0.,0.,0.,0.,0.,0.};
  TString  lab[11]={"BP","PXL1","PXL2","PXL3","support", "TIB1","TIB2","TIB3","TIB4","support","TOB1"};
  
  for (Int_t i = 0; i < n; i++){
    
    Float_t par = parlo+paroff[i];
    if ( i % 2 ) par = parup+paroff[i];

    prel = new TLatex(binR[i]/10.,(max-min)*par,lab[i]);
    prel->SetTextAlign(22);
    prel->SetTextSize(0.08);
    prel->Draw();
    
  }

}

void dataPlot(TH1D * H, Int_t colI){

  H->SetLineColor(colI);
  H->SetStats(kFALSE);
  H->SetLineWidth(1);
  H->SetMarkerStyle(20);
  H->SetMarkerSize(.8);
  H->SetMarkerColor(colI);
  H->Draw("Psame");

}

void mcPlot(TH1D * H, Int_t colI){

  H->SetLineWidth(2);
  H->SetLineColor(colI+1);
  H->SetFillColor(colI+6);
  //  H->SetFillColor(kWhite);
  H->SetMarkerStyle(20);
  H->SetMarkerSize(.6);
  H->SetMarkerColor(colI+1);
  H->SetStats(kFALSE);
  H->Draw("PE2same");
 
}

void plotR_zones_NI(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  //Agguanta gli istogrammi
  //
  //Raw
  TFile * rMCF = new TFile("rMC_R.root");
  TH1D * rMCH = (TH1D*) rMCF->Get("rMC_R");;
  //
  TFile * rSimF = new TFile("rSim_R.root");
  TH1D * rSimH = (TH1D*) rSimF->Get("rSim_R");;
  //
  TFile * rDataF = new TFile("rData_R.root");
  TH1D * rDataH = (TH1D*) rDataF->Get("rData_R");;
  //
  TFile * rFakeF = new TFile("rFake_R.root");
  TH1D * rFakeH = (TH1D*) rFakeF->Get("rFake_R");;
  //
  //Scaled
  TFile * SimF = new TFile("Sim_R.root");
  TH1D * SimH = (TH1D*) SimF->Get("Sim_R");;
  //
  TFile * MCF = new TFile("MC_R.root");
  TH1D * MCH = (TH1D*) MCF->Get("MC_R");;
  //
  TFile * MCFakeSubF = new TFile("MCFakeSub_R.root");
  TH1D * MCFakeSubH = (TH1D*) MCFakeSubF->Get("MCFakeSub_R");;
  //
  TFile * FakeF = new TFile("Fake_R.root");
  TH1D * FakeH = (TH1D*) FakeF->Get("Fake_R");;
  //
  TFile * DataF = new TFile("Data_R.root");
  TH1D * DataH = (TH1D*) DataF->Get("Data_R");;
  //
  TFile * DataFakeSubF = new TFile("DataFakeSub_R.root");
  TH1D * DataFakeSubH = (TH1D*) DataFakeSubF->Get("DataFakeSub_R");;
  //
  //Unfold Stuff
  TFile * UfBay1DataFakeSubF = new TFile("rUfBay1DataFakeSub.root");
  TH1D * UfBay1DataFakeSubH = (TH1D*) UfBay1DataFakeSubF->Get("rUfBay1DataFakeSub");;
  //
  TFile * UfBay2DataFakeSubF = new TFile("rUfBay2DataFakeSub.root");
  TH1D * UfBay2DataFakeSubH = (TH1D*) UfBay2DataFakeSubF->Get("rUfBay2DataFakeSub");;
  //
  TFile * UfBay3DataFakeSubF = new TFile("rUfBay3DataFakeSub.root");
  TH1D * UfBay3DataFakeSubH = (TH1D*) UfBay3DataFakeSubF->Get("rUfBay3DataFakeSub");;
  //
  TFile * UfBay4DataFakeSubF = new TFile("rUfBay4DataFakeSub.root");
  TH1D * UfBay4DataFakeSubH = (TH1D*) UfBay4DataFakeSubF->Get("rUfBay4DataFakeSub");;
  //
  TFile * UfBay1MCFakeSubF = new TFile("rUfBay1MCFakeSub.root");
  TH1D * UfBay1MCFakeSubH = (TH1D*) UfBay1MCFakeSubF->Get("rUfBay1MCFakeSub");;
  //
  TFile * UfBay2MCFakeSubF = new TFile("rUfBay2MCFakeSub.root");
  TH1D * UfBay2MCFakeSubH = (TH1D*) UfBay2MCFakeSubF->Get("rUfBay2MCFakeSub");;
  //
  TFile * UfBay3MCFakeSubF = new TFile("rUfBay3MCFakeSub.root");
  TH1D * UfBay3MCFakeSubH = (TH1D*) UfBay3MCFakeSubF->Get("rUfBay3MCFakeSub");;
  //
  TFile * UfBay4MCFakeSubF = new TFile("rUfBay4MCFakeSub.root");
  TH1D * UfBay4MCFakeSubH = (TH1D*) UfBay4MCFakeSubF->Get("rUfBay4MCFakeSub");;
  //
  
  //
  Float_t nData = rDataH->GetEntries();
  Float_t nMC = rMCH->GetEntries();
  Float_t nFake = rFakeH->GetEntries();
  Float_t nDataE = sqrt(nData);
  Float_t nMCE =  sqrt(nMC);
  Float_t nFakeE = sqrt(nFake);
  Float_t nMCFakeSubE = sqrt(nMC-nFake);
  Float_t iData = rDataH->Integral();
  Float_t iMC = rMCH->Integral();
  Float_t iFake = rFakeH->Integral();
  Float_t fData = iData/nData;
  Float_t fMC = iMC/nMC;
  Float_t iDataE = nDataE*fData;
  Float_t iMCE = nMCE*fMC;
  Float_t iFakeE = nFakeE*fMC;
  //
  std::cout << " d " << nData << "+-" << nDataE << " mc " << nMC << "+-" << nMCE << " f " << nFake << "+-" << nFakeE << std::endl;  std::cout << " d " << iData << "+-" << iDataE << " mc " << iMC << "+-" << iMCE << " f " << iFake << "+-" << iFakeE << std::endl;
  char sData[30];
  sprintf(sData,"%.0f#pm%.0f",iData,iDataE);
  TString TSsData(sData);
  //
  char sDataFakeSub[30];
  sprintf(sDataFakeSub,"%.0f#pm%.0f",iData-iFake,sqrt(iDataE*iDataE+iFakeE*iFakeE));
  TString TSsDataFakeSub(sDataFakeSub);
  //
  char sMC[30];
  sprintf(sMC,"%.0f#pm%.0f",iMC,iMCE);
  TString TSsMC(sMC);
  char sMCFakeSub[30];
  sprintf(sMCFakeSub,"%.0f#pm%.0f",iMC-iFake,nMCFakeSubE*fMC);
  TString TSsMCFakeSub(sMCFakeSub);
  //
  //
  char sFake[30];
  sprintf(sFake,"%.0f#pm%.0f",iFake,iFakeE);
  TString TSsFake(sFake);
  //
  Float_t scaleFact = 1./23876690; //2010-06-04
    //    Float_t scaleFact = 1./10809460; //2010-05-10
    //  Float_t scaleFact = 1.;
  rDataH->Scale(scaleFact);
  rMCH->Scale(scaleFact);
  rFakeH->Scale(scaleFact);
  //
  SimH->Scale(scaleFact); 
  MCH->Scale(scaleFact);
  MCFakeSubH->Scale(scaleFact);
  FakeH->Scale(scaleFact);
  DataH->Scale(scaleFact);
  DataFakeSubH->Scale(scaleFact);
  UfBay1DataFakeSubH->Scale(scaleFact);
  UfBay2DataFakeSubH->Scale(scaleFact);
  UfBay3DataFakeSubH->Scale(scaleFact);
  UfBay4DataFakeSubH->Scale(scaleFact);
  UfBay1MCFakeSubH->Scale(scaleFact);
  UfBay2MCFakeSubH->Scale(scaleFact);
  UfBay3MCFakeSubH->Scale(scaleFact);
  UfBay4MCFakeSubH->Scale(scaleFact);

  //Canvas
  TCanvas *rCan1 = new TCanvas("rCan1","rCan1",1200,600);
  TCanvas *rCan2 = new TCanvas("rCan2","rCan2",1200,600);
  TCanvas *rCan3 = new TCanvas("rCan3","rCan3",1200,600);
  rCan1->SetRightMargin(0.05);
  rCan2->SetRightMargin(0.05);
  rCan3->SetRightMargin(0.05);

  //Frame1
  Double_t Min = 0.00008*DataH->GetBinContent(DataH->GetMaximumBin());;
  Double_t Max = 1.*DataH->GetBinContent(DataH->GetMaximumBin());
  //  Double_t Min = 0.000000009;
  //  Double_t Max = 0.02;
  TH2F *frame1= new TH2F("frame1",";R (cm);Countings/event/dR",1.,rDataH->GetXaxis()->GetBinLowEdge(1),rDataH->GetXaxis()->GetBinUpEdge(rDataH->GetNbinsX()),1.,Min,Max);
  frame1->SetStats(kFALSE);
  frame1->SetTitleOffset(1.1,"y");

  //Plot on the first canvas
  //
  //Raw stuff
  rCan1->cd();
  rCan1->SetLogy();
  frame1->Draw();

  rMCH->SetLineColor(kBlue+1);
  rMCH->SetFillColor(kAzure+6);
  rMCH->SetMarkerStyle(20);
  rMCH->SetMarkerSize(.8);
  rMCH->SetMarkerColor(kBlue+1);
  rMCH->SetStats(kFALSE);
  rMCH->Draw("PE2same");

  rFakeH->SetLineColor(kRed+1);
  rFakeH->SetFillColor(kRed-6);
  rFakeH->SetMarkerStyle(20);
  rFakeH->SetMarkerSize(.8);
  rFakeH->SetMarkerColor(kRed+1);
  rFakeH->SetStats(kFALSE);
  rFakeH->Draw("PE2same");

  rDataH->SetLineColor(kBlack);
  rDataH->SetStats(kFALSE);
  rDataH->SetLineWidth(1);
  rDataH->SetMarkerStyle(20);
  rDataH->SetMarkerSize(.8);
  rDataH->SetMarkerColor(kBlack);
  //  rDataH->Draw("PE0X0same");
  rDataH->Draw("Psame");


  //creating a legend
  TLegend *legend = 0;
  legend = new TLegend(0.3,0.77,0.89,0.89);
  legend->SetNColumns(2);
  legend->SetFillColor(0);
  legend->SetBorderSize(0.);
  TH1D *rMCHLeg   = new TH1D("rMCHLeg"  ,"rMCH", 1, 0, 1.); 
  rMCHLeg->SetLineColor(kBlue+1);
  rMCHLeg->SetFillColor(kAzure+6);
  rMCHLeg->SetMarkerStyle(20);
  rMCHLeg->SetMarkerSize(.8);
  rMCHLeg->SetMarkerColor(kBlue+1);
  TH1D *rFakeHLeg = new TH1D("rFakeHLeg","rFakeH", 1, 0, 1.);
  rFakeHLeg->SetLineColor(kRed+1);
  rFakeHLeg->SetFillColor(kRed-6);
  rFakeHLeg->SetMarkerStyle(20);
  rFakeHLeg->SetMarkerSize(.8);
  rFakeHLeg->SetMarkerColor(kRed+1);
  //  legend->AddEntry(rSimH, "MC truth", "f");
  legend->AddEntry(rMCHLeg, "MC pseudo-data, "+TSsMC, "fp");
  legend->AddEntry(rFakeHLeg, "Fakes, "+TSsFake, "fp");
  legend->AddEntry(DataH, "Data 7TeV, "+TSsData, "pl");
  
  legend->Draw();
  writePrel(frame1, 2.);


  //////////////////////////////////////////////////////////////////////////////
  //Frame2
  //  Double_t Min = 0.1;
  Min = 0.;
  Max = 1.8*SimH->GetBinContent(SimH->GetMaximumBin());
  TH2F *frame2= new TH2F("frame2",";R (cm);#propto N_{hadron} #times P/#lambda_{I}/event",1.,DataH->GetXaxis()->GetBinLowEdge(1),DataH->GetXaxis()->GetBinUpEdge(DataH->GetNbinsX()),1.,Min,Max);
  //  TH2F *frame2= new TH2F("frame2",";R (cm);X_{0}^{-1} (arbitrary units)",1.,DataH->GetXaxis()->GetBinLowEdge(1),DataH->GetXaxis()->GetBinUpEdge(DataH->GetNbinsX()),1.,Min,Max);
  frame2->SetStats(kFALSE);
  frame2->SetTitleOffset(1.1,"y");

  //
  //Plot on the second canvas
  rCan2->cd();
  frame2->Draw();

  writeLabels(frame2);

  SimH->SetLineColor(kGreen+2);
  SimH->SetFillColor(kSpring+1);
  SimH->SetLineWidth(1);
  SimH->SetStats(kFALSE);
  SimH->Draw("HISTsame");

  Int_t colI = kBlue;
  //  mcPlot(UfBay1MCFakeSubH, colI);
  //  mcPlot(UfBay2MCFakeSubH, colI+1);
  mcPlot(UfBay3MCFakeSubH, colI+2);
  //  mcPlot(UfBay4MCFakeSubH, colI+3);

  Int_t colI = kGray;
  //  dataPlot(UfBay1DataFakeSubH, colI);
  // dataPlot(UfBay2DataFakeSubH, colI+1);
  //  dataPlot(UfBay3DataFakeSubH, colI+2);
  //  dataPlot(UfBay4DataFakeSubH, colI+3);
  dataPlot(UfBay3DataFakeSubH, kBlack);

  TLegend *legend = 0;
  legend = new TLegend(0.3,0.7,0.89,0.89);
  legend->SetNColumns(1);
  legend->SetFillColor(0);
  legend->SetBorderSize(0.);
  legend->AddEntry(SimH, "MC truth", "f");
  //  legend->AddEntry(UfBay1DataFakeSubH, "Data 7TeV, fake sub., unfolded Bayes 0 iterations", "pl");
  //  legend->AddEntry(UfBay2DataFakeSubH, "Data 7TeV, fake sub., unfolded Bayes 1 iterations", "pl");
  legend->AddEntry(UfBay3DataFakeSubH, "Data 7TeV, fake sub., unfolded "+TSsDataFakeSub, "pl");
  //  legend->AddEntry(UfBay4DataFakeSubH, "Data 7TeV, fake sub., unfolded Bayes 5 iterations", "pl");
  //  legend->AddEntry(UfBay1MCFakeSubH, "MC reco, fake sub., unfolded Bayes 0 iterations", "pl");
  //  legend->AddEntry(UfBay2MCFakeSubH, "MC reco, fake sub., unfolded Bayes 1 iterations", "pl");
  legend->AddEntry(UfBay3MCFakeSubH, "MC reco, fake sub., unfolded "+TSsMCFakeSub, "fp");
  //  legend->AddEntry(UfBay4MCFakeSubH, "MC reco, fake sub., unfolded Bayes 5 iterations", "pl");

  legend->Draw();
  writePrel(frame2);

  /////////////////////

  //////////////////////////////////////////////////////////////////////////////
  //Frame2
  //  Double_t Min = 0.1;
  Min = 0.00001;
  Max = 0.07*DataH->GetBinContent(DataH->GetMaximumBin());
  TH2F *frame3= new TH2F("frame3",";R (cm);#propto N_{had} #times 1/#lambda_{I}/event",1.,DataH->GetXaxis()->GetBinLowEdge(1),DataH->GetXaxis()->GetBinUpEdge(DataH->GetNbinsX()),1.,Min,Max);
  //  TH2F *frame2= new TH2F("frame2",";R (cm);X_{0}^{-1} (arbitrary units)",1.,DataH->GetXaxis()->GetBinLowEdge(1),DataH->GetXaxis()->GetBinUpEdge(DataH->GetNbinsX()),1.,Min,Max);
  frame3->SetStats(kFALSE);

  //
  //Plot on the second canvas
  rCan3->Divide(1,2,0.,0.);
  rCan3->cd(1);
  rCan3->cd(1)->SetTopMargin(0.13);
  rCan3->cd(1)->SetBottomMargin(0.004);
  rCan3->cd(1)->SetRightMargin(0.01);
  rCan3->cd(1)->SetTicks(0,0);
  //  gPad->SetTickx(0);
  frame3->GetYaxis()->SetNdivisions(505);
  //  TGaxis::SetMaxDigits(2);
  frame3->SetLabelSize(0.1,"xy");
  frame3->SetTitleSize(0.1,"xy");
  frame3->SetTitleOffset(0.4,"y");
  frame3->Draw();

  Int_t colI = kAzure;
  //  mcPlot(UfBay1MCFakeSubH, colI);
  //  mcPlot(UfBay2MCFakeSubH, colI+1);
  mcPlot(MCFakeSubH, colI);
  //  mcPlot(UfBay4MCFakeSubH, colI+3);

  Int_t colI = kGray;
  //  dataPlot(UfBay1DataFakeSubH, colI);
  // dataPlot(UfBay2DataFakeSubH, colI+1);
  //  dataPlot(UfBay3DataFakeSubH, colI+2);
  //  dataPlot(UfBay4DataFakeSubH, colI+3);
  dataPlot(DataFakeSubH, kBlack);

  TLegend *legend = 0;
  Float_t legx1 = 0.45;
  Float_t legx2 = 0.95;
  Float_t legy2 = 0.85;
  Float_t yoff = 0.15;

  legend = new TLegend(legx1,legy2-2.*yoff,legx2,legy2);
  legend->SetNColumns(1);
  legend->SetFillColor(0);
  legend->SetBorderSize(0.);
  //  legend->AddEntry(SimH, "MC truth", "f");
  //  legend->AddEntry(UfBay1DataFakeSubH, "Data 7TeV, fake sub., unfolded Bayes 0 iterations", "pl");
  //  legend->AddEntry(UfBay2DataFakeSubH, "Data 7TeV, fake sub., unfolded Bayes 1 iterations", "pl");
  legend->AddEntry(DataFakeSubH, "Data 7TeV, fake sub., "+TSsDataFakeSub, "pl");
  //  legend->AddEntry(UfBay4DataFakeSubH, "Data 7TeV, fake sub., unfolded Bayes 5 iterations", "pl");
  //  legend->AddEntry(UfBay1MCFakeSubH, "MC reco, fake sub., unfolded Bayes 0 iterations", "pl");
  //  legend->AddEntry(UfBay2MCFakeSubH, "MC reco, fake sub., unfolded Bayes 1 iterations", "pl");
  legend->AddEntry(MCFakeSubH, "MC reco, fake sub., "+TSsMCFakeSub, "fp");
  //  legend->AddEntry(UfBay4MCFakeSubH, "MC reco, fake sub., unfolded Bayes 5 iterations", "pl");

  legend->Draw();
  writePrel(frame3);

  //Frame3
  //  Double_t Min = 0.1;

  //
  // Beam pipe x/X_0 = 0.0800cm/35.28cm
  // x/lI = 0.0800cm/42.10cm
  cout << SimH->GetBinContent(1) << " "<< SimH->GetBinContent(2) << " "<< SimH->GetBinContent(3) << " "<< SimH->GetBinContent(4) << " " << SimH->GetBinContent(5) << endl;
  Double_t scale = (0.08/42.10)/(1.*SimH->GetBinContent(4));

  cout << " scale " << scale << endl; 
  
  SimH->Scale(100.*scale);

  Min = 0.0;
  Max = 1.5*SimH->GetBinContent(SimH->GetMaximumBin());
  TH2F *frame4= new TH2F("frame4",";R (cm);x/#lambda_{I} [%]",1.,DataH->GetXaxis()->GetBinLowEdge(1),DataH->GetXaxis()->GetBinUpEdge(DataH->GetNbinsX()),1.,Min,Max);
  //  TH2F *frame2= new TH2F("frame2",";R (cm);X_{0}^{-1} (arbitrary units)",1.,DataH->GetXaxis()->GetBinLowEdge(1),DataH->GetXaxis()->GetBinUpEdge(DataH->GetNbinsX()),1.,Min,Max);
  frame4->SetStats(kFALSE);
  frame4->SetTitleOffset(1.1,"y");

  rCan3->cd(2);
  rCan3->cd(2)->SetBottomMargin(0.18);
  rCan3->cd(2)->SetTopMargin(0.004);
  rCan3->cd(2)->SetRightMargin(0.01);
  gPad->SetTickx(0);

  frame4->Draw();
  frame4->GetYaxis()->SetNdivisions(505);
  frame4->SetLabelSize(0.09,"xy");
  frame4->SetTitleSize(0.09,"xy");
  frame4->SetTitleOffset(0.5,"y");
  //  frame4->GetYaxis()->SetNdivisions(505);

  legend = new TLegend(legx1,legy2-0.*yoff-0.02,legx2,legy2+yoff-0.02);
  legend->SetNColumns(1);
  legend->SetFillColor(0);
  legend->SetBorderSize(0.);
  legend->AddEntry(SimH, "Material in simulation", "f");

  legend->Draw();
  writeLabels(frame4);


  SimH->SetLineColor(kGreen+2);
  SimH->SetFillColor(kSpring+1);
  SimH->SetLineWidth(1);
  SimH->SetStats(kFALSE);
  SimH->Draw("HISTsame");

  /////////////////////

  rCan1->SaveAs("rCan1.png");
  rCan2->SaveAs("rCan2.png");
  rCan3->SaveAs("rCan3.png");

}
