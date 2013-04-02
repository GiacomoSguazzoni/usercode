#include "elPlot.h"
#include <TVector.h>

elPlot::elPlot(const char* name, Double_t umin, Double_t umax, int uN, Double_t vmin, Double_t vmax, int vN)
{

  std::cout << " Istanzio 2D... " << name << std::endl;

  uNBin = uN;
  vNBin = vN;
  uMin = umin;
  uMax = umax;
  vMin = vmin;
  vMax = vmax;
  uBinW = (uMax - uMin)/uNBin;
  vBinW = (vMax - vMin)/vNBin;

  TSName.Append(name);

  int iHisto = 0;

  //
  for (int iv=0; iv<vNBin; iv++){
    
    double thisBinVMin = vMin + iv*vBinW;
    double vBinC = vMin + (iv+0.5)*vBinW;
    double thisBinVMax = vMin + (iv+1)*vBinW; 

    for (int iu=0; iu<uNBin; iu++){

      double thisBinUMin = uMin + iu*uBinW;
      double uBinC = uMin + (iu+0.5)*uBinW;
      double thisBinUMax = uMin + (iu+1)*uBinW; 

      char binHistoId[100];
      sprintf(binHistoId,"%s_%03d",name, iHisto);

      char binHistoInfo[100];
      sprintf(binHistoInfo,"u:% 4.2f,% 4.2f v:% 4.2f,% 4.2f",thisBinUMin,thisBinUMax,thisBinVMin,thisBinVMax);

      binHistos * thisBH = new binHistos(binHistoId, binHistoInfo, uBinC, vBinC);

      bHVector.push_back(thisBH);

      iHisto++;

      std::cout << " Created binHistos object " << binHistoId << " for range " << binHistoInfo << std::endl;

    }
  }
  //

  //Root file for writing
  theFile = new TFile("./plots/"+TSName+".root","RECREATE");

  //Store the number of bins in the root file
  TVector v(1);
  v[0] = iHisto;
  nBins = iHisto;
  v.Write("nbins");

  otherInits();

}

#include <fstream>

elPlot::elPlot(const char* name)
{

  TSName.Append(name);

  //Root file for reading
  theFile = new TFile("./plots/"+TSName+".root","READ");

  std::cout << " Reading file " << theFile->GetPath() << std::endl;

  TVector *nbinv = (TVector*)theFile->Get("nbins");


  nBins = (*nbinv)[0];

  //
  for (int i=0; i<nBins; i++){

      char binHistoId[100];
      sprintf(binHistoId,"%s_%03d",name, i);

      binHistos * thisBH = new binHistos(binHistoId, theFile);

      bHVector.push_back(thisBH);

      std::cout << " Created binHistos object " << binHistoId << " to be read from file " << std::endl;
  
  }

  otherInits();

}

#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
 
elPlot::~elPlot()
{

  TCanvas can("measVsSim","",1000,1000);
  can.SetRightMargin(0.05);
  can.SetLeftMargin(0.15);
  can.SetTopMargin(0.05);
  can.SetBottomMargin(0.15);

  // Draw frame
  TH2F *frame = new TH2F("frame",";(dP/dx)_{sim} MPV [MeV/cm];(dP/dx)_{meas} MPV [MeV/cm]",1,0.13,0.345,1,0.05,0.54); 
  frame->SetStats(0);
  frame->SetTitle("");
  //  frame->GetYaxis()->SetDecimals(1);
  frame->GetXaxis()->SetNdivisions(505); //Default 510
  frame->GetYaxis()->SetNdivisions(505); //Default 510
  //  frame->GetXaxis()->SetMoreLogLabels();
  //  frame->GetXaxis()->SetNoExponent();
  frame->GetYaxis()->SetTitleOffset(1.6);
  frame->GetXaxis()->SetTitleOffset(1.4);
  frame->SetLabelSize(0.04,"xy");
  frame->SetTitleSize(0.04,"xy");
  frame->Draw();

  /*
  dPdxMeasVsSimGraph.SetMarkerStyle(20);
  dPdxMeasVsSimGraph.SetMarkerSize(0.5);
  dPdxMeasVsSimGraph.SetMarkerColor(kBlue);
  dPdxMeasVsSimGraph.Draw("psame");

  dPdxGausVsSimGraph.SetMarkerStyle(20);
  dPdxGausVsSimGraph.SetMarkerSize(0.5);
  dPdxGausVsSimGraph.SetMarkerColor(kRed);
  dPdxGausVsSimGraph.Draw("psame");
*/


  TLegend * leg = new TLegend(0.3,0.7,0.55,0.95);

  //  meanTrunGraph->draw(leg);
  medianGraph->draw(leg);
  meanGraph->draw(leg);
  gausGraph->draw(leg);
  //  measGraph->draw(leg);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw();

  can.SaveAs("./plots/measVsSim.png");

  //
  // vs Eta and vs Phi
  //

  TCanvas canVsEta("vsEta","",1200,900);
  canVsEta.SetRightMargin(0.05);
  canVsEta.SetLeftMargin(0.15);
  canVsEta.SetTopMargin(0.05);
  canVsEta.SetBottomMargin(0.15);

  // Draw frame
  TH2F *frameVsEta = new TH2F("frameVsEta",";#eta;(dP/dx) [MeV/cm]",1,-2.1,2.1,1,0,0.7); 
  frameVsEta->SetStats(0);
  frameVsEta->SetTitle("");
  //  frameVsEta->GetYaxis()->SetDecimals(1);
  frameVsEta->GetXaxis()->SetNdivisions(505); //Default 510
  frameVsEta->GetYaxis()->SetNdivisions(505); //Default 510
  //  frameVsEta->GetXaxis()->SetMoreLogLabels();
  //  frameVsEta->GetXaxis()->SetNoExponent();
  frameVsEta->GetYaxis()->SetTitleOffset(1.6);
  frameVsEta->GetXaxis()->SetTitleOffset(1.4);
  frameVsEta->SetLabelSize(0.04,"xy");
  frameVsEta->SetTitleSize(0.04,"xy");
  frameVsEta->Draw();

  TLegend * legVsEta = new TLegend(0.2,0.7,0.45,0.95);

  gausVsEtaGraph->draw(legVsEta);
  meanVsEtaGraph->draw(legVsEta);
  medianVsEtaGraph->draw(legVsEta);
  simVsEtaGraph->draw(legVsEta);
  legVsEta->SetFillStyle(0);
  legVsEta->SetBorderSize(0);
  legVsEta->Draw();

  canVsEta.SaveAs("./plots/dpdxVsEta.png");

  //
  // Close file
  //

  std::cout << " Closing the file..." << std::endl;

  theFile->Close();

}

void elPlot::otherInits(){

  //Graphs
  measGraph = new dpdxGraph("Landau Fit",kBlue);
  gausGraph = new dpdxGraph("Gaus Fit",kAzure+3);
  meanGraph = new dpdxGraph("Mean",kSpring-7);
  medianGraph = new dpdxGraph("Median",kOrange-3);
  meanTrunGraph = new dpdxGraph("Trunc Mean",kRed);

  gausVsEtaGraph = new dpdxGraph("Gaus Fit",kAzure+3);
  meanVsEtaGraph = new dpdxGraph("Mean",kSpring-7);
  medianVsEtaGraph = new dpdxGraph("Median",kOrange-3);
  simVsEtaGraph = new dpdxGraph("Sim MPV",kBlue);

  gausVsPhiGraph = new dpdxGraph("Gaus Fit",kAzure+3);
  simVsPhiGraph = new dpdxGraph("Sim MPV",kBlue);

}

void elPlot::Fill(double valSim, double valMeas, double valMeasErr, double u, double v, double wei)
{

  int theBin = FindBinHisto(u, v);

  //  std::cout << " u: " << u << " v: " << v << " The bin is number #:" << theBin << std::endl;

  if ( theBin > -1 ) (bHVector.at(theBin))->Fill(wei, valSim, valMeas, valMeasErr);
  
}

int elPlot::FindBinHisto(double u, double v)
{

  if ( u < uMin ) return -1;
  if ( u > uMax ) return -1;
  if ( v < vMin ) return -1;
  if ( v > vMax ) return -1;
  
  int nubin = (int)floor((u - uMin)/uBinW);
  int nvbin = (int)floor((v - vMin)/vBinW);
  
  return nvbin*uNBin+nubin;

}

void elPlot::PlotAndWriteAll()
{

  for (std::vector<binHistos*>::iterator it=bHVector.begin() ; it < bHVector.end(); it++ ){
    std::cout << " Plotting histos of bin " << (*it)->getName() << std::endl;
    (*it)->PlotAll();
    (*it)->WriteAll(theFile);
  }

}

void elPlot::FitAll()
{

  //
  // Open file for writing
  ofstream myfile;
  myfile.open ("plots/ploss.dat");
  
  for (std::vector<binHistos*>::iterator it=bHVector.begin() ; it < bHVector.end(); it++ ){
    std::cout << " Fitting histos of bin " << (*it)->getName() << std::endl;

    double mpvSim, mpvMeas, mpvGaus, mpvMean, mpvMedian, mpvTruncMean;
    (*it)->FitAll(mpvSim, mpvMeas, mpvGaus, mpvMean, mpvMedian, mpvTruncMean);

    double eta = (*it)->GetBinUCenter();
    double phi = (*it)->GetBinVCenter();

    //
    // Media
    // Gaussiana (prima fit, poi refit tra +-1sigma)
    // Mediana
    // Media con troncamento 0.25 0.75 
    //

    std::cout << "############################ mpvSim: " << mpvSim << " mpvMeas: " << mpvMeas << std::endl; 
    //    dPdxMeasInBinVect.push_back(mpvMeas);
    //    dPdxSimInBinVect.push_back(mpvSim);
    //    dPdxGausInBinVect.push_back(mpvGaus);
    measGraph->addPoint(mpvSim,mpvMeas);
    gausGraph->addPoint(mpvSim,mpvGaus);
    meanGraph->addPoint(mpvSim,mpvMean);
    medianGraph->addPoint(mpvSim,mpvMedian);
    meanTrunGraph->addPoint(mpvSim,mpvTruncMean);
    
    gausVsEtaGraph->addPoint(eta,mpvGaus);
    float shift = 0.03;
    double etashift = eta - shift;
    meanVsEtaGraph->addPoint(etashift,mpvMean);
    etashift = eta + shift;
    medianVsEtaGraph->addPoint(etashift,mpvMedian);
    simVsEtaGraph->addPoint(eta,mpvSim);
    
    gausVsPhiGraph->addPoint(phi,mpvGaus);
    simVsPhiGraph->addPoint(phi,mpvSim);

    // Write values here
    myfile << eta << " " << phi << " " << mpvSim << " " << mpvGaus << " " << mpvMean << " " << mpvMedian << std::endl;


  }

  myfile.close();

}

