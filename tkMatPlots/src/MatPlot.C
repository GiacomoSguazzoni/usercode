#define MatPlot_cxx
#include "MatPlot.h"

void MatPlot::FillData(convR2SforMatPlot* dataR2S, Int_t minE, Int_t maxE)
{
  
#include "MatPlotFillData.cxx"

}

void MatPlot::FillData(niR2SforMatPlot* dataR2S, Int_t minE, Int_t maxE)
{
  
#include "MatPlotFillData.cxx"

}

void MatPlot::FillMC(convR2SforMatPlot* mcRecoR2S, Int_t minE, Int_t maxE)
{

#include "MatPlotFillMC.cxx"

}

void MatPlot::FillMC(niR2SforMatPlot* mcRecoR2S, Int_t minE, Int_t maxE)
{

#include "MatPlotFillMC.cxx"

}

void MatPlot::FillSim(convS2RforMatPlot* mcSimS2R, Int_t minE, Int_t maxE)
{

#include "MatPlotFillSim.cxx"
  
}

void MatPlot::FillSim(niS2RforMatPlot* mcSimS2R, Int_t minE, Int_t maxE)
{

#include "MatPlotFillSim.cxx"
  
}

void MatPlot::Fill(convR2SforMatPlot* r2s, TH2* hist, TH2* fake)
{
  
  r2s->SetGeoCuts(geoCuts);
  r2s->SetEffRadius(effRadius);
  r2s->LoopForFill(hist, fake);

}

void MatPlot::Fill(convR2SforMatPlot* r2s, TH1* hist, TH1* fake)
{

  r2s->SetGeoCuts(geoCuts);
  r2s->SetEffRadius(effRadius);
  r2s->LoopForFill(hist, fake);

}

void MatPlot::Fill(convS2RforMatPlot* s2r, TH2* hist)
{

  s2r->SetGeoCuts(geoCuts);
  s2r->SetEffRadius(effRadius);
  s2r->LoopForFill(hist, 0);

}

void MatPlot::Fill(convS2RforMatPlot* s2r, TH1* hist)
{


  s2r->SetGeoCuts(geoCuts);
  s2r->SetEffRadius(effRadius);
  s2r->LoopForFill(hist, 0);

}

//

void MatPlot::Fill(niR2SforMatPlot* r2s, TH2* hist, TH2* fake)
{

  r2s->SetGeoCuts(geoCuts);
  r2s->SetEffRadius(effRadius);
  r2s->LoopForFill(hist, fake);

}

void MatPlot::Fill(niR2SforMatPlot* r2s, TH1* hist, TH1* fake)
{

  r2s->SetGeoCuts(geoCuts);
  r2s->SetEffRadius(effRadius);
  r2s->LoopForFill(hist, fake);

}

void MatPlot::Fill(niS2RforMatPlot* s2r, TH2* hist)
{

  s2r->SetGeoCuts(geoCuts);
  s2r->SetEffRadius(effRadius);
  s2r->LoopForFill(hist, 0);

}

void MatPlot::Fill(niS2RforMatPlot* s2r, TH1* hist)
{

  s2r->SetGeoCuts(geoCuts);
  s2r->SetEffRadius(effRadius);
  s2r->LoopForFill(hist, 0);

}

void MatPlot::PlotAll()
{

  if( rawSim2DH ) {
    Plot(rawSim2DH);
    Plot(rawMC2DH);
    Plot(rawData2DH);
    Plot(rawFake2DH);
    //
    WBScale2D(Sim2DH);
    Plot(Sim2DH);
    WBScale2D(MC2DH);
    EffScale2D(MC2DH);
    Plot(MC2DH);
    WBScale2D(Data2DH);
    EffScale2D(Data2DH);
    Plot(Data2DH);
    WBScale2D(Fake2DH);
    EffScale2D(Fake2DH);
    Plot(Fake2DH);
    std::cout << ">>>>Data entries " << rawData2DH->GetEntries() << " integral " << rawData2DH->Integral() << std::endl;
    std::cout << ">>>>MC entries " << rawMC2DH->GetEntries() << " integral " << rawMC2DH->Integral() << std::endl;
    std::cout << ">>>>Fake entries " << rawFake2DH->GetEntries() << " integral " << rawFake2DH->Integral() << std::endl;
  }

  if( rawSim1DH ) {
    Plot(rawSim1DH);
    Plot(rawMC1DH);
    Plot(rawMCFs1DH);
    Plot(rawData1DH);
    Plot(rawDataFs1DH);
    Plot(rawFake1DH);
    //
    WBScale1D(Sim1DH);
    Plot(Sim1DH);
    WBScale1D(MC1DH);
    EffScale1D(MC1DH);
    Plot(MC1DH);
    WBScale1D(MCFs1DH);
    EffScale1D(MCFs1DH);
    Plot(MCFs1DH);
    WBScale1D(Data1DH);
    EffScale1D(Data1DH);
    Plot(Data1DH);
    WBScale1D(DataFs1DH);
    EffScale1D(DataFs1DH);
    Plot(DataFs1DH);
    WBScale1D(Fake1DH);
    EffScale1D(Fake1DH);
    Plot(Fake1DH);
    //
#ifdef UNFOLD
    if ( unfoldBay1MCFs1DH ) {
      WBScale1D(unfoldBay1MCFs1DH);
      Plot(unfoldBay1MCFs1DH);
      WBScale1D(unfoldBay2MCFs1DH);
      Plot(unfoldBay2MCFs1DH);
      WBScale1D(unfoldBay3MCFs1DH);
      Plot(unfoldBay3MCFs1DH);
      WBScale1D(unfoldBay4MCFs1DH);
      Plot(unfoldBay4MCFs1DH);
      //
      WBScale1D(unfoldBay1DataFs1DH);
      Plot(unfoldBay1DataFs1DH);
      WBScale1D(unfoldBay2DataFs1DH);
      Plot(unfoldBay2DataFs1DH);
      WBScale1D(unfoldBay3DataFs1DH);
      Plot(unfoldBay3DataFs1DH);
      WBScale1D(unfoldBay4DataFs1DH);
      Plot(unfoldBay4DataFs1DH);
    }
#endif

  }    
}

void MatPlot::WBScale1D(TH1* rWei1DH)
{

  if ( uIndex == 4 && vIndex==0 ) WBScaleR(rWei1DH);
  if ( uIndex == 6 && vIndex==0 ) WBScalePhi(rWei1DH);

}

void MatPlot::EffScale1D(TH1* rWei1DH)
{

  if ( uIndex == 4 && vIndex==0 ) EffScaleR(rWei1DH);

}

void MatPlot::WBScale2D(TH2* rWei2DH)
{

  if ( uIndex == 1 && vIndex==2 ) WBScaleXY(rWei2DH);
  if ( uIndex == 3 && vIndex==4 ) WBScaleRZ(rWei2DH);

}

void MatPlot::EffScale2D(TH2* rWei2DH)
{

  if ( uIndex == 1 && vIndex==2 ) EffScaleXY(rWei2DH);
  if ( uIndex == 3 && vIndex==4 ) EffScaleRZ(rWei2DH);

}

void MatPlot::EffScaleXY(TH2* rWei2DH)
{

  for (Int_t iBin = 1; iBin<rWei2DH->GetNbinsX()+1; iBin++){
    Double_t x = rWei2DH->GetXaxis()->GetBinCenter(iBin); 
    for (Int_t jBin = 1; jBin<rWei2DH->GetNbinsY()+1; jBin++){
      Double_t y = rWei2DH->GetYaxis()->GetBinCenter(jBin); 

      Double_t rAve = sqrt(x*x+y*y);
      Double_t eff = effRadius->EffForRadius(rAve);

      Double_t wei = 1.;
      if ( eff ) wei = 1./eff;
      
      rWei2DH->SetBinError(iBin,jBin,rWei2DH->GetBinError(iBin,jBin)*wei);
      rWei2DH->SetBinContent(iBin,jBin,rWei2DH->GetBinContent(iBin,jBin)*wei);

    }
  }
}

void MatPlot::WBScaleXY(TH2* rWei2DH)
{

  //GeoCuts iterator
  std::vector<GeoCut>::iterator geoIt;
  
  std::cout << " Computing rewei histos XY" << std::endl;

  for (Int_t iBin = 1; iBin<rWei2DH->GetNbinsX()+1; iBin++){
    Double_t x1 = rWei2DH->GetXaxis()->GetBinLowEdge(iBin); 
    Double_t x2 = rWei2DH->GetXaxis()->GetBinUpEdge(iBin); 
    for (Int_t jBin = 1; jBin<rWei2DH->GetNbinsY()+1; jBin++){
      Double_t y1 = rWei2DH->GetYaxis()->GetBinLowEdge(jBin); 
      Double_t y2 = rWei2DH->GetYaxis()->GetBinUpEdge(jBin);
      Double_t zMin = +10000.; 
      Double_t zMax = -10000.;

      Double_t rMax = FindRMax(x1, x2, y1, y2);

      Int_t iValidRange = 0;
     
      for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ )
	{
	  if ( (*geoIt).GetRMin()<rMax && (*geoIt).GetRMax()>rMax )
	    { 
	      iValidRange = 1;
	      if ( (*geoIt).GetZMin() < zMin ) {zMin=(*geoIt).GetZMin();}; 
	      if ( (*geoIt).GetZMax() > zMax ) {zMax=(*geoIt).GetZMax();}; 
	    }
	}

      //
      // FIXME!!! Add Dz segments instead of finding extrema!
      //

      
      Double_t wei = 1.;
      if ( iValidRange ) wei = WBWeiCompDxDyDz(x1, x2, y1, y2, zMin, zMax);
      
      rWei2DH->SetBinError(iBin,jBin,rWei2DH->GetBinError(iBin,jBin)/wei);
      rWei2DH->SetBinContent(iBin,jBin,rWei2DH->GetBinContent(iBin,jBin)/wei);
    }
  }

}

void MatPlot::EffScaleRZ(TH2* rWei2DH)
{

  for (Int_t iBin = 1; iBin<rWei2DH->GetNbinsX()+1; iBin++){
    //    Double_t z = rWei2DH->GetXaxis()->GetBinCenter(iBin); 
    for (Int_t jBin = 1; jBin<rWei2DH->GetNbinsY()+1; jBin++){
      Double_t r = rWei2DH->GetYaxis()->GetBinCenter(jBin); 

      Double_t eff = effRadius->EffForRadius(r);

      Double_t wei = 1.;
      if ( eff ) wei = 1./eff;
      
      rWei2DH->SetBinError(iBin,jBin,rWei2DH->GetBinError(iBin,jBin)*wei);
      rWei2DH->SetBinContent(iBin,jBin,rWei2DH->GetBinContent(iBin,jBin)*wei);

    }
  }
}

void MatPlot::WBScaleRZ(TH2* rWei2DH)
{

  //GeoCuts iterator
  std::vector<GeoCut>::iterator geoIt;
  
  std::cout << " Computing rewei histos RZ" << std::endl;

  for (Int_t iBin = 1; iBin<rWei2DH->GetNbinsX()+1; iBin++){
    Double_t z1 = rWei2DH->GetXaxis()->GetBinLowEdge(iBin); 
    Double_t z2 = rWei2DH->GetXaxis()->GetBinUpEdge(iBin); 
    for (Int_t jBin = 1; jBin<rWei2DH->GetNbinsY()+1; jBin++){
      Double_t r1 = rWei2DH->GetYaxis()->GetBinLowEdge(jBin); 
      Double_t r2 = rWei2DH->GetYaxis()->GetBinUpEdge(jBin);

      Int_t iValidRange = 0;
     
      for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ )
	{
	  if ( (*geoIt).GetRMin()<r1 && (*geoIt).GetRMax()>r2 && (*geoIt).GetZMin() < z1 && (*geoIt).GetZMax() > z2 ) iValidRange = 1;
	}
      
      Double_t wei = 1.;
      if ( iValidRange ) wei = WBWeiCompDrDzDphi(r1, r2, z1, z2);

      rWei2DH->SetBinError(iBin,jBin,rWei2DH->GetBinError(iBin,jBin)/wei);
      rWei2DH->SetBinContent(iBin,jBin,rWei2DH->GetBinContent(iBin,jBin)/wei);
    }
  }

}

void MatPlot::EffScaleR(TH1* rWei1DH)
{

  TH1D * effH = new TH1D("effH","effH",rWei1DH->GetNbinsX(),rWei1DH->GetBinLowEdge(1),rWei1DH->GetXaxis()->GetBinUpEdge(rWei1DH->GetNbinsX()));
  
  for (Int_t iBin = 1; iBin<rWei1DH->GetNbinsX()+1; iBin++){
    Double_t r = rWei1DH->GetXaxis()->GetBinCenter(iBin); 
    Double_t eff = effRadius->EffForRadius(r);
    effH->SetBinContent(iBin,eff);
  }

  effH->Smooth();

  for (Int_t iBin = 1; iBin<rWei1DH->GetNbinsX()+1; iBin++){

      Double_t wei = 1.;
      Double_t eff = effH->GetBinContent(iBin);
      if ( eff ) wei = 1./eff;
      
      rWei1DH->SetBinError(iBin,rWei1DH->GetBinError(iBin)*wei);
      rWei1DH->SetBinContent(iBin,rWei1DH->GetBinContent(iBin)*wei);

  }

  delete effH;

}

void MatPlot::WBScaleR(TH1* rWei1DH)
{
  
  //GeoCuts iterator
  std::vector<GeoCut>::iterator geoIt;
  
  std::cout << " Computing rewei histos R" << std::endl;
  
  for (Int_t iBin = 1; iBin<rWei1DH->GetNbinsX()+1; iBin++){
    Double_t r1 = rWei1DH->GetXaxis()->GetBinLowEdge(iBin); 
    Double_t r2 = rWei1DH->GetXaxis()->GetBinUpEdge(iBin); 
    
    Double_t zMin = +10000.; 
    Double_t zMax = -10000.;
    
    Int_t iValidRange = 0;
    
    for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ )
      {
	if ( (*geoIt).GetRMin()<=r1 && (*geoIt).GetRMax()>=r2 )
	  { 
	    iValidRange = 1;
	    if ( (*geoIt).GetZMin() < zMin ) {zMin=(*geoIt).GetZMin();}; 
	    if ( (*geoIt).GetZMax() > zMax ) {zMax=(*geoIt).GetZMax();}; 
	  }
      }
    
    Double_t wei = 1.;
    if ( iValidRange ) wei = WBWeiCompDrDzDphi(r1, r2, zMin, zMax);
    
    rWei1DH->SetBinError(iBin,rWei1DH->GetBinError(iBin)/wei);
    rWei1DH->SetBinContent(iBin,rWei1DH->GetBinContent(iBin)/wei);
  }
}

void MatPlot::WBScalePhi(TH1* rWei1DH)
{
  
  //GeoCuts iterator
  std::vector<GeoCut>::iterator geoIt;
  
  std::cout << " Computing rewei histos Phi" << std::endl;
  
  for (Int_t iBin = 1; iBin<rWei1DH->GetNbinsX()+1; iBin++){
    Double_t phi1 = rWei1DH->GetXaxis()->GetBinLowEdge(iBin); 
    Double_t phi2 = rWei1DH->GetXaxis()->GetBinUpEdge(iBin); 

    //
    // Simple implementation only with a single geocut...
    //    

    Double_t rMin = 0.; 
    Double_t rMax = 0.;
    Double_t zMin = 0.; 
    Double_t zMax = 0.;

    Int_t iValidRange = 0;
    for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ )
      {
	iValidRange = 1;
	rMin = (*geoIt).GetRMin(); 
	rMax = (*geoIt).GetRMax();
	zMin = (*geoIt).GetZMin(); 
	zMax = (*geoIt).GetZMax();
      }
    
    Double_t wei = 1.;
    if ( iValidRange ) wei = WBWeiCompDrDzDphi(rMin, rMax, zMin, zMax, phi1, phi2);
    
    rWei1DH->SetBinError(iBin,rWei1DH->GetBinError(iBin)/wei);
    rWei1DH->SetBinContent(iBin,rWei1DH->GetBinContent(iBin)/wei);
  }
}

Double_t MatPlot::WBWeiCompDxDyDz(Double_t x1, Double_t x2, Double_t y1, Double_t y2, Double_t z1, Double_t z2){
  
  Double_t x = 0.5*(x1+x2);
  Double_t y = 0.5*(y1+y2);
  Double_t rAve = sqrt(x*x+y*y); 

  /*
    Double_t part1 = log(sqrt(rAve*rAve+z2*z2)+z2)-log(sqrt(rAve*rAve+z1*z1)+z1);
    Double_t wei = (x2-x1)*(y2-y1)*part1/rAve;
  */

  Double_t wei = (x2-x1)*(y2-y1)*(z2-z1)/rAve/rAve;
  
  return wei;
  
}

Double_t MatPlot::WBWeiCompDrDzDphi(Double_t r1, Double_t r2, Double_t z1, Double_t z2, Double_t phi1, Double_t phi2){
  
  //  Double_t PI=3.141592653589793;
  
  /*
    Double_t part1 = r2*log(sqrt(r2*r2+z2*z2)+z2)-r2*log(sqrt(r2*r2+z1*z1)+z1);
    Double_t part2 = r1*log(sqrt(r1*r1+z2*z2)+z2)-r1*log(sqrt(r1*r1+z1*z1)+z1);
    Double_t part3 = z2*log(sqrt(r2*r2+z2*z2)+r2)-z2*log(sqrt(r1*r1+z2*z2)+r1);
    Double_t part4 = z1*log(sqrt(r2*r2+z1*z1)+r2)-z1*log(sqrt(r1*r1+z1*z1)+r1);
    Double_t wei = 2*PI*(part1-part2+part3-part4);
  */
  
  Double_t wei = (phi2-phi1)*(z2-z1)*log(r2/r1);
  
  return wei;
    
}

#include <TCanvas.h>

void MatPlot::Plot(TH2* hist)
{

  std::cout << " E ora plotto " << hist->GetName() << std::endl;

  TCanvas * can = new TCanvas(hist->GetName(),"",1000,1000);

  hist->Draw();

  can->SaveAs(TString(hist->GetName())+".png");
  hist->SaveAs(TString(hist->GetName())+".root");
  
}

void MatPlot::Plot(TH1* hist)
{

  std::cout << " E ora plotto " << hist->GetName() << std::endl;

  TCanvas * can = new TCanvas(hist->GetName(),"",1000,1000);

  hist->Draw();

  can->SaveAs(TString(hist->GetName())+".png");
  hist->SaveAs(TString(hist->GetName())+".root");
  
}

#ifdef UNFOLD
void MatPlot::responseTrain(convS2RforMatPlot* s2r, convR2SforMatPlot* r2s, Int_t minE, Int_t maxE)
{

#include "MatPlotResponseTrain.cxx"

}

void MatPlot::responseTrain(niS2RforMatPlot* s2r, niR2SforMatPlot* r2s, Int_t minE, Int_t maxE)
{

#include "MatPlotResponseTrain.cxx"

}

void MatPlot::doUnfold()
{
  
  unfoldBay1DataFs1DH = doUnfoldHist("rUfBay1DataFs", 0, rawDataFs1DH);
  unfoldBay2DataFs1DH = doUnfoldHist("rUfBay2DataFs", 1, rawDataFs1DH);
  unfoldBay3DataFs1DH = doUnfoldHist("rUfBay3DataFs", 3, rawDataFs1DH);
  unfoldBay4DataFs1DH = doUnfoldHist("rUfBay4DataFs", 5, rawDataFs1DH);
  
  std::cout << " rUfBay1DataFs " << unfoldBay1DataFs1DH->GetName() << " " << unfoldBay1DataFs1DH->GetEntries() <<   std::endl;
  std::cout << " rUfBay2DataFs " << unfoldBay2DataFs1DH->GetName() << " " << unfoldBay2DataFs1DH->GetEntries() <<   std::endl;
  std::cout << " rUfBay3DataFs " << unfoldBay3DataFs1DH->GetName() << " " << unfoldBay3DataFs1DH->GetEntries() <<   std::endl;
  std::cout << " rUfBay4DataFs " << unfoldBay4DataFs1DH->GetName() << " " << unfoldBay4DataFs1DH->GetEntries() <<   std::endl;

  unfoldBay1MCFs1DH = doUnfoldHist("rUfBay1MCFs", 0, rawMCFs1DH);
  unfoldBay2MCFs1DH = doUnfoldHist("rUfBay2MCFs", 1, rawMCFs1DH);
  unfoldBay3MCFs1DH = doUnfoldHist("rUfBay3MCFs", 3, rawMCFs1DH);
  unfoldBay4MCFs1DH = doUnfoldHist("rUfBay4MCFs", 5, rawMCFs1DH);

  std::cout << " rUfBay1MCFs " << unfoldBay1MCFs1DH->GetName() << " " << unfoldBay1MCFs1DH->GetEntries() <<   std::endl;
  std::cout << " rUfBay2MCFs " << unfoldBay2MCFs1DH->GetName() << " " << unfoldBay2MCFs1DH->GetEntries() <<   std::endl;
  std::cout << " rUfBay3MCFs " << unfoldBay3MCFs1DH->GetName() << " " << unfoldBay3MCFs1DH->GetEntries() <<   std::endl;
  std::cout << " rUfBay4MCFs " << unfoldBay4MCFs1DH->GetName() << " " << unfoldBay4MCFs1DH->GetEntries() <<   std::endl;

}

TH1D * MatPlot::doUnfoldHist(const char * name, Int_t nIter, TH1D * in)
{

  unfoldBay = new RooUnfoldBayes(response, in, nIter);
  TH1D* out = (TH1D*) unfoldBay->Hreco();
  out->SetName(name);

  return out;

}

TH2D * MatPlot::doUnfoldHist(const char * name, Int_t nIter, TH2D * in)
{

  unfoldBay = new RooUnfoldBayes(response, in, nIter);
  TH2D* out = (TH2D*) unfoldBay->Hreco();
  out->SetName(name);

  return out;

}
#endif

void MatPlot::test(Double_t val)
{

  /*
  effRadius->Succeed(val, val+1.);
  effRadius->Fail(val);
  std::cout << " The r bin for r=" << val << " is " << effRadius->radiusBin(val) << " and eff " << effRadius->EffForRadius(val) << std::endl;
  */
  effRadius->DumpEfficiency();
  effRadius->DumpRawNumbers();

}
