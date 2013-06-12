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

void MatPlot::Fill(convR2SforMatPlot* r2s, TH2* hist, TH2* fake, TH2* histTrue)
{
  
  r2s->SetGeoCuts(geoCuts, effUV);
  r2s->LoopForFill(hist, fake, histTrue);

}

void MatPlot::Fill(convR2SforMatPlot* r2s, TH1* hist, TH1* fake, TH1* histTrue)
{

  r2s->SetGeoCuts(geoCuts, effUV);
  r2s->LoopForFill(hist, fake, histTrue);

}

void MatPlot::Fill(convS2RforMatPlot* s2r, TH2* hist)
{

  s2r->SetGeoCuts(geoCuts, effUV);
  s2r->LoopForFill(hist);

}

void MatPlot::Fill(convS2RforMatPlot* s2r, TH1* hist)
{

  s2r->SetGeoCuts(geoCuts, effUV);
  s2r->LoopForFill(hist);

}

//

void MatPlot::Fill(niR2SforMatPlot* r2s, TH2* hist, TH2* fake, TH2* histTrue)
{

  r2s->SetGeoCuts(geoCuts, effUV);
  r2s->LoopForFill(hist, fake, histTrue);

}

void MatPlot::Fill(niR2SforMatPlot* r2s, TH1* hist, TH1* fake, TH1* histTrue)
{

  r2s->SetGeoCuts(geoCuts, effUV);
  r2s->LoopForFill(hist, fake, histTrue);

}

void MatPlot::Fill(niS2RforMatPlot* s2r, TH2* hist)
{

  s2r->SetGeoCuts(geoCuts, effUV);
  s2r->LoopForFill(hist);

}

void MatPlot::Fill(niS2RforMatPlot* s2r, TH1* hist)
{

  s2r->SetGeoCuts(geoCuts, effUV);
  s2r->LoopForFill(hist);

}

void MatPlot::PlotAll()
{

  if( rawSim2DH ) {
    Plot(rawSim2DH);
    Plot(rawMC2DH);
    Plot(rawMCFakeSub2DH);
    Plot(rawMCFakeSubWrtSimPosition2DH);
    Plot(rawData2DH);
    Plot(rawMCFake2DH);
    //
    WBScale2D(Sim2DH);
    Plot(Sim2DH);
    //
    WBScale2D(MC2DH);
    EffScale2D(MC2DH);
    Plot(MC2DH);
    //
    WBScale2D(Data2DH);
    EffScale2D(Data2DH);
    Plot(Data2DH);
    //
    WBScale2D(MCFake2DH);
    EffScale2D(MCFake2DH);
    Plot(MCFake2DH);
    std::cout << ">>>>Data entries " << rawData2DH->GetEntries() << " integral " << rawData2DH->Integral() << std::endl;
    std::cout << ">>>>MC entries " << rawMC2DH->GetEntries() << " integral " << rawMC2DH->Integral() << std::endl;
    std::cout << ">>>>Fake entries " << rawMCFake2DH->GetEntries() << " integral " << rawMCFake2DH->Integral() << std::endl;
  }

  if( rawSim1DH ) {
    Plot(rawSim1DH);
    Plot(rawMC1DH);
    Plot(rawMCFakeSub1DH);
    Plot(rawMCFakeSubWrtSimPosition1DH);
    Plot(rawData1DH);
    Plot(rawDataFakeSub1DH);
    Plot(rawMCFake1DH);
    //
    WBScale1D(Sim1DH);
    Plot(Sim1DH);
    WBScale1D(MC1DH);
    EffScale1D(MC1DH);
    Plot(MC1DH);
    WBScale1D(MCFakeSub1DH);
    EffScale1D(MCFakeSub1DH);
    Plot(MCFakeSub1DH);
    WBScale1D(Data1DH);
    EffScale1D(Data1DH);
    Plot(Data1DH);
    WBScale1D(DataFakeSub1DH);
    EffScale1D(DataFakeSub1DH);
    Plot(DataFakeSub1DH);
    WBScale1D(MCFake1DH);
    EffScale1D(MCFake1DH);
    Plot(MCFake1DH);
    //
#ifdef UNFOLD
    if ( unfoldBay1MCFakeSub1DH ) {
      WBScale1D(unfoldBay1MCFakeSub1DH);
      Plot(unfoldBay1MCFakeSub1DH);
      WBScale1D(unfoldBay2MCFakeSub1DH);
      Plot(unfoldBay2MCFakeSub1DH);
      WBScale1D(unfoldBay3MCFakeSub1DH);
      Plot(unfoldBay3MCFakeSub1DH);
      WBScale1D(unfoldBay4MCFakeSub1DH);
      Plot(unfoldBay4MCFakeSub1DH);
      //
      WBScale1D(unfoldBay1DataFakeSub1DH);
      Plot(unfoldBay1DataFakeSub1DH);
      WBScale1D(unfoldBay2DataFakeSub1DH);
      Plot(unfoldBay2DataFakeSub1DH);
      WBScale1D(unfoldBay3DataFakeSub1DH);
      Plot(unfoldBay3DataFakeSub1DH);
      WBScale1D(unfoldBay4DataFakeSub1DH);
      Plot(unfoldBay4DataFakeSub1DH);
    }
#endif

  }    
}

void MatPlot::WBScale1D(TH1* rWei1DH)
{

  if ( uIndex == 4 && vIndex==0 ) WBScaleR(rWei1DH);
  if ( uIndex == 6 && vIndex==0 ) WBScalePhi(rWei1DH);
  if ( uIndex == 5 && vIndex==0 ) WBScaleTheta(rWei1DH);

}

void MatPlot::EffScale1D(TH1* rWei1DH)
{

  EffScaleUV(rWei1DH);

}

void MatPlot::WBScale2D(TH2* rWei2DH)
{

  if ( uIndex == 1 && vIndex==2 ) WBScaleXY(rWei2DH);
  if ( uIndex == 3 && vIndex==4 ) WBScaleRZ(rWei2DH);

}

void MatPlot::EffScale2D(TH2* rWei2DH)
{

  EffScaleUV(rWei2DH);

}

void MatPlot::EffScaleUV(TH2* rWei2DH)
{

  for (Int_t iBin = 1; iBin<rWei2DH->GetNbinsX()+1; iBin++){
    Double_t x = rWei2DH->GetXaxis()->GetBinCenter(iBin); 
    for (Int_t jBin = 1; jBin<rWei2DH->GetNbinsY()+1; jBin++){
      Double_t y = rWei2DH->GetYaxis()->GetBinCenter(jBin); 

      Double_t eff = effUV->EffForUV(x,y);

      Double_t wei = 1.;
      if ( eff ) wei = 1./eff;
      
      rWei2DH->SetBinError(iBin,jBin,rWei2DH->GetBinError(iBin,jBin)*wei);
      rWei2DH->SetBinContent(iBin,jBin,rWei2DH->GetBinContent(iBin,jBin)*wei);

    }
  }
}

void MatPlot::EffScaleUV(TH1* rWei1DH)
{
  
  for (Int_t iBin = 1; iBin<rWei1DH->GetNbinsX()+1; iBin++){
    Double_t x = rWei1DH->GetXaxis()->GetBinCenter(iBin); 
    
    Double_t eff = effUV->EffForUV(x,0.5);
    
    Double_t wei = 1.;
    if ( eff ) wei = 1./eff;
    
    rWei1DH->SetBinError(iBin,rWei1DH->GetBinError(iBin)*wei);
    rWei1DH->SetBinContent(iBin,rWei1DH->GetBinContent(iBin)*wei);
    
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
	  if ( (*geoIt).GetVMin()<rMax && (*geoIt).GetVMax()>rMax )
	    { 
	      iValidRange = 1;
	      if ( (*geoIt).GetUMin() < zMin ) {zMin=(*geoIt).GetUMin();}; 
	      if ( (*geoIt).GetUMax() > zMax ) {zMax=(*geoIt).GetUMax();}; 
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
	  if ( (*geoIt).GetVMin()<r1 && (*geoIt).GetVMax()>r2 && (*geoIt).GetUMin() < z1 && (*geoIt).GetUMax() > z2 ) iValidRange = 1;
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

  //

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
	if ( (*geoIt).GetVMin()<=r1 && (*geoIt).GetVMax()>=r2 )
	  { 
	    iValidRange = 1;
	    if ( (*geoIt).GetUMin() < zMin ) {zMin=(*geoIt).GetUMin();}; 
	    if ( (*geoIt).GetUMax() > zMax ) {zMax=(*geoIt).GetUMax();}; 
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
	rMin = (*geoIt).GetVMin(); 
	rMax = (*geoIt).GetVMax();
	zMin = (*geoIt).GetUMin(); 
	zMax = (*geoIt).GetUMax();
      }
    
    Double_t wei = 1.;
    if ( iValidRange ) wei = WBWeiCompDrDzDphi(rMin, rMax, zMin, zMax, phi1, phi2);
    
    rWei1DH->SetBinError(iBin,rWei1DH->GetBinError(iBin)/wei);
    rWei1DH->SetBinContent(iBin,rWei1DH->GetBinContent(iBin)/wei);
  }
}

void MatPlot::WBScaleTheta(TH1* rWei1DH)
{
  
  //GeoCuts iterator
  std::vector<GeoCut>::iterator geoIt;
  
  std::cout << " Computing rewei histos Phi" << std::endl;
  
  for (Int_t iBin = 1; iBin<rWei1DH->GetNbinsX()+1; iBin++){
    Double_t t1 = rWei1DH->GetXaxis()->GetBinLowEdge(iBin); 
    Double_t t2 = rWei1DH->GetXaxis()->GetBinUpEdge(iBin); 

    //
    // Simple implementation only with a single geocut...
    //    

    Double_t rMin = 0.; 
    Double_t rMax = 0.;
    Double_t phiMin = 0.; 
    Double_t phiMax = 0.;

    Int_t iValidRange = 0;
    for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ )
      {
	iValidRange = 1;
	rMin = (*geoIt).GetVMin(); 
	rMax = (*geoIt).GetVMax();
	phiMin = (*geoIt).GetUMin(); 
	phiMax = (*geoIt).GetUMax();
      }
    
    Double_t wei = 1.;
    if ( iValidRange ) wei = WBWeiCompDrDtethaDphi(rMin, rMax, t1, t2, phiMin, phiMax);
    
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

Double_t MatPlot::WBWeiCompDrDtethaDphi(Double_t r1, Double_t r2, Double_t t1, Double_t t2, Double_t phi1, Double_t phi2){
  
  //  Double_t PI=3.141592653589793;
  
  /*
    Double_t part1 = r2*log(sqrt(r2*r2+z2*z2)+z2)-r2*log(sqrt(r2*r2+z1*z1)+z1);
    Double_t part2 = r1*log(sqrt(r1*r1+z2*z2)+z2)-r1*log(sqrt(r1*r1+z1*z1)+z1);
    Double_t part3 = z2*log(sqrt(r2*r2+z2*z2)+r2)-z2*log(sqrt(r1*r1+z2*z2)+r1);
    Double_t part4 = z1*log(sqrt(r2*r2+z1*z1)+r2)-z1*log(sqrt(r1*r1+z1*z1)+r1);
    Double_t wei = 2*PI*(part1-part2+part3-part4);
  */
  

  Double_t wei = (phi2-phi1)*(r2-r1)*(1./tan(t1)-1./tan(t2));
  
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
  
  unfoldBay1DataFakeSub1DH = doUnfoldHist("rUfBay1DataFakeSub", 0, rawDataFakeSub1DH);
  unfoldBay2DataFakeSub1DH = doUnfoldHist("rUfBay2DataFakeSub", 1, rawDataFakeSub1DH);
  unfoldBay3DataFakeSub1DH = doUnfoldHist("rUfBay3DataFakeSub", 3, rawDataFakeSub1DH);
  unfoldBay4DataFakeSub1DH = doUnfoldHist("rUfBay4DataFakeSub", 5, rawDataFakeSub1DH);
  
  std::cout << " rUfBay1DataFakeSub " << unfoldBay1DataFakeSub1DH->GetName() << " " << unfoldBay1DataFakeSub1DH->GetEntries() <<   std::endl;
  std::cout << " rUfBay2DataFakeSub " << unfoldBay2DataFakeSub1DH->GetName() << " " << unfoldBay2DataFakeSub1DH->GetEntries() <<   std::endl;
  std::cout << " rUfBay3DataFakeSub " << unfoldBay3DataFakeSub1DH->GetName() << " " << unfoldBay3DataFakeSub1DH->GetEntries() <<   std::endl;
  std::cout << " rUfBay4DataFakeSub " << unfoldBay4DataFakeSub1DH->GetName() << " " << unfoldBay4DataFakeSub1DH->GetEntries() <<   std::endl;

  unfoldBay1MCFakeSub1DH = doUnfoldHist("rUfBay1MCFakeSub", 0, rawMCFakeSub1DH);
  unfoldBay2MCFakeSub1DH = doUnfoldHist("rUfBay2MCFakeSub", 1, rawMCFakeSub1DH);
  unfoldBay3MCFakeSub1DH = doUnfoldHist("rUfBay3MCFakeSub", 3, rawMCFakeSub1DH);
  unfoldBay4MCFakeSub1DH = doUnfoldHist("rUfBay4MCFakeSub", 5, rawMCFakeSub1DH);

  std::cout << " rUfBay1MCFakeSub " << unfoldBay1MCFakeSub1DH->GetName() << " " << unfoldBay1MCFakeSub1DH->GetEntries() <<   std::endl;
  std::cout << " rUfBay2MCFakeSub " << unfoldBay2MCFakeSub1DH->GetName() << " " << unfoldBay2MCFakeSub1DH->GetEntries() <<   std::endl;
  std::cout << " rUfBay3MCFakeSub " << unfoldBay3MCFakeSub1DH->GetName() << " " << unfoldBay3MCFakeSub1DH->GetEntries() <<   std::endl;
  std::cout << " rUfBay4MCFakeSub " << unfoldBay4MCFakeSub1DH->GetName() << " " << unfoldBay4MCFakeSub1DH->GetEntries() <<   std::endl;

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

}
