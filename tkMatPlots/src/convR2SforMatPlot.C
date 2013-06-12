#define convR2SforMatPlot_cxx
#include "convR2SforMatPlot.h"

void convR2SforMatPlot::LoopForFill(TH1* hist, TH1* hist2 = 0, TH1* hist3 = 0)
{

  std::cout << " Metodo LoopForFill 1D..." << std::endl;

  std::cout << " minE " << evRangeMin << " maxE " << evRangeMax << std::endl;
  std::cout << " x0 " << x0 << " y0 " << y0 << " z0 " << z0 << std::endl;
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;

  //Number of accepted events for debug
  Long64_t naccev = 0;
  
  //GeoCuts iterator
  std::vector<GeoCut>::iterator geoIt;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    //    if ( ! (event >= evRangeMin && event < evRangeMax) ) continue;
    
    Double_t X = x-x0;
    Double_t Y = y-y0;
    Double_t Z = z-z0;

    Double_t XTrue = x+deltax-x0;
    Double_t YTrue = y+deltay-y0;
    Double_t ZTrue = z+deltaz-z0;
    
    Int_t myAssoc = isAssoc;

#include "VarArray.cxx"
    
#include "LoopForFillR2S.cxx"

  }

  std::cout << " Number af accepted candidates in the step " << naccev << std::endl;

}

void convR2SforMatPlot::LoopForFill(TH2* hist, TH2* hist2 = 0, TH2* hist3 = 0)
{
  
  std::cout << " R2S: Metodo LoopForFill 2D..." << std::endl;
  
  std::cout << " minE " << evRangeMin << " maxE " << evRangeMax << std::endl;
  std::cout << " x0 " << x0 << " y0 " << y0 << " z0 " << z0 << std::endl;

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  
  //Number of accepted events for debug
  Long64_t naccev = 0;
  
  //GeoCuts iterator
  std::vector<GeoCut>::iterator geoIt;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    // Verify event range
    // if ( ! (event >= evRangeMin && event < evRangeMax) ) continue;
    
    Double_t X = x-x0;
    Double_t Y = y-y0;
    Double_t Z = z-z0;

    Double_t XTrue = x+deltax-x0;
    Double_t YTrue = y+deltay-y0;
    Double_t ZTrue = z+deltaz-z0;
    
    Int_t myAssoc = isAssoc;

#include "VarArray.cxx"
    
#include "LoopForFillR2S.cxx"

  }

  std::cout << " Number af accepted candidates in the step " << naccev << std::endl;

}

#ifdef UNFOLD
void convR2SforMatPlot::LoopForTrain(RooUnfoldResponse* response)
{
  
  std::cout << " Metodo LoopForTrain in R2S..." << std::endl;
  
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  
  //GeoCuts iterator
  std::vector<GeoCut>::iterator geoIt;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    //    if ( ! (event >= evRangeMin && event < evRangeMax) ) continue;

    Double_t X = x-x0;
    Double_t Y = y-y0;
    Double_t Z = z-z0;

#include "counter.cxx"

    if ( isAssoc ){
      Double_t xTrue = x+deltax;
      Double_t yTrue = y+deltay;
      Double_t radiusMeas = sqrt(X*X+Y*Y);
      Double_t radiusTrue = sqrt((x+deltax)*(x+deltax)+(y+deltay)*(y+deltay));
      
      Int_t iGeo = 0;
      for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ ){
	iGeo += (*geoIt).GeoCutOk(Z, radiusMeas);
      }
      
      if ( iGeo && QualityCut() ) {
	//
	if ( response->GetDimensionMeasured() == 2) { // 1 for 1dim; 2 for 2dim
	  response->Fill(X, Y, xTrue, yTrue);
	} else {
	  response->Fill(radiusMeas, radiusTrue);
	}
	//
      } else {  
	//
	if ( response->GetDimensionMeasured() == 2) { // 1 for 1dim; 2 for 2dim
	  response->Miss(xTrue, yTrue);
	} else {
	  response->Miss(radiusTrue);
	}
	//
      }
    }

  }

}
#endif //#ifdef UNFOLD

