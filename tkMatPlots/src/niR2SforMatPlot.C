#define niR2SforMatPlot_cxx
#include "niR2SforMatPlot.h"

void niR2SforMatPlot::LoopForFill(TH1* hist, TH1* hist2 = 0, TH1* hist3 = 0)
{

  std::cout << " Metodo LoopForFill 1D..." << std::endl;
  
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

    for ( unsigned int i = 0; i < numberOfPFDV; i++ )
      {

	Double_t X = PFDV_x->at(i)-x0;
	Double_t Y = PFDV_y->at(i)-y0;
	Double_t Z = PFDV_y->at(i)-z0;
	
	Double_t XTrue = 0.;
	Double_t YTrue = 0.;
	Double_t ZTrue = 0.;
	Int_t myAssoc = 0;

	if ( PFDV_isAssociatedMC->at(i) ){
	  //&& fChain->GetBranchStatus( "PFDV_associationMC_TrkVIdx" )

	  int indexAss = PFDV_associationMC_TrkVIdx->at(i); 
	  std::cout << " we have simulation " << indexAss;
	  MC_TrkV_x->at(indexAss);
	  MC_TrkV_y->at(indexAss);
	  MC_TrkV_z->at(indexAss);
	  if ( PFDV_isAssociatedMC->at(i) && MC_TrkV_isNuclearInteraction->at(indexAss) ) myAssoc = 1;

	}

#include "VarArray.cxx"
	
#include "counter.cxx"

    //Build geometrical cut
    Int_t iGeo = 0;
    for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ ){
      iGeo += (*geoIt).GeoCutOk(ACut, BCut);
    }

  //Verify cuts
    if ( iGeo && QualityCut(i) ) { // Within cuts

      hist->Fill(A,B);
      naccev += 1;
      if ( ! hist2 ) effUV->Cand(A,B); // hist2 not there means data

      if ( myAssoc ) { // True or fake?

	//	Double_t radiusTrue = sqrt((x+deltax)*(x+deltax)+(y+deltay)*(y+deltay));

	effUV->Success(ATrue, BTrue, A, B);

	if ( hist3 ) {
	  hist3->Fill(ATrue,BTrue); // for Reco vs. True Position
	}
      } else { //Not associated

	if ( hist2 ) {

	  hist2->Fill(A,B); // for Fakes
	  effUV->Fake(A,B);

	}
      }

    } else { // Not satisfying cuts
      if ( myAssoc ) { // If true count as inefficiency
	//	Double_t radiusTrue = sqrt((x+deltax)*(x+deltax)+(y+deltay)*(y+deltay));
	effUV->Fail(ATrue, BTrue);
      }
    }
    //

	
	
      }
  }  

  std::cout << " Number af accepted candidates in the step " << naccev << std::endl;

}

#ifdef UNFOLD
void niR2SforMatPlot::LoopForTrain(RooUnfoldResponse* response)
{
  
  std::cout << " Metodo LoopForTrain in R2S..." << std::endl;
  
  std::cout << " minE " << evRangeMin << " maxE " << evRangeMax << std::endl;
  std::cout << " x0 " << x0 << " y0 " << y0 << " z0 " << z0 << std::endl;
  
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

#include "counter.cxx"

    Double_t X = x-x0;
    Double_t Y = y-y0;
    Double_t Z = z-z0;

    if ( isAssoc && isNuclSim ){
      Double_t xTrue = x+deltax;
      Double_t yTrue = y+deltay;
      Double_t radiusMeas = sqrt(X*X+Y*Y);
      Double_t radiusTrue = sqrt((x+deltax)*(x+deltax)+(y+deltay)*(y+deltay));
      
      Int_t iGeo = 0;
      for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ ){
	iGeo += (*geoIt).GeoCutOk(Z, radiusMeas);
      }
      
      if ( iGeo && QualityCut(i) ) {
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

