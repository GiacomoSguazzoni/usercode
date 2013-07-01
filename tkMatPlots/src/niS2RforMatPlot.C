#define niS2RforMatPlot_cxx
#include "niS2RforMatPlot.h"

void niS2RforMatPlot::LoopForFill(TH1* hist)
{
  
  std::cout << " Metodo LoopForFill 1D..." << std::endl;
  
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  
  //GeoCuts iterator
  std::vector<GeoCut>::iterator geoIt;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    for ( unsigned int i = 0; i < numberOfMC_TrkV; i++ )
      {
	//    if ( ! (event >= evRangeMin && event < evRangeMax) ) continue;
	
	Double_t XTrue = MC_TrkV_x->at(i);
	Double_t YTrue = MC_TrkV_y->at(i);
	Double_t ZTrue = MC_TrkV_z->at(i);
	
	Double_t X = XTrue;
	Double_t Y = YTrue;
	Double_t Z = ZTrue;
	
#include "VarArray.cxx"
    
#include "counter.cxx"

    //Build geometrical cut
    Int_t iGeo = 0;
    for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ ){
      iGeo += (*geoIt).GeoCutOk(ACut, BCut);
    }

    // Verify cut
    if ( iGeo && QualityCut(i) ) {
      hist->Fill(A,B); // Ok
      if ( ! MC_TrkV_isAssociatedPF->at(i) ) effUV->NotReco(A,B);// Only if not associated count for (in)efficiency
    }
    //
    
      }

  }
  
}

#ifdef UNFOLD
void niS2RforMatPlot::LoopForTrain(RooUnfoldResponse* response)
{
  
  std::cout << " Metodo LoopForTrain in S2R..." << std::endl;
  
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

    if ( isAssoc == 0 ){
      Double_t radius = sqrt(x*x+y*y);
      
      Int_t iGeo = 0;
      //Geocut enlarged!
      for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ ){
        if ( ((*geoIt).GetVMin()-1.<radius) && ((*geoIt).GetVMax()+1.>radius) && ((*geoIt).GetUMin()-1.<z) && ((*geoIt).GetUMax()+1.>z) ) iGeo = 1;
      }
      if ( !iGeo ) continue;
      if ( !QualityCut() ) continue;
  
      //      std::cout << " pippo " << response->GetDimensionMeasured() << std::endl;
      
      if ( response->GetDimensionMeasured() == 2) { // 1 for 1dim; 2 for 2dim
	response->Miss(x,y);
      } else {
	response->Miss(radius);
      }

    }

  }

}
#endif //#ifdef UNFOLD

