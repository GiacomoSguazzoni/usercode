#include "counter.cxx"

    //Build geometrical cut
    Int_t iGeo = 0;
    for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ ){
      iGeo += (*geoIt).GeoCutOk(ACut, BCut);
    }

  //Verify cuts
  if ( iGeo && QualityCut() ) { // Within cuts

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
