#include "counter.cxx"

    //Build geometrical cut
    Int_t iGeo = 0;
    for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ ){
      iGeo += (*geoIt).GeoCutOk(radius, Z);
    }

    //Verify cuts
    if ( iGeo && QualityCut() ) { // Both verified
      hist->Fill(A,B);
      naccev += 1;
      if ( ! hist2 ) effR->Cand(radius); // hist2 not there means data
      if ( myAssoc ) { // True or fake?
	Double_t radiusTrue = sqrt((x+deltax)*(x+deltax)+(y+deltay)*(y+deltay));
	effR->Succeed(radiusTrue, radius);
      } else {
	if ( hist2 ) {
	  hist2->Fill(A,B); // for Fakes
	  effR->Fake(radius);
	}
      }
    } else { // Not satisfying cuts
      if ( myAssoc ) { // If true count as inefficiency
	Double_t radiusTrue = sqrt((x+deltax)*(x+deltax)+(y+deltay)*(y+deltay));
	effR->Fail(radiusTrue);
      }
    }
    //
