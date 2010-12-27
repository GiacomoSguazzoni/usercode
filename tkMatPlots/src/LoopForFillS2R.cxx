#include "counter.cxx"

    //Build geometrical cut
    Int_t iGeo = 0;
    for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ ){
      iGeo += (*geoIt).GeoCutOk(radius, Z);
    }

    // Verify cut
    if ( iGeo && QualityCut() ) {
      hist->Fill(A,B); // Ok
      if ( ! isAssoc ) effR->Fail(radius);// Only if not associated count for (in)efficiency
    }
    //
