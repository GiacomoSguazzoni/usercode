#include "counter.cxx"

    //Build geometrical cut
    Int_t iGeo = 0;
    for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ ){
      iGeo += (*geoIt).GeoCutOk(ACut, BCut);
    }

    // Verify cut
    if ( iGeo && QualityCut() ) {
      hist->Fill(A,B); // Ok
      if ( ! isAssoc ) effUV->NotReco(A,B);// Only if not associated count for (in)efficiency
    }
    //
