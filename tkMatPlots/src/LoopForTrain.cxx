if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast();

Long64_t nbytes = 0, nb = 0;

//GeoCuts iterator
std::vector<GeoCut>::iterator geoIt;

for (Long64_t jentry=0; jentry<nentries;jentry++) {
  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) break;
  nb = fChain->GetEntry(jentry);   nbytes += nb;
  
  Double_t radiusTrue = sqrt(x*x+y*y);
  Double_t xMeas = -1.;
  Double_t yMeas = -1.;
  Double_t radiusMeas = -1.;

  if ( isAssoc == 1 ){
    
    xMeas = x - dx;
    yMeas = y - dy;
    radiusMeas = sqrt(xMeas*xMeas + yMeas*yMeas);
    
    Int_t iGeo = 0;
    for ( geoIt=geoCuts->begin() ; geoIt < geoCuts->end(); geoIt++ ){
      if ( ((*geoIt).GetVMin()<radius) && ((*geoIt).GetVMax()>radius) && ((*geoIt).GetUMin()<z) && ((*geoIt).GetUMax()>z) ){ iGeo = 1; };
    }
    
    if ( !iGeo ) continue;
    if ( !QualityCut() ) continue;
  
  /*
if (measurement_ok)
  response.Fill (x_measured, x_true);
else
  */

  } else {
    response->Miss(radius);
  }

  
 }

