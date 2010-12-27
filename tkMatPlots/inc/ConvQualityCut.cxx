//chi2prob>0.0005
//distOfMinApp <1.
//dPhi <0.2
//Nhits1>4 and Nhits2>4 ==> Nhits1>5 and Nhits2>5
//Possibly track pt>0.5 GeV. 
if ( chi2prob>0.0005 && minapp<1. && nHits1>5 && nHits2>5 && pt1>0.5 && pt2>0.5) {iCut = 1;};
//if ( chi2prob>0.0005 && minapp<1. && nHits1>5 && nHits2>5 && deltaphi<0.2 && pt1>0.5 && pt2>0.5) {iCut = 1;};
//if ( chi2prob>0.0005 && minapp<1. && nHits1>4 && nHits2>4 && deltaphi<0.2 && pt1>0.5 && pt2>0.5) {iCut = 1;};
