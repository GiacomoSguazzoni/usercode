if(
   (isNucl && (nOut-nOutTkStep67Poor-nOutTkStep67Good>0) && (mOut>0.6) && (pt>0.8||pt<0.1) && (ptOut>0.5) && (angle<15))
   ||
   (isNuclLoose && (nOut-nOutTkStep67Poor-nOutTkStep67Good==2) && (mOut>0.6) && (ptOut>0.5) && (angle<15))
   ) 
  {
    iCut = 1;
  };
