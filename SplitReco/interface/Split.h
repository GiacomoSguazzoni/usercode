#ifndef Split_H
#define Split_H

struct myTrack {

  double pAve;
  double pAveErr;
  int nHit;
  int nEffHit;
  double T;
  double rIn;
  double rOut;
  double zIn;
  double zOut;

};

struct Split {
  double p;
  double pErr;
  double T;
  double TErr;
};

struct energyLoss {
  double pLoss;
  double pLossErr;
  double p0;
  double p1;
  double chi2;
  double freePar;
  int nsplits;
  double TIn;
  double TOut;
};

#endif


