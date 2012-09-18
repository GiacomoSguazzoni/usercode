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
  double xIn;
  double yIn;
  double zIn;
  double xOut;
  double yOut;
  double zOut;
  double maxhitchi2;
};

struct Split {
  double pt;
  double ptErr;
  double curv;
  double curvErr;
  double curvz;
  double curvzErr;
  double curvt;
  double curvtErr;
  double T;
  double TErr;
};

struct energyLoss {
  double pLoss;
  double pLossErr;
  double dpdx;
  double dpdxErr;
  double p0;
  double chi2;
  double freePar;
  int nsplits;
  double TIn;
  double TOut;
};

#endif


