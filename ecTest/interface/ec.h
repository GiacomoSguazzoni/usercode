#ifndef ec_H
#define ec_H

struct tsosParams {
  double p;
  double pErr;
  double curv;
  double curvErr;
  double pt;
  double ptErr;
  double pz;
  double pzErr;
  double phi;
  double phiErr;
  double theta;
  double thetaErr;
  uint32_t detid;
  
  void zero(){
    p=0.;
    pErr=0.;
    curv=0.;
    curvErr=0.;
    pt=0.;
    ptErr=0.;
    pz=0.;
    pzErr=0.;
    phi=0.;
    phiErr=0.;
    theta=0.;
    thetaErr=0.;
    detid=0;
  };
  
};

#endif


