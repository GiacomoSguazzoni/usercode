#ifndef ec_H
#define ec_H

struct hitSelector {
  int minhits;
  bool tib;
  bool tid;
  bool tob;
  bool tec;
  bool stereo;
  bool all;

};

struct tsosParams {
  int iok;
  int nhit;
  double p;
  double pErr;
  double curv;
  double curvErr;
  double curvt;
  double curvtErr;
  double pt;
  double ptErr;
  double pz;
  double pzErr;
  double phi;
  double phiErr;
  double theta;
  double thetaErr;
  uint32_t detidTl;
  uint32_t detidIs;
  
  void zero(){
    nhit=0;
    iok=0; 
    p=0.;
    pErr=0.;
    curv=0.;
    curvErr=0.;
    curvt=0.;
    curvtErr=0.;
    pt=0.;
    ptErr=0.;
    pz=0.;
    pzErr=0.;
    phi=0.;
    phiErr=0.;
    theta=0.;
    thetaErr=0.;
    detidTl=0;
    detidIs=0;
  };
  
};

#endif


