#ifndef SplitRefit_h
#define SplitRefit_h

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "Tests/SplitReco/interface/Split.h"

#include <TH1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>

//
// class declaration
//

class SplitRefit{
public:
  SplitRefit(const TrackerGeometry *, const MagneticField *, const TrajectoryFitter *, const TransientTrackingRecHitBuilder *, bool);
  ~SplitRefit();

  myTrack initializeWithTrack(const reco::Track);
  std::vector<Split> doSplitRefits(int, int, double, bool);
  std::vector<energyLoss> energyLossFitOnSplits(std::vector<Split> &);
  energyLoss energyLossFitByPtOnSplits(std::vector<Split> &);
  energyLoss energyLossFitByCurvOnSplits(std::vector<Split> &, int);

 private:
  void setMaterialToKFactor(double);
  double pErrorAtTSOS(TrajectoryStateOnSurface &);
  double ptErrorAtTSOS(TrajectoryStateOnSurface &);
  double pzErrorAtTSOS(TrajectoryStateOnSurface &);
  double lambdaErrorAtTSOS(TrajectoryStateOnSurface &);
  void dumpTSOSInfo(TrajectoryStateOnSurface &);
  void GetPAveFromMeasurements(std::vector<TrajectoryMeasurement>, double&, double&, double&, double&, double&, double&);
  AlgebraicSymMatrix55 rescaleErrorOfTransComponent(TrajectoryStateOnSurface&, double, int, int);
  AlgebraicSymMatrix55 rescaleErrorOfLongComponent(TrajectoryStateOnSurface&, double);
  std::vector<Trajectory> doGenericRefit(TransientTrackingRecHit::RecHitContainer, TrajectoryStateOnSurface);

  void dumpModuleInfo(DetId);


  //Debug switch
  bool myDebug_;

  //Track stuff
  std::vector<TrajectoryMeasurement> theTrajectoryMeasurements;
  TrajectoryStateOnSurface theInitialStateForRefitting;
  reco::Track theTrack;
  double theTrackTheta, theTrackInvSinTheta, theTrackThetaErr;
  double theTrackP, theTrackPSquared, theTrackPt, theTrackPtSquared;
  double BInTesla, theTrackLambda, theTrackInvSinLambda, theTrackInvCosLambda, theTrackCharge;
  myTrack mytrack;

  //Handles
  edm::ESHandle<TrackerGeometry> theG;
  edm::ESHandle<MagneticField> theMF;
  edm::ESHandle<TrajectoryFitter> theF;
  edm::ESHandle<TransientTrackingRecHitBuilder> theB;

  //
  // Quantities needed for split refits
  //
  TransientTrackingRecHit::RecHitContainer hitsAll;       
  std::vector<double> vecT;
  std::vector<int> vecEffHitN;
  std::vector<int> vecEffHitBegin;
  std::vector<int> vecEffHitEnd;
  std::vector<double> vecPerpR;
  std::vector<double> vecX;
  std::vector<double> vecY;
  std::vector<double> vecZ;
  std::vector<double> vecDPerpT;

  // For TSOS error handling
  //
  AlgebraicMatrix55 jacoInvPToInvPt, jacoInvPtToInvP;
  AlgebraicMatrix55 jacoInvPToInvPtTransposed, jacoInvPtToInvPTransposed;

  //
  //Various quantities
  double pAve_;   
  double pAveErr_;
  double nEffHit_;
  double T_;      
  double rIn_;    
  double rOut_;   
  double zIn_;    
  double zOut_;   


};

#endif
