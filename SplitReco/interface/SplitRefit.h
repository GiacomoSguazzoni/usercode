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
  SplitRefit(const TrackerGeometry *, const MagneticField *, const TrajectoryFitter *, const TransientTrackingRecHitBuilder *);
  ~SplitRefit();

  myTrack initializeWithTrack(const reco::Track);
  std::vector<Split> doSplitRefits(int, int, double, double);
  energyLoss energyLossFitOnSplits(std::vector<Split> &);

 private:
  void setMaterialToKFactor(double);
  double pErrorAtTSOS(TrajectoryStateOnSurface &);
  double ptErrorAtTSOS(TrajectoryStateOnSurface &);
  std::vector<Trajectory> doGenericRefit(TransientTrackingRecHit::RecHitContainer, TrajectoryStateOnSurface);


  //Track stuff
  std::vector<TrajectoryMeasurement> theTrajectoryMeasurements;
  TrajectoryStateOnSurface theInitialStateForRefitting;
  reco::Track theTrack;
  double theTrackTheta, theTrackInvSinTheta, theTrackThetaErr;


  bool myDebug_;

  edm::ESHandle<TrackerGeometry> theG;
  edm::ESHandle<MagneticField> theMF;
  edm::ESHandle<TrajectoryFitter> theF;
  edm::ESHandle<TransientTrackingRecHitBuilder> theB;

  //
  // Quantities needed for split refits
  //
  TransientTrackingRecHit::RecHitContainer hitsAll;       

  myTrack mytrack;

  std::vector<double> vecT;
  std::vector<int> vecEffHitN;
  std::vector<int> vecEffHitBegin;
  std::vector<int> vecEffHitEnd;
  std::vector<double> vecPerpR;
  std::vector<double> vecX;
  std::vector<double> vecY;
  std::vector<double> vecZ;
  std::vector<double> vecDPerpT;

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
