#ifndef ecRefit_h
#define ecRefit_h

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

#include "Tests/ecTest/interface/ec.h"

#include <TH1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>

//
// class declaration
//

class ecRefit{
public:
  ecRefit(const TrackerGeometry *, const MagneticField *, const TrajectoryFitter *, const TransientTrackingRecHitBuilder *, bool);
  ~ecRefit();

  tsosParams doWithTrack(const reco::Track, hitSelector);

 private:
  void setMaterialToKFactor(double);
  double pErrorAtTSOS(TrajectoryStateOnSurface &);
  double ptErrorAtTSOS(TrajectoryStateOnSurface &);
  double pzErrorAtTSOS(TrajectoryStateOnSurface &);
  double paramErrorAtTSOS(TrajectoryStateOnSurface &, int);
  void dumpTSOSInfo(TrajectoryStateOnSurface &);
  tsosParams GetTSOSParams(TrajectoryStateOnSurface &);
  int buildHitsVector(const reco::Track, hitSelector, uint32_t&, uint32_t&, int&);
  TrajectoryStateOnSurface buildInitialStateForRefit(const reco::Track, TransientTrackingRecHit::RecHitContainer, const TrackerGeometry *, const MagneticField *);
  TrajectoryStateOnSurface buildInitialStateForTlRefit(TrajectoryStateOnSurface &, const TrackerGeometry *, const MagneticField *);
  std::vector<Trajectory> doGenericRefit(const reco::Track, TransientTrackingRecHit::RecHitContainer, const TrackerGeometry *, const MagneticField *);
  std::vector<Trajectory> doGenericRefitWithTSOS(TrajectoryStateOnSurface &, TransientTrackingRecHit::RecHitContainer, const TrackerGeometry *, const MagneticField *);

  void dumpModuleInfo(DetId);

  //Debug switch
  bool myDebug_;

  //Track stuff
  //  std::vector<TrajectoryMeasurement> theTrajectoryMeasurements;
  reco::Track theTrack;

  //Handles
  edm::ESHandle<TrackerGeometry> theG;
  edm::ESHandle<MagneticField> theMF;
  edm::ESHandle<TrajectoryFitter> theF;
  edm::ESHandle<TransientTrackingRecHitBuilder> theB;

  //
  // Quantities needed for ec refits
  //
  TransientTrackingRecHit::RecHitContainer hitsAll;       
  TransientTrackingRecHit::RecHitContainer hitsTl;       

};

#endif
