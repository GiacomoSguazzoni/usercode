#ifndef ecReco_h
#define ecReco_h

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

#include "Tests/ecTest/interface/ec.h"

#include <TH1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>

//
// class declaration
//

class ecReco : public edm::EDAnalyzer {
public:
  explicit ecReco(const edm::ParameterSet&);
  ~ecReco();
  
  
private:
  virtual void beginJob(const edm::EventSetup & );
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  bool trackPreSelection(const reco::Track &);
  void trackAction(const reco::Track &, uint32_t&, uint32_t&);
  int trackingParticleAction(TrackingParticleRef &, uint32_t, uint32_t);

  double pErrorOfTrack(const reco::Track &);
  double pzErrorOfTrack(const reco::Track &);
  
  void dumpModuleInfo(DetId);

  // ----------member data ---------------------------
  edm::InputTag tracksTag_; 
  edm::InputTag tpTag_; 
  bool myDebug_;
  bool isMC_;
  double ptMinCut_;
  double ptMaxCut_;
  double etaMaxCut_;
  double etaMinCut_;
  int nHitMinCut_;
  int nHitTecMinCut_;
  int nHitSteTecMinCut_;
  int nHitTidMinCut_;
  int nHitSteTidMinCut_;
  int nHitTibMinCut_;
  int nHitSteTibMinCut_;
  int nHitTobMinCut_;
  int nHitSteTobMinCut_;
  hitSelector hs;

  std::string fitterName_;
  std::string associatorName_;
  std::string builderName_;

  edm::ESHandle<TrackerGeometry> theG;
  edm::ESHandle<MagneticField> theMF;
  edm::ESHandle<TrajectoryFitter> theFitter;
  edm::ESHandle<TrackAssociatorBase> theAssociator;
  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;

  reco::RecoToSimCollection RtSC;

  TFile * file;
  int itrack;
  int nTracksInSample;
  //
  // Tree stuff
  TTree * ecTree;
  //
  // run/event/track
  int ecT_nRun;
  int ecT_nEvent;
  int ecT_iEvent;
  int ecT_NTrack;
  int ecT_iTrack;
  //
  // track
  float ecT_p;
  float ecT_pErr;
  float ecT_curv;
  float ecT_curvErr;
  float ecT_curvt;
  float ecT_curvtErr;
  float ecT_pt;
  float ecT_ptErr;
  float ecT_pz;
  float ecT_pzErr;
  float ecT_the;
  float ecT_theErr;
  float ecT_phi;
  float ecT_phiErr;
  float ecT_q;
  float ecT_eta;     
  float ecT_dxy;     
  float ecT_dz;    
  float ecT_vx;     
  float ecT_vy;    
  float ecT_vz;    
  //
  // inner state
  float ecT_pIs;
  float ecT_pErrIs;
  float ecT_curvIs;
  float ecT_curvErrIs;
  float ecT_curvtIs;
  float ecT_curvtErrIs;
  float ecT_ptIs;
  float ecT_ptErrIs;
  float ecT_pzIs;
  float ecT_pzErrIs;
  float ecT_theIs;
  float ecT_theErrIs;
  float ecT_phiIs;
  float ecT_phiErrIs;
  //
  // 
  int   ecT_nHitVal;
  int   ecT_nHitTl;
  int   ecT_nHitPxl;
  int   ecT_nHitTib;
  int   ecT_nHitSteTib;
  int   ecT_nHitTid;
  int   ecT_nHitSteTid;
  int   ecT_nHitTob;
  int   ecT_nHitSteTob;
  int   ecT_nHitTec;
  int   ecT_nHitSteTec;
  int   ecT_nHitMis;
  int   ecT_nHitIna;
  int   ecT_nHitBad;
  float ecT_chi2;
  //
  // simtrack
  int   ecT_iTrackSim;
  float ecT_pSim;
  float ecT_ptSim;
  float ecT_theSim;
  float ecT_phiSim;
  int   ecT_nHitSim;    
  float ecT_hitFrac;    
  //
  // tracklet
  int   ecT_iokTl;
  float ecT_pTl;
  float ecT_pErrTl;
  float ecT_curvTl;
  float ecT_curvErrTl;
  float ecT_curvtTl;
  float ecT_curvtErrTl;
  float ecT_ptTl;
  float ecT_ptErrTl;
  float ecT_pzTl;
  float ecT_pzErrTl;
  float ecT_theTl;
  float ecT_theErrTl;
  float ecT_phiTl;
  float ecT_phiErrTl;
  //
  // simtracklet
  float ecT_pSimTl;
  float ecT_ptSimTl;
  float ecT_theSimTl;
  float ecT_phiSimTl;
  //
  // sim at inner state
  float ecT_pSimIs;
  float ecT_ptSimIs;
  float ecT_theSimIs;
  float ecT_phiSimIs;

};

#endif
