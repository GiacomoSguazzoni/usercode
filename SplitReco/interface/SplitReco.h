#ifndef SplitReco_h
#define SplitReco_h

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


#include <TH1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>

//
// class declaration
//

class SplitReco : public edm::EDAnalyzer {
public:
  explicit SplitReco(const edm::ParameterSet&);
  ~SplitReco();
  
  
private:
  virtual void beginJob(const edm::EventSetup & );
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  bool trackPreSelection(const reco::Track &);
  int trackAction(const reco::Track &);
  int trackingParticleAction(TrackingParticleRef &);

  double pErrorOfTrack(const reco::Track &);
  
  // ----------member data ---------------------------
  edm::InputTag tracksTag_; 
  edm::InputTag tpTag_; 
  edm::InputTag simtracksTag_; 
  bool myDebug_;
  bool isMC_;
  int splitTrackEffHits_;
  int minSplits_;
  double kFactor_;
  double pullCut_;

  double ptMinCut_;
  double ptMaxCut_;
  double etaMaxCut_;
  int nHitMinCut_;

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
  TH1F * deltaPHist;
  TH1F * deltaPSimHist;
  TH1F * pullHist;
  TGraphErrors* graphForFit;
  int itrack;

  //
  // Tree stuff
  TTree * splitTree;
  // run/event/track
  int sT_nRun;
  int sT_nEvent;
  int sT_iEvent;
  int sT_NTrack;
  int sT_iTrack;
  // track
  float sT_p;
  float sT_pErr;
  float sT_pAve;
  float sT_pAveErr;
  float sT_pt;
  float sT_ptErr;
  float sT_eta;
  float sT_the;
  float sT_phi;
  float sT_dxy;     
  float sT_dz;    
  float sT_vx;     
  float sT_vy;    
  float sT_vz;    
  int   sT_nHit;
  int   sT_nEffHit;
  int   sT_q;
  float sT_chi2;
  float sT_T;
  float sT_rIn;
  float sT_rOut;
  float sT_zIn;
  float sT_zOut;
  // simtrack
  int   sT_iTrackSim;
  float sT_pSim;
  float sT_ptSim;
  float sT_etaSim;
  float sT_theSim;
  float sT_phiSim;
  float sT_pLossSim;
  int   sT_nHitSim;    
  float sT_TSim;
  float sT_rInSim; //r first hit
  float sT_rOutSim; //r last hit
  float sT_zInSim; //r first hit
  float sT_zOutSim; //r last hit
  // split track method 0 (not overlapping splits)
  float sT_pLossSplit;
  float sT_pLossErrSplit;
  float sT_p0Split;
  float sT_p1Split;
  float sT_chi2Split;
  int sT_freeParSplit;
  int sT_NSplit;
  float sT_TInSplit;
  float sT_TOutSplit;
  // split track method 1 (super split, overlapping splits)
  float sT_pLossSSplit;
  float sT_pLossErrSSplit;
  float sT_p0SSplit;
  float sT_p1SSplit;
  float sT_chi2SSplit;
  int sT_freeParSSplit;
  int sT_NSSplit;
  float sT_TInSSplit;
  float sT_TOutSSplit;

  //
  // Quantities needed for split refits
  //
  TransientTrackingRecHit::RecHitContainer hitsAll;       
  int effHitN;
  
  std::vector<double> vecT;
  std::vector<int> vecEffHitN;
  std::vector<int> vecEffHitBegin;
  std::vector<int> vecEffHitEnd;
  std::vector<double> vecPerpR;
  std::vector<double> vecZ;
  std::vector<double> vecDPerpT;

};

#endif
