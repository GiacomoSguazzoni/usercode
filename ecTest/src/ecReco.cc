//

#include "Tests/ecTest/interface/ecReco.h"
#include "Tests/ecTest/interface/ecRefit.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h" 
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include <TF1.h>

#include <memory>
#include <iostream>
#include <string>
#include <map>
#include <set>

//
// constructors and destructor
//
ecReco::ecReco(const edm::ParameterSet& iConfig)
{
  //Miscellanea
  std::string FileName_ = iConfig.getUntrackedParameter<std::string>("outputFileName", "ecTrack.root");
  myDebug_ = iConfig.getUntrackedParameter<bool>("myDebug", false);
  isMC_ = iConfig.getUntrackedParameter<bool>("isMC", true);

  ptMinCut_ = iConfig.getUntrackedParameter<double>("ptMinCut", 0.8);
  ptMaxCut_ = iConfig.getUntrackedParameter<double>("ptMaxCut", 5.2);
  nHitMinCut_ = iConfig.getUntrackedParameter<int>("nHitMinCut", 6);
  nHitTecMinCut_ = iConfig.getUntrackedParameter<int>("nHitTecMinCut", 6);
  etaMaxCut_ = iConfig.getUntrackedParameter<double>("etaMaxCut", 2.0);
  etaMinCut_ = iConfig.getUntrackedParameter<double>("etaMinCut", -2.0);
  
  fitterName_ = iConfig.getUntrackedParameter<std::string>("Fitter","KFFittingSmootherWithOutliersRejectionAndRK");
  associatorName_ = iConfig.getUntrackedParameter<std::string>("Associator","TrackAssociatorByHits");
  builderName_ = iConfig.getUntrackedParameter<std::string>("Builder","WithAngleAndTemplate");   

  //Tags
  tracksTag_ = iConfig.getUntrackedParameter<edm::InputTag>("tracks");
  if ( isMC_ ) {
    tpTag_ = iConfig.getUntrackedParameter< edm::InputTag >("tp");
  } 
  
  //now do what ever initialization is needed
  ecT_iEvent=0;
  itrack=0;

  file = new TFile(FileName_.c_str(),"recreate");
  
  //
  // Root-upla
  int bs = 64000;

  ecTree = new TTree("ec","ec");
  // event/track
  ecTree->Branch("nRun",       &ecT_nRun,       "nRun/I",       bs);  
  ecTree->Branch("nEvent",     &ecT_nEvent,     "nEvent/I",     bs);  
  ecTree->Branch("iEvent",     &ecT_iEvent,     "iEvent/I",     bs);  
  ecTree->Branch("NTrack",     &ecT_NTrack,     "NTrack/I",     bs);  
  ecTree->Branch("iTrack",     &ecT_iTrack,     "iTrack/I",     bs);  
  // track	     		   		   	       
  ecTree->Branch("p",          &ecT_p,          "p/F",          bs);       
  ecTree->Branch("pErr",       &ecT_pErr,       "pErr/F",       bs);       
  ecTree->Branch("c",          &ecT_curv,       "c/F",          bs);       
  ecTree->Branch("cErr",       &ecT_curvErr,    "cErr/F",       bs);       
  ecTree->Branch("pt",         &ecT_pt,         "pt/F",         bs);      
  ecTree->Branch("ptErr",      &ecT_ptErr,      "ptErr/F",      bs);       
  ecTree->Branch("pz",         &ecT_pz,         "pz/F",         bs);      
  ecTree->Branch("pzErr",      &ecT_pzErr,      "pzErr/F",      bs);       
  ecTree->Branch("the",        &ecT_the,        "the/F",        bs);     
  ecTree->Branch("theErr",     &ecT_theErr,     "theErr/F",     bs);     
  ecTree->Branch("phi",        &ecT_phi,        "phi/F",        bs);     
  ecTree->Branch("phiErr",     &ecT_phiErr,     "phiErr/F",     bs);     
  ecTree->Branch("nHitVal",    &ecT_nHitVal,    "nHitVal/I",    bs);    
  ecTree->Branch("nHitPxl",    &ecT_nHitPxl,    "nHitPxl/I",    bs);    
  ecTree->Branch("nHitTib",    &ecT_nHitTib,    "nHitTib/I",    bs);    
  ecTree->Branch("nHitTob",    &ecT_nHitTob,    "nHitTob/I",    bs);    
  ecTree->Branch("nHitTec",    &ecT_nHitTec,    "nHitTec/I",    bs);    
  ecTree->Branch("nHitTid",    &ecT_nHitTid,    "nHitTid/I",    bs);    
  ecTree->Branch("nHitMis",    &ecT_nHitMis,    "nHitMis/I",    bs);    
  ecTree->Branch("nHitIna",    &ecT_nHitIna,    "nHitIna/I",    bs);    
  ecTree->Branch("nHitBad",    &ecT_nHitBad,    "nHitBad/I",    bs);    
  ecTree->Branch("q",          &ecT_q,          "q/I",          bs);    
  ecTree->Branch("eta",        &ecT_eta,        "eta/F",        bs);     
  ecTree->Branch("dxy",        &ecT_dxy,        "dxy/F",        bs);     
  ecTree->Branch("dz",         &ecT_dz,         "dz/F",         bs);    
  ecTree->Branch("vx",         &ecT_vx,         "vx/F",         bs);     
  ecTree->Branch("vy",         &ecT_vy,         "vy/F",         bs);    
  ecTree->Branch("vz",         &ecT_vz,         "vz/F",         bs);    
  ecTree->Branch("chi2",       &ecT_chi2,       "chi2/F",       bs);    
  // simtrack	     		   		   	       
  if ( isMC_ ){
    ecTree->Branch("iTrackSim",  &ecT_iTrackSim,  "iTrackSim/I",  bs);    
    ecTree->Branch("pSim",       &ecT_pSim,       "pSim/F",       bs);    
    ecTree->Branch("ptSim",      &ecT_ptSim,      "ptSim/F",      bs);   
    ecTree->Branch("theSim",     &ecT_theSim,     "theSim/F",     bs);  
    ecTree->Branch("phiSim",     &ecT_phiSim,     "phiSim/F",     bs);  
    ecTree->Branch("nHitSim",    &ecT_nHitSim,    "nHitSim/I",    bs);    
    ecTree->Branch("hitFrac",    &ecT_hitFrac,    "hitFrac/F",    bs);    
  }
  // trackLet	     		   		   	       
  ecTree->Branch("pTl",          &ecT_pTl,        "pTl/F",          bs);       
  ecTree->Branch("pErrTl",       &ecT_pErrTl,     "pErrTl/F",       bs);       
  ecTree->Branch("cTl",          &ecT_curvTl,     "cTl/F",          bs);       
  ecTree->Branch("cErrTl",       &ecT_curvErrTl,  "cErrTl/F",       bs);       
  ecTree->Branch("ptTl",         &ecT_ptTl,       "ptTl/F",         bs);      
  ecTree->Branch("ptErrTl",      &ecT_ptErrTl,    "ptErrTl/F",      bs);       
  ecTree->Branch("pzTl",         &ecT_pzTl,       "pzTl/F",         bs);      
  ecTree->Branch("pzErrTl",      &ecT_pzErrTl,    "pzErrTl/F",      bs);       
  ecTree->Branch("theTl",        &ecT_theTl,      "theTl/F",        bs);     
  ecTree->Branch("theErrTl",     &ecT_theErrTl,   "theErrTl/F",     bs);     
  ecTree->Branch("phiTl",        &ecT_phiTl,      "phiTl/F",        bs);     
  ecTree->Branch("phiErrTl",     &ecT_phiErrTl,   "phiErrTl/F",     bs);     
  // simtrack	     		   		   	       
  if ( isMC_ ){
    ecTree->Branch("pSimTl",    &ecT_pSimTl,    "pSimTl/F",  bs);    
    ecTree->Branch("ptSimTl",   &ecT_ptSimTl,   "ptSimTl/F",  bs);    
    ecTree->Branch("theSimTl",  &ecT_theSimTl,  "theSimTl/F",  bs);   
    ecTree->Branch("phiSimTl",  &ecT_phiSimTl,  "phiSimTl/F",  bs);  
  }

  /*
  ecTree->Branch("nHitValTl",    &ecT_nHitValTl,  "nHitValTl/I",    bs);    
  ecTree->Branch("nHitMisTl",    &ecT_nHitMisTl,  "nHitMisTl/I",    bs);    
  ecTree->Branch("nHitInaTl",    &ecT_nHitInaTl,  "nHitInaTl/I",    bs);    
  ecTree->Branch("nHitBadTl",    &ecT_nHitBadTl,  "nHitBadTl/I",    bs);    
  */

}


ecReco::~ecReco()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  file->Write();
  file->Close();
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ecReco::analyze(const edm::Event& iEvent, const edm::EventSetup& setup)
{
  using namespace edm;
  using namespace std;

  using reco::Track;
  using reco::TrackCollection;

  // Store run number
  //
  ecT_nRun = iEvent.run();

  setup.get<TrackerDigiGeometryRecord>().get(theG);
  setup.get<IdealMagneticFieldRecord>().get(theMF);  
  setup.get<TrajectoryFitter::Record>().get(fitterName_,theFitter);
  setup.get<TrackAssociatorRecord>().get(associatorName_,theAssociator);
  setup.get<TransientRecHitRecord>().get(builderName_,theBuilder);

  Handle<View<Track> > trackCollectionH;
  iEvent.getByLabel(tracksTag_,trackCollectionH);
  const View<Track>  tC = *(trackCollectionH.product()); 
  
  if ( isMC_ ) {
    edm::Handle<TrackingParticleCollection>  TPCollectionH;
    iEvent.getByLabel(tpTag_,TPCollectionH);
    const TrackingParticleCollection tPC = *(TPCollectionH.product());
    RtSC = theAssociator->associateRecoToSim(trackCollectionH,TPCollectionH,&iEvent);
  }


  ecT_nEvent = iEvent.id().event();
  ecT_iEvent++;

  //
  // Loop on event tracks
  ecT_iTrack=0;
  ecT_NTrack = tC.size();
  for(View<Track>::size_type i=0; i<tC.size(); ++i) {
    edm::RefToBase<reco::Track> itTrack(trackCollectionH, i);

   //
    // Track preselection
    if ( ! trackPreSelection(*itTrack) ) continue;

    ecT_iTrack++;
    
    ecT_iTrackSim=0;
    ecT_pSim=0.;     
    ecT_ptSim=0.;    
    ecT_theSim=0.;   
    ecT_phiSim=0.;   
    ecT_hitFrac=0.;

    uint32_t detid = trackAction(*itTrack);

    if ( isMC_ ) {

      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = RtSC[itTrack];
	if ( myDebug_ ) std::cout << "Reco Track pT: " << itTrack->pt() 
				  << " matched to " << tp.size() << " MC Tracks" << std::endl;
	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  
	  TrackingParticleRef tpr = it->first;
	  ecT_hitFrac = it->second;
	  if ( myDebug_ ) std::cout << "\t\tMCTrack " << tpr.index() << " pT: " << tpr->pt() << 
	    " NShared: " << ecT_hitFrac << std::endl;
	  trackingParticleAction(tpr, detid);
	  
	}
      } catch (cms::Exception event) {
	if ( myDebug_ ) std::cout << "->   Track pT: " << itTrack->pt() 
				  <<  " matched to 0  MC Tracks" << std::endl;
      }

    }



  }
}


// ------------ method called once each job just before starting event loop  ------------
void ecReco::beginJob(const edm::EventSetup & setup)
{

}

// ------------ method called once each job just after ending the event loop  ------------
void ecReco::endJob() {

}

int ecReco::trackingParticleAction(TrackingParticleRef & tpr, uint32_t iDetId){

  TrackingParticle* tp = const_cast<TrackingParticle*>(tpr.get());
	   
  // If more than one associated sim track only the latter enters in the ntupla
  ecT_iTrackSim=tpr.index(); 
  ecT_pSim=tpr->p();      
  ecT_ptSim=tpr->pt();     
  ecT_theSim=tpr->theta();    
  ecT_phiSim=tpr->phi();    
  ecT_nHitSim=tpr->trackPSimHit().size();

  if ( myDebug_ ) std::cout << " SimTrack p:" << tp->p() << " eta:" << tp->momentum().eta() << std::endl;

  //
  // Store state corresponding to the first hit detid

  for(std::vector<PSimHit>::const_iterator TPhit = tp->pSimHit_begin(); TPhit != tp->pSimHit_end(); TPhit++){ //sguazz
    //
    // get the DetUnit via the DetUnitId and cast it to a StripGeomDetUnit
    //
    // Example code from:
    // http://cmslxr.fnal.gov/lxr/source/Validation/TrackerHits/src/TrackerHitAnalyzer.cc
    DetId detid=DetId(TPhit->detUnitId());

    if (detid == iDetId){ 
      
      const GeomDetUnit * det=(const GeomDetUnit*)theG->idToDetUnit( detid );
      GlobalVector gMomentumAtHit = det->toGlobal(TPhit->momentumAtEntry());
      
      ecT_pSimTl=gMomentumAtHit.mag();
      ecT_ptSimTl=gMomentumAtHit.perp();
      ecT_theSimTl=gMomentumAtHit.theta();
      ecT_phiSimTl=gMomentumAtHit.phi();
      
    }
  }

  return 0;

}

bool ecReco::trackPreSelection(const reco::Track & itTrack){

  ecT_p       = itTrack.p();
  ecT_pErr    = pErrorOfTrack(itTrack);
  ecT_curv    = itTrack.qoverp();
  ecT_curvErr = itTrack.qoverpError();
  ecT_pt      = itTrack.pt();
  ecT_ptErr   = itTrack.ptError();
  ecT_pz      = itTrack.pz();
  ecT_pzErr   = pzErrorOfTrack(itTrack);
  ecT_the     = itTrack.theta();
  ecT_theErr  = itTrack.thetaError();
  ecT_phi     = itTrack.phi();
  ecT_phiErr  = itTrack.phiError();
  ecT_eta     = itTrack.eta();     
  ecT_dxy     = itTrack.dxy();     
  ecT_dz      = itTrack.dz();    
  ecT_vx      = itTrack.vx();     
  ecT_vy      = itTrack.vy();     
  ecT_vz      = itTrack.vz();    
  ecT_q       = itTrack.charge();
  ecT_chi2    = itTrack.normalizedChi2();
  ecT_nHitVal = itTrack.numberOfValidHits();

  //
  // hit pattern of the track
  const reco::HitPattern& hp = itTrack.hitPattern();

  //
  // Number of endcap hits
  ecT_nHitPxl = hp.numberOfValidPixelHits();
  ecT_nHitTib = hp.numberOfValidStripTIBHits();
  ecT_nHitTob = hp.numberOfValidStripTOBHits();
  ecT_nHitTid = hp.numberOfValidStripTIDHits();
  ecT_nHitTec = hp.numberOfValidStripTECHits();

  //
  // Track preselection
  if (
      ecT_eta<etaMinCut_ || 
      ecT_eta>etaMaxCut_ || 
      ecT_nHitVal<nHitMinCut_ ||
      ecT_pt<ptMinCut_ ||
      ecT_pt>ptMaxCut_ ||
      ecT_nHitTec<nHitTecMinCut_
      ) return false;

  //
  // Hit quantities
  ecT_nHitMis = hp.numberOfLostHits();
  ecT_nHitIna = hp.numberOfInactiveHits();
  ecT_nHitBad = hp.numberOfBadHits();

  return true;

}

uint32_t ecReco::trackAction(const reco::Track & itTrack){

  ecRefit thisEcRefit(theG.product(), theMF.product(), theFitter.product(), theBuilder.product(), myDebug_);
  tsosParams myPars = thisEcRefit.doWithTrack(itTrack);

  //
  // Tree stuff
  ecT_pTl =        myPars.p;
  ecT_pErrTl =     myPars.pErr;
  ecT_curvTl =     myPars.curv;
  ecT_curvErrTl =  myPars.curvErr;
  ecT_ptTl =       myPars.pt;
  ecT_ptErrTl =    myPars.ptErr;
  ecT_pzTl =       myPars.pz;
  ecT_pzErrTl =    myPars.pzErr;
  ecT_phiTl =      myPars.phi;
  ecT_phiErrTl =   myPars.phiErr;
  ecT_theTl =      myPars.theta;
  ecT_theErrTl =   myPars.thetaErr;

  //
  // Fill ecTree
  ecTree->Fill();
  
  //
  // return the DetId of the first tracklet hit
  return myPars.detid;

}

double ecReco::pErrorOfTrack(const reco::Track & track){
  
  double perror = track.charge()*track.p()*track.qoverpError()/track.qoverp();
  
  return perror;
  
}

double ecReco::pzErrorOfTrack(const reco::Track & track){
  
  //  AlgebraicSymMatrix66 errors = track.cartesianError().matrix();
  //  return sqrt(errors(5,5))

  float lambda = track.lambda();
  float sinl = sin(lambda);
  float cosl = cos(lambda);

  return track.p()*sqrt(track.covariance(0,0)*sinl*sinl+track.covariance(1,1)*cosl*cosl-2.*track.covariance(0,1)*cosl*sinl/track.qoverp());

}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ecReco);
