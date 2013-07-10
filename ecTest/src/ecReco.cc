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
  nHitMinCut_ = iConfig.getUntrackedParameter<int>("nHitMinCut", 8);
  nHitTibMinCut_ = iConfig.getUntrackedParameter<int>("nHitTibMinCut", 0);
  nHitSteTibMinCut_ = iConfig.getUntrackedParameter<int>("nHitSteTibMinCut", 0);
  nHitTobMinCut_ = iConfig.getUntrackedParameter<int>("nHitTobMinCut", 0);
  nHitSteTobMinCut_ = iConfig.getUntrackedParameter<int>("nHitSteTobMinCut", 0);
  nHitTecMinCut_ = iConfig.getUntrackedParameter<int>("nHitTecMinCut", 0);
  nHitSteTecMinCut_ = iConfig.getUntrackedParameter<int>("nHitSteTecMinCut", 0);
  nHitTidMinCut_ = iConfig.getUntrackedParameter<int>("nHitTidMinCut", 0);
  nHitSteTidMinCut_ = iConfig.getUntrackedParameter<int>("nHitSteTidMinCut", 0);
  etaMaxCut_ = iConfig.getUntrackedParameter<double>("etaMaxCut", 2.0);
  etaMinCut_ = iConfig.getUntrackedParameter<double>("etaMinCut", -2.0);

  
  hs.minhits = iConfig.getUntrackedParameter<int>("nHitTlMinCut", 6);
  hs.tid = iConfig.getUntrackedParameter<bool>("useTID", true);
  hs.tec = iConfig.getUntrackedParameter<bool>("useTEC", true);
  hs.tib = iConfig.getUntrackedParameter<bool>("useTIB", true);
  hs.tob = iConfig.getUntrackedParameter<bool>("useTOB", true);
  hs.mono = iConfig.getUntrackedParameter<bool>("useMono", true);
  hs.stereo = iConfig.getUntrackedParameter<bool>("useStereo", true);
  hs.alt = iConfig.getUntrackedParameter<bool>("useStereoIfGlued", false);
  hs.all = iConfig.getUntrackedParameter<bool>("useAll", false);

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

  //Inizialize track counter
  nTracksInSample = 0;
  std::cout << " >>>>>>>> Total number of track is " << nTracksInSample << std::endl;

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
  ecTree->Branch("ct",         &ecT_curvt,      "ct/F",         bs);       
  ecTree->Branch("ctErr",      &ecT_curvtErr,   "ctErr/F",      bs);       
  ecTree->Branch("pt",         &ecT_pt,         "pt/F",         bs);      
  ecTree->Branch("ptErr",      &ecT_ptErr,      "ptErr/F",      bs);       
  ecTree->Branch("pz",         &ecT_pz,         "pz/F",         bs);      
  ecTree->Branch("pzErr",      &ecT_pzErr,      "pzErr/F",      bs);       
  ecTree->Branch("the",        &ecT_the,        "the/F",        bs);     
  ecTree->Branch("theErr",     &ecT_theErr,     "theErr/F",     bs);     
  ecTree->Branch("phi",        &ecT_phi,        "phi/F",        bs);     
  ecTree->Branch("phiErr",     &ecT_phiErr,     "phiErr/F",     bs);     
  // track inner state
  ecTree->Branch("pIs",          &ecT_pIs,        "pIs/F",          bs);       
  ecTree->Branch("pErrIs",       &ecT_pErrIs,     "pErrIs/F",       bs);       
  ecTree->Branch("cIs",          &ecT_curvIs,     "cIs/F",          bs);       
  ecTree->Branch("cErrIs",       &ecT_curvErrIs,  "cErrIs/F",       bs);       
  ecTree->Branch("ctIs",         &ecT_curvtIs,    "ctIs/F",         bs);       
  ecTree->Branch("ctErrIs",      &ecT_curvtErrIs, "ctErrIs/F",      bs);       
  ecTree->Branch("ptIs",         &ecT_ptIs,       "ptIs/F",         bs);      
  ecTree->Branch("ptErrIs",      &ecT_ptErrIs,    "ptErrIs/F",      bs);       
  ecTree->Branch("pzIs",         &ecT_pzIs,       "pzIs/F",         bs);      
  ecTree->Branch("pzErrIs",      &ecT_pzErrIs,    "pzErrIs/F",      bs);       
  ecTree->Branch("theIs",        &ecT_theIs,      "theIs/F",        bs);     
  ecTree->Branch("theErrIs",     &ecT_theErrIs,   "theErrIs/F",     bs);     
  ecTree->Branch("phiIs",        &ecT_phiIs,      "phiIs/F",        bs);     
  ecTree->Branch("phiErrIs",     &ecT_phiErrIs,   "phiErrIs/F",     bs);     
  //
  ecTree->Branch("nHitVal",    &ecT_nHitVal,    "nHitVal/I",    bs);    
  ecTree->Branch("nHitTl",     &ecT_nHitTl,     "nHitTl/I",     bs);    
  ecTree->Branch("nHitPxl",    &ecT_nHitPxl,    "nHitPxl/I",    bs);    
  ecTree->Branch("nHitTib",    &ecT_nHitTib,    "nHitTib/I",    bs);    
  ecTree->Branch("nHitSteTib", &ecT_nHitSteTib, "nHitSteTib/I", bs);     
  ecTree->Branch("nHitTob",    &ecT_nHitTob,    "nHitTob/I",    bs);    
  ecTree->Branch("nHitSteTob", &ecT_nHitSteTob, "nHitSteTob/I", bs);     
  ecTree->Branch("nHitTid",    &ecT_nHitTid,    "nHitTid/I",    bs);    
  ecTree->Branch("nHitSteTid", &ecT_nHitSteTid, "nHitSteTid/I", bs);     
  ecTree->Branch("nHitTec",    &ecT_nHitTec,    "nHitTec/I",    bs);    
  ecTree->Branch("nHitSteTec", &ecT_nHitSteTec, "nHitSteTec/I", bs);     
  ecTree->Branch("nHitMis",    &ecT_nHitMis,    "nHitMis/I",    bs);    
  ecTree->Branch("nHitIna",    &ecT_nHitIna,    "nHitIna/I",    bs);    
  ecTree->Branch("nHitBad",    &ecT_nHitBad,    "nHitBad/I",    bs);    
  ecTree->Branch("q",          &ecT_q,          "q/F",          bs);    
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
  ecTree->Branch("iok",          &ecT_iokTl,      "iok/I",          bs);       
  ecTree->Branch("pTl",          &ecT_pTl,        "pTl/F",          bs);       
  ecTree->Branch("pErrTl",       &ecT_pErrTl,     "pErrTl/F",       bs);       
  ecTree->Branch("cTl",          &ecT_curvTl,     "cTl/F",          bs);       
  ecTree->Branch("cErrTl",       &ecT_curvErrTl,  "cErrTl/F",       bs);       
  ecTree->Branch("ctTl",         &ecT_curvtTl,    "ctTl/F",         bs);       
  ecTree->Branch("ctErrTl",      &ecT_curvtErrTl, "ctErrTl/F",      bs);       
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
    ecTree->Branch("pSimIs",    &ecT_pSimIs,    "pSimIs/F",  bs);    
    ecTree->Branch("ptSimIs",   &ecT_ptSimIs,   "ptSimIs/F",  bs);    
    ecTree->Branch("theSimIs",  &ecT_theSimIs,  "theSimIs/F",  bs);   
    ecTree->Branch("phiSimIs",  &ecT_phiSimIs,  "phiSimIs/F",  bs);  
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

  std::cout << " >>>>>>>> This sample contains " << nTracksInSample << " tracks in " << ecT_iEvent << " events. Ciao." << std::endl;
 
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
  ecT_iTrack = 0;
  ecT_NTrack = tC.size();
  nTracksInSample += tC.size();
  if ( myDebug_ ) std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  if ( myDebug_ ) std::cout << ">>>>>> This event contains " << ecT_NTrack << " reconstructed tracks and the total is " << nTracksInSample << "." << std::endl;

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
    
    if ( myDebug_ ) std::cout << ">>>>>>---------------------------------------------------------------------------------" << std::endl;
    uint32_t detidTl=0;
    uint32_t detidIs=0;
    trackAction(*itTrack, detidTl, detidIs);
    
    if ( isMC_ ) {
      
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = RtSC[itTrack];
	if ( myDebug_ ) std::cout << "   Reco Track pT: " << itTrack->pt() 
				  << " matched to " << tp.size() << " MC Tracks" << std::endl;
	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  
	  TrackingParticleRef tpr = it->first;
	  ecT_hitFrac = it->second;
	  if ( myDebug_ ) std::cout << "\t\tMCTrack " << tpr.index() << " pT: " << tpr->pt() << 
	    " NShared: " << ecT_hitFrac << std::endl;
	  int itpaout = trackingParticleAction(tpr, detidTl, detidIs);
	  if ( itpaout ) {
	    if ( myDebug_ ) std::cout << " Track rejected due to psimhit double match " << std::endl; 
	    ecT_iokTl=0;
	  }

	}
      } catch (cms::Exception event) {
	if ( myDebug_ ) std::cout << "->   Track pT: " << itTrack->pt() 
				  <<  " matched to 0 MC Tracks" << std::endl;
      }
      
    }
    
    //
    // Fill ecTree
    if ( ecT_iokTl ) ecTree->Fill();

  }

}


// ------------ method called once each job just before starting event loop  ------------
void ecReco::beginJob(const edm::EventSetup & setup) {

}

// ------------ method called once each job just after ending the event loop  ------------
void ecReco::endJob() {

}

int ecReco::trackingParticleAction(TrackingParticleRef & tpr, uint32_t iDetIdTl, uint32_t iDetIdIs){

  TrackingParticle* tp = const_cast<TrackingParticle*>(tpr.get());
	   
  // If more than one associated sim track only the latter enters in the ntupla
  ecT_iTrackSim=tpr.index(); 
  ecT_pSim=tpr->p();      
  ecT_ptSim=tpr->pt();     
  ecT_theSim=tpr->theta();    
  ecT_phiSim=tpr->phi();    
  ecT_nHitSim=tpr->trackPSimHit().size();

  if ( myDebug_ ) std::cout << "   SimTrack p:" << ecT_pSim  
			    << " pt:" << ecT_ptSim 
			    << " the:" << ecT_theSim
			    << " phi:" << ecT_phiSim
			    << std::endl;

  //
  // Store state corresponding to the first hit detid

  int iPSHit = 0;
  int iFirstDetIdMatch = 0;
  int iDoubleMatch = 0;

  for(std::vector<PSimHit>::const_iterator TPhit = tp->pSimHit_begin(); TPhit != tp->pSimHit_end(); TPhit++){ //sguazz
    //
    // get the DetUnit via the DetUnitId and cast it to a StripGeomDetUnit
    //
    // Example code from:
    // http://cmslxr.fnal.gov/lxr/source/Validation/TrackerHits/src/TrackerHitAnalyzer.cc
    DetId detid=DetId(TPhit->detUnitId());


    if ( ! (detid.det() == DetId::Tracker) ) {
      if ( myDebug_ ) std::cout << "   PSimHit #" << iPSHit << " in #det " << detid.det() << "; not in Tracker detector " << std::endl;
      iPSHit++;
      continue;
    }

    const GeomDetUnit * det=(const GeomDetUnit*)theG->idToDetUnit( detid );
    GlobalVector gMomentumAtHit = det->toGlobal(TPhit->momentumAtEntry());

    if ( myDebug_ ) {
      std::cout << "   PSimHit #" << iPSHit;
      dumpModuleInfo(detid);
      std::cout << " p:" << gMomentumAtHit.mag() << 
	" pt:" << gMomentumAtHit.perp() <<
	" the:" << gMomentumAtHit.theta() <<
	" phi:" << gMomentumAtHit.phi() <<
	" ptype:" << TPhit->particleType() << " tkid:" << TPhit->trackId() << " proctyp:" << TPhit->processType();
    }


    if (detid == iDetIdIs){ 
      
      ecT_pSimIs=gMomentumAtHit.mag();
      ecT_ptSimIs=gMomentumAtHit.perp();
      ecT_theSimIs=gMomentumAtHit.theta();
      ecT_phiSimIs=gMomentumAtHit.phi();

      if ( myDebug_ ) std::cout << "   Found detId match Is " << detid << " == " << iDetIdIs; 

    }

    if (detid == iDetIdTl){ 
      
      ecT_pSimTl=gMomentumAtHit.mag();
      ecT_ptSimTl=gMomentumAtHit.perp();
      ecT_theSimTl=gMomentumAtHit.theta();
      ecT_phiSimTl=gMomentumAtHit.phi();

      if ( myDebug_ ) std::cout << "   Found detId match Tl " << detid << " == " << iDetIdTl; 
      if ( iFirstDetIdMatch ) { 
	std::cout << " THIS IS A DOUBLE MATCH!!!" << std::endl;
	iDoubleMatch = 1;
      }

      iFirstDetIdMatch = 1;

    }

    if ( myDebug_ ) std::cout << std::endl;
    iPSHit++;

  }

  if ( iDoubleMatch ) return 1;
  return 0;

}

bool ecReco::trackPreSelection(const reco::Track & itTrack){

  ecT_p        = itTrack.p();
  ecT_pErr     = pErrorOfTrack(itTrack);
  ecT_curv     = itTrack.qoverp();
  ecT_curvErr  = itTrack.qoverpError();
  ecT_pt       = itTrack.pt();
  ecT_ptErr    = itTrack.ptError();
  ecT_pz       = itTrack.pz();
  ecT_pzErr    = pzErrorOfTrack(itTrack);
  ecT_the      = itTrack.theta();
  ecT_theErr   = itTrack.thetaError();
  ecT_phi      = itTrack.phi();
  ecT_phiErr   = itTrack.phiError();
  ecT_eta      = itTrack.eta();     
  ecT_dxy      = itTrack.dxy();     
  ecT_dz       = itTrack.dz();    
  ecT_vx       = itTrack.vx();     
  ecT_vy       = itTrack.vy();     
  ecT_vz       = itTrack.vz();    
  ecT_q        = itTrack.charge();
  ecT_chi2     = itTrack.normalizedChi2();
  ecT_nHitVal  = itTrack.numberOfValidHits();
  ecT_curvt    = ecT_q/ecT_pt;
  ecT_curvtErr = ecT_ptErr*ecT_curvt*ecT_curvt;

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
  ecT_nHitSteTec = hp.getTrackerMonoStereo(6,1);
  ecT_nHitSteTid = hp.getTrackerMonoStereo(4,1);
  ecT_nHitSteTib = hp.getTrackerMonoStereo(3,1);
  ecT_nHitSteTob = hp.getTrackerMonoStereo(5,1);

  if ( myDebug_ ) std::cout << ">> This track: " << " pT:" << ecT_pt << " p:" << ecT_p << " eta:" << ecT_eta << " the:" << ecT_the << " phi:" << ecT_phi << " nValHits:" << ecT_nHitVal << " nTecHits:" << ecT_nHitTec << " nTecSteHits:" << ecT_nHitSteTec << std::endl;  

  //
  // Track preselection
  if (
      fabs(ecT_eta)<etaMinCut_ || 
      fabs(ecT_eta)>etaMaxCut_ || 
      ecT_nHitVal<nHitMinCut_ ||
      ecT_pt<ptMinCut_ ||
      ecT_pt>ptMaxCut_ ||
      ecT_nHitTib<nHitTibMinCut_ ||
      ecT_nHitSteTib<nHitSteTibMinCut_ ||
      ecT_nHitTob<nHitTobMinCut_ ||
      ecT_nHitSteTob<nHitSteTobMinCut_ ||
      ecT_nHitTec<nHitTecMinCut_ ||
      ecT_nHitSteTec<nHitSteTecMinCut_ ||
      ecT_nHitTid<nHitTidMinCut_ ||
      ecT_nHitSteTid<nHitSteTidMinCut_
      ) return false;

  if ( myDebug_ ) std::cout << ">> This track is selected. " << std::endl;  

  //
  // Inner State
  AlgebraicSymMatrix55 cov = itTrack.innerStateCovariance();
  ecT_pIs = itTrack.innerMomentum().R();
  ecT_curvIs = ecT_q/ecT_pIs;
  ecT_curvErrIs = sqrt(cov(0,0));
  ecT_pErrIs = ecT_pIs*ecT_curvErrIs/ecT_curvIs;
  ecT_ptIs = itTrack.innerMomentum().Rho();
  ecT_curvtIs = ecT_q/ecT_ptIs;   
  ecT_pzIs = itTrack.innerMomentum().z();      
  double covterm = 2.*ecT_ptIs*ecT_pzIs*cov(0,1)/ecT_pIs;
  double lamterm = cov(1,1)*ecT_curvIs*ecT_curvIs;
  ecT_ptErrIs = ecT_pIs*sqrt(ecT_ptIs*ecT_ptIs*cov(0,0)+ecT_pzIs*ecT_pzIs*lamterm+2.*covterm);   
  ecT_pzErrIs = ecT_pIs*sqrt(ecT_pzIs*ecT_pzIs*cov(0,0)+ecT_ptIs*ecT_ptIs*lamterm-2.*covterm);
  ecT_curvtErrIs = ecT_ptErrIs*ecT_curvtIs*ecT_curvtIs;
  ecT_theIs = itTrack.innerMomentum().theta();     
  ecT_theErrIs = sqrt(cov(1,1));  
  ecT_phiIs = itTrack.innerMomentum().phi();     
  ecT_phiErrIs = sqrt(cov(2,2));  

  //
  // Hit quantities
  ecT_nHitMis = hp.numberOfLostHits();
  ecT_nHitIna = hp.numberOfInactiveHits();
  ecT_nHitBad = hp.numberOfBadHits();

  return true;

}

void ecReco::trackAction(const reco::Track & itTrack, uint32_t & detidTl, uint32_t & detidIs){

  ecRefit thisEcRefit(theG.product(), theMF.product(), theFitter.product(), theBuilder.product(), myDebug_);
  tsosParams myPars = thisEcRefit.doWithTrack(itTrack, hs);

  //
  // Tree stuff
  ecT_iokTl =      myPars.iok;
  ecT_nHitTl =     myPars.nhit;
  ecT_pTl =        myPars.p;
  ecT_pErrTl =     myPars.pErr;
  ecT_curvTl =     myPars.curv;
  ecT_curvErrTl =  myPars.curvErr;
  ecT_curvtTl =    myPars.curvt;
  ecT_curvtErrTl = myPars.curvtErr;
  ecT_ptTl =       myPars.pt;
  ecT_ptErrTl =    myPars.ptErr;
  ecT_pzTl =       myPars.pz;
  ecT_pzErrTl =    myPars.pzErr;
  ecT_phiTl =      myPars.phi;
  ecT_phiErrTl =   myPars.phiErr;
  ecT_theTl =      myPars.theta;
  ecT_theErrTl =   myPars.thetaErr;

  //
  // return the DetId of the first tracklet hit and the inner state
  detidTl = myPars.detidTl;
  detidIs = myPars.detidIs;

  return;

}

double ecReco::pErrorOfTrack(const reco::Track & track){
  
  return track.charge()*track.p()*track.qoverpError()/track.qoverp();
  
}

double ecReco::pzErrorOfTrack(const reco::Track & track){
  
  //  AlgebraicSymMatrix66 errors = track.cartesianError().matrix();
  //  return sqrt(errors(5,5))

  double mag = track.p();
  double cur = track.qoverp();
  double per = track.pt();
  double lon = track.pz();
  double covterm = 2.*per*lon*track.covariance(0,1)/mag;
  double lamterm = track.covariance(1,1)*cur*cur;
  return mag*sqrt(lon*lon*track.covariance(0,0)+per*per*lamterm-2.*covterm);

}

#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

void ecReco::dumpModuleInfo(DetId hitId){
  
  if (hitId.subdetId() == StripSubdetector::TIB )  
    std::cout  << " TIB " << TIBDetId(hitId).layer();
  else if (hitId.subdetId() == StripSubdetector::TOB ) 
    std::cout  << " TOB " << TOBDetId(hitId).layer();
  else if (hitId.subdetId() == StripSubdetector::TEC ) 
    std::cout  << " TEC " << TECDetId(hitId).wheel() << " ring " << TECDetId(hitId).ring();
  else if (hitId.subdetId() == StripSubdetector::TID ) 
    std::cout  << " TID " << TIDDetId(hitId).wheel() << " ring " << TIDDetId(hitId).ring();
  else if (hitId.subdetId() == (int) PixelSubdetector::PixelBarrel ) 
    std::cout  << " PXB " << PXBDetId(hitId).layer();
  else if (hitId.subdetId() == (int) PixelSubdetector::PixelEndcap )
    std::cout  << " PXF " << PXFDetId(hitId).disk();
  else 
    std::cout  << " UKN ";

}



//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ecReco);
