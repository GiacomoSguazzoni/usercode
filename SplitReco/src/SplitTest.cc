// -*- C++ -*-
//
// Package:    SplitTest
// Class:      SplitTest
// 
/**\class SplitTest SplitTest.cc Tests/RefitTest/src/SplitTest.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Giuseppe Cerati
//         Created:  Fri Aug  7 15:10:58 CEST 2009
// $Id$
//
//
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/SimTracker/TrackAssociation/test/testTrackAssociator.cc?revision=1.17&view=markup&pathrev=CMSSW_2_2_10
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/PatternTools/interface/TrajectoryFitter.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h" 
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include <TH1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TList.h>
#include <TF1.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <map>
#include <set>


//
// class declaration
//

class TrackAssociatorBase;

class SplitTest : public edm::EDAnalyzer {
public:
  explicit SplitTest(const edm::ParameterSet&);
  ~SplitTest();
  
  
private:
  virtual void beginJob(const edm::EventSetup & ) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  double pErrorOfTrack(const reco::Track &);
  double pErrorAtTSOS(TrajectoryStateOnSurface &);
  double ptErrorAtTSOS(TrajectoryStateOnSurface &);
  void setMaterialToKFactor(const TrackerGeometry *, TransientTrackingRecHit::RecHitContainer& , double);
  
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
  std::vector<unsigned int> goodruns_;

  TrackAssociatorBase * associatorByHits;

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
  int   sT_nHit;
  int   sT_nEffHit;
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
  // split track
  float sT_pLossSplit;
  float sT_p0;
  float sT_p1;
  float sT_chi2Split;
  int sT_freeParSplit;
  int sT_NSplit;
  float sT_TInSplit;
  float sT_TOutSplit;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SplitTest::SplitTest(const edm::ParameterSet& iConfig)
{
  //Miscellanea
  std::string FileName_ = iConfig.getUntrackedParameter<std::string>("outputFileName", "SplitTrack.root");
  myDebug_ = iConfig.getUntrackedParameter<bool>("myDebug", false);
  isMC_ = iConfig.getUntrackedParameter<bool>("isMC", true);
  splitTrackEffHits_ = iConfig.getUntrackedParameter<int>("splitTrackEffHits", 4);
  minSplits_ = iConfig.getUntrackedParameter<int>("minSplits", 4);
  kFactor_ = iConfig.getUntrackedParameter<double>("kFactor", 1.);
  pullCut_ = iConfig.getUntrackedParameter<double>("pullCut", 7.);
  std::vector<unsigned int> empty_goodruns;
  goodruns_ = iConfig.getUntrackedParameter<std::vector<unsigned int> >("goodruns", empty_goodruns);
  
  //Tags
  if ( isMC_ ) {
    tracksTag_ = iConfig.getUntrackedParameter<edm::InputTag>("tracks");
    tpTag_ = iConfig.getUntrackedParameter< edm::InputTag >("tp");
    simtracksTag_ = iConfig.getUntrackedParameter< edm::InputTag >("simtracks");
  } else {
    tracksTag_ = iConfig.getUntrackedParameter<edm::InputTag>("cosmics");
  }
  
  //now do what ever initialization is needed
  sT_iEvent=0;
  itrack=0;

  file = new TFile(FileName_.c_str(),"recreate");
  
  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);

  //
  // Summary histos
  deltaPHist = new TH1F("deltaP","deltaP",200,-0.5,0.5);
  deltaPSimHist = new TH1F("deltaPSim","deltaPSim",200,-0.5,0.5);
  pullHist = new TH1F("pull","pull",200,-20.,20.);
  TH1::AddDirectory(oldAddDir);

  //
  // Root-upla
  int bs = 64000;

  splitTree = new TTree("sT","sT");
  // event/track
  splitTree->Branch("nRun",       &sT_nRun,       "nRun/I",       bs);  
  splitTree->Branch("nEvent",     &sT_nEvent,     "nEvent/I",     bs);  
  splitTree->Branch("iEvent",     &sT_iEvent,     "iEvent/I",     bs);  
  splitTree->Branch("NTrack",     &sT_NTrack,     "NTrack/I",     bs);  
  splitTree->Branch("iTrack",     &sT_iTrack,     "iTrack/I",     bs);  
  // track	     		   		   	       
  splitTree->Branch("p",          &sT_p,          "p/F",          bs);       
  splitTree->Branch("pErr",       &sT_pErr,       "pErr/F",       bs);       
  splitTree->Branch("pAve",       &sT_pAve,       "pAve/F",       bs);       
  splitTree->Branch("pAveErr",    &sT_pAveErr,    "pAveErr/F",    bs);       
  splitTree->Branch("pt",         &sT_pt,         "pt/F",         bs);      
  splitTree->Branch("ptErr",      &sT_ptErr,      "ptErr/F",      bs);       
  splitTree->Branch("eta",        &sT_eta,        "eta/F",        bs);     
  splitTree->Branch("the",        &sT_the,        "the/F",        bs);     
  splitTree->Branch("phi",        &sT_phi,        "phi/F",        bs);     
  splitTree->Branch("nHit",       &sT_nHit,       "nHit/I",       bs);    
  splitTree->Branch("nEffHit",    &sT_nEffHit,    "nEffHit/I",    bs);    
  splitTree->Branch("chi2",       &sT_chi2,       "chi2/F",       bs);    
  splitTree->Branch("T",          &sT_T,          "T/F",          bs);       
  splitTree->Branch("rIn",        &sT_rIn,        "rIn/F",        bs); //r first hit
  splitTree->Branch("rOut",       &sT_rOut,       "rOut/F",       bs); //r last hit
  splitTree->Branch("zIn",        &sT_zIn,        "zIn/F",        bs); //z first hit
  splitTree->Branch("zOut",       &sT_zOut,       "zOut/F",       bs); //z last hit
  // simtrack	     		   		   	       
  splitTree->Branch("iTrackSim",  &sT_iTrackSim,  "iTrackSim/I",  bs);    
  splitTree->Branch("pSim",       &sT_pSim,       "pSim/F",       bs);    
  splitTree->Branch("ptSim",      &sT_ptSim,      "ptSim/F",      bs);   
  splitTree->Branch("etaSim",     &sT_etaSim,     "etaSim/F",     bs);  
  splitTree->Branch("theSim",     &sT_theSim,     "theSim/F",     bs);  
  splitTree->Branch("phiSim",     &sT_phiSim,     "phiSim/F",     bs);  
  splitTree->Branch("pLossSim",   &sT_pLossSim,   "pLossSim/F",   bs);
  splitTree->Branch("nHitSim",    &sT_nHitSim,    "nHitSim/I",    bs);    
  splitTree->Branch("TSim",       &sT_TSim,       "TSim/F",       bs); //r first hit
  splitTree->Branch("rInSim",     &sT_rInSim,     "rInSim/F",     bs); //r first hit
  splitTree->Branch("rOutSim",    &sT_rOutSim,    "rOutSim/F",    bs); //r last hit
  splitTree->Branch("zInSim",     &sT_zInSim,     "zInSim/F",     bs); //z first hit
  splitTree->Branch("zOutSim",    &sT_zOutSim,    "zOutSim/F",    bs); //z last hit
  // split track     		   		   	       
  splitTree->Branch("pLossSplit", &sT_pLossSplit, "pLossSplit/F", bs);
  splitTree->Branch("p0Split",    &sT_p0,         "p0Split/F", bs);
  splitTree->Branch("p1Split",    &sT_p1,         "p1Split/F", bs);
  splitTree->Branch("chi2Split",  &sT_chi2Split,  "chi2Split/F",  bs);
  splitTree->Branch("freeParSplit", &sT_freeParSplit,     "freeParSplit/I",     bs);  
  splitTree->Branch("NSplit",     &sT_NSplit,     "NSplit/I",     bs);  
  splitTree->Branch("TInSplit",   &sT_TInSplit,   "TInSplit/F",   bs);
  splitTree->Branch("TOutSplit",  &sT_TOutSplit,  "TOutSplit/F",  bs);
}


SplitTest::~SplitTest()
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
SplitTest::analyze(const edm::Event& iEvent, const edm::EventSetup& setup)
{
   using namespace edm;
   using namespace std;

   using reco::Track;
   using reco::TrackCollection;

   // Store run number
   //
   sT_nRun = iEvent.run();

   // Not a good run? Skipping...
   //
   if (goodruns_.size()>0) {
     if (find(goodruns_.begin(),goodruns_.end(),(unsigned int)sT_nRun)==goodruns_.end())
       {
	 cout << sT_nRun << " not a good run, skipping..."<< endl ;
	 return;
       }
   }

   std::string fitterName = "KFFittingSmootherWithOutliersRejectionAndRK";   
   std::string propagatorName = "RungeKuttaTrackerPropagator";   
   std::string builderName = "WithAngleAndTemplate";   

   edm::ESHandle<TrackerGeometry> theG;
   edm::ESHandle<MagneticField> theMF;
   edm::ESHandle<TrajectoryFitter> theFitter;
   edm::ESHandle<Propagator> thePropagator;
   edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
   setup.get<TrackerDigiGeometryRecord>().get(theG);
   setup.get<IdealMagneticFieldRecord>().get(theMF);  
   setup.get<TrackingComponentsRecord>().get(fitterName,theFitter);
   setup.get<TrackingComponentsRecord>().get(propagatorName,thePropagator);
   setup.get<TransientRecHitRecord>().get(builderName,theBuilder);
 
   Handle<View<Track> > trackCollectionH;
   iEvent.getByLabel(tracksTag_,trackCollectionH);
   const View<Track>  tC = *(trackCollectionH.product()); 
  
   reco::RecoToSimCollection RtSC;
   if ( isMC_ ) {
     Handle<SimTrackContainer> simTrackCollection;
     iEvent.getByLabel(simtracksTag_, simTrackCollection);
     const SimTrackContainer simTC = *(simTrackCollection.product());
     edm::Handle<TrackingParticleCollection>  TPCollectionH;
     iEvent.getByLabel(tpTag_,TPCollectionH);
     const TrackingParticleCollection tPC = *(TPCollectionH.product());
     RtSC = associatorByHits->associateRecoToSim (trackCollectionH,TPCollectionH,&iEvent);
   }


   sT_nEvent = iEvent.id().event();
   sT_iEvent++;

   //
   // Loop on event tracks
   sT_iTrack=0;
   sT_NTrack = tC.size();
   for(View<Track>::size_type i=0; i<tC.size(); ++i) {
     RefToBase<Track> itTrack(trackCollectionH, i);

     sT_iTrack++;
     
     sT_p    = itTrack->p();
     sT_pErr = pErrorOfTrack(*itTrack);
     sT_pt   = itTrack->pt();
     sT_ptErr  = itTrack->ptError();
     sT_eta  = itTrack->eta();
     sT_the  = itTrack->theta();
     sT_phi  = itTrack->phi();
     sT_chi2 = itTrack->normalizedChi2();
     sT_nHit = itTrack->numberOfValidHits();

     // For now no cuts
     //   if (fabs(itTrack->eta())>1. || itTrack->numberOfValidHits()<12) continue;
     if (fabs(itTrack->eta())>2.2 || itTrack->numberOfValidHits()<12) continue;
     
     //////////////////////////////////////////////////////////
     //
     // Before anything else... get the sim track and study it!
     if ( isMC_ ){
       try{ 
	 std::vector<std::pair<TrackingParticleRef, double> > tp = RtSC[itTrack];
	 if ( myDebug_ ) cout << "Reco Track pT: " << itTrack->pt() 
			      << " matched to " << tp.size() << " MC Tracks" << std::endl;
	 for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	      it != tp.end(); ++it) {
	   TrackingParticleRef tpr = it->first;
	   TrackingParticle*   tp=const_cast<TrackingParticle*>(tpr.get());
	   double assocChi2 = it->second;
	   if ( myDebug_ ) cout << "\t\tMCTrack " << tpr.index() << " pT: " << tpr->pt() << 
			     " NShared: " << assocChi2 << endl;
	   
	   // If more than one associated sim track only the latter enters in the ntupla
	   sT_iTrackSim=tpr.index(); 
	   sT_pSim=tpr->p();      
	   sT_ptSim=tpr->pt();     
	   sT_etaSim=tpr->eta();    
	   sT_theSim=tpr->theta();    
	   sT_phiSim=tpr->phi();    
	   //
	   // PLoss of the associated track
	   // We assume only one simtrack per track //FIXME
	   int itphit=0; 
	   double pSimOut=0.; 
	   double x0=0.;
	   double y0=0.;
	   double z0=0.;
	   double T = 0.;
	   std::vector<double> vecTPHitPerpR;
	   std::vector<double> vecTPHitZ;
	   std::vector<double> vecTPHitT;
	   for(std::vector<PSimHit>::const_iterator TPhit = tp->pSimHit_begin(); TPhit != tp->pSimHit_end(); TPhit++){ //sguazz
	     //
	     // get the DetUnit via the DetUnitId and cast it to a StripGeomDetUnit
	     //
	     // Example code from:
	     // http://cmslxr.fnal.gov/lxr/source/Validation/TrackerHits/src/TrackerHitAnalyzer.cc
	     DetId detid=DetId(TPhit->detUnitId());
	     if (detid.det() == DetId::Tracker){ 
	       const GeomDetUnit * det=(const GeomDetUnit*)theG->idToDetUnit( detid );
	       if ( det ) {
		 itphit++;
		 GlobalPoint gpos=det->toGlobal(TPhit->localPosition());
		 
		 double x = gpos.x();
		 double y = gpos.y();
		 double z = gpos.z();
		 
		 double dT=0.;
		 if ( ! x0==0. ) dT = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0));
		 T+=dT;
		 
		 vecTPHitT.push_back(T);
		 vecTPHitPerpR.push_back(gpos.perp());
		 vecTPHitZ.push_back(z);
		 if ( myDebug_ ) cout << itphit << " p=" << TPhit->pabs() << " r=" << gpos.perp() << endl;
		 if ( TPhit->pabs() > 0.5 ) { pSimOut=TPhit->pabs() ;}
		 
		 x0=x;
		 y0=y;
		 z0=z;
		 
	       } else {
		 if ( myDebug_ ) std::cout << " Problem in sim track " << std::endl;
	       }
	     }
	   }
	   sT_pLossSim = tp->p() - pSimOut;
	   sT_nHitSim = itphit;
	   sT_TSim = vecTPHitT.at(itphit-1);
	   sT_rInSim = vecTPHitPerpR.at(0);
	   sT_rOutSim = vecTPHitPerpR.at(itphit-1);
	   sT_zInSim = vecTPHitZ.at(0);
	   sT_zOutSim = vecTPHitZ.at(itphit-1);
	   deltaPSimHist->Fill(sT_pLossSim);
	   if ( myDebug_ ) cout << " SimTrack p:" << tp->p() << " eta:" << tp->momentum().eta() << " pout:" << pSimOut << " pLoss:" << sT_pLossSim << endl;
	   //
	 }
       } catch (Exception event) {
	 cout << "->   Track pT: " << itTrack->pt() 
	      <<  " matched to 0  MC Tracks" << endl;
	 sT_iTrackSim=0;
	 sT_pSim=0.;     
	 sT_ptSim=0.;    
	 sT_etaSim=0.;   
	 sT_theSim=0.;   
	 sT_phiSim=0.;   
	 sT_pLossSim=0.; 
       }
     }
     // End of simtrack section
     //
     /////////////////////////////////////////////////////////

     TrajectoryStateTransform transformer;
     TrajectoryStateOnSurface innerStateFromTrack=transformer.innerStateOnSurface(*itTrack,*theG,theMF.product());
     TrajectoryStateOnSurface outerStateFromTrack=transformer.outerStateOnSurface(*itTrack,*theG,theMF.product());

     //First refit the track with all the hits to have all intermediate TSOS
     //
     // Direction
     PropagationDirection seedDir = itTrack->seedDirection();
     // 
     // Hits
     TransientTrackingRecHit::RecHitContainer hitsAll;       
     for (trackingRecHit_iterator i=itTrack->recHitsBegin(); i!=itTrack->recHitsEnd(); i++){
       hitsAll.push_back(theBuilder->build(&**i ));
     }
     //
     // Initial state
     TrajectoryStateOnSurface theInitialStateForRefitting;
     TrajectoryStateOnSurface initialStateFromTrack = 
       ( (innerStateFromTrack.globalPosition()-hitsAll.front()->det()->position()).mag2() <
	 (outerStateFromTrack.globalPosition()-hitsAll.front()->det()->position()).mag2() ) ? 
       innerStateFromTrack: outerStateFromTrack;       
     initialStateFromTrack.rescaleError(100);
     theInitialStateForRefitting = TrajectoryStateOnSurface(initialStateFromTrack.localParameters(),
							      initialStateFromTrack.localError(),                  
							      initialStateFromTrack.surface(),
							      theMF.product());  
     //
     // Seed
     const TrajectorySeed seed = TrajectorySeed(PTrajectoryStateOnDet(), TrajectorySeed::recHitContainer(), seedDir);
     //
     // Fit!
     std::vector<Trajectory> trajVec = theFitter->fit(seed, hitsAll, theInitialStateForRefitting);

     //
     // Store quantities needed for split refit
     std::vector<double> vecT;
     std::vector<int> vecEffHitN;
     std::vector<int> vecEffHitBegin;
     std::vector<int> vecEffHitEnd;
     std::vector<double> vecPerpR;
     std::vector<double> vecZ;
     std::vector<double> vecDPerpT;
     if (trajVec.size()>0) {

       std::vector<TrajectoryMeasurement> TMeas=trajVec.begin()->measurements();
       std::vector<TrajectoryMeasurement>::iterator itm;
       
       int iHitN=0;
       int iEffHitN=0;
       double x0=(TMeas.end()-1)->updatedState().globalPosition().x();
       double y0=(TMeas.end()-1)->updatedState().globalPosition().y();
       double z0=(TMeas.end()-1)->updatedState().globalPosition().z();
       double T = 0.;
       double pAveTrack=0.;
       double pErrAveTrack=0.;
       for (itm=TMeas.end()-1;itm!=TMeas.begin()-1;itm--){ //loop on hits
	 iHitN++;
	 TrajectoryStateOnSurface upTSOS = itm->updatedState();

	 double x=upTSOS.globalPosition().x();
	 double y=upTSOS.globalPosition().y();
	 double z=upTSOS.globalPosition().z();
	 double perpR=upTSOS.globalPosition().perp();
	 vecPerpR.push_back(perpR);
	 vecZ.push_back(z);
	 pAveTrack+=upTSOS.globalMomentum().mag();
	 pErrAveTrack+=pErrorAtTSOS(upTSOS);
	 double theta=upTSOS.globalMomentum().theta();
	 double dPerpT=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
	 vecDPerpT.push_back(dPerpT);
	 if ( dPerpT > 1. || iHitN == 1 ) {
	   iEffHitN++;
	   vecEffHitBegin.push_back(iHitN);
	   vecEffHitEnd.push_back(iHitN);
	 } else {
	   vecEffHitEnd.at(iEffHitN-1)=iHitN;
	 }
	 vecEffHitN.push_back(iEffHitN);
	 double dT=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0));
	 //	 double dT=dPerpT/sin(theta);
	 T+=dT;
	 vecT.push_back(T);
	 if ( myDebug_ ) cout
	   << " >>>>TSOS>>" << iHitN
	   << " ups=" << upTSOS.globalMomentum().mag() 
	   << " th=" << theta
	   << " sinth=" << sin(theta)
	   << " r=" << perpR
	   << " x=" << x
	   << " y=" << y
	   << " z=" << z
	   << " dx=" << x-x0
	   << " dy=" << y-y0
	   << " dPerpT=" << dPerpT
	   << " effN=" << iEffHitN
	   << " beg:" << vecEffHitBegin.at(iEffHitN-1)
	   << " end:" << vecEffHitEnd.at(iEffHitN-1)
	   << " dT=" << dT
	   << " T=" << T
	   << endl;
	 x0=x;
	 y0=y;
	 z0=z;
       }
       // track average p
       // track average pError
       pAveTrack=pAveTrack/iHitN;
       pErrAveTrack=pErrAveTrack/iHitN;
       //
       // Stuff needed below
       int effHitN=iEffHitN;
       itm=TMeas.end();
       //
       // Tree stuff
       sT_pAve=pAveTrack;
       sT_pAveErr=pErrAveTrack;
       sT_nEffHit=effHitN;
       sT_T=T;
       sT_rIn=vecPerpR.at(0);
       sT_rOut=vecPerpR.at(iHitN-1);
       sT_zIn=vecZ.at(0);
       sT_zOut=vecZ.at(iHitN-1);

       //
       // Change Material!
       double kFactor=kFactor_;
       setMaterialToKFactor(theG.product(), hitsAll, kFactor);

       //
       // Fit of split tracks
       //
       // Minimum number of effective hits
       int nMinHitSplit = splitTrackEffHits_;
       TransientTrackingRecHit::RecHitContainer hitsSplit;       
       int iEffFirst=1;
       int iEffLast=iEffFirst+nMinHitSplit-1;
       //
       // Arrays for Split Track Fit Results
       double arrPSplit[50];
       double arrPErrSplit[50];
       double arrTSplit[50];
       double arrTErrSplit[50];
       int iarr=0;
       while ( iEffLast<effHitN+1 ) {
	 if ( myDebug_ ) cout << " NEXT iEffFirst=" << iEffFirst << " " << " iEffLast=" << iEffLast << " " << endl;
	 if ( myDebug_ ) cout << " Now using hits"; 
	 int iFirst=vecEffHitBegin.at(iEffFirst-1);
	 int iCurrent=iFirst;
	 //
	 // Build the rec hit container for the Split Track
	 while ( iCurrent<vecEffHitEnd.at(iEffLast-1)+1 ) { 
	   // Address first hit itTrack->recHitsBegin()
	   trackingRecHit_iterator i=itTrack->recHitsBegin()+iCurrent-1;
	   if ( myDebug_ ) cout << " #" << iCurrent << " dPerpT=" << vecDPerpT.at(iCurrent-1);
	   iCurrent++;
	   hitsSplit.push_back(theBuilder->build(&**i ));
	   i++;
	 }
	 //
	 // Set the proper TSOS
	 TrajectoryStateOnSurface theInitialStateForSplitRefitting;
	 if ( iFirst == 1 ) {
	   theInitialStateForSplitRefitting=theInitialStateForRefitting;
	   if ( myDebug_ ) cout << " iFirst=" << iFirst << " TSOS=" << theInitialStateForSplitRefitting.globalMomentum().mag(); 
	 } else {
	   TrajectoryStateOnSurface theStateFromTraj = (itm-iFirst)->forwardPredictedState();
	   theStateFromTraj.rescaleError(100);
	   theInitialStateForSplitRefitting = TrajectoryStateOnSurface(theStateFromTraj.localParameters(),
								       theStateFromTraj.localError(),                  
								       theStateFromTraj.surface(),
								       theMF.product());  
	   if ( myDebug_ ) cout << " iFirst=" << iFirst << " TSOS=" << theInitialStateForSplitRefitting.globalMomentum().mag(); 
	 }
	 
	 //
	 // Fit!
	 std::vector<Trajectory> trajVecSplit = theFitter->fit(seed, hitsSplit, theInitialStateForSplitRefitting);
	 if (trajVecSplit.size()>0) {
	   TrajectoryStateOnSurface firstSplitTSOS = trajVecSplit.begin()->lastMeasurement().updatedState();
	   double p = firstSplitTSOS.globalMomentum().mag();
	   //	   double p = trajVecSplit.begin()->firstMeasurement().updatedState().globalMomentum().mag();
	   double perror = pErrorAtTSOS(firstSplitTSOS);
	   if ( myDebug_ ) cout << " p=" << p << "+-" << perror << endl;
	   //
	   //Loop on TSOS of the Split track	   
	   std::vector<TrajectoryMeasurement> TMeas=trajVecSplit.begin()->measurements();
	   std::vector<TrajectoryMeasurement>::iterator itm;
	   int ista=0;
	   double pAveSplit=0.;
	   double pErrAveSplit=0.;
	   double TAveSplit=0.;
	   int nMeas=TMeas.size();
	   for (itm=TMeas.end()-1;itm!=TMeas.begin()-1;itm--){ //loop on hits
	     ista++;
	     TrajectoryStateOnSurface uptsos=itm->updatedState();
	     if ( myDebug_ ) cout << ">>>>SplitTrackTSOS>> " << ista << " p=" << uptsos.globalMomentum().mag() << "+-" << pErrorAtTSOS(uptsos) << " r=" << uptsos.globalPosition().perp() << " T=" << vecT.at(iFirst-1+ista-1) << endl;   
	     //
	     // For now: 
	     // p of the split = plain average 
	     // perr  of the split = plain average of the errors
	     // T of the split = plain average 
	     //
	     pAveSplit += uptsos.globalMomentum().mag()/nMeas; //FIXME now ok... but then use perp and theta? Exclude 1st and last hit?
	     pErrAveSplit += pErrorAtTSOS(uptsos)/nMeas;
	     TAveSplit = vecT.at(iFirst-1+ista-1);
	   }  
	   //
	   // Trajectory lenght average
	   TAveSplit = 0.5*(TAveSplit+vecT.at(iFirst-1));
	   if ( myDebug_ ) cout << ">>>>SplitTrackAve>> pAveSplit="<< pAveSplit << "+-" << pErrAveSplit << " TAveSplit=" << TAveSplit << endl; 
	   //
	   // Pull
	   double pull=(pAveSplit - pAveTrack)/sqrt(pErrAveSplit*pErrAveSplit+pErrAveTrack*pErrAveTrack);
	   //
	   // Fill split-wise summary histos
	   pullHist->Fill(pull);
	   //
	   // For the energy loss fit...
	   if ( abs(pull)<pullCut_ ){ 
	     arrPSplit[iarr]=pAveSplit;
	     arrPErrSplit[iarr]=pErrAveSplit;
	     arrTSplit[iarr]=TAveSplit;
	     arrTErrSplit[iarr]=0.05; //half millimiter error (?)
	     iarr++;
	   }
	 } else {
	   if ( myDebug_ ) cout << " split fit failed" << endl;
	 }
	 
	 iEffFirst++;
	 iEffLast++;
	 //
	 //Empty hits container
	 hitsSplit.clear();
	 //
	 // Put the material back to original values!
	 kFactor=1./kFactor;
	 setMaterialToKFactor(theG.product(), hitsAll, kFactor);
	 
       }
       //
       // def value in case of failure
       sT_pLossSplit=-1.;
       sT_p0=-1.;
       sT_p1=-1.;
       sT_chi2Split=-1.; 
       sT_NSplit=-1;    
       sT_freeParSplit=-1; 
       sT_TInSplit=-1.;  
       sT_TOutSplit=-1.; 
       //
       if ( iarr > minSplits_-1 ) {
	 graphForFit = new TGraphErrors(iarr,arrTSplit,arrPSplit,arrTErrSplit,arrPErrSplit);  
	 int fitResult = graphForFit->Fit("pol1");
	 TF1 *fit = graphForFit->GetFunction("pol1");
	 if ( fitResult ) {
	   if ( myDebug_ ) cout << " Linear fit of splits failed with code: " << fitResult << endl;
	 } else {
	   //
	   // Count tracks with valid fit
	   itrack++;
	   double p0 = fit->GetParameter(0);
	   double p1 = fit->GetParameter(1);
	   double chi2 = fit->GetChisquare();
	   int freePar = fit->GetNumberFreeParameters();
	   if ( myDebug_ ) cout << " linear EL fit p0=" << p0 << " p1=" << p1 << " chi2/freePar=" << chi2 << "/" << freePar << endl; 
	   double TIn  = *(vecT.begin());
	   double TOut = *(vecT.end()-1);
	   double Ttot= TOut - TIn; 
	   double pIn  = p0+p1*TIn;
	   double pOut = p0+p1*TOut;
	   double pLoss = pIn - pOut;
	   cout << " TIn=" << TIn << " pIn=" << pIn << "; TOut=" << TOut << " pOut=" << pOut << " Ttot=" << Ttot << " pLoss=" << pLoss << endl; 
	   
	   //
	   // Fill track-wise summary histos
	   deltaPHist->Fill(pLoss);
	   
	   //
	   // root-ple stuff
	   sT_pLossSplit=pLoss;
	   sT_p0=p0;
	   sT_p1=p1;
	   sT_pLossSplit=pLoss;
	   sT_chi2Split=chi2; 
	   sT_freeParSplit=freePar; 
	   sT_NSplit=iarr;    
	   sT_TInSplit=TIn;  
	   sT_TOutSplit=TOut; 
	   //
	   // Few fit examples for valid fits
	   if ( itrack < 21 ) {
	     char grname[50];
	     sprintf(grname,"track%02d",itrack);
	     graphForFit->SetName(grname);
	     file->cd();
	     graphForFit->Write();
	   }
	 }
       }

       //
       // Fill splitTree
       splitTree->Fill();
       
     }
   }
}


// ------------ method called once each job just before starting event loop  ------------
void SplitTest::beginJob(const edm::EventSetup & setup)
{
  edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
  setup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theHitsAssociator);
  associatorByHits = (TrackAssociatorBase *) theHitsAssociator.product();
}

// ------------ method called once each job just after ending the event loop  ------------
void SplitTest::endJob() {

}


double SplitTest::pErrorOfTrack(const reco::Track & track){
  
  //  float     px    = TSOS.globalMomentum().x();
  //  float     py    = TSOS.globalMomentum().y();
  //  float     pz    = TSOS.globalMomentum().z();
  //  float     pt    = TSOS.globalMomentum().perp();
  //  float     phi   = TSOS.globalMomentum().phi();
  //  float     theta = TSOS.globalMomentum().theta();
  //  float     eta   = TSOS.globalMomentum().eta();
  
  //get the error of the kinimatic parameters
  //  AlgebraicSymMatrix66 errors = TSOS.cartesianError().matrix();
  double partialPterror = track.covariance(3,3)*pow(track.px(),2) + track.covariance(4,4)*pow(track.py(),2);
  //  float     pterror  = sqrt(partialPterror)/TSOS.globalMomentum().perp();
  //  float     pxerror  = sqrt(errors(3,3))/TSOS.globalMomentum().x();
  //  float     pyerror  = sqrt(errors(4,4))/TSOS.globalMomentum().y();
  //  float     pzerror  = sqrt(errors(5,5))/TSOS.globalMomentum().z();
  double     perror   = sqrt(partialPterror+track.covariance(5,5)*pow(track.pz(),2))/track.p();
  //  float     phierror = sqrt(TSOS.curvilinearError().matrix()(2,2));
  //  float     etaerror = sqrt(TSOS.curvilinearError().matrix()(1,1))*fabs(sin(TSOS.globalMomentum().theta()));
  
  return perror;
  
}

double SplitTest::pErrorAtTSOS(TrajectoryStateOnSurface & TSOS){
  
  //  float     px    = TSOS.globalMomentum().x();
  //  float     py    = TSOS.globalMomentum().y();
  //  float     pz    = TSOS.globalMomentum().z();
  //  float     pt    = TSOS.globalMomentum().perp();
  //  float     phi   = TSOS.globalMomentum().phi();
  //  float     theta = TSOS.globalMomentum().theta();
  //  float     eta   = TSOS.globalMomentum().eta();
  
  //get the error of the kinimatic parameters
  AlgebraicSymMatrix66 errors = TSOS.cartesianError().matrix();
  double partialPterror = errors(3,3)*pow(TSOS.globalMomentum().x(),2) + errors(4,4)*pow(TSOS.globalMomentum().y(),2);
  //  float     pterror  = sqrt(partialPterror)/TSOS.globalMomentum().perp();
  //  float     pxerror  = sqrt(errors(3,3))/TSOS.globalMomentum().x();
  //  float     pyerror  = sqrt(errors(4,4))/TSOS.globalMomentum().y();
  //  float     pzerror  = sqrt(errors(5,5))/TSOS.globalMomentum().z();
  double     perror   = sqrt(partialPterror+errors(5,5)*pow(TSOS.globalMomentum().z(),2))/TSOS.globalMomentum().mag();
  //  float     phierror = sqrt(TSOS.curvilinearError().matrix()(2,2));
  //  float     etaerror = sqrt(TSOS.curvilinearError().matrix()(1,1))*fabs(sin(TSOS.globalMomentum().theta()));
  
  return perror;
  
}

double SplitTest::ptErrorAtTSOS(TrajectoryStateOnSurface & TSOS){
  
  //  float     px    = TSOS.globalMomentum().x();
  //  float     py    = TSOS.globalMomentum().y();
  //  float     pz    = TSOS.globalMomentum().z();
  //  float     pt    = TSOS.globalMomentum().perp();
  //  float     phi   = TSOS.globalMomentum().phi();
  //  float     theta = TSOS.globalMomentum().theta();
  //  float     eta   = TSOS.globalMomentum().eta();
  
  //get the error of the kinimatic parameters
  AlgebraicSymMatrix66 errors = TSOS.cartesianError().matrix();
  double partialPterror = errors(3,3)*pow(TSOS.globalMomentum().x(),2) + errors(4,4)*pow(TSOS.globalMomentum().y(),2);
  double     pterror  = sqrt(partialPterror)/TSOS.globalMomentum().perp();
  //  float     pxerror  = sqrt(errors(3,3))/TSOS.globalMomentum().x();
  //  float     pyerror  = sqrt(errors(4,4))/TSOS.globalMomentum().y();
  //  float     pzerror  = sqrt(errors(5,5))/TSOS.globalMomentum().z();
  //  double     perror   = sqrt(partialPterror+errors(5,5)*pow(TSOS.globalMomentum().z(),2))/TSOS.globalMomentum().mag();
  //  float     phierror = sqrt(TSOS.curvilinearError().matrix()(2,2));
  //  float     etaerror = sqrt(TSOS.curvilinearError().matrix()(1,1))*fabs(sin(TSOS.globalMomentum().theta()));
  
  return pterror;
  
}

void SplitTest::setMaterialToKFactor(const TrackerGeometry * theGeometry, TransientTrackingRecHit::RecHitContainer& hits, double kFactor){
  
  std::string outputBlock;
  char buffer [300];
  int nhit=0;

  for (TransientTrackingRecHit::RecHitContainer::const_iterator it=hits.begin(); it!=hits.end();it++){
    nhit++;

    const GeomDet* geomDet = theGeometry->idToDet((*it)->geographicalId());
      if ( geomDet ) {
	double theXi     = geomDet->surface().mediumProperties()->xi();
	double theRadLen = geomDet->surface().mediumProperties()->radLen();
	const MediumProperties theNewMP = MediumProperties(theRadLen, theXi*kFactor);
	BoundPlane* mySurface = (BoundPlane*)&(geomDet->surface());
	mySurface->setMediumProperties(theNewMP); 
	double theNewXi     = geomDet->surface().mediumProperties()->xi();
	double theNewRadLen = geomDet->surface().mediumProperties()->radLen();
	sprintf(buffer, ">>>>setMaterialToKFactor>> nhit= %d kFactor= %f xiOld= %f xiNew= %f radLenOld= %f radLenNew= %f \n", nhit, kFactor, theXi, theNewXi, theRadLen, theNewRadLen);
	if ( myDebug_ ) std::cout << buffer;
      }
  }
  
  //  edm::LogInfo("ELFTrackProducer::setMaterialToKFactor") << outputBlock.data();
}
  


//define this as a plug-in
DEFINE_FWK_MODULE(SplitTest);
