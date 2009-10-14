// -*- C++ -*-
//
// Package:    NHitTest
// Class:      NHitTest
// 
/**\class NHitTest NHitTest.cc Tests/RefitTest/src/NHitTest.cc

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

class NHitTest : public edm::EDAnalyzer {
public:
  explicit NHitTest(const edm::ParameterSet&);
  ~NHitTest();
  
  
private:
  virtual void beginJob(const edm::EventSetup & ) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  double pErrorAtTSOS(TrajectoryStateOnSurface &);
  double ptErrorAtTSOS(TrajectoryStateOnSurface &);
  
  // ----------member data ---------------------------
  edm::InputTag tracksTag_; 
  edm::InputTag tpTag_; 
  edm::InputTag simtracksTag_; 
  bool myDebug_;

  TrackAssociatorBase * associatorByHits;

  TFile * file;
  TH1F * deltaPHist;
  TH1F * deltaPSimHist;
  TH1F * pullHist;
  TGraphErrors* graphForFit;
  int itrack;

  //NHitTest histograms
  std::vector<TH1D*> vecPtResH;
  std::vector<TH1D*> vecPtPullH;
  std::vector<TH1D*> vecPlusPtResH;
  std::vector<TH1D*> vecPlusPtPullH;
  std::vector<TH1D*> vecMinusPtResH;
  std::vector<TH1D*> vecMinusPtPullH;

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
NHitTest::NHitTest(const edm::ParameterSet& iConfig)
{
  //Tags
  tracksTag_ = iConfig.getUntrackedParameter<edm::InputTag>("tracks");
  tpTag_ = iConfig.getUntrackedParameter< edm::InputTag >("tp");
  simtracksTag_ = iConfig.getUntrackedParameter< edm::InputTag >("simtracks");

  //Miscellanea
  std::string FileName_ = iConfig.getUntrackedParameter<std::string>("outputFileName", "SplitTrack.root");
  myDebug_ = iConfig.getUntrackedParameter<bool>("myDebug", false);

  //now do what ever initialization is needed
  itrack=0;

  file = new TFile(FileName_.c_str(),"recreate");
  
  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);

  //
  // Summary histos
  deltaPHist = new TH1F("deltaP","deltaP",200,-0.5,0.5);
  deltaPSimHist = new TH1F("deltaPSim","deltaPSim",200,-0.5,0.5);
  pullHist = new TH1F("pull","pull",200,-20.,20.);

  for(int i=0; i<21; i++){

    char name[50];
    char sign[50];
    
    sprintf(name,"ptres%02d",i);
    TH1D * hist = new TH1D(name, name, 400, -1., 1.);
    vecPtResH.push_back(hist);
    sprintf(name,"ptpull%02d",i);
    hist = new TH1D(name, name, 200, -5., 5.);
    vecPtPullH.push_back(hist);
    
    for ( int isign = 0; isign<2; isign++) {
      sprintf(sign,"Minus");
      if ( isign ) sprintf(sign,"Plus");
      sprintf(name,"ptres%s%02d",sign,i);
      TH1D * hist = new TH1D(name, name, 400, -1., 1.);
      if ( isign ) {
	vecPlusPtResH.push_back(hist);
      } else {
	vecMinusPtResH.push_back(hist);
      }
      sprintf(name,"ptpull%s%02d",sign,i);
      hist = new TH1D(name, name, 200, -5., 5.);
      if ( isign ) {
	vecPlusPtPullH.push_back(hist);
      } else {
	vecMinusPtPullH.push_back(hist);
      }
      
    }
    
  }
  

  TH1::AddDirectory(oldAddDir);


}


NHitTest::~NHitTest()
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
NHitTest::analyze(const edm::Event& iEvent, const edm::EventSetup& setup)
{
   using namespace edm;
   using namespace std;

   using reco::Track;
   using reco::TrackCollection;

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
  
   Handle<SimTrackContainer> simTrackCollection;
   iEvent.getByLabel(simtracksTag_, simTrackCollection);
   const SimTrackContainer simTC = *(simTrackCollection.product());
  
   edm::Handle<TrackingParticleCollection>  TPCollectionH ;
   iEvent.getByLabel(tpTag_,TPCollectionH);
   const TrackingParticleCollection tPC   = *(TPCollectionH.product());

   reco::RecoToSimCollection RtSC = 
     associatorByHits->associateRecoToSim (trackCollectionH,TPCollectionH,&iEvent );

   //
   // Loop on event tracks
   for(View<Track>::size_type i=0; i<tC.size(); ++i) {
     RefToBase<Track> itTrack(trackCollectionH, i);
     
     if (fabs(itTrack->eta())>0.5 || itTrack->numberOfValidHits()<12) continue;

     double pt=itTrack->pt();
     double pterr=itTrack->ptError();
     int charge=itTrack->charge();
     

     //////////////////////////////////////////////////////////
     //
     // Get the sim track!
     double ptSim = -1.;
     try{ 
       std::vector<std::pair<TrackingParticleRef, double> > tp = RtSC[itTrack];
       if ( myDebug_ ) cout << "Reco Track pT: " << itTrack->pt() 
			    << " matched to " << tp.size() << " MC Tracks" << std::endl;
       for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	    it != tp.end(); ++it) {
	 TrackingParticleRef tpr = it->first;
	 double assocChi2 = it->second;
	 ptSim = tpr->pt();
	 if ( myDebug_ ) cout << "\t\tMCTrack " << tpr.index() << " pT: " << tpr->pt() << 
			   " NShared: " << assocChi2 << endl;
       }

       vecPtResH[0]->Fill((pt-ptSim)/ptSim);
       vecPtPullH[0]->Fill((pt-ptSim)/pterr);
       if ( charge > 0 ) {
	 vecPlusPtResH[0]->Fill((pt-ptSim)/ptSim);
	 vecPlusPtPullH[0]->Fill((pt-ptSim)/pterr);
       } else {
	 vecMinusPtResH[0]->Fill((pt-ptSim)/ptSim);
	 vecMinusPtPullH[0]->Fill((pt-ptSim)/pterr);
       }


     } catch (Exception event) {
       cout << "->   Track pT: " << itTrack->pt() 
	    <<  " matched to 0  MC Tracks" << endl;
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
     std::vector<double> vecDPerpT;
     if (trajVec.size()>0) {

       std::vector<TrajectoryMeasurement> TMeas=trajVec.begin()->measurements();
       std::vector<TrajectoryMeasurement>::iterator itm;
       
       int iHitN=0;
       int iEffHitN=0;
       double x0=(TMeas.end()-1)->updatedState().globalPosition().x();
       double y0=(TMeas.end()-1)->updatedState().globalPosition().y();
       double T = 0.;
       double pAveTrack=0.;
       double pErrAveTrack=0.;
       for (itm=TMeas.end()-1;itm!=TMeas.begin()-1;itm--){ //loop on hits
	 iHitN++;
	 TrajectoryStateOnSurface upTSOS = itm->updatedState();

	 double x=upTSOS.globalPosition().x();
	 double y=upTSOS.globalPosition().y();
	 double perpR=upTSOS.globalPosition().perp();
	 vecPerpR.push_back(perpR);
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
	 double dT=dPerpT/sin(theta);
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
       // Fit of split tracks
       //
       // Minimum number of effective hits
       int nMinHitSplit = 3;
       TransientTrackingRecHit::RecHitContainer hitsSplit;       
       int iEffFirst=1;
       int iEffLast=iEffFirst+nMinHitSplit-1;

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
	   double pt = firstSplitTSOS.globalMomentum().perp();
	   //	   double p = trajVecSplit.begin()->firstMeasurement().updatedState().globalMomentum().mag();
	   double pterr = ptErrorAtTSOS(firstSplitTSOS);
	   if ( myDebug_ ) cout << " iEffLast=" << iEffLast << " pt=" << pt << "+-" << pterr << endl;
	   if ( iEffLast<21 && ptSim>0. ) {
	     vecPtResH[iEffLast]->Fill((pt-ptSim)/ptSim);
	     vecPtPullH[iEffLast]->Fill((pt-ptSim)/pterr);
	     if ( charge > 0 ) {
	       vecPlusPtResH[iEffLast]->Fill((pt-ptSim)/ptSim);
	       vecPlusPtPullH[iEffLast]->Fill((pt-ptSim)/pterr);
	     } else {
	       vecMinusPtResH[iEffLast]->Fill((pt-ptSim)/ptSim);
	       vecMinusPtPullH[iEffLast]->Fill((pt-ptSim)/pterr);
	     }
	     
	   }
	   
	 } else {
	   if ( myDebug_ ) cout << " split fit failed" << endl;
	 }
	 
	 //	 iEffFirst++;
	 iEffLast++;
	 //
	 //Empty hits container
	 hitsSplit.clear();
	 
       }
       //
       //TGraphErrors
       itrack++;

     }
   }
}


// ------------ method called once each job just before starting event loop  ------------
void NHitTest::beginJob(const edm::EventSetup & setup)
{
  edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
  setup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theHitsAssociator);
  associatorByHits = (TrackAssociatorBase *) theHitsAssociator.product();
}

// ------------ method called once each job just after ending the event loop  ------------
void NHitTest::endJob() {

}


double NHitTest::pErrorAtTSOS(TrajectoryStateOnSurface & TSOS){
  
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

double NHitTest::ptErrorAtTSOS(TrajectoryStateOnSurface & TSOS){
  
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

//define this as a plug-in
DEFINE_FWK_MODULE(NHitTest);
