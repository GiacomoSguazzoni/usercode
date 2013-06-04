//
// Original Author:  Giuseppe Cerati
//         Created:  Fri Aug  7 15:10:58 CEST 2009
// $Id: ecRefit.cc,v 1.3 2013/05/14 10:14:17 sguazz Exp $
//
//
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/SimTracker/TrackAssociation/test/testTrackAssociator.cc?revision=1.17&view=markup&pathrev=CMSSW_2_2_10
//

#include "Tests/ecTest/interface/ecRefit.h"
#include "Tests/ecTest/interface/ec.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/InvalidTransientRecHit.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h" 
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"

#include <TF1.h>

#include <memory>
#include <iostream>
#include <string>
#include <map>
#include <set>

#define cmsswREMOVEFOR44Version44x


//
// constructors and destructor
//
ecRefit::ecRefit(const TrackerGeometry * theGeometry, const MagneticField * theMagneticField, const TrajectoryFitter * theFitter, const TransientTrackingRecHitBuilder * theBuilder, bool debug)
{

  myDebug_ = debug;

  theG = theGeometry;
  theMF = theMagneticField;
  theF = theFitter;
  theB = theBuilder;

}

int ecRefit::buildHitsVector(const reco::Track track, hitSelector hs, uint32_t& detidIs, uint32_t& detidTl){

  ////////////////////////////////////////////////////////////////////////////////

  if ( myDebug_ ) std::cout << "   Selecting hits for refit" << std::endl;
  if ( myDebug_ ) std::cout << "   Selector: tib:" << hs.tib << " tid:" << hs.tid << " tob:" << hs.tob << " tec:" << hs.tec << " stereo:" << hs.stereo << " all:" << hs.all << std::endl;

  // 
  // Hits
  int iFirstIs = 1;
  int iFirstTl = 1;
  int iCount = 0;
  int iCountSelected = 0;

  //  TransientTrackingRecHit::RecHitContainer hitsTl;       

  for (trackingRecHit_iterator i=theTrack.recHitsBegin(); i!=theTrack.recHitsEnd(); i++){

    DetId id=DetId((theB->build(&**i))->hit()->geographicalId());
    const GeomDet* hitDet = (theB->build(&**i))->det();

    // Store detid corresponding to the first hit of the track
    //
    if ( iFirstIs ) {
      detidIs = id;
      iFirstIs = 0;
    }
    
    if ( myDebug_ ) {
      std::cout << "   #" << iCount;
      dumpModuleInfo(id);
      if ( typeid((theB->build(&**i))->hit()) == typeid(SiStripMatchedRecHit2D) ) std::cout << " Matched ";
      if ( SiStripDetId(id).stereo() ) std::cout << " Stereo ";
      if ( ! SiStripDetId(id).stereo() ) std::cout << " Mono ";
    }
    
    if ( 
	(( 
         ( id.subdetId() == StripSubdetector::TIB && hs.tib )
         ||
         ( id.subdetId() == StripSubdetector::TEC && hs.tec )
         ||
         ( id.subdetId() == StripSubdetector::TID && hs.tid )
         ||
         ( id.subdetId() == StripSubdetector::TOB && hs.tob )
	 )
	&&
	( 
	 hs.stereo 
	 || 
	 ( ! SiStripDetId(id).stereo() )
	 ))
	|| hs.all
	)
      {
	
	//if ( 1 ) {
	hitsTl.push_back(theB->build(&**i)); //Selected hit vector
	iCountSelected++;
	
	if ( myDebug_ ) std::cout << " <--Added to Tracklet";
	
	// Store detid corresponding to the first hit of the tracklet
	//
	if ( iFirstTl ) {
	  detidTl = id;
	}
	
	//
	iFirstTl = 0;
	
      } else if( ! iFirstTl ) {
      
      // Add here the invalid hit
      
      hitsTl.push_back(InvalidTransientRecHit::build(hitDet, TrackingRecHit::missing));
      if ( myDebug_ ) std::cout << " <--Invalid hit added to Tracklet";

    }
    
    if ( myDebug_ ) std::cout << std::endl;
    iCount++;


  }

  if ( myDebug_ ) std::cout << " Number of valid hits selected for refit: " << iCountSelected << std::endl;
  if ( myDebug_ ) std::cout << " Number of hits in the vector: " << hitsTl.size() << std::endl;

  if ( ! iCountSelected ) return 0;

  //Reverse re-loop on hits to drop trailing invalid hits 
  bool iFirstValidNotYetFound = true;
  if ( myDebug_ ) std::cout << " Reverse loop to remove trailing invalid hits" << std::endl;
  int iRecount = hitsTl.size()-1;
  for (TransientTrackingRecHit::RecHitContainer::const_iterator i=hitsTl.end()-1; i!=hitsTl.begin(); i--){

    if ( myDebug_ ) std::cout << " The hit #" << iRecount;
    if ( myDebug_ ) std::cout << " isValid:" << (*i)->isValid();
    if ( (*i)->isValid() ) iFirstValidNotYetFound = false; 
    if (! (*i)->isValid() && iFirstValidNotYetFound ){
      hitsTl.pop_back();
      if ( myDebug_ ) std::cout << " dropped";
    }
    if ( myDebug_ ) std::cout << std::endl;
  
    iRecount--;

  }

  if ( myDebug_ ) std::cout << " Number of hits in the vector after cleaning: " << hitsTl.size() << std::endl;

  return iCountSelected;

}

tsosParams ecRefit::doWithTrack(const reco::Track track, hitSelector hs){

  theTrack = track;

  tsosParams myParams;
  myParams.zero();
  
  //
  // Build hits vector  
  uint32_t detidTl;
  uint32_t detidIs;
  myParams.nhit = buildHitsVector(track, hs, detidIs, detidTl);
  myParams.detidIs = detidIs;
  myParams.detidTl = detidTl;

  if ( myParams.nhit <= hs.minhits ) return myParams;

  //
  // Refit
  std::vector<Trajectory> trajVec = doGenericRefit(theTrack, hitsTl, theG.product(), theMF.product());
    
  if ( myDebug_ ) std::cout << "   Tracklet refit initialization with..." << std::endl;

  if ( ! trajVec.size()>0 ) return myParams;
  
  TrajectoryStateOnSurface upTSOS = trajVec.begin()->lastMeasurement().updatedState();
  myParams = GetTSOSParams(upTSOS);
  myParams.iok = 1;
  
  if ( myDebug_ ) std::cout << "   Tracklet refit done. " << std::endl;
  if ( myDebug_ ) dumpTSOSInfo(upTSOS);
  
  return myParams;

}

ecRefit::~ecRefit()
{
 
}

//
// member functions
//

std::vector<Trajectory> ecRefit::doGenericRefit(const reco::Track t, TransientTrackingRecHit::RecHitContainer h, const TrackerGeometry * tG, const MagneticField * tM){
  
  reco::Track track = t;
  edm::ESHandle<TrackerGeometry> tGeo = tG;
  edm::ESHandle<MagneticField> tMF = tM;

  //
  // Initial state
  TrajectoryStateOnSurface initialStateFromTrack = buildInitialStateForEcRefit(track, h, tGeo.product(), tMF.product());
  
  //
  // Direction
  PropagationDirection seedDir = track.seedDirection();
  //
  // Seed
  const TrajectorySeed seed = TrajectorySeed(PTrajectoryStateOnDet(), TrajectorySeed::recHitContainer(), seedDir);
  //
  // Fit!
  return theF->fit(seed, h, initialStateFromTrack);
  
}

TrajectoryStateOnSurface ecRefit::buildInitialStateForEcRefit(const reco::Track track, TransientTrackingRecHit::RecHitContainer hits, const TrackerGeometry * theGeo, const MagneticField * theMagField){

  reco::Track theTrack = track;
  edm::ESHandle<TrackerGeometry> theG = theGeo;
  edm::ESHandle<MagneticField> theMF = theMagField;

  //
  // Initial state
  const TrajectoryStateOnSurface innerStateFromTrack(trajectoryStateTransform::innerStateOnSurface(track,*theG,theMF.product()));
  const TrajectoryStateOnSurface outerStateFromTrack(trajectoryStateTransform::outerStateOnSurface(track,*theG,theMF.product()));

  TrajectoryStateOnSurface initialStateFromTrack = 
    ( (innerStateFromTrack.globalPosition()-hits.front()->det()->position()).mag2() <
      (outerStateFromTrack.globalPosition()-hits.front()->det()->position()).mag2() ) ? 
    innerStateFromTrack: outerStateFromTrack;       
  initialStateFromTrack.rescaleError(100);
  theInitialStateForRefitting = TrajectoryStateOnSurface(initialStateFromTrack.localParameters(),
							 initialStateFromTrack.localError(),                  
							 initialStateFromTrack.surface(),
							 theMF.product());  
  
  return theInitialStateForRefitting;

}

double ecRefit::pErrorAtTSOS(TrajectoryStateOnSurface & TSOS){
  
  //Relative error on p is the same as error on curvature, i.e. the 

  double p = TSOS.globalMomentum().mag();
  double qoverperror = sqrt(TSOS.curvilinearError().matrix()(0,0));

  return p*p*qoverperror;

}


double ecRefit::ptErrorAtTSOS(TrajectoryStateOnSurface& TSOS){
  
  double px = TSOS.globalMomentum().x();
  double py = TSOS.globalMomentum().y();
  double pt = TSOS.globalMomentum().perp();

  AlgebraicSymMatrix66 errors = TSOS.cartesianError().matrix();

  double dpx2 = errors(3,3);
  double dpy2 = errors(4,4);
  double dpxdpy = errors(3,4);
  
  return sqrt(px*px*dpx2 + py*py*dpy2 + 2.*px*py*dpxdpy)/pt;
  
}

double ecRefit::pzErrorAtTSOS(TrajectoryStateOnSurface& TSOS){

  AlgebraicSymMatrix66 errors = TSOS.cartesianError().matrix();
  
  return sqrt(errors(5,5));
  
}

double ecRefit::paramErrorAtTSOS(TrajectoryStateOnSurface& TSOS, int i){
  
  AlgebraicSymMatrix55 errors = TSOS.curvilinearError().matrix();

  return sqrt(errors(i,i));
  
}

#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

void ecRefit::dumpTSOSInfo(TrajectoryStateOnSurface& TSOS){

  std::cout << "  >>TSOS>> p=" << TSOS.globalMomentum().mag() << "+-" << pErrorAtTSOS(TSOS) <<
    " pt=" << TSOS.globalMomentum().perp() << "+-" << ptErrorAtTSOS(TSOS) <<
    " pz=" << TSOS.globalMomentum().z() << "+-" << pzErrorAtTSOS(TSOS) <<
    " l=" << 0.5*M_PI-TSOS.globalMomentum().theta() << "+-" << paramErrorAtTSOS(TSOS, 1) <<
    " the=" << TSOS.globalMomentum().theta() << "+-" << paramErrorAtTSOS(TSOS, 1) << 
    " phi=" << TSOS.globalMomentum().phi() << "+-" << paramErrorAtTSOS(TSOS, 2) << 
    "; Position: rho:" << 
    TSOS.globalPosition().transverse() << " phi:" << TSOS.globalPosition().phi() << " z:" << 
    TSOS.globalPosition().z() << std::endl;

}

tsosParams ecRefit::GetTSOSParams(TrajectoryStateOnSurface& TSOS){

  tsosParams tp;
  
  tp.p = TSOS.globalMomentum().mag();
  tp.pErr = pErrorAtTSOS(TSOS);
  tp.curv = TSOS.signedInverseMomentum();
  tp.curvErr = paramErrorAtTSOS(TSOS, 0);
  tp.pt = TSOS.globalMomentum().perp();
  tp.ptErr = ptErrorAtTSOS(TSOS);
  tp.curvt = TSOS.charge()/tp.pt;
  tp.curvtErr = tp.ptErr*tp.curvt*tp.curvt;
  tp.pz = TSOS.globalMomentum().z();
  tp.pzErr = pzErrorAtTSOS(TSOS);
  tp.phi = TSOS.globalMomentum().phi();
  tp.phiErr = paramErrorAtTSOS(TSOS, 2);
  tp.theta = TSOS.globalMomentum().theta();
  tp.thetaErr = paramErrorAtTSOS(TSOS, 1);

  return tp;

}

void ecRefit::dumpModuleInfo(DetId hitId){
  
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

