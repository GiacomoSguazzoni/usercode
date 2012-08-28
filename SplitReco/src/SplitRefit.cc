//
// Original Author:  Giuseppe Cerati
//         Created:  Fri Aug  7 15:10:58 CEST 2009
// $Id: SplitRefit.cc,v 1.1 2012/08/02 14:21:37 sguazz Exp $
//
//
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/SimTracker/TrackAssociation/test/testTrackAssociator.cc?revision=1.17&view=markup&pathrev=CMSSW_2_2_10
//

#include "Tests/SplitReco/interface/SplitRefit.h"
#include "Tests/SplitReco/interface/Split.h"

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

#define cmsswVersion44x


//
// constructors and destructor
//
SplitRefit::SplitRefit(const TrackerGeometry * theGeometry, const MagneticField * theMagneticField, const TrajectoryFitter * theFitter, const TransientTrackingRecHitBuilder * theBuilder)
{

  myDebug_ = false;

  theG = theGeometry;
  theMF = theMagneticField;
  theF = theFitter;
  theB = theBuilder;

}

myTrack SplitRefit::initializeWithTrack(const reco::Track track){

  theTrack = track;

  //Store Quantities needed later
  theTrackTheta = theTrack.theta();
  theTrackThetaErr = theTrack.thetaError();
  theTrackInvSinTheta = 1./sin(theTrackTheta);

#ifdef cmsswVersion44x
  // in CMMSW 44x
  TrajectoryStateTransform transformer;
  const TrajectoryStateOnSurface innerStateFromTrack = transformer.innerStateOnSurface(theTrack,*theG,theMF.product());
  const TrajectoryStateOnSurface outerStateFromTrack = transformer.outerStateOnSurface(theTrack,*theG,theMF.product());
#else
  const TrajectoryStateOnSurface innerStateFromTrack(trajectoryStateTransform::innerStateOnSurface(theTrack,*theG,theMF.product()));
  const TrajectoryStateOnSurface outerStateFromTrack(trajectoryStateTransform::outerStateOnSurface(theTrack,*theG,theMF.product()));
#endif
  
  //First refit the track with all the hits to have all intermediate TSOS
  // 
  // Hits
  for (trackingRecHit_iterator i=theTrack.recHitsBegin(); i!=theTrack.recHitsEnd(); i++){
    hitsAll.push_back(theB->build(&**i ));
  }
  //
  // Initial state
  //  TrajectoryStateOnSurface theInitialStateForRefitting;
  TrajectoryStateOnSurface initialStateFromTrack = 
    ( (innerStateFromTrack.globalPosition()-hitsAll.front()->det()->position()).mag2() <
      (outerStateFromTrack.globalPosition()-hitsAll.front()->det()->position()).mag2() ) ? 
    innerStateFromTrack: outerStateFromTrack;       
  initialStateFromTrack.rescaleError(100);
  theInitialStateForRefitting = TrajectoryStateOnSurface(initialStateFromTrack.localParameters(),
							 initialStateFromTrack.localError(),                  
							 initialStateFromTrack.surface(),
							 theMF.product());  
  
  std::vector<Trajectory> trajVec = doGenericRefit(hitsAll, theInitialStateForRefitting);
    
  //
  // Re-initialize quantities needed for split refit
  vecT.clear();
  vecEffHitN.clear();
  vecEffHitBegin.clear();
  vecEffHitEnd.clear();
  vecPerpR.clear();
  vecX.clear();
  vecY.clear();
  vecZ.clear();
  vecDPerpT.clear();

  mytrack.nEffHit=-1;

  if ( myDebug_ ) std::cout << " ++++++++++++++++++++++++++++++++++++++++" << std::endl;
  if ( myDebug_ ) std::cout << " Split refit initialization with track..." << std::endl;

  if ( ! trajVec.size()>0 ) return mytrack;

    theTrajectoryMeasurements=trajVec.begin()->measurements();
    std::vector<TrajectoryMeasurement>::iterator itm;
       
    int iHitN=0;
    int iEffHitN=0;
    double x0=(theTrajectoryMeasurements.end()-1)->updatedState().globalPosition().x();
    double y0=(theTrajectoryMeasurements.end()-1)->updatedState().globalPosition().y();
    double z0=(theTrajectoryMeasurements.end()-1)->updatedState().globalPosition().z();
    double T = 0.;
    double pAveTrack=0.;
    double pErrAveTrack=0.;
    for (itm=theTrajectoryMeasurements.end()-1;itm!=theTrajectoryMeasurements.begin()-1;itm--){ //loop on hits
      iHitN++;
      TrajectoryStateOnSurface upTSOS = itm->updatedState();

      double x=upTSOS.globalPosition().x();
      double y=upTSOS.globalPosition().y();
      double z=upTSOS.globalPosition().z();
      double perpR=upTSOS.globalPosition().perp();
      vecPerpR.push_back(perpR);
      vecX.push_back(x);
      vecY.push_back(y);
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
      if ( myDebug_ ) std::cout
	<< " >>>>TSOS>> nH" << iHitN
	<< " nEffH" << iEffHitN
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
	<< std::endl;
      x0=x;
      y0=y;
      z0=z;
    }

    mytrack.pAve=pAveTrack/iHitN;
    mytrack.pAveErr=pErrAveTrack/iHitN;
    mytrack.nHit=iHitN;
    mytrack.nEffHit=iEffHitN;
    mytrack.T=T;
    mytrack.rIn=vecPerpR.at(0);
    mytrack.rOut=vecPerpR.at(iHitN-1);
    mytrack.xIn=vecX.at(0);
    mytrack.yIn=vecY.at(0);
    mytrack.zIn=vecZ.at(0);
    mytrack.xOut=vecX.at(iHitN-1);
    mytrack.yOut=vecY.at(iHitN-1);
    mytrack.zOut=vecZ.at(iHitN-1);

    if ( myDebug_ ) std::cout << " Initialization successful; track refit done. " << std::endl;

    return mytrack;

}


SplitRefit::~SplitRefit()
{
 
}

//
// member functions
//

std::vector<Trajectory> SplitRefit::doGenericRefit(TransientTrackingRecHit::RecHitContainer hits, TrajectoryStateOnSurface iTSOS){

    //
    // Direction
    PropagationDirection seedDir = theTrack.seedDirection();
    //
    // Seed
    const TrajectorySeed seed = TrajectorySeed(PTrajectoryStateOnDet(), TrajectorySeed::recHitContainer(), seedDir);
    //
    // Fit!
    return theF->fit(seed, hits, iTSOS);

}

std::vector<Split> SplitRefit::doSplitRefits(int nMinHitSplit, int nHitIncrement, double kFactor, double pullCut){

  setMaterialToKFactor(kFactor);

  std::vector<Split> splits;

  TransientTrackingRecHit::RecHitContainer hitsSplit;       
  int iEffFirst=1;
  int iEffLast=iEffFirst+nMinHitSplit-1;

  while ( iEffLast<mytrack.nEffHit+1 ) {

    if ( mytrack.nEffHit - iEffLast < nHitIncrement ) iEffLast = mytrack.nEffHit; //no place for another split

    if ( myDebug_ ) std::cout << " NEXT iEffFirst=" << iEffFirst << " " << " iEffLast=" << iEffLast << " " << std::endl;
    if ( myDebug_ ) std::cout << " Now using hits"; 
    int iFirst=vecEffHitBegin.at(iEffFirst-1);
    int iLast=vecEffHitBegin.at(iEffLast-1);
    int iCurrent=iFirst;
    //
    // Build the rec hit container for the Split Track
    while ( iCurrent<vecEffHitEnd.at(iEffLast-1)+1 ) { 
      // Address first hit theTrack.recHitsBegin()
      trackingRecHit_iterator i=theTrack.recHitsBegin()+iCurrent-1;
      if ( myDebug_ ) std::cout << " #" << iCurrent << " dPerpT=" << vecDPerpT.at(iCurrent-1);
      iCurrent++;
      hitsSplit.push_back(theB->build(&**i ));
      i++;
    }

    //
    // Set the proper TSOS
    TrajectoryStateOnSurface theInitialStateForSplitRefitting;
    if ( iFirst == 1 ) {
      theInitialStateForSplitRefitting=theInitialStateForRefitting;
      if ( myDebug_ ) std::cout << " iFirst=" << iFirst << " TSOS=" << theInitialStateForSplitRefitting.globalMomentum().mag(); 
    } else {
      TrajectoryStateOnSurface theStateFromTraj = (theTrajectoryMeasurements.end()-iFirst)->forwardPredictedState();
      theStateFromTraj.rescaleError(100);
      theInitialStateForSplitRefitting = TrajectoryStateOnSurface(theStateFromTraj.localParameters(),
								  theStateFromTraj.localError(),                  
								  theStateFromTraj.surface(),
								  theMF.product());  
      if ( myDebug_ ) std::cout << " iFirst=" << iFirst << " TSOS=" << theInitialStateForSplitRefitting.globalMomentum().mag(); 
    }
	 
    std::vector<Trajectory> trajVecSplit = doGenericRefit(hitsSplit, theInitialStateForSplitRefitting);

    //
    //
    if (trajVecSplit.size()>0) {

      //
      // For the energy loss fit...
      Split split;

      TrajectoryStateOnSurface firstSplitTSOS = trajVecSplit.begin()->lastMeasurement().updatedState();
      split.pt=firstSplitTSOS.globalMomentum().perp();;
      split.ptErr=ptErrorAtTSOS(firstSplitTSOS);;
      split.T=0.5*(vecT.at(iLast-1)+vecT.at(iFirst-1));
      split.TErr=0.05;  //half millimiter error (?)
      splits.push_back(split);

      if ( myDebug_ ) std::cout << " Split pt=" << split.pt << "+-" << split.ptErr << std::endl;

    } else {

      if ( myDebug_ ) std::cout << " split fit failed" << std::endl;

    }
	 
    iEffFirst+=nHitIncrement;
    iEffLast+=nHitIncrement;
    //
    //Empty hits container
    hitsSplit.clear();
  }

  //
  // Put the material back to original values!
  setMaterialToKFactor(1./kFactor);
    
  return splits;

}

energyLoss SplitRefit::energyLossFitOnSplits(std::vector<Split> &splits){

  int minSplits = 3; //protection
  int nsplits = splits.size();  

  energyLoss eLoss;
  eLoss.chi2 = -1.;
  eLoss.nsplits = nsplits;

  if ( nsplits < minSplits ) return eLoss;
  
  //
  // Fill arrays for graph

#ifdef cmsswVersion44x
  // in CMMSW 44x fixed lenght arrays required
  const int nsize = 50;
#else
  int nsize = nsplits;
#endif

  double arrPtSplit[nsize];
  double arrPtErrSplit[nsize];
  double arrTSplit[nsize];
  double arrTErrSplit[nsize];
  
  int iarr = 0;

  for (std::vector<Split>::iterator it=splits.begin() ; it < splits.end(); it++ ){

    arrPtSplit[iarr]=(*it).pt;
    arrPtErrSplit[iarr]=(*it).ptErr;
    arrTSplit[iarr]=(*it).T;
    arrTErrSplit[iarr]=(*it).TErr;
    iarr++;

  }  

  TGraphErrors graphForFit(nsplits,arrTSplit,arrPtSplit,arrTErrSplit,arrPtErrSplit);  
  int fitResult = graphForFit.Fit("pol1","Q");
  TF1 *fit = graphForFit.GetFunction("pol1");
  if ( fitResult ) {
    if ( myDebug_ ) std::cout << " Linear fit of splits failed with code: " << fitResult << std::endl;
    return eLoss;
  } 

  //
  // In case of valid fit
  eLoss.p0 = fit->GetParameter(0);
  eLoss.dpdx = -1.*fit->GetParameter(1)*theTrackInvSinTheta; //Rescale for theta and revert sign
  eLoss.dpdxErr = fit->GetParError(1)*theTrackInvSinTheta; //Rescale for theta
  eLoss.chi2 = fit->GetChisquare();
  eLoss.freePar = fit->GetNumberFreeParameters();
  if ( myDebug_ ) std::cout << " linear EL fit p0=" << eLoss.p0 << " p1=" << eLoss.dpdx << " chi2/freePar=" << eLoss.chi2 << "/" << eLoss.freePar << std::endl; 
  eLoss.TIn  = (*(splits.begin())).T;
  double TInErr = (*(splits.begin())).TErr;
  eLoss.TOut = (*(splits.end()-1)).T;
  double TOutErr = (*(splits.end()-1)).TErr;
  double TTot = eLoss.TOut - eLoss.TIn;
  eLoss.pLoss = eLoss.dpdx*(eLoss.TOut - eLoss.TIn);
  eLoss.pLossErr = sqrt(TTot*TTot*eLoss.dpdxErr*eLoss.dpdxErr+eLoss.dpdx*eLoss.dpdx*(TInErr*TInErr+TOutErr*TOutErr));
  if ( myDebug_ ) {
    double pIn  = eLoss.p0+eLoss.dpdx*eLoss.TIn;
    double pOut = eLoss.p0+eLoss.dpdx*eLoss.TOut;
    std::cout << " TIn=" << eLoss.TIn << " pIn=" << pIn << "; TOut=" << eLoss.TOut << " pOut=" << pOut << " Ttot=" << TTot << " pLoss=" << eLoss.pLoss << " pLossErr=" << eLoss.pLossErr  << std::endl; 
  }

  return eLoss;

}

double SplitRefit::pErrorAtTSOS(TrajectoryStateOnSurface & TSOS){
  
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

double SplitRefit::ptErrorAtTSOS(TrajectoryStateOnSurface & TSOS){
  
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

//void SplitRefit::setMaterialToKFactor(const TrackerGeometry * theGeometry, TransientTrackingRecHit::RecHitContainer& hits, double kFactor){
void SplitRefit::setMaterialToKFactor(double kFactor){
  
  const TrackerGeometry * theGeometry = theG.product();

  std::string outputBlock;
  char buffer [300];
  int nhit=0;

  for (TransientTrackingRecHit::RecHitContainer::const_iterator it=hitsAll.begin(); it!=hitsAll.end();it++){
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
  
