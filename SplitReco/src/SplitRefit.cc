//
// Original Author:  Giuseppe Cerati
//         Created:  Fri Aug  7 15:10:58 CEST 2009
// $Id: SplitRefit.cc,v 1.3 2012/09/18 15:17:31 sguazz Exp $
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
SplitRefit::SplitRefit(const TrackerGeometry * theGeometry, const MagneticField * theMagneticField, const TrajectoryFitter * theFitter, const TransientTrackingRecHitBuilder * theBuilder, bool debug)
{

  myDebug_ = debug;

  theG = theGeometry;
  theMF = theMagneticField;
  theF = theFitter;
  theB = theBuilder;

}

myTrack SplitRefit::initializeWithTrack(const reco::Track track){

  theTrack = track;

  //
  //Store Quantities needed later
  //
  theTrackEta = theTrack.eta();
  theTrackTheta = theTrack.theta();
  theTrackThetaErr = theTrack.thetaError();
  theTrackInvSinTheta = 1./sin(theTrackTheta);
  theTrackPt = theTrack.pt();
  theTrackPz = theTrack.pz();
  theTrackP = theTrack.p();
  theTrackPtSquared = theTrackPt*theTrackPt;
  theTrackPSquared = theTrackP*theTrackP;
  theTrackLambda = theTrack.lambda();
  theTrackInvSinLambda = 1./sin(theTrackLambda);
  theTrackInvCosLambda = 1./cos(theTrackLambda);
  theTrackCharge = theTrack.charge();
  BInTesla = theMF->inTesla(GlobalPoint(0., 0., 0.)).z(); //Nominal B field

  //
  // Initialize Jacobians
  double sinlambda = sin(theTrackLambda);
  double coslambda = cos(theTrackLambda);

  // To pass from matrix with 1/P (default) to 1/Pt
  jacoInvPToInvPt                = AlgebraicMatrixID();
  jacoInvPToInvPtTransposed      = AlgebraicMatrixID();
  jacoInvPToInvPt(0,0)           = 1./coslambda; 
  jacoInvPToInvPtTransposed(0,0) = jacoInvPToInvPt(0,0);
  jacoInvPToInvPt(0,1)           = sinlambda/coslambda/coslambda/theTrackP; 
  jacoInvPToInvPtTransposed(1,0) = jacoInvPToInvPt(0,1);  

  // To pass from matrix with 1/Pt to 1/P (default)
  jacoInvPtToInvP                = AlgebraicMatrixID();
  jacoInvPtToInvPTransposed      = AlgebraicMatrixID();
  jacoInvPtToInvP(0,0)           = coslambda;
  jacoInvPtToInvPTransposed(0,0) = coslambda;
  jacoInvPtToInvP(0,1)           = -1.*sinlambda/theTrackPt; 
  jacoInvPtToInvPTransposed(1,0) = jacoInvPtToInvP(0,1);

  // To pass from matrix with 1/P (default) to 1/Pz
  jacoInvPToInvPz                = AlgebraicMatrixID();
  jacoInvPToInvPzTransposed      = AlgebraicMatrixID();
  jacoInvPToInvPz(0,0)           = 1./sinlambda; 
  jacoInvPToInvPzTransposed(0,0) = jacoInvPToInvPz(0,0);
  jacoInvPToInvPz(0,1)           = -1.*coslambda/sinlambda/sinlambda/theTrackP; 
  jacoInvPToInvPzTransposed(1,0) = jacoInvPToInvPz(0,1);  

  // To pass from matrix with 1/Pz to 1/P (default)
  jacoInvPzToInvP                = AlgebraicMatrixID();
  jacoInvPzToInvPTransposed      = AlgebraicMatrixID();
  jacoInvPzToInvP(0,0)           = sinlambda;
  jacoInvPzToInvPTransposed(0,0) = sinlambda;
  jacoInvPzToInvP(0,1)           = coslambda/theTrackPz; 
  jacoInvPzToInvPTransposed(1,0) = jacoInvPzToInvP(0,1);

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

    if ( myDebug_ ) std::cout << " ######################################################################################### " << std::endl;
    if ( myDebug_ ) std::cout << " ######################################################################################### " << std::endl;
    if ( myDebug_ ) std::cout << " ######################################################################################### " << std::endl;
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
    double maxhitchi2 = -1.;
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
      //      double theta=upTSOS.globalMomentum().theta();
      double dPerpT=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
      double hitchi2 = itm->estimate();
      if ( hitchi2 > maxhitchi2 ) maxhitchi2 = hitchi2;
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
      //
      if ( myDebug_ ) {
	std::cout
	  << " >>>>Measurement>> nH" << iHitN
	  << " nEffH" << iEffHitN;
	DetId detectorId=DetId(itm->recHit()->geographicalId());
	dumpModuleInfo(detectorId);
	//	<< " sinth=" << sin(theta)
	std::cout << " r=" << perpR
	  //	<< " x=" << x
	  //	<< " y=" << y
		  << " z=" << z
	  //	<< " dx=" << x-x0
	  //	<< " dy=" << y-y0
	  //	<< " dPerpT=" << dPerpT
		  << " effN=" << iEffHitN
		  << " beg:" << vecEffHitBegin.at(iEffHitN-1)
		  << " end:" << vecEffHitEnd.at(iEffHitN-1)
		  << " dT=" << dT
		  << " T=" << T
		  << " chi2=" << hitchi2
		  << " maxchi2=" << maxhitchi2
		  << std::endl;
	TrajectoryStateOnSurface fwTSOS = itm->forwardPredictedState();
	std::cout << ">>>>>fw>>";
	dumpTSOSInfo(fwTSOS);
	//
	TrajectoryStateOnSurface bwTSOS = itm->backwardPredictedState();
	std::cout << ">>>>>bw>>";
	dumpTSOSInfo(bwTSOS);
	//
	std::cout << ">>>>>up>>";
	dumpTSOSInfo(upTSOS);
      }
      
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
    mytrack.maxhitchi2=maxhitchi2;

    if ( myDebug_ ) std::cout << " Initialization successful; track refit done. " << std::endl;
    if ( myDebug_ ) std::cout << " ######################################################################################### " << std::endl;

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

std::vector<Split> SplitRefit::doSplitRefits(int nMinHitSplit, int nHitIncrement, double kFactor, double absEtaBarrelEndcapCut,
					     rescaleParams barrelRPar, rescaleParams endcapRPar){

  setMaterialToKFactor(kFactor);

  std::vector<Split> splits;

  TransientTrackingRecHit::RecHitContainer hitsSplit;       
  int iEffFirst=1;
  int iEffLast=iEffFirst+nMinHitSplit-1;

  while ( iEffLast<mytrack.nEffHit+1 ) {

    if ( mytrack.nEffHit - iEffLast < nHitIncrement ) iEffLast = mytrack.nEffHit; //no place for another split

    if ( myDebug_ ) std::cout << " ######################################################################################### " << std::endl;
    if ( myDebug_ ) std::cout << " ## New split " << std::endl;
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
    if ( myDebug_ ) std::cout << std::endl;

    //
    // Set the proper TSOS
    //    if ( iFirst == 1 ) {
    //      theInitialStateForSplitRefitting=theInitialStateForRefitting;
    //      if ( myDebug_ ) std::cout << " iFirst=" << iFirst << " TSOS=" << theInitialStateForSplitRefitting.globalMomentum().mag(); 
    //    } else {
    

    //Always use the updated state
    TrajectoryStateOnSurface theStateFromTraj = (theTrajectoryMeasurements.end()-iFirst)->updatedState();

    TrajectoryStateOnSurface theInitialStateForSplitRefitting = buildInitialStateForSplitRefit(theStateFromTraj, absEtaBarrelEndcapCut, barrelRPar, endcapRPar);

    if ( myDebug_ ) std::cout << " iFirst=" << iFirst << " TSOS=" << theInitialStateForSplitRefitting.globalMomentum().mag() << std::endl; 
    
    std::vector<Trajectory> trajVecSplit = doGenericRefit(hitsSplit, theInitialStateForSplitRefitting);

    //
    //
    if (trajVecSplit.size()>0) {

      std::vector<TrajectoryMeasurement> theSplitMeasurements;
      theSplitMeasurements=trajVecSplit.begin()->measurements();
 

      if ( myDebug_ ){

	int istate = 0;

	for (std::vector<TrajectoryMeasurement>::iterator itm=theSplitMeasurements.end()-1;itm!=theSplitMeasurements.begin()-1;itm--){ //loop on hits

	  DetId detectorId=DetId(itm->recHit()->geographicalId());
	  
	  std::cout << ">>>>>> Measurement #" << istate;
	  dumpModuleInfo(detectorId);
	  std::cout << " chi2:" << itm->estimate() << std::endl;
	  //
	  TrajectoryStateOnSurface fwTSOS = itm->forwardPredictedState();
	  std::cout << ">>>>>fw>>";
	  dumpTSOSInfo(fwTSOS);
	  //
	  TrajectoryStateOnSurface bwTSOS = itm->backwardPredictedState();
	  std::cout << ">>>>>bw>>";
	  dumpTSOSInfo(bwTSOS);
	  //
	  TrajectoryStateOnSurface upTSOS = itm->updatedState();
	  std::cout << ">>>>>up>>";
	  dumpTSOSInfo(upTSOS);

	  istate++;
	}
      }



      //
      // For the energy loss fit...
      Split split;

      double pAve, pAveErr, pAveT, pAveTErr, pAveZ, pAveZErr;
      GetPAveFromMeasurements(theSplitMeasurements, pAve, pAveErr, pAveT, pAveTErr, pAveZ, pAveZErr);

      TrajectoryStateOnSurface firstSplitTSOS = trajVecSplit.begin()->lastMeasurement().updatedState();
      split.pt=firstSplitTSOS.globalMomentum().perp();
      split.ptErr=ptErrorAtTSOS(firstSplitTSOS);
      // split.curv=abs(firstSplitTSOS.transverseCurvature()); //Apparently this method returns zero...
      split.curv=1./pAve;
      split.curvErr=split.curv*(pAveErr/pAve); //Relative curv error is the same as relative p error
      split.curvt=1./pAveT;
      split.curvtErr=split.curv*(pAveTErr/pAveT); //Relative curv error is the same as relative pt error
      split.curvz=1./pAveZ;
      split.curvzErr=split.curv*(pAveZErr/pAveZ); //Relative curv error is the same as relative pz error
      split.T=0.5*(vecT.at(iLast-1)+vecT.at(iFirst-1));
      split.TErr=0.05;  //half millimiter error (?)
      splits.push_back(split);

      if ( myDebug_ ) {
	double lambda = (M_PI/2.-firstSplitTSOS.globalMomentum().theta());
	double lambdaErr = lambdaErrorAtTSOS(firstSplitTSOS);
	double p = firstSplitTSOS.globalMomentum().mag();
	double pErr = pErrorAtTSOS(firstSplitTSOS);
	double pz = firstSplitTSOS.globalMomentum().z();
	double pzErr = pzErrorAtTSOS(firstSplitTSOS);
	double pFromPt = fabs(split.pt/cos(theTrackLambda)); 
	double pErrFromPt = split.ptErr/cos(theTrackLambda); 
	double pFromPz = pz/sin(theTrackLambda); 
	double pErrFromPz = fabs(pzErr/sin(theTrackLambda));
	double splitChi2 = (trajVecSplit.begin()->chiSquared()/(1.*trajVecSplit.begin()->ndof()));
	std::cout << " ???????????????????????????????????????????????????????????????????????? " << std::endl;
	std::cout << " Split T=" << split.T << " pt=" << split.pt << "+-" << split.ptErr << " pz=" << pz << "+-" << pzErr << " lambda:" << lambda << "+-" << lambdaErr << " chi2/ndof:" << splitChi2 << std::endl;
	std::cout << "  p=" << p << "+-" << pErr << "  pFromPt=" << pFromPt << "+-" << pErrFromPt << " pFromPz=" << pFromPz << "+-" << pErrFromPz << std::endl;
	std::cout << "  pAv=" << pAve << "+-" << pAveErr << "  pAv=" << pAveT << "+-" << pAveTErr << "  pAv=" << pAveZ << "+-" << pAveZErr << std::endl;
	std::cout << " ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; " << std::endl;
	std::cout << firstSplitTSOS.curvilinearError().matrix()  << std::endl;
	std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
      }

    } else {

      if ( myDebug_ ) std::cout << " split refit failed" << std::endl;

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

TrajectoryStateOnSurface SplitRefit::buildInitialStateForSplitRefit(TrajectoryStateOnSurface theStateFromTraj, double absEtaBarrelEndcapCut,
					     rescaleParams barrelRPar, rescaleParams endcapRPar){

  TrajectoryStateOnSurface theInitialStateForSplitRefitting;

  if ( fabs(theTrackEta) < absEtaBarrelEndcapCut ) {
    //
    //Barrel
    if ( barrelRPar.iSpecial ) {
      theStateFromTraj.rescaleError(100);
      AlgebraicSymMatrix55 newErrorMatrix = rescaleErrorOfComponents(theStateFromTraj, 0., barrelRPar);
      // Constructor from global ref system
      theInitialStateForSplitRefitting = TrajectoryStateOnSurface( theStateFromTraj.globalParameters(),
								   CurvilinearTrajectoryError(newErrorMatrix),
								   theStateFromTraj.surface());
    } else {
      theStateFromTraj.rescaleError(100);
      //Constructor from local ref system
      theInitialStateForSplitRefitting = TrajectoryStateOnSurface(theStateFromTraj.localParameters(),
								  theStateFromTraj.localError(),                  
								  theStateFromTraj.surface(),
								  theMF.product());  
    }
    //
  } else {
    //
    //Endcap
    if ( barrelRPar.iSpecial ) {
      theStateFromTraj.rescaleError(100);
      AlgebraicSymMatrix55 newErrorMatrix = rescaleErrorOfComponents(theStateFromTraj, 0., endcapRPar);
      // Constructor from global ref system
      theInitialStateForSplitRefitting = TrajectoryStateOnSurface( theStateFromTraj.globalParameters(),
								   CurvilinearTrajectoryError(newErrorMatrix),
								   theStateFromTraj.surface());
    } else {
      theStateFromTraj.rescaleError(100);
      //Constructor from local ref system
      theInitialStateForSplitRefitting = TrajectoryStateOnSurface(theStateFromTraj.localParameters(),
								  theStateFromTraj.localError(),                  
								  theStateFromTraj.surface(),
								  theMF.product());  
    }
    //
  }

  return theInitialStateForSplitRefitting;

}

std::vector<energyLoss> SplitRefit::energyLossFitOnSplits(std::vector<Split> &splits){

  std::vector<energyLoss> eLosses;

  //  eLosses.push_back(energyLossFitByPtOnSplits(splits));
  eLosses.push_back(energyLossFitByCurvOnSplits(splits,0));
  eLosses.push_back(energyLossFitByCurvOnSplits(splits,1));
  eLosses.push_back(energyLossFitByCurvOnSplits(splits,2));

  return eLosses;

}

energyLoss SplitRefit::energyLossFitByPtOnSplits(std::vector<Split> &splits){

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

  if ( myDebug_ ) std::cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  if ( myDebug_ ) std::cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  if ( myDebug_ ) std::cout << " >>>>>>>>>>>>>>>>>>>>> Starting split value **pt** fit " << std::endl;

  for (std::vector<Split>::iterator it=splits.begin() ; it < splits.end(); it++ ){

    arrPtSplit[iarr]=(*it).pt;
    arrPtErrSplit[iarr]=(*it).ptErr;
    arrTSplit[iarr]=(*it).T;
    arrTErrSplit[iarr]=(*it).TErr;
    iarr++;

  if ( myDebug_ ) std::cout << " Split value for fit #" << iarr << " T:" << (*it).T << " TErr:" << (*it).TErr << " pt:" << (*it).pt << " ptErr:" << (*it).ptErr << std::endl;

  }  

  // Old style fits on pt values
  //
  TGraphErrors graphForFit(nsplits,arrTSplit,arrPtSplit,arrTErrSplit,arrPtErrSplit);  
  int fitResult = graphForFit.Fit("pol1","Q");
  TF1 *fit = graphForFit.GetFunction("pol1");
  if ( fitResult ) {
    if ( myDebug_ ) std::cout << " Linear fit of split pt failed with code: " << fitResult << std::endl;
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

energyLoss SplitRefit::energyLossFitByCurvOnSplits(std::vector<Split> &splits, int iSwitch){

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

  double arrCurvSplit[nsize];
  double arrCurvErrSplit[nsize];
  double arrTSplit[nsize];
  double arrTErrSplit[nsize];
  
  int iarr = 0;

  if ( myDebug_ ) std::cout << " >>>>>>>>>>>>>>>>>>>>> Starting split value **curvature** fit " << std::endl;
  
  for (std::vector<Split>::iterator it=splits.begin() ; it < splits.end(); it++ ){

    if ( iSwitch == 0 ){
      arrCurvSplit[iarr]=(*it).curv;
      arrCurvErrSplit[iarr]=(*it).curvErr;
    }
    if ( iSwitch == 1 ){
      arrCurvSplit[iarr]=(*it).curvt;
      arrCurvErrSplit[iarr]=(*it).curvtErr;
    }
    if ( iSwitch == 2 ){
      arrCurvSplit[iarr]=(*it).curvz;
      arrCurvErrSplit[iarr]=(*it).curvzErr;
    }

    arrTSplit[iarr]=(*it).T;
    arrTErrSplit[iarr]=(*it).TErr;

    if ( myDebug_ ) std::cout << " Split value for fit #" << iarr << " T:" << (*it).T << " TErr:" << (*it).TErr << " curv:" << arrCurvSplit[iarr] << " curvErr:" << arrCurvErrSplit[iarr] << std::endl;

    iarr++;


  }  

  // New style fits on pt values
  //
  try {
    TGraphErrors graphForFit(nsplits,arrTSplit,arrCurvSplit,arrTErrSplit,arrCurvErrSplit);  
    int fitResult = graphForFit.Fit("pol1","Q");
    TF1 *fit = graphForFit.GetFunction("pol1");
    if ( fitResult ) {
      if ( myDebug_ ) std::cout << " Linear fit of split curv failed with code: " << fitResult << std::endl;
      return eLoss;
    } 
    
    //
    // In case of valid fit
    eLoss.p0 = fit->GetParameter(0);
    // 
    // Approximated relation between dCurv/dx and dP/dx
    //
    // dP/dx = -(1/sin theta) (pt*pt/0.3/B) dCurv/dx
    //
    // with 3.33566830114413 = 1./0.29979
    //
    // we fill the vector with Chi = 1./pt . Formula simplifies
    //
    // dP/dx = -(1/sin theta) (pt*pt) dChi/dx
    
    //  double multFactor = theTrackInvSinTheta*theTrackPtSquared*3.33566830114413/BInTesla;
    //  double multFactor = theTrackInvSinTheta*theTrackPtSquared;
    
    double multFactor = theTrackPSquared;
    
    eLoss.dpdx = fit->GetParameter(1)*multFactor;
    eLoss.dpdxErr = fit->GetParError(1)*multFactor; 
    eLoss.chi2 = fit->GetChisquare();
    eLoss.freePar = fit->GetNumberFreeParameters();
    
    
    if ( myDebug_ ) std::cout << " linear EL by Curv fit p0=" << eLoss.p0 << " dpdx=" << eLoss.dpdx << " chi2/freePar=" << eLoss.chi2 << "/" << eLoss.freePar << std::endl; 
    eLoss.TIn  = (*(splits.begin())).T;
    double TInErr = (*(splits.begin())).TErr;
    eLoss.TOut = (*(splits.end()-1)).T;
    double TOutErr = (*(splits.end()-1)).TErr;
    double TTot = eLoss.TOut - eLoss.TIn;
    eLoss.pLoss = eLoss.dpdx*(eLoss.TOut - eLoss.TIn);
    eLoss.pLossErr = sqrt(TTot*TTot*eLoss.dpdxErr*eLoss.dpdxErr+eLoss.dpdx*eLoss.dpdx*(TInErr*TInErr+TOutErr*TOutErr));
    if ( myDebug_ ) {
      std::cout << " TIn=" << eLoss.TIn << "; TOut=" << eLoss.TOut << " Ttot=" << TTot << " pLoss=" << eLoss.pLoss << " pLossErr=" << eLoss.pLossErr  << std::endl; 
    }
    
  } catch (...) {
    
    std::cout << " >>>>> Exception caught in linear fit <<<<<" << std::endl; 
    
  }
  
  return eLoss;

}

double SplitRefit::pErrorAtTSOS(TrajectoryStateOnSurface & TSOS){
  
  //Relative error on p is the same as error on curvature, i.e. the 

  double p = TSOS.globalMomentum().mag();
  double qoverperror = sqrt(TSOS.curvilinearError().matrix()(0,0));

  return p*p*qoverperror;

}

void SplitRefit::GetPAveFromMeasurements(std::vector<TrajectoryMeasurement> theSplitMeasurements, double& pAve, double& pAveErr, double& pAveT, double& pAveTErr, double& pAveZ, double& pAveZErr){

  pAve = 0.;
  pAveErr = 0.;
  pAveT = 0.;
  pAveTErr = 0.;
  pAveZ = 0.;
  pAveZErr = 0.;
  for (std::vector<TrajectoryMeasurement>::iterator itm=theSplitMeasurements.end()-1;itm!=theSplitMeasurements.begin()-1;itm--){ //loop on hits
    TrajectoryStateOnSurface upTSOS = itm->updatedState();
    double pErr2 = pErrorAtTSOS(upTSOS);
    double pTErr2 = ptErrorAtTSOS(upTSOS);
    double pZErr2 = pzErrorAtTSOS(upTSOS)*theTrackInvSinLambda;
    //
    pErr2*=pErr2;
    pAve += upTSOS.globalMomentum().mag()/pErr2;
    pAveErr += 1./pErr2;
    //
    pTErr2*=pTErr2;
    pAveT += upTSOS.globalMomentum().perp()*theTrackInvCosLambda/pTErr2;
    pAveTErr += 1./pTErr2;
    //
    pZErr2*=pZErr2;
    pAveZ += upTSOS.globalMomentum().z()*theTrackInvSinLambda/pZErr2;
    pAveZErr += 1./pZErr2;

  }
  
  pAveErr = 1./pAveErr;
  pAve *= pAveErr;
  pAveErr = sqrt(pAveErr);
  //
  pAveTErr = 1./pAveTErr; 
  pAveT *= pAveTErr;
  pAveTErr = sqrt(pAveTErr);
  //
  pAveZErr = 1./pAveZErr;
  pAveZ *= pAveZErr;
  pAveZErr = sqrt(pAveZErr);

  /*
  // sguazzTemp
  pAveErr = 1./(1./pAveZErr/pAveZErr + 1./pAveTErr/pAveTErr);
  pAve = pAveErr*(pAveZ/pAveZErr/pAveZErr + pAveT/pAveTErr/pAveTErr);
  pAveErr = sqrt(pAveErr);
  */

}

double SplitRefit::ptErrorAtTSOS(TrajectoryStateOnSurface& TSOS){
  
  double px = TSOS.globalMomentum().x();
  double py = TSOS.globalMomentum().y();
  double pt = TSOS.globalMomentum().perp();

  AlgebraicSymMatrix66 errors = TSOS.cartesianError().matrix();

  double dpx2 = errors(3,3);
  double dpy2 = errors(4,4);
  double dpxdpy = errors(3,4);
  
  return sqrt(px*px*dpx2 + py*py*dpy2 + 2.*px*py*dpxdpy)/pt;
  
}

double SplitRefit::pzErrorAtTSOS(TrajectoryStateOnSurface& TSOS){
  
  AlgebraicSymMatrix66 errors = TSOS.cartesianError().matrix();

  return sqrt(errors(5,5));
  
}

double SplitRefit::lambdaErrorAtTSOS(TrajectoryStateOnSurface& TSOS){
  
  AlgebraicSymMatrix55 errors = TSOS.curvilinearError().matrix();

  return sqrt(errors(1,1));
  
}

AlgebraicSymMatrix55 SplitRefit::rescaleErrorOfComponents(TrajectoryStateOnSurface& TSOS, double kFactor, rescaleParams RPar){

  AlgebraicSymMatrix55 oldErrorMatrix = TSOS.curvilinearError().matrix();

  if ( myDebug_ ) {
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << ">> ***[Matrices multiplied by 100. to improve readibility]*** >>>>>>>>>>" << std::endl;
    std::cout << ">>>>>>>>>>>>>>> ErrorMatrix in P before rescaling:" << std::endl;
    std::cout << 100.*oldErrorMatrix << std::endl;
  }
  
  AlgebraicMatrix55 errorMatrixToRescale;
  if ( RPar.iMatrixView == 1 ) {
    errorMatrixToRescale = (jacoInvPToInvPt*oldErrorMatrix)*jacoInvPToInvPtTransposed;
    if ( myDebug_ ) {
      std::cout << ">>>>>>>>>>>>>>> ErrorMatrix in Pt before rescaling:" << std::endl;
      std::cout << 100.*errorMatrixToRescale << std::endl;
    }
  } else if ( RPar.iMatrixView == 2) {
    errorMatrixToRescale = (jacoInvPToInvPz*oldErrorMatrix)*jacoInvPToInvPzTransposed;
    if ( myDebug_ ) {
      std::cout << ">>>>>>>>>>>>>>> ErrorMatrix in Pz before rescaling:" << std::endl;
      std::cout << 100.*errorMatrixToRescale << std::endl;
    }
  } else {
    errorMatrixToRescale = oldErrorMatrix;
  }

  
  //Let's rescale only the phi and dxy component
  if ( RPar.iQp )  { for (int j=0;j<5;j++) { errorMatrixToRescale(0,j)*=kFactor; errorMatrixToRescale(j,0)*=kFactor; };}
  if ( RPar.iLam ) { for (int j=0;j<5;j++) { errorMatrixToRescale(1,j)*=kFactor; errorMatrixToRescale(j,1)*=kFactor; };}
  if ( RPar.iPhi ) { for (int j=0;j<5;j++) { errorMatrixToRescale(2,j)*=kFactor; errorMatrixToRescale(j,2)*=kFactor; };}
  if ( RPar.iDxy ) { for (int j=0;j<5;j++) { errorMatrixToRescale(3,j)*=kFactor; errorMatrixToRescale(j,3)*=kFactor; };}
  if ( RPar.iDsz ) { for (int j=0;j<5;j++) { errorMatrixToRescale(4,j)*=kFactor; errorMatrixToRescale(j,4)*=kFactor; };}

  AlgebraicMatrix55 newErrorMatrix;

  if ( RPar.iMatrixView == 1 ) {
    newErrorMatrix = (jacoInvPtToInvP*errorMatrixToRescale)*jacoInvPtToInvPTransposed;
    if ( myDebug_ ) {
      std::cout << ">>>>>>>>>>>>>>> ErrorMatrix in Pt **AFTER** rescaling:" << std::endl;
      std::cout << 100.*newErrorMatrix << std::endl;
    }
  } else if ( RPar.iMatrixView == 2) {
    newErrorMatrix = (jacoInvPzToInvP*errorMatrixToRescale)*jacoInvPzToInvPTransposed;
    if ( myDebug_ ) {
      std::cout << ">>>>>>>>>>>>>>> ErrorMatrix in Pz **AFTER** rescaling:" << std::endl;
      std::cout << 100.*newErrorMatrix << std::endl;
    }
  } else {
    newErrorMatrix = errorMatrixToRescale;
  }

  if ( myDebug_ ) {
    std::cout << ">>>>>>>>>>>>>>> ErrorMatrix in P **AFTER** rescaling:" << std::endl;
    std::cout << 100.*newErrorMatrix << std::endl;
  }

  AlgebraicSymMatrix55 newSymErrorMatrix = newErrorMatrix.LowerBlock();
  
  if ( myDebug_ ) {
    std::cout << ">>>>>>>>>>>>>>> returned Sym ErrorMatrix in P **AFTER** rescaling (for cross check):" << std::endl;
    std::cout << 100.*newSymErrorMatrix << std::endl;
  }

  return newSymErrorMatrix;

}

/*
AlgebraicSymMatrix55 SplitRefit::rescaleErrorOfLongComponent(TrajectoryStateOnSurface& TSOS, double kFactor){

  AlgebraicSymMatrix55 errorMatrixToRescale = TSOS.curvilinearError().matrix();

  if ( myDebug_ ) {
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << ">> ***[Matrices multiplied by 100. to improve readibility]*** >>>>>>>>>>" << std::endl;
    std::cout << ">>>>>>>>>>>>>>> ErrorMatrix before long rescaling:" << std::endl;
    std::cout << 100.*errorMatrixToRescale << std::endl;
  }
  
  //Let's rescale only the lambda and dsz component
  for (int j=0;j<5;j++) { errorMatrixToRescale(1,j)*=kFactor; errorMatrixToRescale(j,1)*=kFactor; }
  for (int j=0;j<5;j++) { errorMatrixToRescale(4,j)*=kFactor; errorMatrixToRescale(j,4)*=kFactor; }

  if ( myDebug_ ) {
    std::cout << ">>>>>>>>>>>>>>> ErrorMatrix **AFTER** long rescaling:" << std::endl;
    std::cout << 100.*errorMatrixToRescale << std::endl;
  }

  AlgebraicSymMatrix55 newSymErrorMatrix = errorMatrixToRescale.LowerBlock();
  
  if ( myDebug_ ) {
    std::cout << ">>>>>>>>>>>>>>> returned Sym ErrorMatrix **AFTER** long rescaling (for cross check):" << std::endl;
    std::cout << 100.*newSymErrorMatrix << std::endl;
  }

  return newSymErrorMatrix;

}
*/

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
      //      if ( myDebug_ ) std::cout << buffer;
    }
  }
  
  //  edm::LogInfo("ELFTrackProducer::setMaterialToKFactor") << outputBlock.data();
}
  
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

void SplitRefit::dumpModuleInfo(DetId hitId){
  
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

void SplitRefit::dumpTSOSInfo(TrajectoryStateOnSurface& TSOS){

  std::cout << ">>TSOS>> p=" << TSOS.globalMomentum().mag() << "+-" << pErrorAtTSOS(TSOS) <<
    " pt=" << TSOS.globalMomentum().perp() << "+-" << ptErrorAtTSOS(TSOS) <<
    " pz=" << TSOS.globalMomentum().z() << "+-" << pzErrorAtTSOS(TSOS) <<
    " l=" << 0.5*M_PI-TSOS.globalMomentum().theta() << "+-" << lambdaErrorAtTSOS(TSOS) <<
    std::endl;

}
