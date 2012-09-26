//
// Original Author:  Giuseppe Cerati
//         Created:  Fri Aug  7 15:10:58 CEST 2009
// $Id: SplitReco.cc,v 1.3 2012/09/18 15:17:31 sguazz Exp $
//
//
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/SimTracker/TrackAssociation/test/testTrackAssociator.cc?revision=1.17&view=markup&pathrev=CMSSW_2_2_10
//

#include "Tests/SplitReco/interface/SplitReco.h"
#include "Tests/SplitReco/interface/SplitRefit.h"

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
SplitReco::SplitReco(const edm::ParameterSet& iConfig)
{
  //Miscellanea
  std::string FileName_ = iConfig.getUntrackedParameter<std::string>("outputFileName", "SplitTrack.root");
  myDebug_ = iConfig.getUntrackedParameter<bool>("myDebug", false);
  isMC_ = iConfig.getUntrackedParameter<bool>("isMC", true);
  splitTrackEffHits_ = iConfig.getUntrackedParameter<int>("splitTrackEffHits", 4);
  minSplits_ = iConfig.getUntrackedParameter<int>("minSplits", 4);
  kFactor_ = iConfig.getUntrackedParameter<double>("kFactor", 1.);

  superSplit_ = iConfig.getUntrackedParameter<bool>("superSplit", false);

  barrelRPar.iSpecial = iConfig.getUntrackedParameter<bool>("specialRescale", false);
  barrelRPar.iMatrixView = iConfig.getUntrackedParameter<int>("matrixView", 0);
  barrelRPar.iQp = iConfig.getUntrackedParameter<bool>("qOverPComp", false);
  barrelRPar.iLam = iConfig.getUntrackedParameter<bool>("lambdaComp", false);
  barrelRPar.iPhi = iConfig.getUntrackedParameter<bool>("phiComp", false);
  barrelRPar.iDxy = iConfig.getUntrackedParameter<bool>("dxyComp", false);
  barrelRPar.iDsz = iConfig.getUntrackedParameter<bool>("dszComp", false);

  endcapRPar.iSpecial = iConfig.getUntrackedParameter<bool>("EndCapSpecialRescale", false);
  endcapRPar.iMatrixView = iConfig.getUntrackedParameter<int>("EndCapMatrixView", 0);
  endcapRPar.iQp = iConfig.getUntrackedParameter<bool>("EndCapQOverPComp", false);
  endcapRPar.iLam = iConfig.getUntrackedParameter<bool>("EndCapSLambdaComp", false);
  endcapRPar.iPhi = iConfig.getUntrackedParameter<bool>("EndCapPhiComp", false);
  endcapRPar.iDxy = iConfig.getUntrackedParameter<bool>("EndCapDxyComp", false);
  endcapRPar.iDsz = iConfig.getUntrackedParameter<bool>("EndCapDszComp", false);

  absEtaBarrelEndcapCut_ = iConfig.getUntrackedParameter<double>("etaBarrelEndcapCut", 0.9);

  ptMinCut_ = iConfig.getUntrackedParameter<double>("ptMinCut", 0.8);
  ptMaxCut_ = iConfig.getUntrackedParameter<double>("ptMaxCut", 5.2);
  nHitMinCut_ = iConfig.getUntrackedParameter<int>("nHitMinCut", 6);
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
  sT_iEvent=0;
  itrack=0;

  file = new TFile(FileName_.c_str(),"recreate");
  
  /*
  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);

  //
  // Summary histos
  deltaPHist = new TH1F("deltaP","deltaP",200,-0.5,0.5);
  deltaPSimHist = new TH1F("deltaPSim","deltaPSim",200,-0.5,0.5);
  pullHist = new TH1F("pull","pull",200,-20.,20.);
  TH1::AddDirectory(oldAddDir);
  */

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
  splitTree->Branch("nHitVal",    &sT_nHitVal,    "nHitVal/I",    bs);    
  splitTree->Branch("nHitMis",    &sT_nHitMis,    "nHitMis/I",    bs);    
  splitTree->Branch("nHitIna",    &sT_nHitIna,    "nHitIna/I",    bs);    
  splitTree->Branch("nHitBad",    &sT_nHitBad,    "nHitBad/I",    bs);    
  splitTree->Branch("nEffHit",    &sT_nEffHit,    "nEffHit/I",    bs);    
  splitTree->Branch("q",          &sT_q,          "q/I",          bs);    
  splitTree->Branch("dxy",        &sT_dxy,        "dxy/F",        bs);     
  splitTree->Branch("dz",         &sT_dz,         "dz/F",         bs);    
  splitTree->Branch("vx",         &sT_vx,         "vx/F",         bs);     
  splitTree->Branch("vy",         &sT_vy,         "vy/F",         bs);    
  splitTree->Branch("vz",         &sT_vz,         "vz/F",         bs);    
  splitTree->Branch("chi2",       &sT_chi2,       "chi2/F",       bs);    
  splitTree->Branch("maxchi2",    &sT_maxchi2,    "maxchi2/F",    bs);    
  splitTree->Branch("T",          &sT_T,          "T/F",          bs);       
  splitTree->Branch("rIn",        &sT_rIn,        "rIn/F",        bs); //r first hit
  splitTree->Branch("rOut",       &sT_rOut,       "rOut/F",       bs); //r last hit
  splitTree->Branch("xIn",        &sT_xIn,        "xIn/F",        bs); //x first hit
  splitTree->Branch("yIn",        &sT_yIn,        "yIn/F",        bs); //y first hit
  splitTree->Branch("zIn",        &sT_zIn,        "zIn/F",        bs); //z first hit
  splitTree->Branch("xOut",       &sT_xOut,       "xOut/F",       bs); //x last hit
  splitTree->Branch("yOut",       &sT_yOut,       "yOut/F",       bs); //y last hit
  splitTree->Branch("zOut",       &sT_zOut,       "zOut/F",       bs); //z last hit
  // simtrack	     		   		   	       
  if ( isMC_ ){
    splitTree->Branch("iTrackSim",  &sT_iTrackSim,  "iTrackSim/I",  bs);    
    splitTree->Branch("pSim",       &sT_pSim,       "pSim/F",       bs);    
    splitTree->Branch("ptSim",      &sT_ptSim,      "ptSim/F",      bs);   
    splitTree->Branch("etaSim",     &sT_etaSim,     "etaSim/F",     bs);  
    splitTree->Branch("theSim",     &sT_theSim,     "theSim/F",     bs);  
    splitTree->Branch("phiSim",     &sT_phiSim,     "phiSim/F",     bs);  
    splitTree->Branch("pLossSim",   &sT_pLossSim,   "pLossSim/F",   bs);
    splitTree->Branch("dpdxSim",    &sT_dpdxSim,    "dpdxSim/F",    bs);
    splitTree->Branch("nHitSim",    &sT_nHitSim,    "nHitSim/I",    bs);    
    splitTree->Branch("hitFrac",    &sT_hitFrac,    "hitFrac/F",    bs);    
    splitTree->Branch("TSim",       &sT_TSim,       "TSim/F",       bs); //r first hit
    splitTree->Branch("rInSim",     &sT_rInSim,     "rInSim/F",     bs); //r first hit
    splitTree->Branch("rOutSim",    &sT_rOutSim,    "rOutSim/F",    bs); //r last hit
    splitTree->Branch("xInSim",     &sT_xInSim,     "xInSim/F",     bs); //x first hit
    splitTree->Branch("yInSim",     &sT_yInSim,     "yInSim/F",     bs); //y first hit
    splitTree->Branch("zInSim",     &sT_zInSim,     "zInSim/F",     bs); //z first hit
    splitTree->Branch("xOutSim",    &sT_xOutSim,    "xOutSim/F",    bs); //x last hit
    splitTree->Branch("yOutSim",    &sT_yOutSim,    "yOutSim/F",    bs); //y last hit
    splitTree->Branch("zOutSim",    &sT_zOutSim,    "zOutSim/F",    bs); //z last hit
  }
  // split track method 0 (not overlapping splits)
  splitTree->Branch("NSplit",     &sT_NSplit,     "NSplit/I",     bs);  
  splitTree->Branch("TInSplit",   &sT_TInSplit,   "TInSplit/F",   bs);
  splitTree->Branch("TOutSplit",  &sT_TOutSplit,  "TOutSplit/F",  bs);

  splitTree->Branch("dpdxSplit", &sT_dpdxSplit, "dpdxSplit/F", bs);
  splitTree->Branch("dpdxErrSplit", &sT_dpdxErrSplit, "dpdxSplitErr/F", bs);
  splitTree->Branch("chi2Split",  &sT_chi2Split,  "chi2Split/F",  bs);
  splitTree->Branch("freeParSplit", &sT_freeParSplit,     "freeParSplit/I",     bs);  

  splitTree->Branch("dpdxTSplit", &sT_dpdxTSplit, "dpdxTSplit/F", bs);
  splitTree->Branch("dpdxTErrSplit", &sT_dpdxTErrSplit, "dpdxTSplitErr/F", bs);
  splitTree->Branch("chi2TSplit",  &sT_chi2TSplit,  "chi2TSplit/F",  bs);
  splitTree->Branch("freeParTSplit", &sT_freeParTSplit,     "freeParTSplit/I",     bs);  

  splitTree->Branch("dpdxZSplit", &sT_dpdxZSplit, "dpdxZSplit/F", bs);
  splitTree->Branch("dpdxZErrSplit", &sT_dpdxZErrSplit, "dpdxZSplitErr/F", bs);
  splitTree->Branch("chi2ZSplit",  &sT_chi2ZSplit,  "chi2ZSplit/F",  bs);
  splitTree->Branch("freeParZSplit", &sT_freeParZSplit, "freeParZSplit/I",     bs);  

  // split track method 1 (overlapping splits) super split   		   		   	       
  if ( superSplit_ ){
    splitTree->Branch("NSSplit",     &sT_NSSplit,     "NSSplit/I",     bs);  
    
    splitTree->Branch("dpdxSSplit", &sT_dpdxSSplit, "dpdxSSplit/F", bs);
    splitTree->Branch("dpdxErrSSplit", &sT_dpdxErrSSplit, "dpdxErrSSplit/F", bs);
    splitTree->Branch("chi2SSplit",  &sT_chi2SSplit,  "chi2SSplit/F",  bs);
    splitTree->Branch("freeParSSplit", &sT_freeParSSplit,     "freeParSSplit/I",     bs);  
    
    splitTree->Branch("dpdxTSSplit", &sT_dpdxTSSplit, "dpdxTSSplit/F", bs);
    splitTree->Branch("dpdxTErrSSplit", &sT_dpdxTErrSSplit, "dpdxTSSplitErr/F", bs);
    splitTree->Branch("chi2TSSplit",  &sT_chi2TSSplit,  "chi2TSSplit/F",  bs);
    splitTree->Branch("freeParTSSplit", &sT_freeParTSSplit, "freeParTSSplit/I",     bs);  
    
    splitTree->Branch("dpdxZSSplit", &sT_dpdxZSSplit, "dpdxZSSplit/F", bs);
    splitTree->Branch("dpdxZErrSSplit", &sT_dpdxZErrSSplit, "dpdxZSSplitErr/F", bs);
    splitTree->Branch("chi2ZSSplit",  &sT_chi2ZSSplit,  "chi2ZSSplit/F",  bs);
    splitTree->Branch("freeParZSSplit", &sT_freeParZSSplit, "freeParZSSplit/I",     bs);  
  }

#ifdef extra
  splitTree->Branch("pLossSplit", &sT_pLossSplit, "pLossSplit/F", bs);
  splitTree->Branch("pLossErrSplit", &sT_pLossErrSplit, "pLossSplitErr/F", bs);
  splitTree->Branch("p0Split",    &sT_p0Split,         "p0Split/F", bs);
  if ( superSplit_ ){
    splitTree->Branch("pLossSSplit", &sT_pLossSSplit, "pLossSSplit/F", bs);
    splitTree->Branch("pLossErrSSplit", &sT_pLossErrSSplit, "pLossErrSSplit/F", bs);
    splitTree->Branch("p0SSplit",    &sT_p0SSplit,         "p0SSplit/F", bs);
  }
#endif

}


SplitReco::~SplitReco()
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
SplitReco::analyze(const edm::Event& iEvent, const edm::EventSetup& setup)
{
  using namespace edm;
  using namespace std;

  using reco::Track;
  using reco::TrackCollection;

  // Store run number
  //
  sT_nRun = iEvent.run();

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


  sT_nEvent = iEvent.id().event();
  sT_iEvent++;

  //
  // Loop on event tracks
  sT_iTrack=0;
  sT_NTrack = tC.size();
  for(View<Track>::size_type i=0; i<tC.size(); ++i) {
    edm::RefToBase<reco::Track> itTrack(trackCollectionH, i);

   //
    // Track preselection
    if ( ! trackPreSelection(*itTrack) ) continue;

    sT_iTrack++;
    
    sT_iTrackSim=0;
    sT_pSim=0.;     
    sT_ptSim=0.;    
    sT_etaSim=0.;   
    sT_theSim=0.;   
    sT_phiSim=0.;   
    sT_pLossSim=0.; 
    sT_hitFrac=0.;

    if ( isMC_ ) {

      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = RtSC[itTrack];
	if ( myDebug_ ) std::cout << "Reco Track pT: " << itTrack->pt() 
				  << " matched to " << tp.size() << " MC Tracks" << std::endl;
	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  
	  TrackingParticleRef tpr = it->first;
	  sT_hitFrac = it->second;
	  if ( myDebug_ ) std::cout << "\t\tMCTrack " << tpr.index() << " pT: " << tpr->pt() << 
	    " NShared: " << sT_hitFrac << std::endl;
	  trackingParticleAction(tpr);
	  
	}
      } catch (cms::Exception event) {
	if ( myDebug_ ) std::cout << "->   Track pT: " << itTrack->pt() 
				  <<  " matched to 0  MC Tracks" << std::endl;
      }



    }

    trackAction(*itTrack);

  }
}


// ------------ method called once each job just before starting event loop  ------------
void SplitReco::beginJob(const edm::EventSetup & setup)
{

}

// ------------ method called once each job just after ending the event loop  ------------
void SplitReco::endJob() {

}

int SplitReco::trackingParticleAction(TrackingParticleRef & tpr){

  TrackingParticle* tp = const_cast<TrackingParticle*>(tpr.get());
	   
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
  //FIX ME
  double x0=0.;
  double y0=0.;
  double z0=0.;
  double T = 0.;
  std::vector<double> vecTPHitPerpR;
  std::vector<double> vecTPHitX;
  std::vector<double> vecTPHitY;
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
	vecTPHitX.push_back(x);
	vecTPHitY.push_back(y);
	vecTPHitZ.push_back(z);
	if ( myDebug_ ) std::cout << itphit << " TPHit p=" << TPhit->pabs() << " r=" << gpos.perp() << std::endl;
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
  sT_dpdxSim = sT_pLossSim/sT_TSim;
  sT_rInSim = vecTPHitPerpR.at(0);
  sT_rOutSim = vecTPHitPerpR.at(itphit-1);
  sT_xInSim = vecTPHitX.at(0);
  sT_yInSim = vecTPHitY.at(0);
  sT_zInSim = vecTPHitZ.at(0);
  sT_xOutSim = vecTPHitX.at(itphit-1);
  sT_yOutSim = vecTPHitY.at(itphit-1);
  sT_zOutSim = vecTPHitZ.at(itphit-1);
  if ( myDebug_ ) std::cout << " SimTrack p:" << tp->p() << " eta:" << tp->momentum().eta() << " pout:" << pSimOut << " pLoss:" << sT_pLossSim << std::endl;
  //

  return 0;

}

bool SplitReco::trackPreSelection(const reco::Track & itTrack){

  sT_p     = itTrack.p();
  sT_pErr  = pErrorOfTrack(itTrack);
  sT_pt    = itTrack.pt();
  sT_ptErr = itTrack.ptError();
  sT_eta   = itTrack.eta();
  sT_the   = itTrack.theta();
  sT_phi   = itTrack.phi();
  sT_dxy   = itTrack.dxy();     
  sT_dz    = itTrack.dz();    
  sT_vx    = itTrack.vx();     
  sT_vy    = itTrack.vy();     
  sT_vz    = itTrack.vz();    
  sT_q     = itTrack.charge();
  sT_chi2  = itTrack.normalizedChi2();
  sT_nHitVal = itTrack.numberOfValidHits();

  //
  // Track preselection
  if (
      sT_eta<etaMinCut_ || 
      sT_eta>etaMaxCut_ || 
      sT_nHitVal<nHitMinCut_ ||
      sT_pt<ptMinCut_ ||
      sT_pt>ptMaxCut_
      ) return false;

  //
  // hit pattern of the track
  const reco::HitPattern& hp = itTrack.hitPattern();
  sT_nHitMis = hp.numberOfLostHits();
  sT_nHitIna = hp.numberOfInactiveHits();
  sT_nHitBad = hp.numberOfBadHits();

  return true;

}

int SplitReco::trackAction(const reco::Track & itTrack){

  SplitRefit thisSplitRefit(theG.product(), theMF.product(), theFitter.product(), theBuilder.product(), myDebug_);
  myTrack mytrack = thisSplitRefit.initializeWithTrack(itTrack);

  //
  // Tree stuff
  sT_pAve   =mytrack.pAve;
  sT_pAveErr=mytrack.pAveErr;
  sT_nEffHit=mytrack.nEffHit;
  sT_T      =mytrack.T;
  sT_rIn    =mytrack.rIn;
  sT_rOut   =mytrack.rOut;
  sT_xIn    =mytrack.xIn;
  sT_yIn    =mytrack.yIn;
  sT_zIn    =mytrack.zIn;
  sT_xOut   =mytrack.xOut;
  sT_yOut   =mytrack.yOut;
  sT_zOut   =mytrack.zOut;
  sT_maxchi2 = mytrack.maxhitchi2;


  /*
  //'offline' track selection for debugging purposes
  int icut = 0;
  if (
      //ripulire
      sT_rIn < 8. 
      &&
      ((sT_rOut>100.)||(sT_zOut>260.)||(sT_zOut<-260.))
      &&
      (sT_nHitIna+sT_nHitMis+sT_nHitBad)==0
      &&
      //To better identify the direction
      sT_dz<5.
      &&
      sT_dz>-5.
      ) icut = 1;
  
  if ( icut ) return 1;
  */

  //do splits refits
  if ( myDebug_ ) std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%% Split refit " << std::endl; 
  std::vector<Split> splits = thisSplitRefit.doSplitRefits(splitTrackEffHits_, splitTrackEffHits_, kFactor_, absEtaBarrelEndcapCut_, 
							   barrelRPar, endcapRPar);
  std::vector<energyLoss> eLosses = thisSplitRefit.energyLossFitOnSplits(splits);
  if ( myDebug_ ) std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%% Split refit done" << std::endl; 
  
  sT_NSplit=eLosses.at(0).nsplits;
  sT_TInSplit=eLosses.at(0).TIn;
  sT_TOutSplit=eLosses.at(0).TOut;

  sT_dpdxSplit=eLosses.at(0).dpdx;
  sT_dpdxErrSplit=eLosses.at(0).dpdxErr;
  sT_chi2Split=eLosses.at(0).chi2;
  sT_freeParSplit=eLosses.at(0).freePar;

  sT_dpdxTSplit=eLosses.at(1).dpdx;
  sT_dpdxTErrSplit=eLosses.at(1).dpdxErr;
  sT_chi2TSplit=eLosses.at(1).chi2;
  sT_freeParTSplit=eLosses.at(1).freePar;

  sT_dpdxZSplit=eLosses.at(2).dpdx;
  sT_dpdxZErrSplit=eLosses.at(2).dpdxErr;
  sT_chi2ZSplit=eLosses.at(2).chi2;
  sT_freeParZSplit=eLosses.at(2).freePar;

  //do super splits refits
  if ( superSplit_ ){
    
    if ( myDebug_ ) std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%% Super Split refit " << std::endl; 
    std::vector<Split> splits2 = thisSplitRefit.doSplitRefits(splitTrackEffHits_, 1, kFactor_, absEtaBarrelEndcapCut_,
							      barrelRPar, endcapRPar);
    std::vector<energyLoss> eLosses2 = thisSplitRefit.energyLossFitOnSplits(splits2);
    if ( myDebug_ ) std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%% Super Split refit done " << std::endl; 
    
    sT_NSSplit=eLosses2.at(0).nsplits;
    
    sT_dpdxSSplit=eLosses2.at(0).dpdx;
    sT_dpdxErrSSplit=eLosses2.at(0).dpdxErr;
    sT_chi2SSplit=eLosses2.at(0).chi2;
    sT_freeParSSplit=eLosses2.at(0).freePar;
    
    sT_dpdxTSSplit=eLosses2.at(1).dpdx;
    sT_dpdxTErrSSplit=eLosses2.at(1).dpdxErr;
    sT_chi2TSSplit=eLosses2.at(1).chi2;
    sT_freeParTSSplit=eLosses2.at(1).freePar;
    
    sT_dpdxZSSplit=eLosses2.at(2).dpdx;
    sT_dpdxZErrSSplit=eLosses2.at(2).dpdxErr;
    sT_chi2ZSSplit=eLosses2.at(2).chi2;
    sT_freeParZSSplit=eLosses2.at(2).freePar;
  }

  //
  // Extras
#ifdef extra
  sT_pLossSplit=eLosses.at(0).pLoss;
  sT_pLossErrSplit=eLosses.at(0).pLossErr;
  sT_p0Split=eLosses.at(0).p0;
  if ( superSplit_ ){
    sT_pLossSSplit=eLosses2.at(0).pLoss;
    sT_pLossErrSSplit=eLosses2.at(0).pLossErr;
    sT_p0SSplit=eLosses.at(0).p0;
  }
#endif


  //
  // Fill splitTree
  splitTree->Fill();
  
  return 0;

}

double SplitReco::pErrorOfTrack(const reco::Track & track){
  
  double perror = track.charge()*track.p()*track.qoverpError()/track.qoverp();
  
  return perror;
  
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SplitReco);
