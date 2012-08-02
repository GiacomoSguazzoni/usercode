//
// Original Author:  Giuseppe Cerati
//         Created:  Fri Aug  7 15:10:58 CEST 2009
// $Id: SplitReco.cc,v 1.1 2009/10/14 16:02:21 sguazz Exp $
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
  pullCut_ = iConfig.getUntrackedParameter<double>("pullCut", 7.);

  ptMinCut_ = iConfig.getUntrackedParameter<double>("ptMinCut", 0.8);
  ptMaxCut_ = iConfig.getUntrackedParameter<double>("ptMaxCut", 5.2);
  nHitMinCut_ = iConfig.getUntrackedParameter<int>("nHitMinCut", 6);
  etaMaxCut_ = iConfig.getUntrackedParameter<double>("etaMaxCut", 2.0);
  
  fitterName_ = iConfig.getUntrackedParameter<std::string>("Fitter","KFFittingSmootherWithOutliersRejectionAndRK");
  associatorName_ = iConfig.getUntrackedParameter<std::string>("Associator","TrackAssociatorByHits");
  builderName_ = iConfig.getUntrackedParameter<std::string>("Builder","WithAngleAndTemplate");   

  //Tags
  tracksTag_ = iConfig.getUntrackedParameter<edm::InputTag>("tracks");
  if ( isMC_ ) {
    tpTag_ = iConfig.getUntrackedParameter< edm::InputTag >("tp");
    simtracksTag_ = iConfig.getUntrackedParameter< edm::InputTag >("simtracks");
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
  splitTree->Branch("nHit",       &sT_nHit,       "nHit/I",       bs);    
  splitTree->Branch("nEffHit",    &sT_nEffHit,    "nEffHit/I",    bs);    
  splitTree->Branch("q",          &sT_q,          "q/I",          bs);    
  splitTree->Branch("dxy",        &sT_dxy,        "dxy/F",        bs);     
  splitTree->Branch("dz",         &sT_dz,         "dz/F",         bs);    
  splitTree->Branch("vx",         &sT_vx,         "vx/F",         bs);     
  splitTree->Branch("vy",         &sT_vy,         "vy/F",         bs);    
  splitTree->Branch("vz",         &sT_vz,         "vz/F",         bs);    
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
  // split track method 0 (not overlapping splits)
  splitTree->Branch("pLossSplit", &sT_pLossSplit, "pLossSplit/F", bs);
  splitTree->Branch("pLossErrSplit", &sT_pLossErrSplit, "pLossSplitErr/F", bs);
  splitTree->Branch("p0Split",    &sT_p0Split,         "p0Split/F", bs);
  splitTree->Branch("p1Split",    &sT_p1Split,         "p1Split/F", bs);
  splitTree->Branch("chi2Split",  &sT_chi2Split,  "chi2Split/F",  bs);
  splitTree->Branch("freeParSplit", &sT_freeParSplit,     "freeParSplit/I",     bs);  
  splitTree->Branch("NSplit",     &sT_NSplit,     "NSplit/I",     bs);  
  splitTree->Branch("TInSplit",   &sT_TInSplit,   "TInSplit/F",   bs);
  splitTree->Branch("TOutSplit",  &sT_TOutSplit,  "TOutSplit/F",  bs);
  // split track method 1 (overlapping splits) super split   		   		   	       
  splitTree->Branch("pLossSSplit", &sT_pLossSSplit, "pLossSSplit/F", bs);
  splitTree->Branch("pLossErrSSplit", &sT_pLossErrSSplit, "pLossErrSSplit/F", bs);
  splitTree->Branch("p0SSplit",    &sT_p0SSplit,         "p0SSplit/F", bs);
  splitTree->Branch("p1SSplit",    &sT_p1SSplit,         "p1SSplit/F", bs);
  splitTree->Branch("chi2SSplit",  &sT_chi2SSplit,  "chi2SSplit/F",  bs);
  splitTree->Branch("freeParSSplit", &sT_freeParSSplit,     "freeParSSplit/I",     bs);  
  splitTree->Branch("NSSplit",     &sT_NSSplit,     "NSSplit/I",     bs);  
  splitTree->Branch("TInSSplit",   &sT_TInSSplit,   "TInSSplit/F",   bs);
  splitTree->Branch("TOutSSplit",  &sT_TOutSSplit,  "TOutSSplit/F",  bs);
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
    edm::Handle<SimTrackContainer> simTrackCollection;
    iEvent.getByLabel(simtracksTag_, simTrackCollection);
    const SimTrackContainer simTC = *(simTrackCollection.product());
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

    if ( isMC_ ) {

      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = RtSC[itTrack];
	if ( myDebug_ ) std::cout << "Reco Track pT: " << itTrack->pt() 
				  << " matched to " << tp.size() << " MC Tracks" << std::endl;
	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {

	  TrackingParticleRef tpr = it->first;
	  double assocChi2 = it->second;
	  if ( myDebug_ ) std::cout << "\t\tMCTrack " << tpr.index() << " pT: " << tpr->pt() << 
	    " NShared: " << assocChi2 << std::endl;
	  trackingParticleAction(tpr);

	}
      } catch (cms::Exception event) {
	std::cout << "->   Track pT: " << itTrack->pt() 
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

  TrackingParticle*   tp=const_cast<TrackingParticle*>(tpr.get());
	   
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
	if ( myDebug_ ) std::cout << itphit << " p=" << TPhit->pabs() << " r=" << gpos.perp() << std::endl;
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
  sT_nHit  = itTrack.numberOfValidHits();

  //
  // Track preselection
  if (
      fabs(sT_eta)>etaMaxCut_ || 
      sT_nHit<nHitMinCut_ ||
      sT_pt<ptMinCut_ ||
      sT_pt>ptMaxCut_
      ) return false;

  return true;

}

int SplitReco::trackAction(const reco::Track & itTrack){

  SplitRefit thisSplitRefit(theG.product(), theMF.product(), theFitter.product(), theBuilder.product());
  myTrack mytrack = thisSplitRefit.initializeWithTrack(itTrack);

  //
  // Tree stuff
  sT_pAve   =mytrack.pAve;
  sT_pAveErr=mytrack.pAveErr;
  sT_nEffHit=mytrack.nEffHit;
  sT_T      =mytrack.T;
  sT_rIn    =mytrack.rIn;
  sT_rOut   =mytrack.rOut;
  sT_zIn    =mytrack.zIn;
  sT_zOut   =mytrack.zOut;

  //do splits refits
  std::vector<Split> splits = thisSplitRefit.doSplitRefits(splitTrackEffHits_, splitTrackEffHits_, kFactor_, pullCut_);
  energyLoss eLoss = thisSplitRefit.energyLossFitOnSplits(splits);
  
  sT_pLossSplit=eLoss.pLoss;
  sT_pLossErrSplit=eLoss.pLossErr;
  sT_p0Split=eLoss.p0;
  sT_p1Split=eLoss.p1;
  sT_chi2Split=eLoss.chi2;
  sT_freeParSplit=eLoss.freePar;
  sT_NSplit=eLoss.nsplits;
  sT_TInSplit=eLoss.TIn;
  sT_TOutSplit=eLoss.TOut;

  //do super splits refits
  std::vector<Split> splits2 = thisSplitRefit.doSplitRefits(splitTrackEffHits_, 1, kFactor_, pullCut_);
  energyLoss eLoss2 = thisSplitRefit.energyLossFitOnSplits(splits2);
  
  sT_pLossSSplit=eLoss2.pLoss;
  sT_pLossErrSSplit=eLoss2.pLossErr;
  sT_p0SSplit=eLoss2.p0;
  sT_p1SSplit=eLoss2.p1;
  sT_chi2SSplit=eLoss2.chi2;
  sT_freeParSSplit=eLoss2.freePar;
  sT_NSSplit=eLoss2.nsplits;
  sT_TInSSplit=eLoss2.TIn;
  sT_TOutSSplit=eLoss2.TOut;

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