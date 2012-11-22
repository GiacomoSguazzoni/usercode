#include "splitTreeForPlot.h"

splitTreeForPlot::splitTreeForPlot(const char* filename)
{
  
  std::cout << " Istanzio splitTreeForPlot..." << std::endl;
  std::cout << " Cerco di aprire questo file: " << filename << std::endl;

  TFile * f = new TFile(filename);
  
  TTree * tree = (TTree*)f->Get("sT");

  Init(tree);

  SetEvRangeMax(-1);

}

splitTreeForPlot::splitTreeForPlot(TTree *tree)
{
  
  std::cout << " Istanzio splitTreeForPlot tramite un tree..." << std::endl;

  Init(tree);

  SetEvRangeMax(-1);

  iMC = false;

  if ( tree->GetBranch( "pLossSim" )) {
    std::cout << " This is a MC file... " << std::endl;
    iMC = true ;
  }

}

splitTreeForPlot::~splitTreeForPlot()
{

  std::cout << " Distruggo splitTreeForPlot..." << std::endl;

}

int splitTreeForPlot::QualityCut()
{

  //chi2Split is used to check for good fit; failed fits have chi2Split == -1.
  if (
      pt>0.9
      &&
      pt<1.1
      //      pt>1.8
      //	 &&
      //	 pt<2.2
      &&
      chi2<2. 
      &&
      maxchi2<4. 
      &&
      rIn < 8. 
      &&
      ((rOut>100.)||(zOut>260.)||(zOut<-260.))
      ///////////////////////////////      &&
      ///////////////////////////////      pLossSim>0.
      &&
      (nHitIna+nHitMis+nHitBad)==0
      &&
      //To better identify the direction
      dz<5.
      &&
      dz>-5.
      ) return 1;
  
  //      pt>0.9&&pt<1.1&&chi2<2.&&rIn<8.&&((rOut>100.)||(zOut>260.)||(zOut<-260.))&&pLossSim>0.&&(nHitIna+nHitMis+nHitBad)==0&&dz<5.&&dz>-5.

  return 0;

}

#include <cmath>

void splitTreeForPlot::LoopForFill(elPlot* elplot)
{

  std::cout << " Max number of Tracks: " << evRangeMax << std::endl;
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  if ( evRangeMax > 0 ) nentries = evRangeMax; 

  Long64_t nbytes = 0, nb = 0;

  //Number of accepted events for debug
  Long64_t naccev = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if ( ! (jentry % 100000) ) std::cout << ">>>>Entry: " << jentry << " Event:" << iEvent << "\n";
    
    //Verify cuts
    if ( QualityCut() ) { 
      
      // Factor 1000 to go in MeV/cm
      //      double dpdxSim = 0.; //FIXME
      double dPdxTrue = 1000.*dpdxSim;
      
      double dPdx, dPdxErr, dPdxChi2;

      int iSplit = 1; 

      //1 p, split
      //2 pt, splt
      //3 pz, split
      //4 weighed average of dpdx from pz and pt

      //11 p, split
      //22 pt, splt
      //33 pz, split

      float nDof = NSplit - 2.;

      if ( iSplit == 1 ) {
	dPdx = 1000.*dpdxSplit;
        dPdxErr = 1000.*dpdxErrSplit;
	dPdxChi2 = chi2Split/nDof;
      } 
      if ( iSplit == 2 ) {
	dPdx = 1000.*dpdxTSplit;
        dPdxErr = 1000.*dpdxTErrSplit;
	dPdxChi2 = chi2TSplit/nDof;
      } 
      if ( iSplit == 3 ) {
	dPdx = 1000.*dpdxZSplit;
        dPdxErr = 1000.*dpdxZErrSplit;
	dPdxChi2 = chi2ZSplit/nDof;
      } 
      if ( iSplit == 4 ) {
	float invTErr2 = 1./(dpdxTErrSplit*dpdxTErrSplit);
	float invZErr2 = 1./(dpdxZErrSplit*dpdxZErrSplit);
	dPdxErr = (invTErr2+invZErr2);
	dPdx = 1000.*(invZErr2*dpdxZSplit + invTErr2*dpdxTSplit)/dPdxErr;
	dPdxErr = sqrt(dPdxErr);
	dPdxChi2 = max(chi2ZSplit/nDof, chi2TSplit/nDof);
	if ( chi2ZSplit == 0. ) dPdxChi2 = 0.;
	if ( chi2TSplit == 0. ) dPdxChi2 = 0.;
      } 
      //

#define removeforSPLIT
#ifdef SPLIT
      if ( iSplit == 11 ) {
	dPdx = 1000.*dpdxSSplit;
        dPdxErr = 1000.*dpdxErrSSplit;
	dPdxChi2 = chi2SSplit/nDof;
      } 
      if ( iSplit == 22 ) {
	dPdx = 1000.*dpdxTSSplit;
        dPdxErr = 1000.*dpdxTErrSSplit;
	dPdxChi2 = chi2TSSplit/nDof;
      } 
      if ( iSplit == 33 ) {
	dPdx = 1000.*dpdxZSSplit;
        dPdxErr = 1000.*dpdxZErrSSplit;
	dPdxChi2 = chi2ZSSplit/nDof;
      } 
#endif

      if ( dPdxChi2 < 4. && dPdxChi2>0. ){
	//if ( dPdxChi2>0. ){

	//      double phiOut =  sguazzNormalizedPhi(phi-q*asin(0.5*0.01*rOut*0.2998*3.8/pt)); //Bring phi in the range 0 to 2PI
	double phiOut =  sguazzNormalizedPhi(atan2(yOut,xOut)); //Bring phi in the range 0 to 2PI

	elplot->Fill(dPdxTrue, dPdx, dPdxErr, eta, phiOut);

      }

    }

  }

  std::cout << " Number af accepted candidates in the step " << naccev << std::endl;

}

double splitTreeForPlot::normalizedPhi(double phi) {
   static double const TWO_PI = M_PI * 2;
   while ( phi < -M_PI ) phi += TWO_PI;
   while ( phi >  M_PI ) phi -= TWO_PI;
   return phi;
}

double splitTreeForPlot::sguazzNormalizedPhi(double phi) {
   static double const TWO_PI = M_PI * 2;
   static double const myPhiZero = -0.125*M_PI;
   static double const myTwoPi = M_PI * 2 + myPhiZero;
   while ( phi < myPhiZero ) phi += TWO_PI;
   while ( phi > myTwoPi ) phi -= TWO_PI;
   return phi;
}
