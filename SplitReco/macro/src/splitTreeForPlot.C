#include "splitTreeForPlot.h"

splitTreeForPlot::splitTreeForPlot(TTree *tree, int evMax, int evMaxNorm)
{
  
  std::cout << " Istanzio splitTreeForPlot tramite un tree..." << std::endl;

  evRangeMax = evMax;
  evRangeMaxNorm = evMaxNorm;


  Init(tree);

  iMC = false;

  if ( tree->GetBranch( "pLossSim" )) {
    std::cout << " This is a MC file... " << std::endl;
    iMC = true ;
  }

  fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname

   fChain->SetBranchStatus("p", 1);
   fChain->SetBranchStatus("pt", 1);
   fChain->SetBranchStatus("eta", 1);
   fChain->SetBranchStatus("phi", 1);
   fChain->SetBranchStatus("nHitVal", 1);
   fChain->SetBranchStatus("nHitMis", 1);
   fChain->SetBranchStatus("nHitIna", 1);
   fChain->SetBranchStatus("nHitBad", 1);
   fChain->SetBranchStatus("dz", 1);
   fChain->SetBranchStatus("chi2", 1);
   fChain->SetBranchStatus("maxchi2", 1);
   fChain->SetBranchStatus("T", 1);
   fChain->SetBranchStatus("rIn", 1);
   fChain->SetBranchStatus("rOut", 1);
   fChain->SetBranchStatus("xOut", 1);
   fChain->SetBranchStatus("yOut", 1);
   fChain->SetBranchStatus("zOut", 1);
   fChain->SetBranchStatus("pLossSim", 1);
   fChain->SetBranchStatus("dpdxSim", 1);
   fChain->SetBranchStatus("TSim", 1);
   fChain->SetBranchStatus("NSplit", 1);
   fChain->SetBranchStatus("dpdxSplit", 1);
   fChain->SetBranchStatus("dpdxErrSplit", 1);
   fChain->SetBranchStatus("chi2Split", 1);
   fChain->SetBranchStatus("dpdxTSplit", 1);
   fChain->SetBranchStatus("dpdxTErrSplit", 1);
   fChain->SetBranchStatus("chi2TSplit", 1);
   fChain->SetBranchStatus("dpdxZSplit", 1);
   fChain->SetBranchStatus("dpdxZErrSplit", 1);
   fChain->SetBranchStatus("chi2ZSplit", 1);


   ptmin = 0.9;
   ptmax = 1.1;
   fillNormHisto();

}

/*
splitTreeForPlot::splitTreeForPlot(const char* filename)
{
  
  std::cout << " Istanzio splitTreeForPlot..." << std::endl;
  std::cout << " Cerco di aprire questo file: " << filename << std::endl;

  TFile * f = new TFile(filename);
  
  TTree * tree = (TTree*)f->Get("sT");

  Init(tree);


}
*/


splitTreeForPlot::~splitTreeForPlot()
{

  std::cout << " Distruggo splitTreeForPlot..." << std::endl;

  histoPtWei->SaveAs("plots/histoPtWei.root");
  histoPt2->SaveAs("plots/histoPt2.root");

}

void splitTreeForPlot::setDPdxQuantities(){
  
  // Factor 1000 to go in MeV/cm
  //      double dpdxSim = 0.; //FIXME
  //      double dPdxTrue = 1000.*dpdxSim; //FIXME
  //      double dPdxTrue = 1000.*pLossSim/T; //FIXME //due to a bug dpdxSim is wrong; using plossSim/T instead
  dPdxTrue=0;
  if (((fabs(eta)<0.3)&&(TSim<130.))||((fabs(eta)>0.3)&&(TSim<154.545*fabs(eta)+83.636))) dPdxTrue = 1000.*pLossSim/TSim; // temporary TSim cut
  //FIXME //due to a bug dpdxSim is wrong; using plossSim/TSim with a TSim cut instead
  
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
  //Transverse split
  if ( iSplit == 2 ) {
    dPdx = 1000.*dpdxTSplit;
    dPdxErr = 1000.*dpdxTErrSplit;
    dPdxChi2 = chi2TSplit/nDof;
  } 
  //Z split
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
  
}


int splitTreeForPlot::AnalysisCut()
{
  
  if ( QualityCut() ) { 
    
    setDPdxQuantities();
    if ( dPdxChi2 < 4. && dPdxChi2>0. ) return 1;      
    
  }
  
  return 0;

}

int splitTreeForPlot::QualityCut()
{

  //chi2Split is used to check for good fit; failed fits have chi2Split == -1.
  if (
      pt>ptmin
      &&
      pt<ptmax
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

void splitTreeForPlot::fillNormHisto(){

  std::cout << " Setting up normalization infrastructure... " << std::endl;

  int normBins = 50;

  normHisto = new TH1F("norm","norm",normBins,ptmin,ptmax);
  histoPtWei = new TH1F("ptWei","ptWei",normBins,ptmin,ptmax);
  histoPt2 = new TH1F("pt2","pt2",normBins,ptmin,ptmax);
  TH1F * histoPt = new TH1F("pt","pt",normBins,ptmin,ptmax);

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  if ( evRangeMaxNorm > 0 ) nentries = evRangeMaxNorm; 

  //  nentries = 20000; 

  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if ( ! (jentry % 1000000) ) std::cout << ">>>>Normalization loop. Entry: " << jentry << "\n";
    if ( AnalysisCut() ) histoPt->Fill(pt);
  
  }

  histoPt->Scale(1./(histoPt->Integral()));
  histoPt->SaveAs("plots/histoPt.root");


  float factor = 1./normBins;

  for (int i=1; i<normBins+1; i++){

    normHisto->SetBinContent(i, factor/histoPt->GetBinContent(i));

  }


  normHisto->SaveAs("plots/normHisto.root");

  std::cout << " Done. " << std::endl;

}

float splitTreeForPlot::getWeiFromNormHisto(float& val){

  return normHisto->GetBinContent(normHisto->FindBin(val));
 
}

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
    
    if ( ! (jentry % 1000000) ) std::cout << ">>>>Entry: " << jentry << "\n";
    
    //Verify cuts
    if ( AnalysisCut() ) { 

      float wei = getWeiFromNormHisto(pt);
      
      histoPtWei->Fill(pt, wei);
      histoPt2->Fill(pt);
      
      double phiOut =  sguazzNormalizedPhi(atan2(yOut,xOut)); //Bring phi in the range 0 to 2PI
      
      elplot->Fill(dPdxTrue, dPdx, dPdxErr, eta, phiOut, wei); //Add wei
      
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
