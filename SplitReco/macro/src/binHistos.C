#include "splitRecoFunctions.h"
#include "binHistos.h"

binHistos::binHistos(const char* id, const char* binInfo, double uC, double vC){
  
  TSId.Append(id);
  TString TSBinInfo;
  TSBinInfo.Append(binInfo);
  dPdxSim = new TH1D("dPdxSim"+TSId,"dPdxSim_"+TSId+" "+TSBinInfo+";(dP/dx)_{sim}  [MeV/cm]",79,0.01,0.8); //No entries if dpdxsim = 0
  dPdxMeas = new TH1D("dPdxMeas"+TSId,"dPdxMeas_"+TSId+" "+TSBinInfo+";(dP/dx)_{meas} [MeV/cm]",80,-2.,4.);
  dPdxPull = new TH1D("dPdxPull"+TSId,"dPdxPull_"+TSId+" "+TSBinInfo+";(Meas - Sim)/MeasErr",80,-8.,8.);
  dPdxMeasVSdPdxSim = new TH2D("dPdxMeasVSdPdxSim"+TSId,"dPdxMeasVSdPdxSim_"+TSId+" "+TSBinInfo+";(dP/dx)_{Sim} [MeV/cm];(dP/dx)_{Meas} [MeV/cm]",100,0.,0.8,100,-0.98,3.);

  uBinCenter = uC;
  vBinCenter = vC;

}
  
binHistos::binHistos(const char* id, TFile * file){
  
  TSId.Append(id);

  TString simHistoName("dPdxSim"+TSId);
  TString measHistoName("dPdxMeas"+TSId);

  dPdxSim = getHistoFromFile(simHistoName, file);
  dPdxMeas = getHistoFromFile(measHistoName, file);

  float umin, umax, vmin, vmax;
  TString formatString("dPdxSim_"+TSId+" u:% 4.2f,% 4.2f v:% 4.2f,% 4.2f");
  TString title=dPdxSim->GetTitle();
  sscanf(title.Data(),"%*s %*2c%f%*c%f%*3c%f%*c%f", &umin, &umax, &vmin, &vmax);
  std::cout << title.Data() << " " << std::endl;
  std::cout << " Bin center u: " << umin << " " << umax << " " << 0.5*(umin+umax) << std::endl;
  std::cout << " Bin center v: " << vmin << " " << vmax << " " << 0.5*(vmin+vmax) << std::endl;

  uBinCenter = 0.5*(umin+umax);
  vBinCenter = 0.5*(vmin+vmax);

}

binHistos::~binHistos()
{

  std::cout << " Cleaning binHistos object..." << std::endl;

}

void binHistos::Fill(double wei, double sim, double meas, double measErr){

  dPdxMeas->Fill(meas, wei);

  if (sim>0.){
    dPdxSim->Fill(sim, wei);
    dPdxPull->Fill((meas-sim)/measErr, wei);
    dPdxMeasVSdPdxSim->Fill(sim,meas,wei);
  }

  /// Qui posso scegliere che valori usare


}

#include <TCanvas.h>

void binHistos::PlotAll()
{

  Plot(dPdxSim);
  Plot(dPdxMeas);
  Plot(dPdxPull);
  Plot(dPdxMeasVSdPdxSim);

}

void binHistos::WriteAll(TFile *file)
{

  Write(dPdxSim, file);
  Write(dPdxMeas, file);
  Write(dPdxPull, file);
  Write(dPdxMeasVSdPdxSim, file);

}

void binHistos::FitAll(double &mpvSim, double &mpvMeas, double &mpvGaus, double &mpvMean, double &mpvMedian, double &mpvTrunMean)
{

  mpvSim = 0.;
  if ( dPdxSim->Integral() ) mpvSim = FitSim(dPdxSim);
  //  FitRoot(dPdxSim);
  mpvGaus = FitGaus(dPdxMeas);
  //  mpvMeas = FitRoot(dPdxMeas);
  mpvMean = dPdxMeas->GetMean();
 
  int nprob = 3;
  double prob[3] = {0.25, 0.50, 0.75};
  double res[3];

  dPdxMeas->GetQuantiles(nprob,res,prob);

  //  mpvMedian = Median(dPdxMeas);
  mpvMedian = res[1];
  mpvTrunMean = meanInRange(dPdxMeas,res[0],res[2]);

  //  FitMeas(dPdxMeas);

}

void binHistos::Plot(TH2* hist)
{

  std::cout << " E ora plotto " << hist->GetName() << std::endl;

  TCanvas * can = new TCanvas(hist->GetName(),"",1000,1000);

  hist->Draw("zcol");

  can->SaveAs("./plots/"+TString(hist->GetName())+".png");
  //  hist->SaveAs("./plots/"+TString(hist->GetName())+".root");
  
}

void binHistos::Write(TH1* hist, TFile* file)
{

  std::cout << " >>>> Writing histo " << hist->GetName() << " to file " << file->GetPath() << std::endl;
  file->cd();
  hist->Write();
  
}

void binHistos::Plot(TH1* hist)
{

  std::cout << " E ora plotto " << hist->GetName() << std::endl;

  TCanvas * can = new TCanvas(hist->GetName(),"",1000,1000);

  hist->Draw();

  can->SaveAs("./plots/"+TString(hist->GetName())+".png");
  //  hist->SaveAs("./plots/"+TString(hist->GetName())+".root");
  
}

TH1D * binHistos::getHistoFromFile(TString name, TFile* file){

  TH1D * histo = (TH1D*)file->Get(name);

  return histo;

}

#include <RooCmdArg.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooLandau.h>
#include <RooFFTConvPdf.h>
#include <RooPlot.h>
#include <TLatex.h>

double binHistos::FitMeas(TH1* hist){

  // Construct observable
  RooRealVar t("t","t", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 

  setMyRange(hist, &t);
  double mpvEst, rmsEst;
  estMyValues(hist, mpvEst, rmsEst);

  // Construct landau(t,ml,sl) ;
  //  RooRealVar ml("ml","mean landau",0.5,-1.0,1.0);
  //  RooRealVar sl("sl","sigma landau",0.05,0.001,1.0);
  RooRealVar ml("ml","mean landau",1.5*mpvEst,0.5*mpvEst,2.0*mpvEst);
  RooRealVar sl("sl","sigma landau",0.01*rmsEst,0.001*rmsEst,0.02*rmsEst);
  RooLandau landau("lx","lx",t,ml,sl);
  
  // Construct gauss(t,mg,sg)
  RooRealVar mg("mg","mg",0.);
  RooRealVar sg("sg","sg",rmsEst,0.1,2.*rmsEst);
  RooGaussian gauss("gauss","gauss",t,mg,sg);

  // C o n s t r u c t   c o n v o l u t i o n   p d f 
  // ---------------------------------------

  // Set #bins to be used for FFT sampling to 10000
  t.setBins(10000,"cache") ; 

  // Construct landau (x) gauss
  RooFFTConvPdf lxg("lxg","landau (X) gauss",t,landau,gauss) ;

  // F i t   a n d   p l o t   c o n v o l u t e d   p d f 
  // ----------------------------------------------------------------------

  RooDataHist data(hist->GetName(),hist->GetTitle(),t,hist); 
  RooPlot* tframe = t.frame(RooFit::Name(hist->GetName()),RooFit::Title(hist->GetTitle())) ;

  //RooPlot* tframe = new RooPlot(hist->GetName(), hist->GetTitle(), t, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), hist->GetNbinsX());

  // Fit gxlx to data
  lxg.fitTo(data);

  // Plot
  data.plotOn(tframe);
  lxg.plotOn(tframe) ;
  //  landau.plotOn(tframe,RooFit::LineStyle(kSolid),RooFit::LineColor(kRed),RooFit::LineWidth(1)) ;

#define oldGetVal

#ifdef oldGetVal
  double mpv = ml.getVal();
  double fwhm = 4.*sl.getVal();
  double sigma = sg.getVal();
#else
  double mpv = ml.getValV();
  double fwhm = 4.*sl.getValV();
  double sigma = sg.getValV();
#endif

  double mpvHi = ml.getErrorHi();
  //  double mpvLo = ml.getErrorLo();

  double fwhmHi = 4.*sl.getErrorHi();
  //  double fwhmLo = 4.*sl.getErrorLo();

  double sigmaHi = sg.getErrorHi();
  //  double sigmaLo = sg.getErrorLo();

  // Draw frame on canvas
  TCanvas * can = new TCanvas("FitCanvas_"+TString(hist->GetName()),"",1000,1000);
  can->SetRightMargin(0.03);
  can->SetLeftMargin(0.13);
  can->SetTopMargin(0.06);
  can->SetBottomMargin(0.10);

  tframe->Draw();
  tframe->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
  tframe->GetYaxis()->SetTitleOffset(1.5);
  tframe->GetXaxis()->SetTitleOffset(1.2);
  
  int iFont = 132;
  double xLabel = 0.935; 
  double yLabel = 0.90;
  double yLabelStep = 0.06;

  char landauMpvLatexChar[100];
  //  sprintf(landauMpvLatexChar,"MPV = %5.4f^{%+5.4f}_{%+5.4f}",mpv,mpvHi,mpvLo);
  sprintf(landauMpvLatexChar,"#font[%d]{MPV = %5.4f #pm %5.4f}",iFont,mpv,mpvHi);
  TLatex landauMpvLatex(xLabel, yLabel, landauMpvLatexChar);
  landauMpvLatex.SetNDC();
  landauMpvLatex.SetTextAlign(33);
  landauMpvLatex.Draw();

  char landauFWHMLatexChar[100];
  //  sprintf(landauFWHMLatexChar,"FWHM = %5.4f^{%+5.4f}_{%+5.4f}",fwhm,fwhmHi,fwhmLo);
  sprintf(landauFWHMLatexChar,"#font[%d]{FWHM = %5.4f #pm %5.4f}",iFont,fwhm,fwhmHi);
  TLatex landauFWHMLatex(xLabel, yLabel-yLabelStep, landauFWHMLatexChar);
  landauFWHMLatex.SetNDC();
  landauFWHMLatex.SetTextAlign(33);
  landauFWHMLatex.Draw();

  char landauSIGMALatexChar[100];
  //  sprintf(landauSIGMALatexChar,"#sigma = %5.4f^{%+5.4f}_{%+5.4f}",sigma,sigmaHi,sigmaLo);
  sprintf(landauSIGMALatexChar,"#font[%d]{#sigma = %5.4f #pm %5.4f}",iFont,sigma,sigmaHi);
  TLatex landauSIGMALatex(xLabel, yLabel-2.*yLabelStep, landauSIGMALatexChar);
  landauSIGMALatex.SetNDC();
  landauSIGMALatex.SetTextAlign(33);
  landauSIGMALatex.Draw();

  can->SaveAs("./plots/Fit_"+TString(hist->GetName())+".png");

  return mpv;

}

#include <TPaveStats.h>

double binHistos::FitGaus(TH1* hist){

    TPaveStats* stat = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    stat->SetX1NDC(0.55);
    stat->SetX2NDC(0.97);
    stat->SetY1NDC(0.55);
    stat->SetY2NDC(0.93);
    stat->SetOptStat(000002210);
    stat->SetOptFit(111);
    stat->SetBorderSize(0);

    double min, max;
    double mpv = hist->GetBinCenter(hist->GetMaximumBin());
    double rms = hist->GetRMS();
  
    min = mpv - 1.8*rms;
    max = mpv - 1.8*rms;

    TF1 *gfit1 = new TF1("Gaussian1","gaus",min,max); // Create the fit function
    hist->Fit(gfit1,"RQ");

    //    Double_t amp = gfit->GetParameter(0); //value of 0th parameter
    //    Double_t eamp = gfit->GetParError(0); //error on 0th parameter
    Double_t mean = gfit1->GetParameter(1); //value of 1st parameter
    //    Double_t emean = gfit->GetParError(1); //error on 1st parameter
    Double_t sig = gfit1->GetParameter(2); //value of 1st parameter

    TF1 *gfit2 = new TF1("Gaussian2","gaus",mean-1.0*sig,mean+1.0*sig); // Create the fit function
    hist->Fit(gfit2,"RQ");

    // Draw frame on canvas
    TCanvas * can = new TCanvas("FitGausCanvas_"+TString(hist->GetName()),"",1000,1000);
    can->SetRightMargin(0.03);
    can->SetLeftMargin(0.13);
    can->SetTopMargin(0.06);
    can->SetBottomMargin(0.10);

    //    hist->GetXaxis()->SetRange(0,70);

    hist->SetMarkerStyle(20);
    hist->Draw("PE1X0");
    gfit2->SetLineColor(kBlue);
    gfit2->SetLineWidth(2);
    gfit2->Draw("lsame");
    
    can->SaveAs("./plots/FitGaus_"+TString(hist->GetName())+".png");

    return gfit2->GetParameter(1);

}

double binHistos::FitSim(TH1* hist){

  // Construct observable
  RooRealVar t("t","t", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 

  //  setMyRange(hist, &t);
  double mpvEst, rmsEst;
  estMyValues(hist, mpvEst, rmsEst);

  // Construct landau(t,ml,sl) ;
  //  RooRealVar ml("ml","mean landau",0.5,-1.0,1.0);
  //  RooRealVar sl("sl","sigma landau",0.05,0.001,1.0);
  RooRealVar ml("ml","mean landau",mpvEst,0.5*mpvEst,1.5*mpvEst);
  RooRealVar sl("sl","sigma landau",0.01,0.001,0.1);
  RooLandau landau("lx","lx",t,ml,sl);
  
 // Construct gauss(t,mg,sg)
 RooRealVar mg("mg","mg",0.);
 RooRealVar sg("sg","sg",0.03,0.001,0.05);
 RooGaussian gauss("gauss","gauss",t,mg,sg);

 // C o n s t r u c t   c o n v o l u t i o n   p d f 
 // ---------------------------------------

 // Set #bins to be used for FFT sampling to 10000
 t.setBins(10000,"cache") ; 

 // Construct landau (x) gauss
 RooFFTConvPdf lxg("lxg","landau (X) gauss",t,landau,gauss) ;

  // F i t   a n d   p l o t   c o n v o l u t e d   p d f 
  // ----------------------------------------------------------------------

  RooDataHist data(hist->GetName(),hist->GetTitle(),t,hist); 
  RooPlot* tframe = t.frame(RooFit::Name(hist->GetName()),RooFit::Title(hist->GetTitle())) ;

  //RooPlot* tframe = new RooPlot(hist->GetName(), hist->GetTitle(), t, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), hist->GetNbinsX());

  // Fit gxlx to data
  lxg.fitTo(data);

  // Plot
  data.plotOn(tframe);
  lxg.plotOn(tframe) ;
  landau.plotOn(tframe,RooFit::LineStyle(kSolid),RooFit::LineColor(kRed),RooFit::LineWidth(1)) ;

#ifdef oldGetVal
  double mpv = ml.getVal();
  double fwhm = 4.*sl.getVal();
  double sigma = sg.getVal();
#else
  double mpv = ml.getValV();
  double fwhm = 4.*sl.getValV();
  double sigma = sg.getValV();
#endif

  double mpvHi = ml.getErrorHi();
  //  double mpvLo = ml.getErrorLo();

  double fwhmHi = 4.*sl.getErrorHi();
  //  double fwhmLo = 4.*sl.getErrorLo();

  double sigmaHi = sg.getErrorHi();
  //  double sigmaLo = sg.getErrorLo();

  // Draw frame on canvas
  TCanvas * can = new TCanvas("FitCanvas_"+TString(hist->GetName()),"",1000,1000);
  can->SetRightMargin(0.03);
  can->SetLeftMargin(0.13);
  can->SetTopMargin(0.06);
  can->SetBottomMargin(0.10);

  tframe->Draw();
  tframe->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
  tframe->GetYaxis()->SetTitleOffset(1.5);
  tframe->GetXaxis()->SetTitleOffset(1.2);
  
  int iFont = 132;
  double xLabel = 0.935; 
  double yLabel = 0.90;
  double yLabelStep = 0.06;

  char landauMpvLatexChar[100];
  //  sprintf(landauMpvLatexChar,"MPV = %5.4f^{%+5.4f}_{%+5.4f}",mpv,mpvHi,mpvLo);
  sprintf(landauMpvLatexChar,"#font[%d]{MPV = %5.4f #pm %5.4f}",iFont,mpv,mpvHi);
  TLatex landauMpvLatex(xLabel, yLabel, landauMpvLatexChar);
  landauMpvLatex.SetNDC();
  landauMpvLatex.SetTextAlign(33);
  landauMpvLatex.Draw();

  char landauFWHMLatexChar[100];
  //  sprintf(landauFWHMLatexChar,"FWHM = %5.4f^{%+5.4f}_{%+5.4f}",fwhm,fwhmHi,fwhmLo);
  sprintf(landauFWHMLatexChar,"#font[%d]{FWHM = %5.4f #pm %5.4f}",iFont,fwhm,fwhmHi);
  TLatex landauFWHMLatex(xLabel, yLabel-yLabelStep, landauFWHMLatexChar);
  landauFWHMLatex.SetNDC();
  landauFWHMLatex.SetTextAlign(33);
  landauFWHMLatex.Draw();

  char landauSIGMALatexChar[100];
  //  sprintf(landauSIGMALatexChar,"#sigma = %5.4f^{%+5.4f}_{%+5.4f}",sigma,sigmaHi,sigmaLo);
  sprintf(landauSIGMALatexChar,"#font[%d]{#sigma = %5.4f #pm %5.4f}",iFont,sigma,sigmaHi);
  TLatex landauSIGMALatex(xLabel, yLabel-2.*yLabelStep, landauSIGMALatexChar);
  landauSIGMALatex.SetNDC();
  landauSIGMALatex.SetTextAlign(33);
  landauSIGMALatex.Draw();
  
  can->SaveAs("./plots/Fit_"+TString(hist->GetName())+".png");

  return mpv;

}

double binHistos::FitRoot(TH1* hist){

    // Fitting SNR histo
    printf("Fitting...\n");


    TPaveStats* stat = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    stat->SetX1NDC(0.55);
    stat->SetX2NDC(0.97);
    stat->SetY1NDC(0.55);
    stat->SetY2NDC(0.93);
    stat->SetOptStat(000002210);
    stat->SetOptFit(111);
    stat->SetBorderSize(0);


    double min = 0.;
    double max = 0.;

    // Setting fit range and start values
    //
    // Once again, here are the Landau * Gaussian parameters:
    //   par[0]=Width (scale) parameter of Landau density
    //   par[1]=Most Probable (MP, location) parameter of Landau density
    //   par[2]=Total area (integral -inf to inf, normalization constant)
    //   par[3]=Width (sigma) of convoluted Gaussian function
    //
    // Variables for langaufit call:
    //   his             histogram to fit
    //   fitrange[2]     lo and hi boundaries of fit range
    //   startvalues[4]  reasonable start values for the fit
    //   parlimitslo[4]  lower parameter limits
    //   parlimitshi[4]  upper parameter limits
    //   fitparams[4]    returns the final fit parameters
    //   fiterrors[4]    returns the final fit errors
    //   ChiSqr          returns the chi square
    //   NDF             returns ndf
    Double_t fr[2];
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];

    double mpv = hist->GetBinCenter(hist->GetMaximumBin());
    double rms = hist->GetRMS();
    Double_t chisqr;
    Int_t    ndf;

    min = mpv - 1.0*rms;
    max = mpv + 2.0*rms;
    fr[0]=min;
    fr[1]=max;
    pllo[0]=0.01; pllo[1]=0.01*mpv; pllo[2]=10.; pllo[3]=0.01;
    plhi[0]=0.3; plhi[1]=10.*mpv; plhi[2]=3000.; plhi[3]=2.0;
    sv[0]=0.05; sv[1]=mpv; sv[2]=1000.; sv[3]=0.3;

    TF1 *fitsnr = srf::langaufit(hist,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);

    if ( chisqr/(1.*ndf) > 5. ) {

      sv[0]=0.1; sv[3]=0.2;
      fitsnr = srf::langaufit(hist,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
      
    }

    Double_t SNRPeak, SNRFWHM;
    srf::langaupro(fp,SNRPeak,SNRFWHM);

    printf("Fitting done\nPlotting results...\n");

    // Global style settings
    // Draw frame on canvas
    TCanvas * can = new TCanvas("FitCanvas_"+TString(hist->GetName()),"",1000,1000);
    can->SetRightMargin(0.03);
    can->SetLeftMargin(0.13);
    can->SetTopMargin(0.06);
    can->SetBottomMargin(0.10);

    //    hist->GetXaxis()->SetRange(0,70);

    hist->SetMarkerStyle(20);
    hist->Draw("PE1X0");
    fitsnr->SetLineColor(kBlue);
    fitsnr->SetLineWidth(2);
    fitsnr->Draw("lsame");
    
    can->SaveAs("./plots/Fit_"+TString(hist->GetName())+".png");

    return (double)fp[1];

}

double binHistos::setMyRange(TH1* hist, RooRealVar* t){

  TH1D *htemp = (TH1D*)hist->Clone("htemp");

  htemp->Rebin();

  double mpv = htemp->GetBinCenter(htemp->GetMaximumBin());
  double rms = htemp->GetRMS();
  
  double min = mpv - 1.0*rms;
  //  double min = 0.;
  double max = mpv + 2.0*rms;

  t->setRange(min,max);

  delete htemp;

  return mpv;

}

void binHistos::estMyValues(TH1* hist, double& mpv, double& rms){

  TH1D *htemp = (TH1D*)hist->Clone("htemp");

  htemp->Rebin();

  mpv = htemp->GetBinCenter(htemp->GetMaximumBin());
  rms = htemp->GetRMS();
  
  delete htemp;

}

#include "TMath.h"

double binHistos::Median(const TH1D* h1) { 

   int n = h1->GetXaxis()->GetNbins();  
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray(); 
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]); 

}

double binHistos::meanInRange(const TH1D* hist, double& min, double& max) { 

  TH1D *htemp = (TH1D*)hist->Clone("htemp");

  htemp->GetXaxis()->SetRangeUser(min,max);

  return htemp->GetMean(1); 

}
