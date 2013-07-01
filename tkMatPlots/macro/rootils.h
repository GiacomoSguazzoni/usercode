#include <TString.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <cmath>

namespace rootils {

  void SetFancyGrayscalePalette();
  void SetFancyLightGreenToLightRedPalette();
  void SetFancyColorPalette();
  TString generateRandomName();
  TCanvas* makeNiceCanvasByPixMargins(Int_t pixelPerBinX, Int_t pixelPerBinY, Int_t nbinx, Int_t nbiny, Int_t top, Int_t bottom, Int_t left, Int_t right);
  TCanvas* makeNiceCanvasByFracMargins(Int_t pixelPerBin, Int_t nbinx, Int_t nbiny, Double_t top, Double_t bottom, Double_t left, Double_t right) ;
  void drawArrow(Double_t u, Double_t v, Double_t du, Double_t due, Double_t dv, Double_t dve);
  Double_t ComputeRelSize(Int_t pxl);
  Double_t ComputeTitleOffset(Int_t fontpxl, Int_t pxl);
  TH2F* drawFrame(TString nameU, TString nameV, Double_t minU, Double_t maxU, Double_t minV, Double_t maxV);
  void drawEtaValues(Float_t, Float_t, Float_t, Float_t);
  void drawEtaValuesOnThetaPlot(Float_t, Float_t, Float_t, Float_t);
  Double_t max(Double_t, Double_t);
  Double_t min(Double_t, Double_t);

  inline double deltaPhi(double phi1, double phi2) { 
    double result = phi1 - phi2;
    while (result > M_PI) result -= 2*M_PI;
    while (result <= -M_PI) result += 2*M_PI;
    return result;
  }

  inline double deltaPhi(float phi1, double phi2) {
    return deltaPhi(static_cast<double>(phi1), phi2);
  }
  
  inline double deltaPhi(double phi1, float phi2) {
    return deltaPhi(phi1, static_cast<double>(phi2));
  }
  

  inline float deltaPhi(float phi1, float phi2) { 
    float result = phi1 - phi2;
    while (result > float(M_PI)) result -= float(2*M_PI);
    while (result <= -float(M_PI)) result += float(2*M_PI);
    return result;
  }

  template<typename T1, typename T2>
    inline double deltaPhi(T1& t1, T2 & t2) {
    return deltaPhi(t1.phi(), t2.phi());
  }      

  template <typename T> 
    inline T deltaPhi (T phi1, T phi2) { 
    T result = phi1 - phi2;
    while (result > M_PI) result -= 2*M_PI;
    while (result <= -M_PI) result += 2*M_PI;
    return result;
  }

}

