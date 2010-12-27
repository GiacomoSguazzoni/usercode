// Common includes
//
//MC reco
  TChain *mc = new TChain("ntupleR2S");
  mc->Add("../7TeV/ntuple_conversion_minbias7TeV_July12.root");
