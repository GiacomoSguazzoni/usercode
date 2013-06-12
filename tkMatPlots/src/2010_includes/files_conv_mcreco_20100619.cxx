// Common includes
//
//MC reco
  TChain *convmc = new TChain("ntupleR2S");
  convmc->Add("../7TeV/ntuple_conversion_CMSSW358p3_minbias7TeV_2010-06-04.root");
  convmc->Add("../7TeV/ntuple_conversion_CMSSW358p3_minbias7TeV_2010-06-04_v2.root");
