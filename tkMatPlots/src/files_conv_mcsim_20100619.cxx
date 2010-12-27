// Common includes
//
//MC sim
  TChain *convsim = new TChain("ntupleS2R");
  convsim->Add("../7TeV/ntuple_conversion_CMSSW358p3_minbias7TeV_2010-06-04.root");
  convsim->Add("../7TeV/ntuple_conversion_CMSSW358p3_minbias7TeV_2010-06-04_v2.root");
