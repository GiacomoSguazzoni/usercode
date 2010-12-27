// Common includes
//
//MC sim
  TChain *convsim = new TChain("ntupleS2R");
  convsim->Add("../7TeV/ntuple_conversion_minbias7TeV_July12.root");
