// Common includes
//
//MC sim
  TChain *convsim = new TChain("ntupleS2R");
  convsim->Add("../7TeV/362/ntuple_conversion_minbias7TeV_file1.root");
  convsim->Add("../7TeV/362/ntuple_conversion_minbias7TeV_file2.root");
  convsim->Add("../7TeV/362/ntuple_conversion_minbias7TeV_file3.root");
  convsim->Add("../7TeV/362/ntuple_conversion_minbias7TeV_file4.root");
  convsim->Add("../7TeV/362/ntuple_conversion_minbias7TeV_file5.root");
  convsim->Add("../7TeV/362/ntuple_conversion_minbias7TeV_file6.root");
