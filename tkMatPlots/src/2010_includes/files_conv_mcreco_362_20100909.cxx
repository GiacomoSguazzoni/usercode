// Common includes
//
//MC reco
  TChain *convmc = new TChain("ntupleR2S");
  convmc->Add("../7TeV/362/ntuple_conversion_minbias7TeV_file1.root");
  convmc->Add("../7TeV/362/ntuple_conversion_minbias7TeV_file2.root");
  convmc->Add("../7TeV/362/ntuple_conversion_minbias7TeV_file3.root");
  convmc->Add("../7TeV/362/ntuple_conversion_minbias7TeV_file4.root");
  convmc->Add("../7TeV/362/ntuple_conversion_minbias7TeV_file5.root");
  convmc->Add("../7TeV/362/ntuple_conversion_minbias7TeV_file6.root");







