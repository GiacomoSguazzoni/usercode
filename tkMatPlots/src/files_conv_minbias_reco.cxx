// Common includes
//
//MC reco
  TChain *mc = new TChain("ntupleR2S");
  mc->Add("/raid/sguazz/maps/minbias_mc_reco_conv.root");
