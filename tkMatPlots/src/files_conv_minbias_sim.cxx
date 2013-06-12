// Common includes
//
//MC sim
  TChain *sim = new TChain("ntupleS2R");
  sim->Add("/raid/sguazz/maps/minbias_mc_reco_conv.root");
