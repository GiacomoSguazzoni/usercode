// Common includes
//
//Data
  TChain *data = new TChain("ntupleR2S");
  data->Add("/raid/sguazz/maps/minbias_run2012A_aod_conv.root");
  data->Add("/raid/sguazz/maps/minbias_run2012B_aod_conv.root");
  data->Add("/raid/sguazz/maps/minbias_run2012Cpart1_aod_conv.root");
  data->Add("/raid/sguazz/maps/minbias_run2012Cpart2_aod_conv.root");
  data->Add("/raid/sguazz/maps/minbias_run2012D_aod_conv.root");
