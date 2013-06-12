// Common includes
//
//MC sim
  TChain *sim = new TChain("ntupleS2R");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_01.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_02.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_03.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_04.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_05.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_06.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_07.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_08.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_09.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_10.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_11.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_12.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_13.root");
  sim->Add("/data01/sguazz/MinBiasMC_386p1_altTrig/ntuple_nuclint_MinBiasMC_14.root");
