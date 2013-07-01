// Common includes
//
//MC sim
  TChain *sim = new TChain("MyNtupleMaking/NuclearInteractionsTree");
  sim->Add("root://eoscms.cern.ch//eos/cms/store/group/comm_tracker/Material/Ntuple__MinimumBias__Run2012D-PromptReco-v1__MBZBRndTriggers__203768-208686/Ntuple__MinimumBias__Run2012D-PromptReco-v1__MBZBRndTriggers_1000_1_nye.root");   
