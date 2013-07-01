// Common includes
//
//MC reco
  TChain *mc = new TChain("MyNtupleMaking/NuclearInteractionsTree");
  mc->Add("root://eoscms.cern.ch//eos/cms/store/group/comm_tracker/Material/Ntuple__MinimumBias__Run2012D-PromptReco-v1__MBZBRndTriggers__203768-208686/Ntuple__MinimumBias__Run2012D-PromptReco-v1__MBZBRndTriggers_1001_1_qFQ.root");   

