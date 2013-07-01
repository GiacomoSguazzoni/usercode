std::cout << " E ora fillo i dati..." << std::endl;

dataR2S->SetEvRangeMin(minE);
dataR2S->SetEvRangeMax(maxE);
dataR2S->SetUIndex(uIndex, uCutIndex);
dataR2S->SetVIndex(vIndex, vCutIndex);

//2D
  if ( rawData2DH ) {
    /*
    if ( ! rawMCFake2DH->GetEntries() ) {
      std::cout << " >>>>>>> ERROR!!! Run FillMC first " << std::endl;
      exit(1);
    }
    if ( ! rawSim2DH->GetEntries() ) {
      std::cout << " >>>>>>> ERROR!!! Run FillSim first " << std::endl;
      exit(1);
    }
    */
    Fill(dataR2S, rawData2DH);
    rawData2DH->Sumw2();
    rawData2DH->Scale(DataScaleFact);
    Double_t scaleFact = rawData2DH->Integral()/rawMC2DH->Integral();
    if ( MCScaleFact < 0. ) {
      rawMC2DH->Scale(scaleFact);
      rawMCFake2DH->Scale(scaleFact);
      rawMCFakeSub2DH->Scale(scaleFact);
      rawMCFakeSubWrtSimPosition2DH->Scale(scaleFact);
    } else {
      rawMC2DH->Scale(MCScaleFact);
      rawMCFake2DH->Scale(MCScaleFact);
      rawMCFakeSub2DH->Scale(MCScaleFact);
      rawMCFakeSubWrtSimPosition2DH->Scale(MCScaleFact);
    }
    //
    if ( SimScaleFact < 0. ) {
      rawSim2DH->Scale(scaleFact);
    } else {
      rawSim2DH->Scale(SimScaleFact);
    }
    rawDataFakeSub2DH->Add(rawData2DH,rawMCFake2DH,1.,-1.);
    rawDataFakeSub2DH->Sumw2();
    //
    Data2DH = (TH2D*)rawData2DH->Clone("Data_"+TSName);
    DataFakeSub2DH = (TH2D*)rawDataFakeSub2DH->Clone("DataFakeSub_"+TSName);
    MC2DH = (TH2D*)rawMC2DH->Clone("MC_"+TSName);
    MCFakeSub2DH = (TH2D*)rawMCFakeSub2DH->Clone("MCFakeSub_"+TSName);
    MCFake2DH = (TH2D*)rawMCFake2DH->Clone("MCFake_"+TSName);
    Sim2DH = (TH2D*)rawSim2DH->Clone("Sim_"+TSName);
  }
//1D
  if ( rawData1DH ) {
    /*
    if ( ! rawMCFake1DH->GetEntries() ) {
      std::cout << " >>>>>>> ERROR!!! Run FillMC first " << std::endl;
      exit(1);
    }
    if ( ! rawSim1DH->GetEntries() ) {
      std::cout << " >>>>>>> ERROR!!! Run FillSim first " << std::endl;
      exit(1);
    }
    */
    Fill(dataR2S, rawData1DH);
    std::cout << " rawData1DH GetEntries " << rawData1DH->GetEntries() << " Integral " << rawData1DH->Integral() << std::endl; 
    rawData1DH->Sumw2();
    std::cout << " rawData1DH after sumw2 GetEntries " << rawData1DH->GetEntries() << " Integral " << rawData1DH->Integral() << std::endl; 
    rawData1DH->Scale(DataScaleFact);
    std::cout << " rawData1DH after scale GetEntries " << rawData1DH->GetEntries() << " Integral " << rawData1DH->Integral() << std::endl; 
    Double_t scaleFact = rawData1DH->Integral()/rawMC1DH->Integral();
    if ( MCScaleFact < 0. ) {
      rawMC1DH->Scale(scaleFact);
      rawMCFake1DH->Scale(scaleFact);
      rawMCFakeSub1DH->Scale(scaleFact);
      rawMCFakeSubWrtSimPosition1DH->Scale(scaleFact);
    } else {
      rawMC1DH->Scale(MCScaleFact);
      rawMCFake1DH->Scale(MCScaleFact);
      rawMCFakeSub1DH->Scale(MCScaleFact);
      rawMCFakeSubWrtSimPosition1DH->Scale(MCScaleFact);
    }
    if ( SimScaleFact < 0. ) {
      rawSim1DH->Scale(scaleFact);
    } else {
      rawSim1DH->Scale(SimScaleFact);
    }
    rawDataFakeSub1DH->Add(rawData1DH,rawMCFake1DH,1.,-1.);
    rawDataFakeSub1DH->Sumw2();
    //
    Data1DH = (TH1D*)rawData1DH->Clone("Data_"+TSName);
    DataFakeSub1DH = (TH1D*)rawDataFakeSub1DH->Clone("DataFakeSub_"+TSName);
    MC1DH = (TH1D*)rawMC1DH->Clone("MC_"+TSName);
    MCFakeSub1DH = (TH1D*)rawMCFakeSub1DH->Clone("MCFakeSub_"+TSName);
    MCFake1DH = (TH1D*)rawMCFake1DH->Clone("MCFake_"+TSName);
    Sim1DH = (TH1D*)rawSim1DH->Clone("Sim_"+TSName);
  }
