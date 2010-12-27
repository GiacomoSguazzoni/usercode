std::cout << " E ora fillo i dati..." << std::endl;

dataR2S->SetEvRangeMin(minE);
dataR2S->SetEvRangeMax(maxE);
dataR2S->SetUIndex(uIndex);
dataR2S->SetVIndex(vIndex);

//2D
  if ( rawData2DH ) {
    if ( ! rawFake2DH->GetEntries() ) {
      std::cout << " >>>>>>> ERROR!!! Run FillMC first " << std::endl;
      exit(1);
    }
    if ( ! rawSim2DH->GetEntries() ) {
      std::cout << " >>>>>>> ERROR!!! Run FillSim first " << std::endl;
      exit(1);
    }
    Fill(dataR2S, rawData2DH);
    rawData2DH->Sumw2();
    rawData2DH->Scale(DataScaleFact);
    Double_t scaleFact = rawData2DH->Integral()/rawMC2DH->Integral();
    if ( MCScaleFact < 0. ) {
      rawMC2DH->Scale(scaleFact);
      rawFake2DH->Scale(scaleFact);
      rawMCFs2DH->Scale(scaleFact);
    } else {
      rawMC2DH->Scale(MCScaleFact);
      rawFake2DH->Scale(MCScaleFact);
      rawMCFs2DH->Scale(MCScaleFact);
    }
    //
    if ( SimScaleFact < 0. ) {
      rawSim2DH->Scale(scaleFact);
    } else {
      rawSim2DH->Scale(SimScaleFact);
    }
    rawDataFs2DH->Add(rawData2DH,rawFake2DH,1.,-1.);
    rawDataFs2DH->Sumw2();
    //
    Data2DH = (TH2D*)rawData2DH->Clone("Data_"+TSName);
    DataFs2DH = (TH2D*)rawDataFs2DH->Clone("DataFs_"+TSName);
    MC2DH = (TH2D*)rawMC2DH->Clone("MC_"+TSName);
    MCFs2DH = (TH2D*)rawMCFs2DH->Clone("MCFs_"+TSName);
    Fake2DH = (TH2D*)rawFake2DH->Clone("Fake_"+TSName);
    Sim2DH = (TH2D*)rawSim2DH->Clone("Sim_"+TSName);
  }
//1D
  if ( rawData1DH ) {
    if ( ! rawFake1DH->GetEntries() ) {
      std::cout << " >>>>>>> ERROR!!! Run FillMC first " << std::endl;
      exit(1);
    }
    if ( ! rawSim1DH->GetEntries() ) {
      std::cout << " >>>>>>> ERROR!!! Run FillSim first " << std::endl;
      exit(1);
    }
    Fill(dataR2S, rawData1DH);
    std::cout << " rawData1DH GetEntries " << rawData1DH->GetEntries() << " Integral " << rawData1DH->Integral() << std::endl; 
    rawData1DH->Sumw2();
    std::cout << " rawData1DH after sumw2 GetEntries " << rawData1DH->GetEntries() << " Integral " << rawData1DH->Integral() << std::endl; 
    rawData1DH->Scale(DataScaleFact);
    std::cout << " rawData1DH after scale GetEntries " << rawData1DH->GetEntries() << " Integral " << rawData1DH->Integral() << std::endl; 
    Double_t scaleFact = rawData1DH->Integral()/rawMC1DH->Integral();
    if ( MCScaleFact < 0. ) {
      rawMC1DH->Scale(scaleFact);
      rawFake1DH->Scale(scaleFact);
      rawMCFs1DH->Scale(scaleFact);
    } else {
      rawMC1DH->Scale(MCScaleFact);
      rawFake1DH->Scale(MCScaleFact);
      rawMCFs1DH->Scale(MCScaleFact);
    }
    if ( SimScaleFact < 0. ) {
      rawSim1DH->Scale(scaleFact);
    } else {
      rawSim1DH->Scale(SimScaleFact);
    }
    rawDataFs1DH->Add(rawData1DH,rawFake1DH,1.,-1.);
    rawDataFs1DH->Sumw2();
    //
    Data1DH = (TH1D*)rawData1DH->Clone("Data_"+TSName);
    DataFs1DH = (TH1D*)rawDataFs1DH->Clone("DataFs_"+TSName);
    MC1DH = (TH1D*)rawMC1DH->Clone("MC_"+TSName);
    MCFs1DH = (TH1D*)rawMCFs1DH->Clone("MCFs_"+TSName);
    Fake1DH = (TH1D*)rawFake1DH->Clone("Fake_"+TSName);
    Sim1DH = (TH1D*)rawSim1DH->Clone("Sim_"+TSName);
  }
