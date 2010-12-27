  std::cout << " E ora fillo i mcReco..." << std::endl;
  mcRecoR2S->SetEvRangeMin(minE);
  mcRecoR2S->SetEvRangeMax(maxE);
  mcRecoR2S->SetUIndex(uIndex);
  mcRecoR2S->SetVIndex(vIndex);

//2D
  if ( rawMC2DH ){
    Fill(mcRecoR2S, rawMC2DH, rawFake2DH);
    rawMC2DH->Sumw2();
    rawFake2DH->Sumw2();
    rawMCFs2DH->Add(rawMC2DH,rawFake2DH,1.,-1.);
    rawMCFs2DH->Sumw2();
    //
  }
  if ( rawMC1DH ) {
    Fill(mcRecoR2S, rawMC1DH, rawFake1DH);
    rawMC1DH->Sumw2();
    rawFake1DH->Sumw2();
    rawMCFs1DH->Add(rawMC1DH,rawFake1DH,1.,-1.);
    rawMCFs1DH->Sumw2();
    //
  }
