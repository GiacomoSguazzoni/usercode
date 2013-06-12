  std::cout << " E ora fillo i mcReco..." << std::endl;
  mcRecoR2S->SetEvRangeMin(minE);
  mcRecoR2S->SetEvRangeMax(maxE);
  mcRecoR2S->SetUIndex(uIndex, uCutIndex);
  mcRecoR2S->SetVIndex(vIndex, vCutIndex);

//2D
  if ( rawMC2DH ){
    Fill(mcRecoR2S, rawMC2DH, rawMCFake2DH, rawMCFakeSubWrtSimPosition2DH);
    rawMC2DH->Sumw2();
    rawMCFakeSubWrtSimPosition2DH->Sumw2();
    rawMCFake2DH->Sumw2();
    rawMCFakeSub2DH->Add(rawMC2DH,rawMCFake2DH,1.,-1.);
    rawMCFakeSub2DH->Sumw2();
    //
  }
  if ( rawMC1DH ) {
    Fill(mcRecoR2S, rawMC1DH, rawMCFake1DH, rawMCFakeSubWrtSimPosition1DH);
    rawMC1DH->Sumw2();
    rawMCFakeSubWrtSimPosition1DH->Sumw2();
    rawMCFake1DH->Sumw2();
    rawMCFakeSub1DH->Add(rawMC1DH,rawMCFake1DH,1.,-1.);
    rawMCFakeSub1DH->Sumw2();
    //
  }
