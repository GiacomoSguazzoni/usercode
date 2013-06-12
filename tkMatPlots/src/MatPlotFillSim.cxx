  std::cout << " E ora fillo i mcSim..." << std::endl;
  mcSimS2R->SetEvRangeMin(minE);
  mcSimS2R->SetEvRangeMax(maxE);
  mcSimS2R->SetUIndex(uIndex, uCutIndex);
  mcSimS2R->SetVIndex(vIndex, vCutIndex);

  if ( rawSim2DH )
    {
      Fill(mcSimS2R, rawSim2DH);
      rawSim2DH->Sumw2();
    } 
  if ( rawSim1DH ) 
    {
      Fill(mcSimS2R, rawSim1DH);
      rawSim1DH->Sumw2();
    } 
