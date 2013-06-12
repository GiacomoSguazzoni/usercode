  std::cout << "Unfold training... " << std::endl;
  s2r->SetEvRangeMin(minE);
  s2r->SetEvRangeMax(maxE);
  r2s->SetEvRangeMin(minE);
  r2s->SetEvRangeMax(maxE);

s2r->SetGeoCuts(geoCuts, effUV);
  s2r->LoopForTrain(response);
r2s->SetGeoCuts(geoCuts, effUV);
  r2s->LoopForTrain(response);

