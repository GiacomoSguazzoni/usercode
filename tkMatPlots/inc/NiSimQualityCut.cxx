Int_t sguazzCut = 1;
double mOut = MC_TrkV_momentumOut_mass->at(i);
int nOut = MC_TrkV_momentumOut_numberOfParticles->at(i);

if ( mOut>0.104 && mOut<0.106 && nOut==1) sguazzCut = 0;
if ( mOut>0.0005 && mOut<0.00052 && nOut==1) sguazzCut = 0;
if ( mOut>0.994 && mOut<0.996 && nOut==6) sguazzCut = 0;
if ( mOut>1.280 && mOut<1.284 && nOut==4) sguazzCut = 0;
if ( mOut>0.708 && mOut<0.712 && nOut==4) sguazzCut = 0;
// no muon
// no electrons
// no 995.2 with nOut == 6
// no 1.282 with nOut == 4
// no 0.7102 with nOut == 4

if ( MC_TrkV_isNuclearInteraction->at(i) && sguazzCut ) iCut = 1;


