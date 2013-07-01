//double angle = // Angle in degree between direction x,y,z and outcoming momentum //TMath::Pi()*180.0;

TVector3 vertex(PFDV_x->at(i),PFDV_y->at(i),PFDV_z->at(i)); 
TVector3 pout(1.,0.,0.); 
pout.SetPhi(PFDV_momentumOut_phi->at(i));
pout.SetTheta(PFDV_momentumOut_theta->at(i));

double deta = vertex.Eta()-pout.Eta();
double dphi = ROOT::Math::VectorUtil::DeltaPhi(vertex,pout);
double rcone = sqrt(deta*deta+dphi*dphi);

bool itobtec = (vertex.Perp()>56.)||(fabs(vertex.Z())>95.);

if(
   (PFDV_isNuclear->at(i) && ((PFDV_momentumOut_numberOfTracks->at(i)>2||itobtec)

			      //-nOutTkStep67Poor-nOutTkStep67Good //Not available
			      ) && (PFDV_momentumOut_mass->at(i)>0.6) && (PFDV_momentumInc_pt->at(i)>0.8||PFDV_momentumInc_pt->at(i)<0.1) && (PFDV_momentumOut_pt->at(i)>0.5)
    && (rcone<0.25)
    //&& (angle<15)
    )
   ||
   (PFDV_isNuclearLoose->at(i) && ((PFDV_momentumOut_numberOfTracks->at(i)>1||itobtec)
				   //-nOutTkStep67Good==2 //Not available
				   ) && (PFDV_momentumOut_mass->at(i)>0.6) && (PFDV_momentumOut_pt->at(i)>0.5)
    && (rcone<0.25)
    //&& (angle<15)
    )
   ) 
  {
    iCut = 1;
  };
