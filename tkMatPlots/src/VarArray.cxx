Double_t vAr[8];
Double_t vArTrue[8];

vAr[0]=1.;
vAr[1]=X;
vAr[2]=Y;
vAr[3]=Z;
Double_t radius = sqrt(X*X+Y*Y);
vAr[4]=radius;
vAr[5]=atan2(radius,Z); //theta
vAr[6]=atan2(Y,X); //phi
//vAr[7]=eta; //eta

Double_t A = vAr[uIndex];
Double_t B = vAr[vIndex];

Double_t ACut = vAr[uCutIndex];
Double_t BCut = vAr[vCutIndex];

vArTrue[0]=1.;
vArTrue[1]=XTrue;
vArTrue[2]=YTrue;
vArTrue[3]=ZTrue;
Double_t radiusTrue = sqrt(XTrue*XTrue+YTrue*YTrue);
vArTrue[4]=radiusTrue;
vArTrue[5]=atan2(radiusTrue,ZTrue); //theta
vArTrue[6]=atan2(YTrue,XTrue); //phi
//vAr[7]=eta; //eta

Double_t ATrue = vArTrue[uIndex];
Double_t BTrue = vArTrue[vIndex];

Double_t ACutTrue = vArTrue[uCutIndex];
Double_t BCutTrue = vArTrue[vCutIndex];




