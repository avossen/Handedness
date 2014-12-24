#include "EventMixMap.h"
#include "TwoHadAsymsCommons.h"


void EventMixMap::addHadQuadArray(HadronQuadArray& quad, MEvent& event)
{
  float thrustTheta=event.thrustThetaCMS;
  float thrustPhi=event.thrustPhiCMS;
  normalizeAngle(thrustPhi);
  ////----------------
  ///since the event mix does mix hp2 with hp2, we flip the thrust theta so that effectively hp1 and hp2 are combined...

  thrustTheta=TMath::Pi()-thrustTheta;
  thrustPhi=thrustPhi-TMath::Pi();
  normalizeAngle(thrustTheta);
  normalizeAngle(thrustPhi);
  //  cout <<"adding theta: " << thrustTheta <<" ph i: " << thrustPhi <<endl;
  //////////////---------

  bool thrustFlipped=false;

  //  if(thrustTheta>TMath::Pi()/2)
  if(thrustTheta>TMath::Pi())
    {
      thrustTheta=TMath::Pi()-thrustTheta;
      thrustFlipped=true;
    }

  //  int numThetaBin=getHBin(thrustTheta,numThrustThetaBins,TMath::Pi()/2);
  int numThetaBin=getHBin(thrustTheta,numThrustThetaBins,TMath::Pi());
  numThetaBin--;
  int numPhiBin=getHBin(thrustPhi,numThrustPhiBins,TMath::Pi()*2);
  numPhiBin--;
  ////  cout <<" thrustTheta: " << thrustTheta <<" phi: "<< thrustPhi<<endl;
  //  cout <<"filling theta bin: " << numThetaBin <<" phi bin: " << numPhiBin <<endl;
  (*hadQuads[numThetaBin][numPhiBin])=quad;
  //  cout <<"bin filled .." <<endl;
  binFilled[numThetaBin][numPhiBin]++;
  //  cout <<"done" <<endl;
}


HadronQuadArray* EventMixMap::getHadQuadArray(MEvent& event)
{
  float thrustTheta=event.thrustThetaCMS;
  float thrustPhi=event.thrustPhiCMS;
  normalizeAngle(thrustPhi);
  bool thrustFlipped=false;
  //  if(thrustTheta>TMath::Pi()/2)
  if(thrustTheta>TMath::Pi())
    {
      thrustTheta=TMath::Pi()-thrustTheta;
      thrustFlipped=true;
    }
  //  int numThetaBin=getHBin(thrustTheta,numThrustThetaBins,TMath::Pi()/2);
  //  cout <<" retrieving theta: " << thrustTheta <<" ph i: " << thrustPhi <<endl;
  int numThetaBin=getHBin(thrustTheta,numThrustThetaBins,TMath::Pi());
  numThetaBin--;
  int numPhiBin=getHBin(thrustPhi,numThrustPhiBins,TMath::Pi()*2);
  numPhiBin--;

  //   cout <<" thrustTheta: " << thrustTheta <<" phi: "<< thrustPhi<<endl;
  //      cout<<"checking numthetaBin: " << numThetaBin <<" num phi Bin " << numPhiBin <<endl;
  if(binFilled[numThetaBin][numPhiBin]==0)
    {
      //          cout<<" does not exist ..." <<endl;
    return 0;
    }
  //  cout <<"exists.. " <<endl;
  return hadQuads[numThetaBin][numPhiBin];

}
