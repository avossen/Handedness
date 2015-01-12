#include "EventMixMap.h"
#include "TwoHadAsymsCommons.h"


void EventMixMap::addHadQuadArray(HadronQuadArray& quad, MEvent& event)
{

  //use jet axis instead of thrust axis
  float thrustTheta=event.thrustThetaCMS;
  float thrustPhi=event.thrustPhiCMS;

  //we mix hp2 with hp2, so we have to take the second jet, but when we return, we have to match it to the first jet
  float jetTheta=event.jet2Theta;
  float jetPhi=event.jet2Phi;

  normalizeAngle(jetPhi);
  ////----------------
  ///since the event mix does mix hp2 with hp2, we flip the thrust theta so that effectively hp1 and hp2 are combined...
  //  thrustTheta=TMath::Pi()-thrustTheta;
  //  thrustPhi=thrustPhi-TMath::Pi();

  jetTheta=TMath::Pi()-jetTheta;
  jetPhi=jetPhi-TMath::Pi();
  normalizeAngle(jetTheta);
  normalizeAngle(jetPhi);
  //  normalizeAngle(thrustPhi);
  //  cout <<"adding theta: " << thrustTheta <<" ph i: " << thrustPhi <<endl;
  //////////////---------

  //  bool thrustFlipped=false;
  bool jetFlipped=false;
  //  if(thrustTheta>TMath::Pi()/2)
  //  if(thrustTheta>TMath::Pi())
  //    {
  //      thrustTheta=TMath::Pi()-thrustTheta;
  //      thrustFlipped=true;
  //    }

  if(jetTheta>TMath::Pi())
    {
      jetTheta=TMath::Pi()-jetTheta;
      jetFlipped=true;
    }
  normalizeAngle(jetTheta);
  //  int numThetaBin=getHBin(thrustTheta,numThrustThetaBins,TMath::Pi()/2);
     //  int numThetaBin=getHBin(thrustTheta,numThrustThetaBins,TMath::Pi());
     //  numThetaBin--;
     //  int numPhiBin=getHBin(thrustPhi,numThrustPhiBins,TMath::Pi()*2);
     //  numPhiBin--;

     int numThetaBin=getHBin(jetTheta,numThrustThetaBins,TMath::Pi());
     numThetaBin--;
       int numPhiBin=getHBin(jetPhi,numThrustPhiBins,TMath::Pi()*2);
     numPhiBin--;

     //         cout <<" thrustTheta: " << jetTheta <<" phi: "<< jetPhi<<endl;
     //         cout <<"filling theta bin: " << numThetaBin <<" phi bin: " << numPhiBin <<endl;
     int locBinPointer=binPointer[numThetaBin][numPhiBin];
     //          cout <<"binPointer: "<< locBinPointer <<endl;
     (*hadQuads[numThetaBin][numPhiBin][locBinPointer])=quad;

  //  cout <<"bin filled .." <<endl;
     //set the new pointer...(max queueLength -1...)
     binPointer[numThetaBin][numPhiBin]=(++locBinPointer)%(queueLength);
     binFilled[numThetaBin][numPhiBin]=(locBinPointer==0 ? queueLength : locBinPointer);

     //     cout <<" bin pointer is: now " << binPointer[numThetaBin][numPhiBin] <<" num fill: " << binFilled[numThetaBin][numPhiBin] <<endl;
  //  cout <<"done" <<endl;
}


HadronQuadArray** EventMixMap::getHadQuadArray(MEvent& event,int& numQuads)
{
  float jetTheta=event.jet1Theta;
  float jetPhi=event.jet1Phi;

  //  float thrustTheta=event.thrustThetaCMS;
  //  float thrustPhi=event.thrustPhiCMS;
  //  normalizeAngle(thrustPhi);

  normalizeAngle(jetPhi);

  //  bool thrustFlipped=false;
  bool jetFlipped=true;
  //  if(thrustTheta>TMath::Pi()/2)
  //  if(thrustTheta>TMath::Pi())
  //    {
  //      thrustTheta=TMath::Pi()-thrustTheta;
  //      thrustFlipped=true;
  //    }

  if(jetTheta>TMath::Pi())
    {
      jetTheta=TMath::Pi()-jetTheta;
      jetFlipped=true;
    }

  //  int numThetaBin=getHBin(thrustTheta,numThrustThetaBins,TMath::Pi()/2);
  //  cout <<" retrieving theta: " << thrustTheta <<" ph i: " << thrustPhi <<endl;
  int numThetaBin=getHBin(jetTheta,numThrustThetaBins,TMath::Pi());
  numThetaBin--;
  int numPhiBin=getHBin(jetPhi,numThrustPhiBins,TMath::Pi()*2);
  numPhiBin--;

  //   cout <<" thrustTheta: " << thrustTheta <<" phi: "<< thrustPhi<<endl;
  //      cout<<"checking numthetaBin: " << numThetaBin <<" num phi Bin " << numPhiBin <<endl;

  numQuads=binFilled[numThetaBin][numPhiBin];

  if(binFilled[numThetaBin][numPhiBin]==0)
    {
      //          cout<<" does not exist ..." <<endl;
      return 0;
    }
  //  cout <<"exists.. " <<endl;
  //  cout <<" we have "<< numQuads<< " quads in the queue " << endl;
  return hadQuads[numThetaBin][numPhiBin];

}
