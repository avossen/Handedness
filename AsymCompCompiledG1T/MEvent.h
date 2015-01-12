#ifndef M_EVENT_H
#define  M_EVENT_H
#include <string>
#include <vector>
#include "TChain.h"
#include "ReaderBase.h"
#include "TVector3.h"
//structure to hold event information

class MEvent:public ReaderBase
{
 public:

  float Thrust;
  float E_miss;
  float thrustThetaCMS;
  float thrustPhiCMS;
  float thrustThetaLab;
  float thrustPhi;
  float thetaEThrust;
  float transProj;
  float longProj;
  float jet1Phi;
  float jet2Phi;
  float jet1Theta;
  float jet2Theta;
  float jet1E;
  float jet2E;


  bool cutEvent;

  TVector3 jet1;
  TVector3 jet2;

  const float minJetE;
  const float maxMissingEnergy;
  const float lowerThrustThetaCut;
  const float upperThrustThetaCut;
  const float lowerThrustThetaCutCMS;
  const float upperThrustThetaCutCMS;
  const float thrustThetaCMSMaxProj;
  const float jetMaxCMSProj;
  const float jetMinCMSProj;




  //hard cuts 1.2-1.6
  //1 - 2.14 is 0.7-1.76 in cms
  //acos(0.5e) is 1.05 rad --> in lab 0.73 -.1.7
  //changed from thetamax proj 0.3 and lower 1.2 , upper 1.9 (cms) to open
  //the non-cms values where 0.14159 and 2.5
  //1.35r is the old, 1.7
  MEvent(TChain* chain, int mMCFlag=mcFlagNone):ReaderBase(mMCFlag),lowerThrustThetaCut(0),upperThrustThetaCut(3.7), thrustThetaCMSMaxProj(1.3), lowerThrustThetaCutCMS(0.0), upperThrustThetaCutCMS(4.9), maxMissingEnergy(2.0),jet1(1,0,0),jet2(1,0,0), minJetE(2.75),jetMaxCMSProj(100),jetMinCMSProj(0.0)
  {
    myChain=chain;
    if(chain)
      {
	branchPointers.push_back(&thrustThetaCMS);
	if(mMCFlag!=mcFlagWoA)
	  {
	    branchPointers.push_back(&jet1E);
	    branchPointers.push_back(&jet2E);
	    branchPointers.push_back(&jet1Phi);
	    branchPointers.push_back(&jet2Phi);
	    branchPointers.push_back(&jet1Theta);
	    branchPointers.push_back(&jet2Theta);
	    branchPointers.push_back(&Thrust);
	    branchPointers.push_back(&E_miss);
	    branchPointers.push_back(&thrustPhiCMS);
	  }
	branchPointers.push_back(&thrustThetaLab);
	//    branchPointers.push_back(&thrustPhi);
	branchPointers.push_back(&thetaEThrust);
	//the field thrustTheta_mc exists but encodes only the difference, so that leads to all cuts to fail
	//    branchNames.push_back("thrustTheta"+addendum);
	if(mMCFlag==mcFlagWoA)
	  branchNames.push_back("thrustTheta"+addendum);
	else
	  {
	    branchNames.push_back("thrustTheta");

	    branchNames.push_back("jetE1"+addendum);
	    branchNames.push_back("jetE2"+addendum);
	    branchNames.push_back("jet1Phi"+addendum);
	    branchNames.push_back("jet2Phi"+addendum);
	    branchNames.push_back("jet1Theta"+addendum);
	    branchNames.push_back("jet2Theta"+addendum);

	   
	    branchNames.push_back("Thrust"+addendum);
	    branchNames.push_back("E_miss"+addendum);
	    branchNames.push_back("thrustPhi");
	  }
	//    branchNames.push_back("thrustThetaLab"+addendum);
	if(mMCFlag==mcFlagWoA)
	  branchNames.push_back("thrustThetaLab"+addendum);////doesn't exist as mc
	else
	  branchNames.push_back("thrustThetaLab");////doesn't exist as mc
	//    branchNames.push_back("thrustPhi"+addendum);
	branchNames.push_back("thetaEThrust"+addendum);
	doAllBranching();
      }
  }
  //decide if event is good or not
  void afterFill()
  {
    jet1.SetMag(1.0);
    jet2.SetMag(1.0);
    jet1.SetPhi(jet1Phi);
    jet2.SetPhi(jet2Phi);
    jet1.SetTheta(jet1Theta);
    jet2.SetTheta(jet2Theta);

    //    cout <<" thrustThetaLab? " << thrustThetaLab <<" cms : " << thrustThetaCMS << " proj: " << fabs(cos(thrustThetaCMS)) <<endl;
    cutEvent=false;

    if(jet1E< minJetE || jet2E< minJetE)
      cutEvent=true;
    
    //institute a cut against too much reconstructed energy...
    if(E_miss<-1)
      cutEvent=true;
    if(E_miss>maxMissingEnergy)
      cutEvent=true;
    if(thrustThetaLab>upperThrustThetaCut)
      cutEvent=true;
    if(thrustThetaLab<lowerThrustThetaCut)
      cutEvent=true;
    if(thrustThetaCMS>upperThrustThetaCutCMS)
      cutEvent=true;
    if(thrustThetaCMS<lowerThrustThetaCutCMS)
      cutEvent=true;

    if(fabs(cos(thrustThetaCMS))>thrustThetaCMSMaxProj)
      cutEvent=true;
    //    cout <<"mcFlag: "<< mMCFlag<<endl;
    //        cout <<"angle between jets: "<< jet1.Angle(jet2)<<" jet1 phi/theta: " << jet1Phi <<"/ " << jet1Theta <<" jet2: " << jet2Phi <<"/ " << jet2Theta<<endl;
    //	cout <<"jet E1: "<< jet1E <<" jet2 E: " << jet2E<<endl;
    if(jet1.Angle(jet2)>1.1)
      cutEvent=true;
    //was 1.34 & 1.8
         if(jet1.Theta()<1.38 || jet1.Theta()>1.75)
           cutEvent=true;
        if(jet2.Theta()<1.38 || jet2.Theta()>1.75)
          cutEvent=true;

    transProj=sin(thetaEThrust)*sin(thetaEThrust)/(1+cos(thetaEThrust)*cos(thetaEThrust));
    longProj=sqrt(1-transProj*transProj);
  }
}; 

#endif
