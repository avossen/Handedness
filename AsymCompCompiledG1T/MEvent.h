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


  int D0Tag;
  int DStarTag;

  int runNr;
  int evtNr;

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
  const float maxJetCMSTheta;
  const float minJetCMSTheta;

  const bool onlyDTagged;
  const bool onlyDStarTagged;


  //hard cuts 1.2-1.6
  //1 - 2.14 is 0.7-1.76 in cms
  //acos(0.5e) is 1.05 rad --> in lab 0.73 -.1.7
  //changed from thetamax proj 0.3 and lower 1.2 , upper 1.9 (cms) to open
  //the non-cms values where 0.14159 and 2.5
  //1.35r is the old, 1.7
  //max was 1.76, min: 1.38 for jet cms
  //1.52 <-> 175 for thrust theta
  //maxJetCMSTheta war 1.761593
  /////  MEvent(TChain* chain, int mMCFlag=mcFlagNone):ReaderBase(mMCFlag),lowerThrustThetaCut(-1.52),upperThrustThetaCut(10.75), thrustThetaCMSMaxProj(1.3), lowerThrustThetaCutCMS(0.0), upperThrustThetaCutCMS(4.9), maxMissingEnergy(2.0),jet1(1,0,0),jet2(1,0,0), minJetE(3.75),jetMaxCMSProj(100),jetMinCMSProj(0.0),maxJetCMSTheta(1.75), minJetCMSTheta(1.38), onlyDTagged(false),onlyDStarTagged(false)
  MEvent(TChain* chain, int mMCFlag=mcFlagNone):ReaderBase(mMCFlag),lowerThrustThetaCut(-1.52),upperThrustThetaCut(10.75), thrustThetaCMSMaxProj(1.3), lowerThrustThetaCutCMS(0.0), upperThrustThetaCutCMS(4.9), maxMissingEnergy(2.0),jet1(1,0,0),jet2(1,0,0), minJetE(3.75),jetMaxCMSProj(100),jetMinCMSProj(0.0),maxJetCMSTheta(1.7616), minJetCMSTheta(1.38), onlyDTagged(false),onlyDStarTagged(false)
  {
    myChain=chain;
    if(chain)
      {
	branchPointers.push_back(&thrustThetaCMS);
	if(mMCFlag!=mcFlagWoA)
	  {
	    branchPointersI.push_back(&runNr);
	    branchPointersI.push_back(&evtNr);
	    branchPointersI.push_back(&D0Tag);
	    branchPointersI.push_back(&DStarTag);
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
	    branchNamesI.push_back("runNr");
	    branchNamesI.push_back("eventNr");
	    branchNamesI.push_back("D0Tag");
	    branchNamesI.push_back("DStarTag");
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
    //    cout <<"runNr: "<< runNr <<" event nr: "<< evtNr <<endl;
    if(jet1E< minJetE || jet2E< minJetE)
      {
	//	cout <<"cutting event " << evtNr <<" due to jet E : "<< jet1E <<" or " << jet2E <<endl;
      cutEvent=true;
      }
    
    if(onlyDTagged && D0Tag!=1)
      cutEvent=true;
    if(onlyDStarTagged && DStarTag!=1)
      cutEvent=true;
    //    if(!cutEvent)
      //      cout <<"found DStar" <<endl;
    //    cout <<"onlyDSTarTagged: "<< onlyDStarTagged <<" DStarTag: " << DStarTag <<" cut event" << cutEvent <<endl;
    //institute a cut against too much reconstructed energy...
    if(E_miss<-1)
      {
	//cout <<" event " <<evtNr<< " emiss cut " <<endl;
      cutEvent=true;
      }
    if(E_miss>maxMissingEnergy)
      {
	//		cout <<"event " << evtNr<< " second emiss cut " << E_miss<<endl;
      cutEvent=true;
      }
    if(thrustThetaLab>upperThrustThetaCut)
      {
	//		cout <<"thrustThetaLab cut event " << evtNr <<" thrustTHetaLab: " << thrustThetaLab <<endl;
	cutEvent=true;
      }
    if(thrustThetaLab<lowerThrustThetaCut)
      {
	//		cout <<"lower cut event " << evtNr <<" thrustTHetaLab: " << thrustThetaLab <<endl;
	cutEvent=true;
      }
    if(thrustThetaCMS>upperThrustThetaCutCMS)
      {
	//		cout <<"upper theta cms cut event " << evtNr <<" thrustTHetaCMS: " << thrustThetaCMS <<endl;
      cutEvent=true;
      }
    if(thrustThetaCMS<lowerThrustThetaCutCMS)
      {
	//		cout <<"lower theta cms cut event " << evtNr <<" thrustTHetaCMS: " << thrustThetaCMS <<endl;
      cutEvent=true;
      }

    if(fabs(cos(thrustThetaCMS))>thrustThetaCMSMaxProj)
      {
	//		cout <<"theta cms proj cut event " << evtNr <<" thrustTHetaCMSProj: " << fabs(cos(thrustThetaCMS)) <<endl;
	cutEvent=true;
      }
    //    cout <<"mcFlag: "<< mMCFlag<<endl;
    //        cout <<"angle between jets: "<< jet1.Angle(jet2)<<" jet1 phi/theta: " << jet1Phi <<"/ " << jet1Theta <<" jet2: " << jet2Phi <<"/ " << jet2Theta<<endl;
    //	cout <<"jet E1: "<< jet1E <<" jet2 E: " << jet2E<<endl;

    //    if(jet1.Angle(jet2)>0.3)
    //    if(jet1.Angle(jet2)>0.3)
    //      cutEvent=true;
    //was 1.34 & 1.8
    if(jet1.Theta()<minJetCMSTheta || jet1.Theta()>maxJetCMSTheta)
      {

	//	cout <<"event : " << evtNr <<" jet1 cms theta cuts: " << jet1.Theta() <<" (second one: "<< jet2.Theta() <<endl;
      cutEvent=true;

      }
    if(jet2.Theta()<minJetCMSTheta || jet2.Theta()>maxJetCMSTheta)
      {
	//	cout <<"event : " << evtNr <<" jet2 cms theta cuts: " << jet2.Theta() <<" (second one: "<< jet1.Theta() <<endl;
      cutEvent=true;
      }
    transProj=sin(thetaEThrust)*sin(thetaEThrust)/(1+cos(thetaEThrust)*cos(thetaEThrust));
    longProj=sqrt(1-transProj*transProj);
    //    cout <<"cut this event? " << cutEvent <<endl;
  }
}; 

#endif
