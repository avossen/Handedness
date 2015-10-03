#ifndef HADRON_PAIR_ARRAY_H
#define HADRON_PAIR_ARRAY_H
#include "ReaderBase.h"
#include "TMath.h"

//300 is also the number in the Treesaver, so cannot be larger...
#define Max_ArrSize 300

struct HadronPairArray:public ReaderBase
{
  //cuts:
  const float OpeningCut;
  const float maxOpeningCut;
  const float pi0LowCut;
  const float pi0HighCut;
  const float pi0LowBgCut;
  const float pi0HighBgCut;
  const float zCut;
  const float zUpperCut;
  const float singleZCut;
  const float singleZUpperCut;
  const float maxLabCosTheta;
  const float minLabCosTheta;

  int hadPairNum;
  int numPairs;
  float z[Max_ArrSize];
  float z1[Max_ArrSize];
  float z2[Max_ArrSize];
  //sum over all...
  float zRatio[Max_ArrSize];
  float mass[Max_ArrSize];
  float phiR[Max_ArrSize];
  float phiZero[Max_ArrSize];
  float decayTheta[Max_ArrSize];
  float thrustProj1[Max_ArrSize];
  float thrustProj2[Max_ArrSize];
  float theta1[Max_ArrSize];
  float theta2[Max_ArrSize];


  float phi1[Max_ArrSize];
  float phi2[Max_ArrSize];
  float pi0mass1[Max_ArrSize];
  float pi0mass2[Max_ArrSize];

  int chargeType[Max_ArrSize];
  int particleType[Max_ArrSize];

  //only for mc
  int motherGenId1[Max_ArrSize];
  int motherGenId2[Max_ArrSize];

  int cut[Max_ArrSize];
  int pi0Sig[Max_ArrSize];
  int pi0Bg[Max_ArrSize];

  float hadOpening[Max_ArrSize];
  //considering the thrust axis resolution of 0.16+-0.09 rad, a max opening cut of 0.99 is even too large...
  //changed in nov: opening cut 0.2 -->1.2
  HadronPairArray(TChain* chain,int hadNum, int MCFlag=mcFlagNone):ReaderBase(MCFlag),OpeningCut(0.0),maxOpeningCut(1.9999),pi0LowCut(0.12), pi0HighCut(0.15), pi0LowBgCut(0.21), pi0HighBgCut(0.3), zCut(0.2),singleZCut(0.1), zUpperCut(1.4),singleZUpperCut(1.3),maxLabCosTheta(100.9),minLabCosTheta(-100.6)
  {
    //no chain implies standalone. Cannot branch on the same field twice, this would override
    //    cout <<" do we have a chain? " << chain<<endl;
    if(chain)
      {
	//cout <<"yes " <<endl;
	myChain=chain;
	branchPointersI.push_back(&numPairs);
	
	branchPointers.push_back(z);
	branchPointers.push_back(zRatio);
	branchPointers.push_back(mass);

	branchPointers.push_back(thrustProj1);
	branchPointers.push_back(thrustProj2);
	branchPointers.push_back(decayTheta);

	branchPointers.push_back(theta1);
	branchPointers.push_back(theta2);
	branchPointers.push_back(phi1);
	branchPointers.push_back(phi2);
	branchPointers.push_back(pi0mass1);
	branchPointers.push_back(pi0mass2);
	
	branchPointersI.push_back(chargeType);
	branchPointersI.push_back(particleType);

	if(mMCFlag!=mcFlagNone)
	  {
	    branchPointersI.push_back(motherGenId1);
	    branchPointersI.push_back(motherGenId2);
	  }
	
	string sAdd="1";
	if(hadNum==2)
	  {
	    sAdd="2";
	  }
	hadPairNum=hadNum;
	
	branchNamesI.push_back("mass"+sAdd+addendum+"Counter");
	branchNames.push_back("z"+sAdd+addendum);
	branchNames.push_back("z"+sAdd+"Ratio"+addendum);
	branchNames.push_back("mass"+sAdd+addendum);
	branchNames.push_back("thrustProj"+sAdd+"1"+addendum);
	branchNames.push_back("thrustProj"+sAdd+"2"+addendum);
	branchNames.push_back("decayTheta"+sAdd+addendum);
	branchNames.push_back("theta"+sAdd+"1"+addendum);
	branchNames.push_back("theta"+sAdd+"2"+addendum);

	branchNames.push_back("phi"+sAdd+"1"+addendum);
	branchNames.push_back("phi"+sAdd+"2"+addendum);
	//pi0 mass doesn't exist as mc
	//		branchNames.push_back("pi0Mass"+sAdd+"1"+addendum);
		//	branchNames.push_back("pi0Mass"+sAdd+"2"+addendum);
	branchNames.push_back("pi0Mass"+sAdd+"1");
	branchNames.push_back("pi0Mass"+sAdd+"2");
	
	branchNamesI.push_back("chargeType"+sAdd+addendum);
	branchNamesI.push_back("particleType"+sAdd+addendum);

	//unfortunately the naming convention was not kept...
	if(mMCFlag==mcFlagWoA)
	  {
	    //	    branchNamesI.push_back("motherGenID"+sAdd+"1"+addendum);
	    //	    branchNamesI.push_back("motherGenID"+sAdd+"2"+addendum);
	    	    branchNamesI.push_back("motherGenId_"+sAdd+"1"+addendum);
	    	    branchNamesI.push_back("motherGenId_"+sAdd+"2"+addendum);
	  }
	if(mMCFlag==mcFlagMC)
	  {
	    branchNamesI.push_back("motherGenId_"+sAdd+"1"+addendum);
	    branchNamesI.push_back("motherGenId_"+sAdd+"2"+addendum);
	  }
	doAllBranching();
      }
  }

  void afterFill()
  {
    for(int i=0;i<numPairs;i++)
      {
	//	cout <<"hp after fill!" <<endl;
	cut[i]=0;

	//////////-----------test: restrict to hadron from uds
	if(mMCFlag==mcFlagMC)
	  {
	    //	    cout <<"mc mother Gen1: " << motherGenId1[i] <<" mother gen 2 : " << motherGenId2[i]<<endl;
	      //	    if(abs(motherGenId1[i]) >3 || abs(motherGenId2[i]) > 3)
	    //

	    //	    cout <<"testing mothergen " <<endl;
	    if(!(motherGenId1[i]==10022 && motherGenId2[i] == 10022))
	      {
		//		cout <<"mother gen id not uds : " << motherGenId1[i] <<" and " << motherGenId2[i] <<endl;
		//		cut[i]=1;
		//		cout <<"done " <<endl;
	      }
	    else
	      {
		//		cout <<"hadron comes from uds..." <<endl;
	      }
	  }
	//it seems that the WoA has the evtGen codes. So parents don't seem to be quarks but the virutal photon...
	//nope, not true...
	if(mMCFlag==mcFlagWoA)
	  {
	    if(!(motherGenId1[i]==10022 && motherGenId2[i] == 10022))
	      {
		//		cout <<"mother gen id not uds : " << motherGenId1[i] <<" and " << motherGenId2[i] <<endl;
		//				cut[i]=1;
	      }
	    else
	      {
		//		cout <<"hadron comes from uds..." <<endl;
	      }
	  }

	///------------------

	//       if(cos(decayTheta[i])<0)
	//	  	  cut[i]=1;
	if(cos(theta1[i])<minLabCosTheta && cos(theta2[i])<minLabCosTheta)
	  {
	    cut[i]=1;
	    //	    cout <<"cut hadron due to dec theta " << endl;
	  }
	if(cos(theta1[i])>maxLabCosTheta && cos(theta2[i])>maxLabCosTheta)
	  {
	    cut[i]=1;
	    //	    cout <<"max cos theta " << endl;
	  }
	if(z[i]<=0 || z[i] >1.1)
	  {
	        if(particleType[i]==0 && chargeType[i]==0)
		  //	    	      cout <<" cut due to wrong z: " << z[i] <<endl;
	    cut[i]=1;
	  }

	if(z[i] >zUpperCut)
	  {
	    //	    cout <<"cut due to upper z cut.." <<endl;
	    cut[i]=1;
	  }

	if(z[i] <zCut)
	  {
	    if(particleType[i]==0 && chargeType[i]==0)
	      //	      cout <<" zcut... " << z[i] <<endl;
	    cut[i]=1;
	  }
      
	z1[i]=zRatio[i]/z[i];
	if(z1[i]>0)
	  z1[i]=1/z1[i];
	z2[i]=z[i]-z1[i];

	//	cout <<"hp " <<"  z1["<<i<<"]: " << z1[i] <<" z2: " << z2[i] <<endl;


	if(z1[i]>singleZUpperCut || z2[i]>singleZUpperCut)
	  {
	    //	    cout <<"single z cut.. " <<endl;
	    cut[i]=1;
	  }
	if(z1[i]<singleZCut || z2[i]<singleZCut)
	  {
	    //	       if(particleType[i]==0 && chargeType[i]==0)
		 //	    	      cout <<" singleZcut: "<< z1 <<" " << z2<<endl;
	    cut[i]=1;
	  }

	if(mass[i]<0 || mass[i] > 2.0)
	  {
	    //	      cout <<"cut, mass is wrong... " <<endl;
	  cut[i]=1;
	  }
	//	cout <<"getting dec theta: "<< decayTheta[i]<<", z is: "<< z[i] <<endl;;
	hadOpening[i]=fabs(thrustProj1[i]);
	if(fabs(thrustProj2[i])<fabs(thrustProj1[i]))
	  hadOpening[i]=fabs(thrustProj2[i]);


	//	cout <<"hadOpening: " << hadOpening[i]<<endl;

	if(fabs(thrustProj1[i])<OpeningCut || fabs(thrustProj2[i])<OpeningCut)
	  {
	    cut[i]=1;
	    //	      cout <<"opening cut! " << thrustProj1[i] <<" or: " << thrustProj2[i] <<endl;
	  }
	if(fabs(thrustProj1[i])>maxOpeningCut || fabs(thrustProj2[i])>maxOpeningCut)
	  {
	    //	    cout <<" thrust projection " <<endl;
	  cut[i]=1;
	  }
	if(isnan(mass[i]) || isnan(z[i] || isnan(phiR[i])))
	  {
	   cut[i]=1;
	   //	   cout <<"nan!" <<endl;
	  }
	//one of the particles is pion?
	if(pi0mass1[i]>0 || pi0mass2[i]>0)
	  {
	    float pi0Mass=pi0mass1[i];
	    //	    cout <<"found pi0mass " << pi0Mass <<" cut: " << cut[i] <<endl;
	    if(!(pi0Mass>0))
	      pi0Mass=pi0mass2[i];

	    if(pi0Mass<pi0LowCut || pi0Mass> pi0HighBgCut)
	      {
		pi0Sig[i]=0;
		pi0Bg[i]=0;
	      }
	    if(pi0Mass>pi0LowCut && pi0Mass<pi0HighCut)
	      {
		pi0Sig[i]=1;
		pi0Bg[i]=0;
	      }
	    if(pi0Mass > pi0LowBgCut && pi0Mass < pi0HighBgCut)
	      {
		pi0Sig[i]=0;
		pi0Bg[i]=1;
	      }
	  }

      }
  }

  //set cut flag for all events that don't have valid pi0, set flag to true to mask out any events that do not have pi0s
  void selectPi0Sig(bool onlyPi0Ev=false);
  //set cut flag for all events that don't have pi0 Bg
  void selectPi0Bg();

  //clone this guy...
  HadronPairArray& operator =(HadronPairArray rhs);
  //assign +- randomly, i.e. change direction of R randomly to check for false asymmetries
  void mixItUp();
  void setElement(int pairCounter,HadronPairArray& hp,int index);
  void print();
  protected:
  void setSingleElement(int pairCounter,HadronPairArray& hp,int index);


};


#endif
