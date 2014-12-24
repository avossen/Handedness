#ifndef READER_BASE_H
#define READER_BASE_H
#include "TChain.h"
#include <vector>
#include <string>
#include "TMath.h"
#include <iostream>
#include "TwoHadAsymsCommons.h"


using namespace std;

class ReaderBase
{
 public:
  int mMCFlag;
  vector<float*> branchPointers;
  vector<string> branchNames;
  vector<int*> branchPointersI;
  vector<string> branchNamesI;

  

  TChain* myChain;
  ReaderBase(int mcFlag=mcFlagNone)
    {
      mMCFlag=mcFlag;
      if(mcFlagMC==mMCFlag)
	addendum="_mc";
      if(mcFlagWoA==mMCFlag)
	addendum="_mcWoA";
    }


  void doAllBranching()
  {
    for(unsigned int i=0;i<branchPointers.size();i++)
      {
	myChain->SetBranchAddress(branchNames[i].c_str(),branchPointers[i]);
	cout <<" branching: " << branchNames[i]<<endl;
      }
    for(unsigned int i=0;i<branchPointersI.size();i++)
      {
	myChain->SetBranchAddress(branchNamesI[i].c_str(),branchPointersI[i]);
	cout <<" branching: " << branchNamesI[i]<<endl;
      }
    cout <<"done with branching.. " <<endl;
  }

  void afterFill()
  {
  };

  void normalizeAngle(float& angle)
  {

    while(angle>2*TMath::Pi())
      {
	angle-=2*TMath::Pi();
      }
    while(angle<0)
      {
	angle+=2*TMath::Pi();
      }
    //for the cases where we potentially add 2pi...

  }

 protected:
  string addendum;
};

#endif
