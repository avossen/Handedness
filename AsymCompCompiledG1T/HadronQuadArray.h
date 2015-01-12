#ifndef HADRON_QUAD_ARRAY_H
#define HADRON_QUAD_ARRAY_H

#include "ReaderBase.h"
#include "HadronPairArrays.h"
#include "TMath.h"


struct HadronQuadArray: public ReaderBase
{
  HadronPairArray hp1;
  HadronPairArray hp2;
  //float 
  int numHadQuads;
  float phiRSum[Max_ArrSize];
  float phiR1[Max_ArrSize];

  float phiZeroR[Max_ArrSize];
  float phiZero1[Max_ArrSize];

  float phiRDiff[Max_ArrSize];
  float twoPhiRDiff[Max_ArrSize];

  float phiZeroDiff[Max_ArrSize];
  float phiZeroSum[Max_ArrSize];
  float twoPhiZeroDiff[Max_ArrSize];

  float qT[Max_ArrSize];

  float weight[Max_ArrSize];
  float weightZero[Max_ArrSize];

  bool debugPrint;

  int minPairs;
  const float qTCut;
  const float minQt;
  //qt cut used to be 4.7
  HadronQuadArray(TChain* chain, int MCFlag=mcFlagNone):ReaderBase(MCFlag),hp1(chain, 1,MCFlag), hp2(chain, 2,MCFlag),qTCut(2.7), debugPrint(false), minPairs(0), minQt(0.0)
  {
    //        cout <<"had quad array chain: " << chain <<endl;
    myChain=chain;
    if(chain)
      {
	branchPointersI.push_back(&(this->numHadQuads));
	//all should be the same..
	branchNamesI.push_back("z1"+addendum+"Counter");
	branchNames.push_back("phiRSum"+addendum);
	branchNames.push_back("phiR1"+addendum);
	//for some reason was named so..., i.e phiZero1 is different for both..
	if(mMCFlag==mcFlagWoA)
	  {
	    branchNames.push_back("phiZero1"+addendum);
	    branchNames.push_back("twoPhiZero"+addendum);
	  }
	else
	  {
	    branchNames.push_back("phiZeroR"+addendum);
	    branchNames.push_back("phiZero1"+addendum);
	  }

	if(mMCFlag==mcFlagWoA)
	  branchNames.push_back("qT"+addendum);
	else
	  branchNames.push_back("qT");

	branchPointers.push_back(phiRSum);
	branchPointers.push_back(phiR1);
	branchPointers.push_back(phiZeroR);
	branchPointers.push_back(phiZero1);

	branchPointers.push_back(qT);
	doAllBranching();
      }
  }

  void afterFill()
  {

    hp1.afterFill();
    hp2.afterFill();
    //so we can check in the log files...
    if(hp1.numPairs>Max_ArrSize)
      cout<<" too many pairs!!!: " << hp1.numPairs<<endl;

    //    cout <<"after fill " << endl;
    int numPairsHemi1=-1;
    int numPairsHemi2=-1;

    numPairs(numPairsHemi1,numPairsHemi2);
    //	cout <<" event pair counts, hemi1 : " << numPairsHemi1 <<" hemi2: "<< numPairsHemi2 << ", overall: " << hp1.numPairs<<endl;
    if(numPairsHemi1<minPairs || numPairsHemi2<minPairs)
      {
	//cut all
	for(int i=0;i<hp1.numPairs;i++)
	  {
	    hp1.cut[i]=1;
	    hp2.cut[i]=1;
	  }
      }

    for(int i=0;i<hp1.numPairs;i++)
      {
	weight[i]=1.0;
	weightZero[i]=1.0;

	if(qT[i]<minQt)
	  {
	    hp1.cut[i]=1;
	    hp2.cut[i]=1;
	  }

	if(qT[i]>qTCut)
	  {
	    //t	    cout <<" qt cut: " << qT[i] <<endl;
	    hp1.cut[i]=1;
	    hp2.cut[i]=1;
	  }
	if(isnan(phiRSum[i]))
	  {
	    //	    cout <<"isnan phirsum" <<endl;
	    hp1.cut[i]=1;
	    hp2.cut[i]=1;
	  }
	if(isnan(phiZeroR[i]))
	  {
	    //	    cout <<"isnan phirzero" <<endl;
	    hp1.cut[i]=1;
	    hp2.cut[i]=1;
	  }
	hp1.phiZero[i]=phiZeroR[i];
	hp2.phiZero[i]=phiZero1[i];
	phiZeroDiff[i]=hp1.phiZero[i]-hp2.phiZero[i];
	phiZeroSum[i]=hp1.phiZero[i]+hp2.phiZero[i];
	twoPhiZeroDiff[i]=2*(hp1.phiZero[i]-hp2.phiZero[i]);
	//	cout <<"we got charge 1: " << hp1.chargeType[i] <<" and second charge: " << hp2.chargeType[i] <<endl;
	normalizeAngle(phiZeroR[i]);
	normalizeAngle(phiZero1[i]);

	hp1.phiR[i]=phiR1[i];
	hp2.phiR[i]=phiRSum[i]-phiR1[i];

	//		cout <<" set phir1: " << hp1.phiR[i] << " phiR2: " << hp2.phiR[i] << " i: " << i <<endl;
	normalizeAngle(hp1.phiR[i]);
	normalizeAngle(hp2.phiR[i]);
	normalizeAngle(hp1.phiZero[i]);
	normalizeAngle(hp2.phiZero[i]);
	phiRDiff[i]=hp1.phiR[i]-hp2.phiR[i];
	twoPhiRDiff[i]=2*(hp1.phiR[i]-hp2.phiR[i]);

	normalizeAngle(phiRSum[i]);
	normalizeAngle(phiRDiff[i]);
	normalizeAngle(phiZeroSum[i]);
	normalizeAngle(phiZeroDiff[i]);
	normalizeAngle(twoPhiRDiff[i]);
	normalizeAngle(twoPhiZeroDiff[i]);

      }
  }

  //set cut flag for all events that don't have valid pi0, set flag to true to mask out any events that do not have pi0s
  void selectPi0Sig(bool onlyPi0Ev=false);
  //set cut flag for all events that don't have pi0 Bg
  void selectPi0Bg();

  //////put weighting based on MC here...
  void doWeighting(HadronQuadArray& hadQuad, float a1, float a2, float a3)
  {
    if(hadQuad.hp1.numPairs!=this->hp1.numPairs)
      {
	cout <<"weighting not possible!, different number of pairs... " << hadQuad.hp1.numPairs <<" vs " << this->hp1.numPairs<<endl;
	return;
	//	exit(0);
      }

    for(int i=0;i<hp1.numPairs;i++)
      {
	if(hadQuad.hp1.cut[i] || hadQuad.hp2.cut[i])
	  {
	    if(debugPrint)
	      {
		cout <<"base quad for weighting is cut... " <<endl;
		if(!hp1.cut[i]&&!hp2.cut[i] ) 
		  cout<<"but this quad is not cut!!" <<endl;
	      }
	    if(this->hp1.chargeType[i]==0 && this->hp2.chargeType[i]==0 && this->hp1.particleType[i]==0 && this->hp2.particleType[i]==0)
	      {
		if(!(hp1.cut[i] || hp2.cut[i] ))
		  {
		    //		    cout <<"weighee is cut!, qt weighee: " << hadQuad.qT[i] <<" qT weighted: "<< qT[i]<<endl;

		  }
	      }

	    ///let's not do this because the misid is so rare that the weighee is usually only cut if the weigted is cut. And there was a bug in thrustProj_mc leading to more stuff cut for the weighee...
	    this->hp1.cut[i]=1;
	    this->hp2.cut[i]=1;
	  }
	//	float kinFact=0.02+0.5*hadQuad.hp1.z[i];
	float kinFact=hadQuad.hp1.z[i];
	if(debugPrint)
	  {
	    cout <<"old weights: " << weight[i] <<" and " << weightZero[i] <<endl;
	    cout <<"weighing with " << (1+kinFact*(a1*cos(hadQuad.phiRSum[i])+a2*cos(hadQuad.phiRDiff[i])+a3*cos(hadQuad.twoPhiRDiff[i]))) << " and " << (1+kinFact*(a1*cos(hadQuad.phiZeroSum[i])+a2*cos(hadQuad.phiZeroDiff[i])+a3*cos(hadQuad.twoPhiZeroDiff[i]))) <<endl;
	  }

	//	weight[i]=1.0;
	//		weightZero[i]=1.0;
		weight[i]*=(1+kinFact*(a1*cos(hadQuad.phiRSum[i])+a2*cos(hadQuad.phiRDiff[i])+a3*cos(hadQuad.twoPhiRDiff[i])));
		weightZero[i]*=(1+kinFact*(a1*cos(hadQuad.phiZeroSum[i])+a2*cos(hadQuad.phiZeroDiff[i])+a3*cos(hadQuad.twoPhiZeroDiff[i])));
	if(debugPrint)
	  cout <<"weight: " << weight[i] <<" weightZero: "<< weightZero[i]<<endl;



	if(isnan(weight[i]) || isnan(weightZero[i]))
	  {
	    cout <<"had quad phirsum: " << hadQuad.phiRSum[i] <<" diff: " <<  hadQuad.phiRDiff[i] << " twoDiff: "<< hadQuad.twoPhiRDiff[i]<<endl;
	    cout <<"z: "<< hadQuad.hp1.z[i] <<endl;
	    cout <<"phir1 : " << hadQuad.hp1.phiR[i] << " phir2: "<< hadQuad.hp2.phiR[i]<<endl;
	    cout <<" weights are not numbers: " << weight[i] << " " << weightZero[i]<<endl;
	    cout <<"this phir: " << hp1.phiR[i] <<" and phir2: " << hp2.phiR[i] <<endl;
	    cout <<"this phirZero: " << hp1.phiZero[i] <<" and phir2: " << hp2.phiZero[i] <<endl;
	    cout <<"kinFact: " << kinFact <<" phiZeroSum[i] " << phiZeroSum[i] << " diff: "<< phiZeroDiff[i] <<endl;

	    cout <<"other phirZero: " << hadQuad.hp1.phiZero[i] <<" and phir2: " << hadQuad.hp2.phiZero[i] <<endl;
	    cout <<"kinFact: " << kinFact <<" other phiZeroSum[i] " << hadQuad.phiZeroSum[i] << " diff: "<< hadQuad.phiZeroDiff[i] <<endl;
	    weight[i]=1.0;
	    weightZero[i]=1.0;
	  }
      }
  }

//add pi to the angles
  void mixItUp();
void print();
//mix hp1 from one with hp2 of the other
  void mixEvent(HadronQuadArray& hadQuad);
  void numPairs(int& firstHemi, int& secondHemi);
  HadronQuadArray& operator=(HadronQuadArray rhs);

};



inline void HadronQuadArray::selectPi0Sig(bool onlyPi0Ev)
{
  hp1.selectPi0Sig(onlyPi0Ev);
  hp2.selectPi0Sig(onlyPi0Ev);
};

inline void HadronQuadArray::selectPi0Bg()
{
  hp1.selectPi0Bg();
  hp2.selectPi0Bg();
};
#endif
