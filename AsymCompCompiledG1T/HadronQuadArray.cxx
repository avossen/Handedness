#include "HadronQuadArray.h"


HadronQuadArray& HadronQuadArray::operator=(HadronQuadArray rhs)
{
  //hopefully using operator= of the hadron pair
  hp1=rhs.hp1;
  hp2=rhs.hp2;

  numHadQuads=rhs.numHadQuads;


  for(int i=0;i<numHadQuads;i++)
    {
      phiRSum[i]=rhs.phiRSum[i];
      phiR1[i]=rhs.phiR1[i];
      phiZeroR[i]=rhs.phiZeroR[i];
      phiZero1[i]=rhs.phiZero1[i];
      phiRDiff[i]=rhs.phiRDiff[i];

      twoPhiRDiff[i]=rhs.twoPhiRDiff[i];
      phiZeroDiff[i]=rhs.phiZeroDiff[i];

      phiZeroSum[i]=rhs.phiZeroSum[i];

      twoPhiZeroDiff[i]=rhs.phiZeroDiff[i];
      qT[i]=rhs.qT[i];
      weight[i]=rhs.weight[i];
      weightZero[i]=rhs.weightZero[i];

    }
  return *this;
};


//add pi to the angles
void HadronQuadArray::mixItUp()
{
  hp1.mixItUp();
  hp2.mixItUp();
  for(int i=0;i<hp1.numPairs;i++)
    {
      phiRDiff[i]=hp1.phiR[i]-hp2.phiR[i];
      twoPhiRDiff[i]=2*(hp1.phiR[i]-hp2.phiR[i]);
      phiRSum[i]=hp1.phiR[i]+hp2.phiR[i];
      normalizeAngle(hp1.phiR[i]);
      normalizeAngle(hp2.phiR[i]);
      normalizeAngle(phiRDiff[i]);
      normalizeAngle(twoPhiRDiff[i]);
      normalizeAngle(phiRSum[i]);
    }
};


void HadronQuadArray::numPairs(int& hemi1, int& hemi2)
{

  int pairCounter=0;
  int numPairs=this->hp1.numPairs;
  float phiR1_1=-100;
  for(int i=0;i<numPairs;i++)
    {
      //check if we saw all pairs
      if(phiR1_1!=this->phiR1[i] && i!=0)
	break;
      phiR1_1=this->phiR1[i];
      pairCounter++;
    }
  hemi1=pairCounter;
  if(hemi1>0)
    hemi2=numPairs/hemi1;

  //  cout <<" overall we have " << numPairs << " combinations " << hemi1 <<" * " << hemi2 <<endl;
}


//mix hp1 from one with hp2 of the other
void HadronQuadArray::mixEvent(HadronQuadArray& hadQuad)
{
  //  cout <<"mix hp1 : " << endl;
  //  this->hp2.print();
  //  cout <<" with: "<<endl;
  //  hadQuad.hp2.print();
  int numPairs1=this->hp1.numPairs;
  int numPairs2=hadQuad.hp1.numPairs;

  //needed because otherwise we start replacing fields that we mix with lateron
  HadronQuadArray tmpQuadArray(0,this->mMCFlag);

  ////since the second hadron is the inner loop when combining hadrons we can check when phiR1 is changing to get all the combinations...
  float phiR1_1=-100;
  float phiR1_2=-100;

  int pairCounter=0;
  //
  //    for(int i=0;i<numPairs1;i++)
  //    {
  //      cout <<" first pair phi :"<< this->hp1.phiR[i] << " and phi2: " << this->hp2.phiR[i] <<endl;
  //      cout <<"z : " << this->hp1.z[i]<<endl;
  //    }
  //
  //      for(int j=0;j<numPairs2;j++)
  //	{
  //	  cout <<"second pair phi: " << hadQuad.hp1.phiR[j] <<" and phi2: " << hadQuad.hp2.phiR[j]<<endl;
  //	  cout <<"z : " << hadQuad.hp2.z[j]<<endl;
  //	}
  //
  for(int i=0;i<numPairs1;i++)
    {
      //check if we saw all pairs
      if(phiR1_1!=this->hp1.phiR[i] && i!=0)
	break;
      phiR1_1=this->hp1.phiR[i];

      for(int j=0;j<numPairs2;j++)
	{
	  if(phiR1_2!=hadQuad.hp1.phiR[j]&& j!=0)
	    break;
	  phiR1_2=hadQuad.hp1.phiR[j];

	  //ZZ has many cominations and is as per now useless (would be eta/pi0 but in reality has a lot of pi0/pi0 combinations...
	  if(this->hp2.chargeType[i]==ZZ ||hadQuad.hp2.chargeType[j]==ZZ)
	    continue;
	  //take the second hadron pair (inner loop when building quads...)
	  tmpQuadArray.hp1.setElement(pairCounter,this->hp2,i);
	  tmpQuadArray.hp2.setElement(pairCounter,hadQuad.hp2,j);
	  

	  ///since both of these are 'R1's we have to flip the sign of one. In particular the one that we use as R1
	  //not so sure if that is the same for phiZero, after all this is an asymmetric system, but let's just do it...
	  float phiR1=this->hp2.phiR[i];
	  phiR1=(-1)*phiR1;
	  normalizeAngle(phiR1);
	  float phiR2=hadQuad.hp2.phiR[j];
	  //	  cout <<"mixing " << phiR1 << " with: " << phiR2 <<endl;
	  float phiZero1=this->hp2.phiZero[i];
	  phiZero1=(-1)*phiZero1;
	  normalizeAngle(phiZero1);
	  float phiZero2=hadQuad.hp2.phiZero[j];
	  //	  cout <<"mixing phir1 :  " <<phiR1<< " and " <<phiR2<<endl;
	  tmpQuadArray.hp1.phiR[pairCounter]=phiR1;
	  tmpQuadArray.hp2.phiR[pairCounter]=phiR2;
	  tmpQuadArray.phiR1[pairCounter]=phiR1;
	  tmpQuadArray.phiRSum[pairCounter]=phiR1+phiR2;
	  tmpQuadArray.phiRDiff[pairCounter]=phiR1-phiR2;
	  tmpQuadArray.twoPhiRDiff[pairCounter]=2*(phiR1-phiR2);
	  //	  cout <<"phiR1: " << phiR1 << " phiR2: " << phiR2 <<" sum: " << phiR1+phiR2 <<" diff: "<< phiR1-phiR2 <<" two diff: " << 2*(phiR1-phiR2)<<endl;
	  tmpQuadArray.phiZero1[pairCounter]=phiZero1;
	  //shouldn't be used for fits anyways
	  tmpQuadArray.phiZeroR[pairCounter]=this->phiZeroR[i];
	  tmpQuadArray.phiZeroSum[pairCounter]=phiZero1+phiZero2;
	  tmpQuadArray.hp1.phiZero[pairCounter]=phiZero1;
	  tmpQuadArray.hp2.phiZero[pairCounter]=phiZero2;
	  tmpQuadArray.phiZeroDiff[pairCounter]=phiZero1-phiZero2;
	  tmpQuadArray.twoPhiZeroDiff[pairCounter]=2*(phiZero1-phiZero2);

	  normalizeAngle(tmpQuadArray.phiRSum[pairCounter]);
	  normalizeAngle(tmpQuadArray.phiRDiff[pairCounter]);
	  normalizeAngle(tmpQuadArray.twoPhiRDiff[pairCounter]);

	  normalizeAngle(tmpQuadArray.phiZeroSum[pairCounter]);
	  normalizeAngle(tmpQuadArray.phiZeroDiff[pairCounter]);
	  normalizeAngle(tmpQuadArray.twoPhiZeroDiff[pairCounter]);

	  if(pairCounter<this->hp1.numPairs)
	    {
	      tmpQuadArray.qT[pairCounter]=this->qT[i];
	      tmpQuadArray.weight[pairCounter]=this->weight[i];
	      tmpQuadArray.weightZero[pairCounter]=this->weightZero[i];
	    }
	  else
	    {
	      if(pairCounter<hadQuad.hp1.numPairs)
		{
		  tmpQuadArray.qT[pairCounter]=hadQuad.qT[j];
		  tmpQuadArray.weight[pairCounter]=hadQuad.weight[j];
		  tmpQuadArray.weightZero[pairCounter]=hadQuad.weightZero[j];
		}
	      else
		{
		  break;
		}
	    }

	  pairCounter++;
	}
    }


  tmpQuadArray.hp1.numPairs=pairCounter;
  tmpQuadArray.hp2.numPairs=pairCounter;

  //not necessary to set number, already done in constructor...
  tmpQuadArray.numHadQuads=pairCounter;

  //  cout<<endl<<endl;
  //  cout <<"numHadQuads: " << tmpQuadArray.numHadQuads<<endl;
  (*this)=tmpQuadArray;
    //  this->print();
  //  cout <<"setting num comb to " << pairCounter << endl;
};


void HadronQuadArray::print()
{
  cout <<"---------------------------------Hadron Quadruple -------------------------"<<endl;
  cout <<"first Hadron: " << endl;
  hp1.print();
  cout <<"second Hadron: " << endl;
  hp2.print();
  cout <<" phiRSum | phiR1 | phiZeroR | phiZero1 | phiRDiff | twoPhirRDiff | phiZeroDiff | phiZeroSum | twoPhiZeroDiff | qT | weight | weightZero | " <<endl;
  for(int i=0;i<hp1.numPairs;i++)
    {
      cout << phiRSum[i] << " | " << phiR1[i] << " | " << phiZeroR[i] << " | " << phiZero1[i] << " | " << phiRDiff[i] << " | " << twoPhiRDiff[i] << " | " << phiZeroDiff[i] << " | " << phiZeroSum[i] << " | " << twoPhiZeroDiff[i] << " | " << qT[i] << " | " << weight[i] << " | " << weightZero[i] << " | " << endl;
    }

}
