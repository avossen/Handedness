#include "HadronPairArrays.h"

HadronPairArray& HadronPairArray::operator=(HadronPairArray rhs)
  {
    hadPairNum=rhs.hadPairNum;
    numPairs=rhs.numPairs;
    mMCFlag=rhs.mMCFlag;
    addendum=rhs.addendum;
    for(int i=0;i<rhs.numPairs;i++)
      {
	//use this function so that all and single elements change the same...
	setSingleElement(i,rhs,i);
	/*(	cut[i]=rhs.cut[i];
	z[i]=rhs.z[i];
	mass[i]=rhs.mass[i];
	phiR[i]=rhs.phiR[i];
	phiZero[i]=rhs.phiZero[i];
	thrustProj1[i]=rhs.thrustProj1[i];
	thrustProj2[i]=rhs.thrustProj2[i];
	theta1[i]=rhs.theta1[i];
	theta2[i]=rhs.theta2[i];
	chargeType[i]=rhs.chargeType[i];
	particleType[i]=rhs.particleType[i];
	pi0Sig[i]=rhs.pi0Sig[i];
	pi0Bg[i]=rhs.pi0Bg[i];
	pi0mass1[i]=rhs.pi0mass1[i];
	pi0mass2[i]=rhs.pi0mass2[i];*/
      }
    return *this;
  };



//the the element 'pairCounter' of this pair to the element 'index' of the argument hadpair. This is used for event mixing...

void HadronPairArray::setSingleElement(int pairCounter,HadronPairArray& hp,int index)
{

  z[pairCounter]=hp.z[index];
  zRatio[pairCounter]=hp.zRatio[index];
  z1[pairCounter]=hp.z1[index];
  z2[pairCounter]=hp.z2[index];
  mass[pairCounter]=hp.mass[index];
  phiR[pairCounter]=hp.phiR[index];
  phiZero[pairCounter]=hp.phiZero[index];
  thrustProj1[pairCounter]=hp.thrustProj1[index];
  thrustProj2[pairCounter]=hp.thrustProj2[index];

  decayTheta[pairCounter]=hp.decayTheta[pairCounter];

  theta1[pairCounter]=hp.theta1[index];
  theta2[pairCounter]=hp.theta2[index];

  phi1[pairCounter]=hp.phi1[index];
  phi2[pairCounter]=hp.phi2[index];

  pi0mass1[pairCounter]=hp.pi0mass1[index];
  pi0mass2[pairCounter]=hp.pi0mass2[index];
  chargeType[pairCounter]=hp.chargeType[index];
  particleType[pairCounter]=hp.particleType[index];

  cut[pairCounter]=hp.cut[index];
  pi0Sig[pairCounter]=hp.pi0Sig[index];
  pi0Bg[pairCounter]=hp.pi0Bg[index];
  hadOpening[pairCounter]=hp.hadOpening[index];

  motherGenId1[pairCounter]=hp.motherGenId1[index];
  motherGenId2[pairCounter]=hp.motherGenId2[index];

};


void HadronPairArray::setElement(int pairCounter,HadronPairArray& hp,int index)
{
    mMCFlag=hp.mMCFlag;
    addendum=hp.addendum;
    setSingleElement(pairCounter,hp,index);
};

void HadronPairArray::print()
{

  cout <<"Hadron Index | z |  zRatio | mass | phiR |phiZero | thrustProj1 | thrustProj2 | theta1 | theta2 | chargeType | particleType | cut |  hadOpening " <<endl;
  for(int i=0;i<numPairs;i++)
    {
      cout <<i<<" | " << z[i] << " | " << zRatio[i] << " | " << mass[i] << " | " << phiR[i] << " | " << phiZero[i] << " | " << thrustProj1[i] << " | " <<thrustProj2[i] << " | " <<theta1[i] << " | " <<theta2[i] << " | " <<chargeType[i] << " | " <<particleType[i] << " | " <<cut[i] << " | " <<hadOpening[i] << " | " <<endl;
    }



}
//effect of changing charges is effectively reversing the direction of R which equals a shift by Pi
void HadronPairArray::mixItUp()
{
  for(int i=0;i<numPairs;i++)
    {
      //      cout <<"mixing pair  " << i <<endl;
      //shift by pi...
      if(rand() % 100 <50)
	{
	  //	  cout << " shifting r: " << phiR[i] << endl;
	  phiR[i]+=TMath::Pi();
	  //	  cout << " now r: " << phiR[i] << endl;
	}
    }


};


void HadronPairArray::selectPi0Sig(bool onlyPi0Ev)
{
  for(int i=0;i<numPairs;i++)
    {
      if(!pi0Sig[i] && ( onlyPi0Ev || (pi0mass1[i] > 0 || pi0mass2[i] >0)))
	cut[i]=1;
    }
};


void HadronPairArray::selectPi0Bg()
{
  for(int i=0;i<numPairs;i++)
    {
      if(!pi0Bg[i])
	cut[i]=1;
    }
};


