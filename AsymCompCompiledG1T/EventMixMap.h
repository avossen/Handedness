#ifndef EVENT_MIX_MAP_HH
#define EVENT_MIX_MAP_HH

#include "HadronQuadArray.h"
#include "MEvent.h"

/*
provide interface to get the right hadron quads for event mixing so that the thrust axis' are somewhat aligned

inherits from ReaderBase for normalizeAngle
*/
class EventMixMap:public ReaderBase 
{
 public:
  EventMixMap(int mNumThrustThetaBins, int mNumThrustPhiBins, int dataMCFlag):queueLength(1)
    {
      numThrustThetaBins=mNumThrustThetaBins;
      numThrustPhiBins=mNumThrustPhiBins;
      hadQuads=new HadronQuadArray***[numThrustThetaBins];
      for(int iT=0;iT<numThrustThetaBins;iT++)
	{
	  hadQuads[iT]=new HadronQuadArray**[numThrustPhiBins];
	  for(int iP=0;iP<numThrustPhiBins;iP++)
	    {
	      hadQuads[iT][iP]=new HadronQuadArray*[queueLength];
	      for(int iB=0;iB<queueLength;iB++)
		{
		  hadQuads[iT][iP][iB]=new HadronQuadArray(0,dataMCFlag);
		}
	    }
	}
      binPointer=allocateArray<int>(numThrustThetaBins,numThrustPhiBins);
      binFilled=allocateArray<int>(numThrustThetaBins,numThrustPhiBins);
    };
  ~EventMixMap()
    {
      for(int iT=0;iT<numThrustThetaBins;iT++)
	{
	  for(int iP=0;iP<numThrustPhiBins;iP++)
	    {
	      for(int iB=0;iB<queueLength;iB++)
		{
		  delete hadQuads[iT][iP][iB];
		}
	      delete[] hadQuads[iT][iP];
	    }
	  delete[] hadQuads[iT];
	}
      delete hadQuads;
    }

  void addHadQuadArray(HadronQuadArray& quad, MEvent& event);
  HadronQuadArray** getHadQuadArray(MEvent& event, int& numQuads);
 protected:
  int numThrustThetaBins;
  int numThrustPhiBins;
  HadronQuadArray**** hadQuads;
  int** binFilled;
  int** binPointer;
  int queueLength;

};

#endif
