#ifndef ACCEPTANCE_MAP_H
#define ACCEPTANCE_MAP_H
#include "NamedExp.h"
#include "MEvent.h"
#include <sstream>
#include "TH2D.h"
#include "TwoHadAsymsCommons.h"
#include "MultiFitter.h"
class AcceptanceMap: public NamedExp,ReaderBase
{
 public:
  AcceptanceMap(const char* filenameBase,string nameAdd, int exNr, bool onRes, bool uds, bool charm,bool mc, MultiFitter* mfBinning):NamedExp(filenameBase,nameAdd,exNr,onRes,uds,charm,mc),m_weightingType(MZ), m_mfBinning(mfBinning)
    {
      maxOpenAngle=acos(0.8);
      maxZValue=1.1;
      maxMValue=2.0;

      stringstream ss;

      if(mc)
	{
	  ss<<"Mc_";
	}
      else
	{
	  ss<<"_data_";
	}
      if(uds)
	ss<<"_uds_";
      else
	{
	  if(charm)
	    ss<<"_charm_";
	}
      cout <<"constructing histos.. " <<endl;
      //theta,phi map
      rAcceptanceMapsPosDir=new TH2D*[numThetaThrustBins];
      rAcceptanceMapsNegDir=new TH2D*[numThetaThrustBins];
      rAcceptanceMaps_MZ1=new TH2D*[numMBins];
      rAcceptanceMaps_MZ2=new TH2D*[numMBins];


      binningSingleZ.push_back(0.15);
      binningSingleZ.push_back(0.2);
      binningSingleZ.push_back(0.3);
      binningSingleZ.push_back(0.4);
      binningSingleZ.push_back(0.5);
      binningSingleZ.push_back(0.6);
      binningSingleZ.push_back(1.5);

      acceptanceMap=new TH2D*[binningSingleZ.size()];

      for(int i=1;i<20;i++)
	{
	  thrustThetaBinning.push_back(i*TMath::Pi()/(float)16);
	}
      for(int i=1;i<20;i++)
	{
	  thrustPhiBinning.push_back(2*TMath::Pi()/(float)16*i);
	}

      char buffer[100];

      arrAccMap_MZ1=allocateArray<double>(thrustThetaBinning.size(),thrustPhiBinning.size(),m_mfBinning->binningM[0].size(),m_mfBinning->binningZ[0].size(), numPhiRBins);
      arrAccMap_MZ2=allocateArray<double>(thrustThetaBinning.size(),thrustPhiBinning.size(),m_mfBinning->binningM[0].size(),m_mfBinning->binningZ[0].size(), numPhiRBins);


      for(int iM=0;iM<numMBins;iM++)
	{
	  sprintf(buffer,"%d",iM);
	  rAcceptanceMaps_MZ1[iM]=new TH2D((string("rAccMapMZ1_")+buffer+nameAdd).c_str(),(string("rAccMapMZ1_")+buffer+nameAdd).c_str(),numPhiRBins,0,2*TMath::Pi(),numZBins,0,maxZValue);
	  rAcceptanceMaps_MZ2[iM]=new TH2D((string("rAccMapMZ2_")+buffer+nameAdd).c_str(),(string("rAccMapMZ2_")+buffer+nameAdd).c_str(),numPhiRBins,0,2*TMath::Pi(),numZBins,0,maxZValue);
	}

      for(int iC=0;iC<numThetaThrustBins;iC++)
	{
	  sprintf(buffer,"%d",iC);
	  rAcceptanceMapsPosDir[iC]=new TH2D((string("rAccMapPosDir_")+buffer+nameAdd).c_str(),(string("rAccMap_")+buffer+nameAdd).c_str(),numPhiRBins,0,2*TMath::Pi(),numThetaOpenBins,0,maxOpenAngle);
	  rAcceptanceMapsNegDir[iC]=new TH2D((string("rAccMapnegDir_")+buffer+nameAdd).c_str(),(string("rAccMap_")+buffer+nameAdd).c_str(),numPhiRBins,0,2*TMath::Pi(),numThetaOpenBins,0,maxOpenAngle);
	}
    
      for(int iZ=0;iZ<binningSingleZ.size();iZ++)
	{	
	  sprintf(buffer,"zBin_%d",iZ);		       
	  acceptanceMap[iZ]=new TH2D((string("accMap")+buffer+nameAdd).c_str(),(string("accMap")+buffer+nameAdd).c_str(),numThetaBins,0,TMath::Pi(),numPhiBins,0,2*TMath::Pi());
	}
      acceptanceMapZCut=new TH2D((string("accMapZCut")+nameAdd).c_str(),(string("accMapZCut")+nameAdd).c_str(),numThetaBins,0,TMath::Pi(),numPhiBins,0,2*TMath::Pi());
      acceptanceMapThrust=new TH2D((string("accMapThrust")+nameAdd).c_str(),(string("accMapThrust")+nameAdd).c_str(),numThetaBins,0,TMath::Pi(),numPhiBins,0,2*TMath::Pi());
      phiR1VsThrustTheta=new TH2D((string("phiR1VsThrustTheta")+nameAdd).c_str(),(string("phiR1VsThrustTheta")+nameAdd).c_str(),numThetaBins,0,TMath::Pi(),numPhiBins,0,2*TMath::Pi());
      phiRDiffVsThrustTheta=new TH2D((string("phiRDiffVsThrustTheta")+nameAdd).c_str(),(string("phiRDiffVsThrustTheta")+nameAdd).c_str(),numThetaBins,0,TMath::Pi(),numPhiBins,0,2*TMath::Pi());
      phiRTwoDiffVsThrustTheta=new TH2D((string("phiRTwoDiffVsThrustTheta")+nameAdd).c_str(),(string("phiRTwoDiffVsThrustTheta")+nameAdd).c_str(),numThetaBins,0,TMath::Pi(),numPhiBins,0,2*TMath::Pi());
      phiRSumVsThrustTheta=new TH2D((string("phiRSumVsThrustTheta")+nameAdd).c_str(),(string("phiRSumVsThrustTheta")+nameAdd).c_str(),numThetaBins,0,TMath::Pi(),numPhiBins,0,2*TMath::Pi());
      cout <<"done " <<endl;
}
    //
    WeightingType m_weightingType;
    void doWeighting(HadronQuadArray& hq, MEvent& event);
    void addHadQuadArray(HadronQuadArray* hq, MEvent& event);
    void normalize();
    void save();
    void saveMZHisto(int mbin, int zbin);

 protected:
    float maxOpenAngle;
    float maxZValue;
    float maxMValue;
    TH2D** rAcceptanceMapsPosDir;
    TH2D** rAcceptanceMapsNegDir;

    TH2D** rAcceptanceMaps_MZ1;
    TH2D** rAcceptanceMaps_MZ2;

    //shouldn't need two, but easier to debug....
    double***** arrAccMap_MZ1;
    double***** arrAccMap_MZ2;

    MultiFitter* m_mfBinning;
    TH2D** acceptanceMap;
    TH2D* acceptanceMapZCut;
    TH2D* acceptanceMapThrust;
    TH2D* phiR1VsThrustTheta;
    TH2D* phiRDiffVsThrustTheta;
    TH2D* phiRDiffVsThrustPhi;
    TH2D* phiRTwoDiffVsThrustTheta;
    TH2D* phiRSumVsThrustTheta;

    static const int numThetaBins;
    static const int numPhiBins;
    static const int numThetaThrustBins;
    static const int numPhiRBins;
    static const int numThetaOpenBins;
    static const int numMBins;
    static const int numZBins;

    vector<float> thrustThetaBinning;
    vector<float> thrustPhiBinning;
    vector<float> binningSingleZ;

};




#endif
