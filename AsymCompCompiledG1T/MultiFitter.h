#ifndef MULTI_FITTER_H
#define MULTI_FITTER_H
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "BinningContainer.h"
#include "NamedExp.h"
#include "MEvent.h"
#include "TwoHadAsymsCommons.h"
#include "TMath.h"
#include "FitResults.h"
#include "TFile.h"
#include "HadronQuadArray.h" 
#include "stdlib.h"
#include <cstdlib>
#include <iostream>
#include <fstream>





/// DO NOT FORGET to change the number of kinematic binnings in the cxx file!!!
/// and if NEW FIELDS ARE ADDED, they also have to be added in the '=' operator implementation of HadronQuad and PairArray..
enum binningType{binType_m_m, binType_z_z,binType_z_m, binType_m_z, binType_labTheta_z, binType_kinFact_z, binType_zOnly, binType_mOnly, binType_labThetaOnly, binType_qTOnly, binType_kinFactOnly, binType_hadOpeningOnly, binType_ThrustThetaPhi,binType_ThrustPhiTheta, binType_multOnly, binType_EmissOnly, binType_ThrustOnly,binType_sinDecThetaOnly,binType_cosDecThetaOnly,binType_end};
enum quadType{quadPN, quadPZ_ZN, quadPN_PZ, quadPN_ZN,quadPZ, quadZN,quadPP_NN, quadUnknownCharge, quadTypeEnd};
enum plotType{plotType_2D, plotType_1D, plotType_DR,plotType_RF,plotType_end};

class MultiFitter: public ReaderBase, NamedExp, BinningContainer//for the normalize angle
{
 public:
  MultiFitter(const char* filenameBase,string nameAdd, int exNr, bool onRes, bool uds, bool charm,bool mc,int numAngBins=16):NamedExp(filenameBase,nameAdd,exNr,onRes,uds,charm,mc),numAngBins(numAngBins), minCounts(12)
    {
      zeroBin=0;
      rFile.mkdir("fitHistos");
      rFile.cd("fitHistos");
      rFile.mkdir("oneDHistos");
      rFile.mkdir("twoDHistos");
      rFile.mkdir("RFitDHistos");
      rFile.cd();
      rFile.mkdir("graphs");
      for(int i=0;i<numAngBins;i++)
	{
	  binningAng.push_back((i+1)*2*TMath::Pi()/((float)numAngBins));
	  cout <<"binning ang: " << (i+1)*2*TMath::Pi()/((float)numAngBins) <<endl;
	}
      loadBinning(binningM, binningZ);

      maxKinBins=binningZ[0].size();
      if(binningM[0].size()>binningZ[0].size())
	maxKinBins=binningM[0].size();
      loadThetaBinnings();
      if(binningLabTheta.size()>maxKinBins)
	maxKinBins=binningLabTheta.size();
      if(binningKinFact.size()>maxKinBins)
	maxKinBins=binningKinFact.size();
      if(binningQt.size()>maxKinBins)
	maxKinBins=binningQt.size();
      if(binningHadOpen.size()>maxKinBins)
	maxKinBins=binningHadOpen.size();
      if(binningCmsThrustTheta.size()>maxKinBins)
	maxKinBins=binningCmsThrustTheta.size();
      if(binningCmsThrustPhi.size()>maxKinBins)
	maxKinBins=binningCmsThrustPhi.size();
      if(binningMult.size()>maxKinBins)
	maxKinBins=binningMult.size();
      if(binningThrust.size()>maxKinBins)
	maxKinBins=binningThrust.size();
      if(binningEmiss.size()>maxKinBins)
	maxKinBins=binningEmiss.size();
      if(binningSinDecTheta.size()>maxKinBins)
	maxKinBins=binningSinDecTheta.size();
      if(binningCosDecTheta.size()>maxKinBins)
	maxKinBins=binningCosDecTheta.size();

      hChi2OverNdf=new TH1D("chi2OverNdf","chi2OverNdf",100,0,6);
      hOneDVsTwoDA1=new TH1D("oneDVs2DA1","oneDVs2DA1",100,-3,3);
      hOneDVsTwoDA2=new TH1D("oneDVs2DA2","oneDVs2DA2",100,-3,3);
      hOneDVsTwoDA3=new TH1D("oneDVs2DA3","oneDVs2DA3",100,-3,3);

      setBinningMap();
      fitResults=new FitResults[numKinematicBinning*NumCharges*maxKinBins*maxKinBins];
      fitResults1D=new FitResults[numKinematicBinning*NumCharges*maxKinBins*maxKinBins];
      fitResultsDR=new FitResults[numKinematicBinning*NumCharges*maxKinBins*maxKinBins];


      fitResVect.push_back(fitResults);
      fitResVect.push_back(fitResults1D);
      fitResVect.push_back(fitResultsDR);
      for(int i=0;i<numKinematicBinning*NumCharges*maxKinBins*maxKinBins;i++)
	{
	  for(vector<FitResults*>::iterator it=fitResVect.begin();it!=fitResVect.end();it++)
	    {
	      (*it)[i].meanKinBin1=0.0;
	      (*it)[i].meanKinBin2=0.0;
	      (*it)[i].chi2=0.0;
	      (*it)[i].chi2OverNdf=0.0;
	      (*it)[i].A1=0.0;
	      (*it)[i].eA1=100000000.0;
	      (*it)[i].eA1=100000000.0;
	      (*it)[i].A2=0.0;
	      (*it)[i].eA2=100000000.0;
	      (*it)[i].A3=0.0;
	      (*it)[i].eA3=100000000.0;
	    }
	}

      cout <<"allocating RooDataSets.." <<endl;
      rooPhi1=new RooRealVar("phi1","phi1",0,2*TMath::Pi());
      rooZ1=new RooRealVar("z1","z1",0,1.0);
      rooM1=new RooRealVar("m1","m1",0,3.0);
      rooZ2=new RooRealVar("z2","z2",0,1.0);
      rooM2=new RooRealVar("m2","m2",0,3.0);

      rooPhi2=new RooRealVar("phi2","phi2",0,2*TMath::Pi());
      rooPhiSum=new RooRealVar("phiSum","phiSum",0,2*TMath::Pi());
      rooTheta1=new RooRealVar("theta1","theta1",0,TMath::Pi());
      rooTheta2=new RooRealVar("theta2","theta2",0,TMath::Pi());

      rooPhi1->setBins(16);
      rooPhi2->setBins(16);
      rooTheta1->setBins(8);
      rooTheta2->setBins(8);
      
      rooBin = new RooCategory("iBin","iBin");
      rooKin1=new RooCategory("iKin1","iKin1");
      rooKin2=new RooCategory("iKin2","iKin2");
    
      for(int i=0;i<maxKinBins;i++)
	{
	  stringstream ss; 
	  ss <<"kin1Bin";
	  ss << i;
	  rooKin1->defineType(ss.str().c_str(),i);
	  cout <<"defining type " << ss.str().c_str()<<endl;
	  stringstream ss2; 
	  ss2 <<"kin2Bin";
	  ss2 << i;
	  rooKin2->defineType(ss2.str().c_str(),i);
	  cout <<"defining type " << ss2.str().c_str()<<endl;
	}
      for(int iBt=0;iBt<numKinematicBinning;iBt++)
	{
	  //try to keep down the rooHist
	  if(iBt>binType_mOnly)
	    continue;
	  stringstream ss; 
	  ss <<"kinBin";
	  ss << iBt;
	  rooBin->defineType(ss.str().c_str(),iBt);
	  cout <<"defining type " << ss.str().c_str()<<" index: "<< rooBin->getIndex()<<endl;
	}

      unbinnedData=new RooDataSet("data","data",RooArgSet(*rooZ1,*rooM1,*rooZ2,*rooM2,*rooTheta1,*rooTheta2,*rooPhi1,*rooPhi2));
      cout <<"constructiong histogram with dimensions " << maxKinBins <<", " << numKinematicBinning <<", "<< rooTheta1->getBins() <<", " << rooTheta2->getBins() <<", " << rooPhi1->getBins()<<", " << rooPhi2->getBins()<<endl;

      //      binnedData=new RooDataHist("binnedData","binnedData",RooArgSet(*rooTheta1,*rooTheta2,*rooPhi1,*rooPhi2));
      //      cout <<"binned data has " << binnedData->numEntries()<<" bins" << endl;
      //      binnedData=new RooDataHist("binnedData","binnedData",RooArgSet(*rooBin,*rooKin1,*rooKin2,*rooTheta1,*rooTheta2,*rooPhi1,*rooPhi2));
      //      cout <<"now binned data has " << binnedData->numEntries()<<" bins" << endl;
      //      binnedData=new RooDataHist("binnedData","binnedData",RooArgSet(*rooTheta1,*rooTheta2,*rooPhi1,*rooPhi2));
      cout <<"done with construction..." <<endl;
      unbinnedDataIff=new RooDataSet("dataIff","dataIff",RooArgSet(*rooBin,*rooKin1,*rooKin2,*rooPhiSum));

      cout <<"allocating " << numKinematicBinning <<" * " << NumCharge << " * " << maxKinBins <<" * " << maxKinBins <<" * " <<numAngBins <<" * " <<numAngBins<<endl;

      counts=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numAngBins,numAngBins);
      countsSumR=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numAngBins);//for comparison
      countsDiffR=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numAngBins);//for comparison
      countsTwoDiffR=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numAngBins);//for comparison

      counts_wSq=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numAngBins,numAngBins);
      countsSumR_wSq=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numAngBins);//for comparison
      countsDiffR_wSq=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numAngBins);//for comparison
      countsTwoDiffR_wSq=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numAngBins);//for comparison


      rawCounts=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numAngBins,numAngBins);
      rawCountsSumR=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numAngBins);//for comparison
      rawCountsDiffR=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numAngBins);//for comparison
      rawCountsTwoDiffR=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numAngBins);//for comparison

      meanValues_kin1=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins);
      meanValues_kin2=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins);


      eventCounts=new TH1D***[numKinematicBinning];
      char buffer[200];
      for(int iB=0;iB<numKinematicBinning;iB++)
	{
	  eventCounts[iB]=new TH1D**[NumCharges];
	  for(int iC=0;iC<=quadPN;iC++)
	    {
	      eventCounts[iB][iC]=new TH1D*[maxKinMap[iB].first+1];
	      for(int iK=0;iK<=maxKinMap[iB].first;iK++)
		{
		  sprintf(buffer,"%s_eventCounts_%s",filenameBase,getBinName(iB,iC,iK,0).c_str());
		  cout <<" creating eventCountHisto: "<< buffer << " bt: "<< iB << " kin: " << iK <<endl;
  		  eventCounts[iB][iC][iK]=new TH1D(buffer,buffer,maxKinMap[iB].second,0,maxKinMap[iB].second);
		}
	    }

	}

      openXCheckFiles();

    };

    //
    //    void setFitReuslt();
    void openXCheckFiles();
    void addHadQuadArray(HadronQuadArray* hq, MEvent& event,bool usePhiZero=false, bool print=false);
    void setBinningMap();
    void doFits(MultiFitter* mfMix=0);
    void doDRFits(double** locCounts, double& As, double& AsErr, string binName, bool iffLike=false);
    void doDRFitsG1T(double** locCounts, double& As, double& AsErr, string binName);
    void savePlot(int binningType,int chargeType, plotType mPlotType=plotType_2D);
    void getIntAsymmetry(float a[3], float ea[3],int binningType,int chargeType, bool save1D=false);
    void saveAsymmetries(plotType mPlotType=plotType_2D);
    void setName(string s);
    string getName();
    FitResults* fitResults;
    FitResults* fitResults1D;
    FitResults* fitResultsDR;
    vector<FitResults*> fitResVect;

    //maps the type of binning to two bin pointers
    vector< pair<int*,int*> > binningMap;
    vector< pair<float*,float*> > meanMap;
    vector< pair<int, int> > maxKinMap;
    inline int getResIdx(int binningType,int chargeType, int firstKinBin, int secondKinBin)
      {
	return binningType*NumCharges*maxKinBins*maxKinBins+chargeType*maxKinBins*maxKinBins+firstKinBin*maxKinBins+secondKinBin;
      }
 protected:
    bool checkMinCounts(double** counts);
    ofstream ****xCheckEventLists;
    ofstream* fullXCheckEventList;


    //to reorder arrays used for fitting such that the x values are ascending and do not wrap around
    //should work because y is moved as well
    void reorder(float* mX, float* mY, float* mYErr, int numBins);
    static const int numKinematicBinning;
    static const int NumCharges;
    static const int numParticles;
    //to differentiate between different fit results...

    TH1D* hChi2OverNdf;
    TH1D* hOneDVsTwoDA1;
    TH1D* hOneDVsTwoDA2;
    TH1D* hOneDVsTwoDA3;


    //give negative values for first or second bin if they should not be part of the name
    string getBinName(int binningType,int chargeType, int firstBin, int secondBin);
    string getBinningName(int binningType, int chargeType);
    string getXAxisName(int binningType);
    void loadThetaBinnings();
 public:
    vector<float> binningAng;

    vector<float> binningLabTheta;
    vector<float> binningKinFact;
    vector<float> binningQt;
    vector<float> binningHadOpen;
    vector<float> binningCmsThrustTheta;
    vector<float> binningCmsThrustPhi;
    vector<float> binningMult;
    vector<float> binningThrust;
    vector<float> binningEmiss;
    vector<float> binningSinDecTheta;
    vector<float> binningCosDecTheta;


 protected:
    int zbin1;
    int zbin2;
    int mbin1;
    int mbin2;
    //always zero
    int zeroBin;
    int labThetaBin;
    int cmsThrustThetaBin;
    int cmsThrustPhiBin;
    int multBin;
    int kinFactorBin;
    int qTBin;
    int hadOpenBin1;
    int hadOpenBin2;
    int thrustBin;
    int eMissBin;
    int sinDecThetaBin1;
    int sinDecThetaBin2;

    int cosDecThetaBin1;
    int cosDecThetaBin2;

    float hadronOpening1;
    float hadronOpening2;


    float z1;
    float z2;
    float m1;
    float m2;
    float labTheta;
    float multiplicity;
    float kinFactor;
    float cmsThrustTheta;
    float cmsThrustPhi;
    float qT;
    float thrust;
    float Emiss;
    float sinDecTheta1;
    float sinDecTheta2;
    float cosDecTheta1;
    float cosDecTheta2;

    float decayTheta1;
    float decayTheta2;

    TH1D**** eventCounts;


    unsigned int minCounts;
    unsigned int numAngBins;
    unsigned int maxKinBins;
    double****** counts;
    double****** counts_wSq;
    double****** rawCounts;
    //the below are just for x-check
    double***** countsSumR;
    double***** countsSumR_wSq;
    double***** countsDiffR;
    double***** countsDiffR_wSq;
    double***** countsTwoDiffR;
    double***** countsTwoDiffR_wSq;
    double***** rawCountsSumR;
    double***** rawCountsDiffR;
    double***** rawCountsTwoDiffR;

    double**** meanValues_kin1;
    double**** meanValues_kin2;
    RooDataSet* unbinnedData;
    RooDataSet* unbinnedDataIff;
    RooDataHist* binnedData;


    RooRealVar* rooPhiSum;
    RooRealVar* rooPhi1;
    RooRealVar* rooPhi2;
    RooRealVar* rooTheta1;
    RooRealVar* rooTheta2;
    RooRealVar* rooZ1;
    RooRealVar* rooM1;
    RooRealVar* rooZ2;
    RooRealVar* rooM2;

    RooCategory* rooKin1;
    RooCategory* rooKin2;
    RooCategory* rooBin;


    RooRealVar* rooCollPhi;

    RooRealVar* rooPhi0;



};

inline string MultiFitter::getName()
{
  return nameAddition;
};

inline void MultiFitter::setName(string s)
{
  nameAddition=s;
};

#endif
