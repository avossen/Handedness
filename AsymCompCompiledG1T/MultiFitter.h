#ifndef MULTI_FITTER_H
#define MULTI_FITTER_H
#include "NamedExp.h"
#include "MEvent.h"
#include "TwoHadAsymsCommons.h"
#include "TMath.h"
#include "FitResults.h"
#include "TFile.h"
#include "HadronQuadArray.h" 
#include "stdlib.h"
#include <cstdlib>

#define NumParticles 7
enum binningType{binType_m_m, binType_z_z,binType_z_m, binType_m_z, binType_labTheta_z, binType_kinFact_z, binType_zOnly, binType_mOnly, binType_labThetaOnly, binType_qTOnly, binType_kinFactOnly, binType_hadOpeningOnly, binType_ThrustThetaPhi,binType_ThrustPhiTheta, binType_multOnly, binType_EmissOnly, binType_ThrustOnly,binType_end};
enum quadType{quadPN, quadPZ_ZN, quadPN_PZ, quadPN_ZN,quadPZ, quadZN,quadPP_NN, quadUnknownCharge, quadTypeEnd};
enum plotType{plotType_2D, plotType_1D, plotType_DR,plotType_end};

class MultiFitter: public ReaderBase, NamedExp//for the normalize angle
{
 public:
  MultiFitter(const char* filenameBase,string nameAdd, int exNr, bool onRes, bool uds, bool charm,bool mc,int numAngBins=16):NamedExp(filenameBase,nameAdd,exNr,onRes,uds,charm,mc),numAngBins(numAngBins), minCounts(20)
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

    };

    //
    //    void setFitReuslt();

    void addHadQuadArray(HadronQuadArray* hq, MEvent& event,bool usePhiZero=false);
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
    string getXAxisName(int binningType);
    void loadThetaBinnings();
 public:
    vector<float> binningAng;
    vector<float> binningM[NumParticles];
    vector<float> binningZ[NumParticles];
    vector<float> binningLabTheta;
    vector<float> binningKinFact;
    vector<float> binningQt;
    vector<float> binningHadOpen;
    vector<float> binningCmsThrustTheta;
    vector<float> binningCmsThrustPhi;
    vector<float> binningMult;
    vector<float> binningThrust;
    vector<float> binningEmiss;


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
