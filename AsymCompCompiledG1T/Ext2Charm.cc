#define KAON_LESS_CLASS
#define ONLY_DIST
//#define ONLY_MASS
#define D0Mass 1.865
//essentially resolution, just a guess
#define D0Width 0.6

//#define MAX_EVENTS 100000
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TStyle.h"
#include "TH2.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include <TROOT.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <time.h>

#include <set>
#include "TMVA/Factory.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "TROOT.h"
#include "TPluginManager.h"

#include "MvaClass.h"
#include "TwoHadAsymsCommons.h"


//#include "TMVA/TMVAGui.C"
// #include "TMVA/Method"

using namespace std;
char c_tmp[100];
vector<float>* binningM;
vector<float>* binningZ;
vector<float>* binningTheta;
#include "Ext2Charm.h"
#define OPENING_CUT 0.8
#define MAX_OPENING_CUT 2.0

void setZero(float* f1, float* f2, float* f3, float* f4, float* f5, float* f6,int num)
{
  for(int i=0;i<num;i++)
    {
      f1[i]=0;
      f2[i]=0;
      f3[i]=0;
      f4[i]=0;
      f5[i]=0;
      f6[i]=0;
    }
}
void setZero(float* f1, float* f2, float* f3, int num)
{
  for(int i=0;i<num;i++)
    {
      f1[i]=0;
      f2[i]=0;
      f3[i]=0;
     }
}
void setZero(float* f1, float* f2, float* f3, float* f4, float* f5, int num)
{
  for(int i=0;i<num;i++)
    {
      f1[i]=0;
      f2[i]=0;
      f3[i]=0;
      f4[i]=0;
      f5[i]=0;
     }
}
int main(int argc, char** argv)
{
  double errorSumUdsZ=0;
  double errorSumUdsM=0;
  double errorSumCharmZ=0;
  double errorSumCharmM=0;
  gROOT->SetStyle("Plain");
  int kinBin1=0;
  int kinBin2=0;
  int binning=0;
  int iCharm=0;
  int iOther=0;
  float meanM=0;
  float meanZ=0;
  float meanMFailCuts=0;
  float meanZFailCuts=0;
  float meanMPassCuts=0;
  float meanZPassCuts=0;

  vector<pair<float ,float> > * asymAllM;
  vector<pair<float ,float> > * asymCharmM;
  vector<pair<float ,float> > * asymUdsM;
  vector<pair<float ,float> > * asymAllZ;
  vector<pair<float ,float> > * asymCharmZ;
  vector<pair<float ,float> > * asymUdsZ;
  vector<pair<float ,float> > * asymAllMTheta;
  vector<pair<float ,float> > * asymCharmMTheta;
  vector<pair<float ,float> > * asymUdsMTheta;
  vector<pair<float ,float> > * asymAllZTheta;
  vector<pair<float ,float> > * asymCharmZTheta;
  vector<pair<float ,float> > * asymUdsZTheta;
  int maxKin=-1;
  vector<string> fieldNamesF;
  vector<string> fieldNamesI;
  vector<pair<float ,float> > * asymsNoCut;
  vector<pair<float ,float> > * asymsWithCut;

  binningM=new vector<float>[NumParticle];
  binningZ=new vector<float>[NumParticle];
  binningTheta=new vector<float>[NumParticle];
  loadBinning(binningM, binningZ);
  if(binningM[PiPi].size()>binningZ[PiPi].size())
      maxKin=binningM[PiPi].size();
  else
      maxKin=binningZ[PiPi].size();

  for(int i=0;i<maxKin;i++)
    {
      binningTheta[PiPi].push_back((i+1)*1.0/maxKin);
      cout <<" in b theta: " << (i+1)*1.0/maxKin <<endl;
    }  int numBins=0;
  int** passAsCharmZ=new int*[maxKin];
  int** passAsCharmM=new int*[maxKin];
  int** noPassAsCharmZ=new int*[maxKin];
  int** noPassAsCharmM=new int*[maxKin];
  int** passAsCharmZUds=new int*[maxKin];
  int** passAsCharmMUds=new int*[maxKin];
  int** noPassAsCharmZUds=new int*[maxKin];
  int** noPassAsCharmMUds=new int*[maxKin];

  int** passAsCharmZTheta=new int*[maxKin];
  int** passAsCharmMTheta=new int*[maxKin];
  int** noPassAsCharmZTheta=new int*[maxKin];
  int** noPassAsCharmMTheta=new int*[maxKin];
  int** passAsCharmZThetaUds=new int*[maxKin];
  int** passAsCharmMThetaUds=new int*[maxKin];
  int** noPassAsCharmZThetaUds=new int*[maxKin];
  int** noPassAsCharmMThetaUds=new int*[maxKin];



  //uds, charm
  float**** meanKin=allocateArray<float>(2,4,2,maxKin);
  float**** meanKinPassCuts=allocateArray<float>(2,4,2,maxKin);
  float**** meanKinFailCuts=allocateArray<float>(2,4,2,maxKin);


  float mvaCut=0.4; //bdtd is pos/neg

  if(argc==3)
    {
      char* c_mvaCut=argv[2];
      mvaCut=atof(c_mvaCut);
    }
  sprintf(c_tmp,"%0.2f",mvaCut);
  replDot(c_tmp);
#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
  TFile mcharmFile((string("charmContribOnlyDist")+c_tmp+".root").c_str());
  cout <<"loading: "<<"charmContribOnlyDist"<<c_tmp<<".root"<<endl;
#else
  TFile mcharmFile((string("charmContribKaonLess")+c_tmp+".root").c_str());
  cout <<"loading: "<<"charmContribKaonLess"<<c_tmp<<".root"<<endl;
#endif
#else
  TFile mcharmFile((string("charmContrib")+c_tmp+".root").c_str());
  cout <<"loading: "<<"charmContrib"<<c_tmp<<".root"<<endl;
#endif
#else
  TFile mcharmFile((string("charmContribOnlyMass")+c_tmp+".root").c_str());
  cout <<"loading: "<<"charmContribOnlyMass"<<c_tmp<<".root"<<endl;
#endif
  ofstream oFile("errorOut",ios_base::app);
  oFile << "running with cut of " << c_tmp<<endl;

  TTree* charmTree=(TTree*)mcharmFile.Get("charmTree");
  TTree* udsTree=(TTree*)mcharmFile.Get("udsTree");


  charmTree->SetBranchAddress("kinBin1",&kinBin1);
  charmTree->SetBranchAddress("kinBin2",&kinBin2);
  charmTree->SetBranchAddress("binning",&binning);
  charmTree->SetBranchAddress("numCorr",&iCharm);
  charmTree->SetBranchAddress("numFalse",&iOther);
  charmTree->SetBranchAddress("meanM",&meanM);
  charmTree->SetBranchAddress("meanZ",&meanZ);
  charmTree->SetBranchAddress("meanMFailCuts",&meanMFailCuts);
  charmTree->SetBranchAddress("meanZFailCuts",&meanZFailCuts);
  charmTree->SetBranchAddress("meanMPassCuts",&meanMPassCuts);
  charmTree->SetBranchAddress("meanZPassCuts",&meanZPassCuts);

  udsTree->SetBranchAddress("kinBin1",&kinBin1);
  udsTree->SetBranchAddress("kinBin2",&kinBin2);
  udsTree->SetBranchAddress("binning",&binning);
  udsTree->SetBranchAddress("numCorr",&iCharm);
  udsTree->SetBranchAddress("numFalse",&iOther);
  udsTree->SetBranchAddress("meanM",&meanM);
  udsTree->SetBranchAddress("meanZ",&meanZ);
  udsTree->SetBranchAddress("meanMFailCuts",&meanMFailCuts);
  udsTree->SetBranchAddress("meanZFailCuts",&meanZFailCuts);
  udsTree->SetBranchAddress("meanMPassCuts",&meanMPassCuts);
  udsTree->SetBranchAddress("meanZPassCuts",&meanZPassCuts);

  float* mLowCuts=new float[NumParticle];
  float* mHighCuts=new float[NumParticle];

  float* zLowCuts=new float[NumParticle];
  float* zHighCuts=new float[NumParticle];

  vector<float> binningAng;
  int numAngBins=16;
  double pi=3.14159265;
  for(int i=0;i<numAngBins;i++)
    {
      binningAng.push_back((i+1)*2*pi/((float)numAngBins)-pi);
      cout <<"ang binning: " << binningAng[i]<<endl;
    }

  mLowCuts[PiPi]=0.25;
  mLowCuts[PiK]=0.5;
  mLowCuts[KPi]=0.5;
  mLowCuts[KK]=0.5;
  mLowCuts[UNKNOWN]=0.25;

  mHighCuts[PiPi]=2.0;
  mHighCuts[PiK]=2.0;
  mHighCuts[KPi]=2.0;
  mHighCuts[KK]=1000.0;
  mHighCuts[UNKNOWN]=1000.0;

  zHighCuts[PiPi]=1.0;
  zHighCuts[PiK]=1.0;
  zHighCuts[KPi]=1.0;
  zHighCuts[KK]=1.0;
  zHighCuts[UNKNOWN]=1.0;

  char mvaMethod[]="MLP method";
  //  float mvaCut=0.9;

  cout <<"maxKin " <<maxKin<<endl;
  for(int i=0;i<maxKin;i++)
    {
      passAsCharmM[i]=new int[maxKin];
      noPassAsCharmM[i]=new int[maxKin];
      passAsCharmZ[i]=new int[maxKin];
      noPassAsCharmZ[i]=new int[maxKin];
      passAsCharmMUds[i]=new int[maxKin];
      noPassAsCharmMUds[i]=new int[maxKin];
      passAsCharmZUds[i]=new int[maxKin];
      noPassAsCharmZUds[i]=new int[maxKin];

      passAsCharmMTheta[i]=new int[maxKin];
      noPassAsCharmMTheta[i]=new int[maxKin];
      passAsCharmZTheta[i]=new int[maxKin];
      noPassAsCharmZTheta[i]=new int[maxKin];
      passAsCharmMThetaUds[i]=new int[maxKin];
      noPassAsCharmMThetaUds[i]=new int[maxKin];
      passAsCharmZThetaUds[i]=new int[maxKin];
      noPassAsCharmZThetaUds[i]=new int[maxKin];

      for(int j=0;j<maxKin;j++)
	{
	  passAsCharmM[i][j]=0;
	  passAsCharmZ[i][j]=0;
	  noPassAsCharmM[i][j]=0;
	  noPassAsCharmZ[i][j]=0;
	  passAsCharmMUds[i][j]=0;
	  passAsCharmZUds[i][j]=0;
	  noPassAsCharmMUds[i][j]=0;
	  noPassAsCharmZUds[i][j]=0;

	  passAsCharmMTheta[i][j]=0;
	  passAsCharmZTheta[i][j]=0;
	  noPassAsCharmMTheta[i][j]=0;
	  noPassAsCharmZTheta[i][j]=0;
	  passAsCharmMThetaUds[i][j]=0;
	  passAsCharmZThetaUds[i][j]=0;
	  noPassAsCharmMThetaUds[i][j]=0;
	  noPassAsCharmZThetaUds[i][j]=0;
	}
    }

  cout <<"branching at uds, we have "<<(Int_t)udsTree->GetEntries() << "entries"  <<endl;
  long iUdsOtherZ=0;
  long iUdsCharmZ=0;
  long iCharmOtherZ=0;
  long iCharmCharmZ=0;

  long iUdsOtherM=0;
  long iUdsCharmM=0;
  long iCharmOtherM=0;
  long iCharmCharmM=0;

  long iUdsOtherZTheta=0;
  long iUdsCharmZTheta=0;
  long iCharmOtherZTheta=0;
  long iCharmCharmZTheta=0;

  long iUdsOtherMTheta=0;
  long iUdsCharmMTheta=0;
  long iCharmOtherMTheta=0;
  long iCharmCharmMTheta=0;

  for(int i=0;i<(Int_t)udsTree->GetEntries();i++)
    {
      udsTree->GetEntry(i);
      meanKin[0][binning][zBinning][kinBin1]=meanZ;
      meanKin[0][binning][mBinning][kinBin1]=meanM;
      cout <<"uds loading 0, binning: " << binning <<" zbin, kinBin1: " << kinBin1 << " mean: " << meanZ <<" meanM: " << meanM << " meanzpass: " << meanZPassCuts <<endl;
      meanKinPassCuts[0][binning][zBinning][kinBin1]=meanZPassCuts;
      meanKinPassCuts[0][binning][mBinning][kinBin1]=meanMPassCuts;
      meanKinFailCuts[0][binning][zBinning][kinBin1]=meanZFailCuts;
      meanKinFailCuts[0][binning][mBinning][kinBin1]=meanMFailCuts;
      switch(binning)
	{
	case zBinning: 
	  cout <<"zbinning: " << iCharm <<" other: " <<iOther<<endl;
	  cout <<"kinbin1: " << kinBin1 <<" kinBin2: " << kinBin2 <<endl;
	  passAsCharmZUds[kinBin1][kinBin2]=iCharm;
	  noPassAsCharmZUds[kinBin1][kinBin2]=iOther;
	  iUdsOtherZ+=iOther;
	  iUdsCharmZ+=iCharm;
	  break;

	case mBinning:
	  cout <<"mbinning: " << iCharm <<" other: " <<iOther<<endl;
	  passAsCharmMUds[kinBin1][kinBin2]=iCharm;
	  noPassAsCharmMUds[kinBin1][kinBin2]=iOther;
	  iUdsOtherM+=iOther;
	  iUdsCharmM+=iCharm;
	  break;
	case multBinning:
	  cout <<"mult bin 1: "<<kinBin1 <<" bin2: "<< kinBin2  <<endl;
	  passAsCharmMThetaUds[kinBin1][kinBin2]=iCharm;
	  noPassAsCharmMThetaUds[kinBin1][kinBin2]=iOther;
	  iUdsOtherMTheta+=iOther;
	  iUdsCharmMTheta+=iCharm;
	  cout <<"mult bin2 " <<endl;
	  break;
	case mixedBinning:
	  cout <<"mixed bin 1: "<<kinBin1 <<" bin2: "<< kinBin2  <<endl;
	  passAsCharmZThetaUds[kinBin1][kinBin2]=iCharm;
	  noPassAsCharmZThetaUds[kinBin1][kinBin2]=iOther;
	  iUdsOtherZTheta+=iOther;
	  iUdsCharmZTheta+=iCharm;
	  cout <<"end mixed bin2 " <<endl;
	  break;
	default:
	  cout<<"wrong binnign 1 uds" << endl;
	}
    }
  cout <<"branching at charm, we have "<<(Int_t)charmTree->GetEntries() << "entries"  <<endl;
  for(int i=0;i<(Int_t)charmTree->GetEntries();i++)
    {
      charmTree->GetEntry(i);
      meanKin[1][binning][zBinning][kinBin1]=meanZ;
      meanKin[1][binning][mBinning][kinBin1]=meanM;
      cout <<"loading charm 0, binning: " << binning <<" zbin, kinBin1: " << kinBin1 << " mean: " << meanZ <<" meanM: " << meanM << endl;
      meanKinPassCuts[1][binning][zBinning][kinBin1]=meanZPassCuts;
      meanKinPassCuts[1][binning][mBinning][kinBin1]=meanMPassCuts;
      meanKinFailCuts[1][binning][zBinning][kinBin1]=meanZFailCuts;
      meanKinFailCuts[1][binning][mBinning][kinBin1]=meanMFailCuts;
      switch(binning)
	{
	case zBinning: 
	  cout <<"zbinning, charm: " << iCharm <<" other: " <<iOther<<endl;
	  passAsCharmZ[kinBin1][kinBin2]=iCharm;
	  noPassAsCharmZ[kinBin1][kinBin2]=iOther;
	  iCharmOtherZ+=iOther;
	  iCharmCharmZ+=iCharm;
	  break;

	case mBinning:
	  cout <<"mbinning, charm: " << iCharm <<" other: " <<iOther<<endl;
	  passAsCharmM[kinBin1][kinBin2]=iCharm;
	  noPassAsCharmM[kinBin1][kinBin2]=iOther;
	  iCharmOtherM+=iOther;
	  iCharmCharmM+=iCharm;
	  break;
	case multBinning:
	  passAsCharmMTheta[kinBin1][kinBin2]=iCharm;
	  noPassAsCharmMTheta[kinBin1][kinBin2]=iOther;
	  iCharmOtherMTheta+=iOther;
	  iCharmCharmMTheta+=iCharm;
	  break;
	case mixedBinning:
	  passAsCharmZTheta[kinBin1][kinBin2]=iCharm;
	  noPassAsCharmZTheta[kinBin1][kinBin2]=iOther;
	  iCharmOtherZTheta+=iOther;
	  iCharmCharmZTheta+=iCharm;
	  break;

	default:
	  cout<<"wrong binnign 1 charm" << endl;
	}
    }

  cout <<"mbinning charm: " << iCharmOtherM << " " << iCharmCharmM <<endl;
  cout <<"zbinning charm: " << iCharmOtherZ << " " << iCharmCharmZ <<endl;
  cout <<"mbinning uds: " << iUdsOtherM << " " << iUdsCharmM <<endl;
  cout <<"zbinning uds: " << iUdsOtherZ << " " << iUdsCharmZ <<endl;
    cout <<"mThetabinning charm: " << iCharmOtherMTheta << " " << iCharmCharmMTheta <<endl;
  cout <<"zThetabinning charm: " << iCharmOtherZTheta << " " << iCharmCharmZTheta <<endl;
  cout <<"mThetabinning uds: " << iUdsOtherMTheta << " " << iUdsCharmMTheta <<endl;
  cout <<"zThetabinning uds: " << iUdsOtherZTheta << " " << iUdsCharmZTheta <<endl;

  //let it run over uds and charm separately, but outputs in one root file, then use that for evaluation
  bool isCharm=false;
  int argc2=argc;
  char** argv2=new char*[argc];
  for(int i=0;i<argc;i++)
    {
      argv2[i]=new char[200];
      strcpy(argv2[i],argv[i]);
    }

  char* rootPath=argv2[1];
  srand(time(NULL));
  cout <<"Root path is: " << rootPath <<endl;
  string sRootPath(rootPath);
  if(sRootPath.find_last_of('/')==sRootPath.length()-1)
    {
      sRootPath.erase(sRootPath.find_last_of('/'));
    }
  size_t found=sRootPath.find_last_of('/');
  string folderName=sRootPath.substr(found+1);
  bool onResonance=false;

  int numPos=folderName.find("ex");
  int expNumber=-1;
  if(numPos!=string::npos)
    {
      char tmpNum[3];
      tmpNum[0]=folderName[numPos+2];
      cout <<"first num: " << tmpNum[0] << ", pos: " << numPos <<endl;
      tmpNum[1]=folderName[numPos+3];
      cout <<"sec num: " << tmpNum[1] << ", pos: " << numPos <<endl;
      tmpNum[2]='\n';
      expNumber=atoi(tmpNum);
    }
  cout<<"looking at exp " << expNumber <<" onres?: " << onResonance <<endl;
  cout <<"hallo, Welt" <<endl;
  vector<float*> memLocsF;
  vector<int*> memLocsI;

  int z1Counter;
  float z1[2000];
  int z2Counter;
  float z2[2000];
  int mass1Counter;
  float mass1[2000];
  int mass2Counter;
  float mass2[2000];
  int phiRCounter;
  float phiRSum[2000];
  float thrustTheta;
  int chargeCounter;
  int chargeType[2000];
  int charge1Counter;
  int chargeType1[2000];
  int particle1Counter;
  int particleType1[2000];
  int charge2Counter;
  int chargeType2[2000];
  int particle2Counter;
  int particleType2[2000];

  float thrustProj11[2000];
  int thrustProj11Counter;
  float thrustProj12[2000];
  int thrustProj12Counter;
  float thrustProj21[2000];
  int thrustProj21Counter;
  float thrustProj22[2000];
  int thrustProj22Counter;

  float thetaEThrust=pi/(float)2; //for now until the new data is generated //correction factor is 1/4 sin^2(theta)
  int z1Loc=0;
  int z2Loc=1;
  int mass1Loc=2;
  int mass2Loc=3;
  int phiRSumLoc=4;

  float minMassDiffToD=0;
  float minMassDiffToD_PN=0;
  float maxZ_PN=0;
  float maxZ=0;
  float numKaons=0;
  float meanDist=0;
  float meanDistK=0;
  float maxCosThetaK=0;
  float maxCosTheta_PN=0;
  float validMva=0;

  long nevents=0;
  //binned in the costheta/sintheta
  int****** countsTheta=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  int****** countsUdsTheta=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  int****** countsCharmTheta=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);

  int****** counts=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  int****** countsUds=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  int****** countsCharm=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);

  float***** xVals=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** xValsNum=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** xValsTheta=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** xValsNumTheta=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);

  float***** xValsUds=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** xValsNumUds=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** xValsUdsTheta=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** xValsNumUdsTheta=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);

  float***** valsOtherKin=allocateArray<float>(2,NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin);
  int***** valsOtherKinNum=allocateArray<int>(2,NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin);
  //  float***** xValsOtherKin=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  //  int***** xValsNumOtherKin=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** yVals=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** yValsNum=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** yValsUds=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** yValsNumUds=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);

  float***** xValsCharm=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** xValsNumCharm=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** yValsCharm=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** yValsNumCharm=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** yValsTheta=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** yValsNumTheta=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** yValsUdsTheta=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** yValsNumUdsTheta=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** xValsCharmTheta=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** xValsNumCharmTheta=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** yValsCharmTheta=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** yValsNumCharmTheta=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);

  int nanEv=0;
  int cuttedEv=0;
  int accEv=0;
  const float zCut=1.0;//is also the max for the z binning. Rought estimate, should follow from resolution of z
  const float mCut=1000;
  const float thetaCut=1000;


  fieldNamesI.push_back("z1Counter");
  fieldNamesF.push_back("z1");

  memLocsI.push_back(&z1Counter);
  memLocsF.push_back(z1);

  fieldNamesI.push_back("z2Counter");
  fieldNamesF.push_back("z2");

  memLocsI.push_back(&z2Counter);
  memLocsF.push_back(z2);

  fieldNamesI.push_back("mass1Counter");
  fieldNamesF.push_back("mass1");
  memLocsI.push_back(&mass1Counter);
  memLocsF.push_back(mass1);

  fieldNamesI.push_back("mass2Counter");
  fieldNamesF.push_back("mass2");

  memLocsI.push_back(&mass2Counter);
  memLocsF.push_back(mass2);

  fieldNamesI.push_back("phiRSumCounter");
  fieldNamesF.push_back("phiRSum");
  memLocsI.push_back(&phiRCounter);
  memLocsF.push_back(phiRSum);

  fieldNamesI.push_back("chargeTypeCounter");
  fieldNamesI.push_back("chargeType");
  memLocsI.push_back(&chargeCounter);
  memLocsI.push_back(chargeType);

  fieldNamesI.push_back("chargeType1Counter");
  fieldNamesI.push_back("chargeType1");
  memLocsI.push_back(&charge1Counter);
  memLocsI.push_back(chargeType1);

  fieldNamesI.push_back("particleType1Counter");
  fieldNamesI.push_back("particleType1");
  memLocsI.push_back(&particle1Counter);
  memLocsI.push_back(particleType1);

  fieldNamesI.push_back("chargeType2Counter");
  fieldNamesI.push_back("chargeType2");
  memLocsI.push_back(&charge2Counter);
  memLocsI.push_back(chargeType2);

  fieldNamesI.push_back("particleType2Counter");
  fieldNamesI.push_back("particleType2");
  memLocsI.push_back(&particle2Counter);
  memLocsI.push_back(particleType2);

  fieldNamesF.push_back("thrustTheta");
  memLocsF.push_back(&thrustTheta);

  fieldNamesF.push_back("thetaEThrust");
  memLocsF.push_back(&thetaEThrust);
  fieldNamesF.push_back("thrustProj11");
  fieldNamesI.push_back("thrustProj11Counter");
  memLocsF.push_back(thrustProj11);
  memLocsI.push_back(&thrustProj11Counter);
  fieldNamesF.push_back("thrustProj12");
  fieldNamesI.push_back("thrustProj12Counter");
  memLocsF.push_back(thrustProj12);
  memLocsI.push_back(&thrustProj12Counter);
  fieldNamesF.push_back("thrustProj21");
  fieldNamesI.push_back("thrustProj21Counter");
  memLocsF.push_back(thrustProj21);
  memLocsI.push_back(&thrustProj21Counter);
  fieldNamesF.push_back("thrustProj22");
  fieldNamesI.push_back("thrustProj22Counter");
  memLocsF.push_back(thrustProj22);
  memLocsI.push_back(&thrustProj22Counter);


  /*
  pTreeSaver->addFieldF("maxIntraVertDistK");
  pTreeSaver->addFieldF("maxExtraVertDistK");
  pTreeSaver->addFieldF("maxIntraVertDistPN");
  pTreeSaver->addFieldF("maxExtraVertDistPN");
  pTreeSaver->addFieldF("meanDistPN");
  pTreeSaver->addFieldF("kaonPionRatio");
  pTreeSaver->addFieldF("dr");
  pTreeSaver->addFieldF("minMassDiffToD");
  pTreeSaver->addFieldF("minMassDiffToD_PN");
  pTreeSaver->addFieldF("maxZ_PN");
  pTreeSaver->addFieldF("maxZ");
  pTreeSaver->addFieldF("numKaons");
  pTreeSaver->addFieldF("meanDist");
  pTreeSaver->addFieldF("meanDistK");
  pTreeSaver->addFieldF("maxCosThetaK");
  pTreeSaver->addFieldF("maxCosTheta_PN");
  pTreeSaver->addFieldF("validMva");
  */
  float maxIntraVertDist=0;
  float maxExtraVertDist=0;
  float maxIntraVertDistPN=0;
  float maxExtraVertDistPN=0;
  float meanDistPN=0;
  float kaonPionRatio=0;
  float dr=0;
  float maxExtraVertDistK=0;
  float maxIntraVertDistK=0;
  int numCharm=0;
  int numUds=0;
  fieldNamesF.push_back("maxIntraVertDistK");
  memLocsF.push_back(&maxIntraVertDistK);
  fieldNamesF.push_back("maxExtraVertDistK");
  memLocsF.push_back(&maxExtraVertDistK);
  fieldNamesF.push_back("maxIntraVertDistPN");
  memLocsF.push_back(&maxIntraVertDistPN);
  fieldNamesF.push_back("maxExtraVertDistPN");
  memLocsF.push_back(&maxExtraVertDistPN);
  fieldNamesF.push_back("meanDistPN");
  memLocsF.push_back(&meanDistPN);
  fieldNamesF.push_back("kaonPionRatio");
  memLocsF.push_back(&kaonPionRatio); //don't forget to transform the variable before the input
  fieldNamesF.push_back("dr");
  memLocsF.push_back(&dr);
  fieldNamesF.push_back("minMassDiffToD");
  memLocsF.push_back(&minMassDiffToD);
  fieldNamesF.push_back("minMassDiffToD_PN");
  memLocsF.push_back(&minMassDiffToD_PN);
  fieldNamesF.push_back("maxZ_PN");
  memLocsF.push_back(&maxZ_PN);
  fieldNamesF.push_back("maxZ");
  memLocsF.push_back(&maxZ);
  fieldNamesF.push_back("numKaons");
  memLocsF.push_back(&numKaons);
  fieldNamesF.push_back("meanDist");
  memLocsF.push_back(&meanDist);
  fieldNamesF.push_back("meanDistK");
  memLocsF.push_back(&meanDistK);
  fieldNamesF.push_back("maxCosThetaK");
  memLocsF.push_back(&maxCosThetaK);
  fieldNamesF.push_back("maxCosTheta_PN");
  memLocsF.push_back(&maxCosTheta_PN);
  fieldNamesF.push_back("validMva");
  memLocsF.push_back(&validMva);

  /* factory->AddVariable( "shiftDr := log(abs(dr))", "expressDr", "mm", 'F' );
   factory->AddVariable( "maxZ", "maxZ", "mm", 'F' );
   factory->AddVariable( "maxZ_PN", "maxZ_PN", "mm", 'F' );
   factory->AddVariable( "thrust", "thrust", "", 'F' );
   factory->AddVariable( "LogmeanDist := log(meanDist)", "LogmeanDist", "mm", 'F' );
   factory->AddVariable( "LogmeanDistPN := log(meanDistPN)", "LogmeanDistPN", "mm", 'F' );
   factory->AddVariable( "LogmeanDistK := log(meanDistK)", "LogmeanDistK", "mm", 'F' );
   factory->AddVariable( "LogmaxExtraVertDistK := log(maxExtraVertDistK)", "LogmaxExtraVertDistK", "mm", 'F' );
   factory->AddVariable( "LogmaxIntraVertDistK := log(maxIntraVertDistK)", "LogmaxIntraVertDistK", "mm", 'F' );
   factory->AddVariable( "LogmaxExtraVertDistPN := log(maxExtraVertDistPN)", "LogmaxExtraVertDistPN", "mm", 'F' );
   factory->AddVariable( "LogmaxIntraVertDistPN := log(maxIntraVertDistPN)", "LogmaxIntraVertDistPN", "mm", 'F' );
   factory->AddVariable( "maxCosThetaK", "maxCosThetaK", "", 'F' );
   factory->AddVariable( "maxCosTheta_PN", "maxCosTheta_PN", "", 'F' );
   factory->AddVariable("visEnergy","visEnergy","GeV",'F');
   factory->AddVariable("visEnergyOnFile","visEnergyOnFile","GeV",'F');
   factory->AddVariable("maxKt","maxKt","GeV",'F');
   //see if we are mainly sensitive to inv. mass
   factory->AddVariable( "LogminMassDiffToD := log(minMassDiffToD)", "LogminMassDiffToD", "mm", 'F' );
   factory->AddVariable( "LogminMassDiffToD_PN := log(minMassDiffToD_PN)", "LogminMassDiffToD_PN", "mm", 'F' );
*/

  //respect order!!!
  vector<string> varNamesF;
  vector<string> varNamesI;
  varNamesF.push_back("shiftDr := log(abs(dr))");
#ifndef KAON_LESS_CLASS
  varNamesF.push_back("maxZ");
#endif
#ifndef ONLY_DIST
  varNamesF.push_back("maxZ_PN");
  varNamesF.push_back("thrust");
#endif
  varNamesF.push_back("LogmeanDist := log(meanDist)");
#ifndef KAON_LESS_CLASS
  varNamesF.push_back("LogmeanDistPN := log(meanDistPN)");
  varNamesF.push_back("LogmeanDistK := log(meanDistK)");
  varNamesF.push_back("LogmaxExtraVertDistK := log(maxExtraVertDistK)");
#endif
  //  varNamesF.push_back("LogmaxIntraVertDistK := log(maxIntraVertDistK)");
  varNamesF.push_back("LogmaxExtraVertDistPN := log(maxExtraVertDistPN)");
  varNamesF.push_back("LogmaxIntraVertDistPN := log(maxIntraVertDistPN)");
#ifndef ONLY_DIST
#ifndef KAON_LESS_CLASS
  varNamesF.push_back("maxCosThetaK");
#endif
  varNamesF.push_back("maxCosTheta_PN");
  varNamesF.push_back("visEnergy");
#ifndef KAON_LESS_CLASS
  varNamesF.push_back("visEnergyOnFile");
#endif
  varNamesF.push_back("maxKt");
#ifndef KAON_LESS_CLASS
  varNamesF.push_back("LogminMassDiffToD := log(minMassDiffToD)");
#endif
  varNamesF.push_back("LogminMassDiffToD_PN := log(minMassDiffToD_PN)");
#endif
  MvaClass myMva;

  myMva.initialize(varNamesF,varNamesI);
  cout <<"done init" <<endl;
  //in principle I also have the maxKt and misEnergy because I know that the counter on the mva tree only advances for validMva=true
  TChain* chAll=new TChain("DataTree");
  TChain* chMvaTree=new TChain("mvaTree"); //has all the vars, but we only need the once in addition that we don't have in the other one
  chAll->Add((string(rootPath)+"/singlerootfiles/*.mr.root").c_str());
  chMvaTree->Add((string(rootPath)+"/singlerootfiles/*.mr.root").c_str());

  float visEnergy;
  float visEnergyOnFile;
  float maxKt;
  float thrust;
  float mva_minMassDiffToD;
  chMvaTree->SetBranchAddress("visEnergy",&visEnergy);
  chMvaTree->SetBranchAddress("visEnergyOnFile",&visEnergyOnFile);
  chMvaTree->SetBranchAddress("maxKt",&maxKt);
  chMvaTree->SetBranchAddress("thrust",&thrust);
  chMvaTree->SetBranchAddress("minMassDiffToD", &mva_minMassDiffToD); //to xcheck

  for(int i=0;i<fieldNamesF.size();i++)
    {
      cout <<"mva trying to branch on " << fieldNamesF[i]<<endl;
      float* memLoc=memLocsF[i];
      chAll->SetBranchAddress(fieldNamesF[i].c_str(),memLoc);
    }
  for(int i=0;i<fieldNamesI.size();i++)
    {
      int* memLoc2=memLocsI[i];
      chAll->SetBranchAddress(fieldNamesI[i].c_str(),memLoc2);
    }
  cout <<"done branching "<<endl;
  nevents=chAll->GetEntries();
  cout <<"we have " << nevents <<endl;
  //   factory->AddVariable( "shiftDr := log(abs(dr))", "expressDr", "mm", 'F' );
  /*   factory->AddVariable( "LogmeanDist := log(meanDist)", "LogmeanDist", "mm", 'F' );
       factory->AddVariable( "LogmeanDistPN := log(meanDistPN)", "LogmeanDistPN", "mm", 'F' );
       factory->AddVariable( "LogmeanDistK := log(meanDistK)", "LogmeanDistK", "mm", 'F' );
       factory->AddVariable( "LogmaxExtraVertDistK := log(maxExtraVertDistK)", "LogmaxExtraVertDistK", "mm", 'F' );
       factory->AddVariable( "LogmaxIntraVertDistK := log(maxIntraVertDistK)", "LogmaxIntraVertDistK", "mm", 'F' );
       factory->AddVariable( "LogmaxExtraVertDistPN := log(maxExtraVertDistPN)", "LogmaxExtraVertDistPN", "mm", 'F' );
       factory->AddVariable( "LogmaxIntraVertDistPN := log(maxIntraVertDistPN)", "LogmaxIntraVertDistPN", "mm", 'F' );
       factory->AddVariable( "LogminMassDiffToD := log(minMassDiffToD)", "LogminMassDiffToD", "mm", 'F' );
       factory->AddVariable( "LogminMassDiffToD_PN := log(minMassDiffToD_PN)", "LogminMassDiffToD_PN", "mm", 'F' );*/
  int mvaCounter=0;
  int nonValidMvaCounter=0;
  long allHads=0;
  long allHadsPiPiPN=0;
  long allHadsThrust=0;
  long allHadsOpening=0;
  long allHadsHighCuts=0;
  long allHadsBeforeFill=0;
  long allHadsBeforeFillM=0;
  for(long i=0;i<nevents;i++)
    {
      if(!(i%100000))
	cout <<"processing event nr " << i << " of " << nevents << "(" << 100*i/(float)nevents<< "% )"<<endl;
#ifdef MAX_EVENTS
      if(i>MAX_EVENTS)
	break;
#endif

      chAll->GetEntry(i);
      if(validMva==1)
	{
	  chMvaTree->GetEntry(mvaCounter);
	  //	  cout <<" comp: " << mva_minMassDiffToD <<" or : " << minMassDiffToD <<endl; //worked
	  mvaCounter++;
	}
      else
	{
	  nonValidMvaCounter++;
#ifndef KAON_LESS_CLASS
	  continue;
#endif
	}
      //      cout <<"valid mva" <<endl;
      float shiftDr=log(fabs(dr));
      float LogmeanDist=log(meanDist);
      float LogmeanDistPN=log(meanDistPN);
      float LogmeanDistK=log(meanDistK);
      float LogmaxExtraVertDistK=log(maxExtraVertDistK);
      //      cout <<"maxIntraVertDistK: " << maxIntraVertDistK <<endl;

      float LogmaxIntraVertDistK=log(maxIntraVertDistK);
#ifndef KAON_LESS_CLASS
      if(maxIntraVertDistK==0)
	{
	  LogmaxIntraVertDistK=0;
	  //	  cout <<"cont wg. intravertdistk=0" <<endl;
	  continue;//be safe
	}
#endif
      float LogmaxExtraVertDistPN=log(maxExtraVertDistPN);
      float LogmaxIntraVertDistPN=log(maxIntraVertDistPN);
      float LogminMassDiffToD=log(minMassDiffToD);
      float LogminMassDiffToD_PN=log(minMassDiffToD_PN);
      vector<float> valsF;
      valsF.push_back(shiftDr);
#ifndef KAON_LESS_CLASS
      valsF.push_back(maxZ);
#endif
#ifndef ONLY_DIST
      valsF.push_back(maxZ_PN);
      valsF.push_back(thrust);
#endif
      valsF.push_back(LogmeanDist);
#ifndef KAON_LESS_CLASS
      valsF.push_back(LogmeanDistPN);
      valsF.push_back(LogmeanDistK);
      valsF.push_back(LogmaxExtraVertDistK);
#endif
      //      valsF.push_back(LogmaxIntraVertDistK);
      valsF.push_back(LogmaxExtraVertDistPN);
      valsF.push_back(LogmaxIntraVertDistPN);
#ifndef ONLY_DIST
#ifndef KAON_LESS_CLASS
      valsF.push_back(maxCosThetaK);
#endif
      valsF.push_back(maxCosTheta_PN);
      valsF.push_back(visEnergy);
#ifndef KAON_LESS_CLASS
      valsF.push_back(visEnergyOnFile);
#endif
      valsF.push_back(maxKt);
#ifndef KAON_LESS_CLASS
      valsF.push_back(LogminMassDiffToD);
#endif
      valsF.push_back(LogminMassDiffToD_PN);
#endif
      //pass as charm?
      vector<int> valsI;
      for(int i=0;i<valsF.size();i++)
	{
	  //	      cout <<"i " << i << " " << valsF[i]<<endl;
	}
      //	  cout <<"mva ouput: " << myMva.evaluate(valsF,valsI,mvaMethod)<<endl;
      //	      cout <<"zbin1: " << zbin1 <<", zbin2: " << zbin2 <<" maxKin: " << maxKin <<endl;
      //	      cout <<"mbin1: " << zbin1 <<", mbin2: " << zbin2 <<" maxKin: " << maxKin <<endl;

#ifdef ONLY_MASS
	  if(fabs(minMassDiffToD)< mvaCut*D0Width)
#else
      if(myMva.evaluate(valsF,valsI,mvaMethod)>mvaCut)	
#endif
	{
	  isCharm=true;
	  numCharm++;
	}
      else
	{
	  numUds++;
	  isCharm=false;
	}
#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
      if(meanDist==0 || dr==0 || maxExtraVertDistPN==0 || maxIntraVertDistPN==0)
#else
      if(meanDist==0 || meanDistPN ==0 || dr==0 || maxZ_PN==0 || thrust ==0 || maxExtraVertDistPN==0 || maxIntraVertDistPN==0 || maxCosTheta_PN==0 || visEnergy==0 ||  maxKt==0 ||  minMassDiffToD_PN==0)
#endif
	{
	  //	  cout << " meanDist: " << meanDist << " meanDistPN " << meanDistPN <<" dr: " << dr <<" maxZ_PN: " << maxZ_PN << " thrust " << thrust <<" maxExtra: "<< maxExtraVertDistPN << " maxIntraVertDist: " << maxIntraVertDistPN <<" maxcosttheta: " << maxCosTheta_PN << " visenergy: " << visEnergy <<  " maxKT: " << maxKt << " minMassDiff " << minMassDiffToD_PN <<endl;
	  //	  cout <<"cont due to some var ==0" <<endl;
	  continue;
	}
#endif
#endif

      for(int iHadsInEv=0;iHadsInEv<z1Counter;iHadsInEv++) 
	{
	  allHads++;
	  if(!((particleType1[iHadsInEv]==PiPi) && (particleType2[iHadsInEv]==PiPi)))
	    continue;
	  if(!((chargeType1[iHadsInEv]==PN) && (chargeType2[iHadsInEv]==PN)))
	    continue;
	  allHadsPiPiPN++;
	  if(fabs(thrustProj11[iHadsInEv])>MAX_OPENING_CUT || fabs(thrustProj12[iHadsInEv])>MAX_OPENING_CUT ||fabs(thrustProj21[iHadsInEv])>MAX_OPENING_CUT ||fabs(thrustProj22[iHadsInEv])>MAX_OPENING_CUT)
	    continue;
	  allHadsThrust++;
	  if(fabs(thrustProj11[iHadsInEv])<OPENING_CUT || fabs(thrustProj12[iHadsInEv])<OPENING_CUT ||fabs(thrustProj21[iHadsInEv])<OPENING_CUT ||fabs(thrustProj22[iHadsInEv])<OPENING_CUT)
	    {
	      continue;
	    }
	  allHadsOpening++;
	  int chargeIndx=ind(chargeType1[iHadsInEv],chargeType2[iHadsInEv],NumCharge);
	  int mIndex=ind(particleType1[iHadsInEv],particleType2[iHadsInEv],NumParticle);	  
	  //	  cout <<"mindex; "<< mIndex <<" pt1: " << particleType1[iHadsInEv] << " pt2: " << particleType2[iHadsInEv] << " num: " << NumParticle <<endl;

	  float zValue1=((float*)(memLocsF[z1Loc]))[iHadsInEv];
	  float mValue1=((float*)(memLocsF[mass1Loc]))[iHadsInEv];
	  float zValue2=((float*)(memLocsF[z2Loc]))[iHadsInEv];
	  float mValue2=((float*)(memLocsF[mass2Loc]))[iHadsInEv];


	  if(mValue1>mHighCuts[particleType1[iHadsInEv]] || mValue2>mHighCuts[particleType1[iHadsInEv]])
	    continue;
	  if(zValue1 >1 || zValue1 <0 || zValue2>1 || zValue2 <0)
	    continue;
	  allHadsHighCuts++;
	  int zbin1=getBin(binningZ[particleType1[iHadsInEv]],zValue1);
	  int zbin2=getBin(binningZ[particleType1[iHadsInEv]],zValue2);
	  int mbin1=getBin(binningM[particleType1[iHadsInEv]],mValue1);
	  int mbin2=getBin(binningM[particleType1[iHadsInEv]],mValue2);
	  //	  if(zbin1<0 || mbin1 <0)
	  //	    cout <<"z: " << zValue1 <<" m: " << mValue1 <<endl;

	  int m_chType1=chargeType1[iHadsInEv];
	  int m_pt1=particleType1[iHadsInEv];

	  int m_chType2=chargeType2[iHadsInEv];
	  int m_pt2=particleType2[iHadsInEv];
	  float value1=((float*)(memLocsF[z1Loc]))[iHadsInEv];

	  float thetaValue=sin(thetaEThrust)*sin(thetaEThrust)/(1+cos(thetaEThrust)*cos(thetaEThrust));
	  //	  cut <<"thetaValue: " << thetaValue<<endl;

	  if(isnan(value1))
	    {
	      cout <<"value 1 nan" <<endl;
	      continue;
	    }
	  //		    cout <<"1+1 " <<endl;
	  float value2=((float*)(memLocsF[z2Loc]))[iHadsInEv];
	  if(isnan(value2))
	    {
	      cout <<"value 2 nan" <<endl;
	      continue;
	    }
	  float value3=((float*)(memLocsF[phiRSumLoc]))[iHadsInEv];

	  if(isnan(value3))
	    {
	      cout <<"value 3 nan" <<endl;
	      continue;
	    }
	  allHadsBeforeFill++;
	  for(int iBin=zBinning;iBin<=mBinning;iBin++)
	    {
	      switch(iBin)
		{
		case zBinning:
		  {
		    if(value2>zHighCuts[particleType1[iHadsInEv]] || value1 > zHighCuts[particleType1[iHadsInEv]])
		      {
			cuttedEv++;
			continue;
		      }
		    //		    cout <<"fill: " << value1 <<", val2: " << value2 <<" val3: " << value3 <<endl;
		    fillBorders(binningZ[particleType1[iHadsInEv]],binningAng,value1,value2,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);
		    fillBordersMix(binningZ[particleType1[iHadsInEv]],binningTheta[particleType1[iHadsInEv]],binningAng,value1,thetaValue,value3,countsTheta[iBin][chargeIndx][mIndex],xValsTheta[iBin][chargeIndx][mIndex],xValsNumTheta[iBin][chargeIndx][mIndex],yValsTheta[iBin][chargeIndx][mIndex],yValsNumTheta[iBin][chargeIndx][mIndex]);
		    if(isCharm)
		      {
			fillBorders(binningZ[particleType1[iHadsInEv]],binningAng,value1,value2,value3,countsCharm[iBin][chargeIndx][mIndex],xValsCharm[iBin][chargeIndx][mIndex],xValsNumCharm[iBin][chargeIndx][mIndex],yValsCharm[iBin][chargeIndx][mIndex],yValsNumCharm[iBin][chargeIndx][mIndex]);
		    fillBordersMix(binningZ[particleType1[iHadsInEv]],binningTheta[particleType1[iHadsInEv]],binningAng,value1,thetaValue,value3,countsCharmTheta[iBin][chargeIndx][mIndex],xValsCharmTheta[iBin][chargeIndx][mIndex],xValsNumCharmTheta[iBin][chargeIndx][mIndex],yValsCharmTheta[iBin][chargeIndx][mIndex],yValsNumCharmTheta[iBin][chargeIndx][mIndex]);
		      }
		    if(!isCharm)
		      {
			fillBorders(binningZ[particleType1[iHadsInEv]],binningAng,value1,value2,value3,countsUds[iBin][chargeIndx][mIndex],xValsUds[iBin][chargeIndx][mIndex],xValsNumUds[iBin][chargeIndx][mIndex],yValsUds[iBin][chargeIndx][mIndex],yValsNumUds[iBin][chargeIndx][mIndex]);
		    fillBordersMix(binningZ[particleType1[iHadsInEv]],binningTheta[particleType1[iHadsInEv]],binningAng,value1,thetaValue,value3,countsUdsTheta[iBin][chargeIndx][mIndex],xValsUdsTheta[iBin][chargeIndx][mIndex],xValsNumUdsTheta[iBin][chargeIndx][mIndex],yValsUdsTheta[iBin][chargeIndx][mIndex],yValsNumUdsTheta[iBin][chargeIndx][mIndex]);
		      }
		    break;
		  }
		case mBinning:
		  {
		    float value1=((float*)(memLocsF[mass1Loc]))[iHadsInEv];
		    //		    cout <<"val1: " << value1 <<endl;
		    if(isnan(value1))
		      {
			cout <<"value 1 m nan" <<endl;
			nanEv++;
			continue;
		      }
		    float value2=((float*)(memLocsF[mass2Loc]))[iHadsInEv];

		    if(isnan(value2))
		      {
			cout <<"value 2 m nan" <<endl;
			nanEv++;
			continue;
		      }

		    if(value2>mCut || value1 > mCut)
		      {
			cout <<"m1 value: " << value1 << " m2 value: " << value2 <<endl;
			cuttedEv++;
			continue;
		      }
		    float value3=((float*)(memLocsF[phiRSumLoc]))[iHadsInEv];

		    if(isnan(value3))
		      {
			cout <<"value 3 m nan" <<endl;
			nanEv++;
			continue;
		      }
		    if(value1<mLowCuts[particleType1[iHadsInEv]] || value2<mLowCuts[particleType1[iHadsInEv]])
		      {
			cuttedEv++;
			continue;
		      }
		    if(value1>mHighCuts[particleType1[iHadsInEv]] || value2>mHighCuts[particleType1[iHadsInEv]])
		      {
			cuttedEv++;
			continue;
		      }
		    allHadsBeforeFillM++;
		    accEv++;
		    //		    cout <<"mindex: "<< mIndex <<" charg: " << chargeIndx <<endl;
		    fillBorders(binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);

		    fillBordersMix(binningM[particleType1[iHadsInEv]],binningTheta[particleType1[iHadsInEv]],binningAng,value1,thetaValue,value3,countsTheta[iBin][chargeIndx][mIndex],xValsTheta[iBin][chargeIndx][mIndex],xValsNumTheta[iBin][chargeIndx][mIndex],yValsTheta[iBin][chargeIndx][mIndex],yValsNumTheta[iBin][chargeIndx][mIndex]);
		    if(isCharm)
		      {
			fillBorders(binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,countsCharm[iBin][chargeIndx][mIndex],xValsCharm[iBin][chargeIndx][mIndex],xValsNumCharm[iBin][chargeIndx][mIndex],yValsCharm[iBin][chargeIndx][mIndex],yValsNumCharm[iBin][chargeIndx][mIndex]);
			fillBordersMix(binningM[particleType1[iHadsInEv]],binningTheta[particleType1[iHadsInEv]],binningAng,value1,thetaValue,value3,countsCharmTheta[iBin][chargeIndx][mIndex],xValsCharmTheta[iBin][chargeIndx][mIndex],xValsNumCharmTheta[iBin][chargeIndx][mIndex],yValsCharmTheta[iBin][chargeIndx][mIndex],yValsNumCharmTheta[iBin][chargeIndx][mIndex]);
		      }
		    if(!isCharm)
		      {
			fillBorders(binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,countsUds[iBin][chargeIndx][mIndex],xValsUds[iBin][chargeIndx][mIndex],xValsNumUds[iBin][chargeIndx][mIndex],yValsUds[iBin][chargeIndx][mIndex],yValsNumUds[iBin][chargeIndx][mIndex]);
			fillBordersMix(binningM[particleType1[iHadsInEv]],binningTheta[particleType1[iHadsInEv]],binningAng,value1,thetaValue,value3,countsUdsTheta[iBin][chargeIndx][mIndex],xValsUdsTheta[iBin][chargeIndx][mIndex],xValsNumUdsTheta[iBin][chargeIndx][mIndex],yValsUdsTheta[iBin][chargeIndx][mIndex],yValsNumUdsTheta[iBin][chargeIndx][mIndex]);
		      }
		    break;
		  }
		default:
		  cout << "binning not recognized1" << endl;
		  exit(0);
		}
	    }
	}
    }
  ofstream* errorVglFile;
  ofstream outputFile("outputFile");
  ofstream xcheckFile("xcheck");
  errorVglFile=new ofstream("errorVgl");

  for(int iBin=zBinning;iBin<=mBinning;iBin++)
    {
      //      for(int iCh=PN;iCh<=ZZ;iCh++)
      cout <<"iBin: " << iBin <<endl;
      int iCh1=PN;
      int iCh2=PN;
      int iPa1=PiPi;
      int iPa2=PiPi;
      vector<float>* t_chi2=new vector<float>;
      vector<int>* t_ndfs=new vector<int>;
      vector<float>* t_chi2Z=new vector<float>;
      vector<int>* t_ndfsZ=new vector<int>;
      vector<float>* t_chi2Theta=new vector<float>;
      vector<int>* t_ndfsTheta=new vector<int>;
      vector<vector<pair<float,float> >* > zAsyms;
      vector<vector<pair<float,float> >* >mAsyms;
      vector<vector<pair<float,float> >* > zAsymsTheta;
      vector<vector<pair<float,float> >* >mAsymsTheta;

      vector<vector<float>* > v_chi2Z;
      vector<vector<int>* > v_ndfsZ;
      vector<vector<float>* > v_chi2M;
      vector<vector<int>* > v_ndfsM;
      vector<vector<float>* > v_chi2ZTheta;
      vector<vector<int>* > v_ndfsZTheta;
      vector<vector<float>* > v_chi2MTheta;
      vector<vector<int>* > v_ndfsMTheta;
      switch(iBin)
	{
	case zBinning:
	  {
	    cout <<"zb: " << endl;
	    numBins=binningZ[iPa1].size();
	    cout <<"zb1: " << endl;
	    asymAllZ=fitTheSh__(counts[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2,t_ndfs);
	    cout <<"zb2: " << endl;
	    asymCharmZ=fitTheSh__(countsCharm[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2,t_ndfs);
	    asymUdsZ=fitTheSh__(countsUds[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2,t_ndfs);
	    cout <<"theta asym Z" <<endl;
	    asymAllZTheta=fitTheSh__(countsTheta[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2Theta,t_ndfsTheta);
	    asymCharmZTheta=fitTheSh__(countsCharmTheta[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2Theta,t_ndfsTheta);
	    asymUdsZTheta=fitTheSh__(countsUdsTheta[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2Theta,t_ndfsTheta);
	    v_chi2Z.push_back(t_chi2);
	    v_ndfsZ.push_back(t_ndfs);
	    break;
	  }
	case mBinning:
	  {
	    cout <<"mb: " << endl;
	    numBins=binningM[iPa1].size();
	    cout <<"mb1: " << endl;
	    asymAllM=fitTheSh__(counts[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2,t_ndfs);
	    cout <<"mb2: " << endl;
	    asymCharmM=fitTheSh__(countsCharm[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2,t_ndfs);
	    asymUdsM=fitTheSh__(countsUds[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2,t_ndfs);
	    cout <<"theta asym M" <<endl;
	    asymAllMTheta=fitTheSh__(countsTheta[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2Theta,t_ndfsTheta);
	    asymCharmMTheta=fitTheSh__(countsCharmTheta[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2Theta,t_ndfsTheta);
	    asymUdsMTheta=fitTheSh__(countsCharmTheta[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2Theta,t_ndfsTheta);
	    cout <<"mb3: " << endl;
	    v_chi2M.push_back(t_chi2);
	    v_ndfsM.push_back(t_ndfs);
	    break;
	  }
	default:
	  cout << "binning not recognized2" << endl;
	  exit(0);
	}
    }

  //danger: assumption for drawthetaasyms: that all binnings have the same size
  drawThetaAsyms(asymAllMTheta,asymAllZTheta,asymCharmMTheta,asymCharmZTheta, xValsTheta,xValsNumTheta,xValsCharmTheta,xValsNumCharmTheta);
  cout <<" z binning: " << endl;
  int iCh1=PN;
  int iCh2=PN;
  int iPa1=PiPi;
  int iPa2=PiPi;

  numBins=binningZ[iPa1].size();
  oFile <<"Zbinning"<<endl;


  float xRevCutInt[numBins];
  float xAllInt[numBins];
  float xCutInt[numBins];
  float xUdsInt[numBins];
  float xCharmInt[numBins];

  float yRevCutInt[numBins];
  float yAllInt[numBins];
  float yCutInt[numBins];
  float yUdsInt[numBins];
  float yCharmInt[numBins];

  float exInt[numBins];

  float eyRevCutInt[numBins];
  float eyAllInt[numBins];
  float eyCutInt[numBins];
  float eyUdsInt[numBins];
  float eyCharmInt[numBins];

  float wAllSum[numBins];
  float wRevCutSum[numBins];
  float wCutSum[numBins];
  float wUdsSum[numBins];
  float wCharmSum[numBins];

  float xAllS[numBins];
  float xCutS[numBins];
  float xRevCutS[numBins];
  float yAllS[numBins];
  float yCutS[numBins];
  float yRevCutS[numBins];

  float eyAllS2[numBins];
  float eyCutS2[numBins];

  float xUdsS[numBins];
  float yUdsS[numBins];
  float eyUdsS2[numBins];
  float xCharmS[numBins];
  float yCharmS[numBins];
  float eyCharmS2[numBins];
  float ex[numBins];

  float meanZinM[numBins];
  float meanMinZ[numBins];

  setZero(xUdsInt,xCharmInt,yUdsInt,exInt,eyUdsInt,eyCharmInt,numBins);
  setZero(wUdsSum,wCharmSum,xUdsS,yUdsS,eyUdsS2,xCharmS,numBins);
  setZero(wAllSum,wCutSum,xAllS,yAllS,xCutS,yCutS,numBins);
  setZero(yCharmS,eyCharmS2,ex,numBins);
  setZero(xRevCutInt,yRevCutInt,yRevCutS,wRevCutSum,xRevCutS,numBins);
  for(int iZ=0;iZ<numBins;iZ++)
    {
      float wUds=0;
      float wCharm=0;
      float wAll=0;
      float wCut=0;
      float wRevCut=0;
      float xUds[numBins];
      float xCharm[numBins];
      float yUds[numBins];
      float yCharm[numBins];

      float eyUds[numBins];
      float eyCharm[numBins];

      float yAll[numBins];
      float yCut[numBins];
      float yRevCut[numBins];
      float xAll[numBins];
      float xCut[numBins];
      float xRevCut[numBins];
      float eyAll[numBins];
      float eyCut[numBins];
      float eyRevCut[numBins];

      oFile <<"bin " << iZ <<endl;
      for(int iZ2=0;iZ2<numBins;iZ2++)
	{
	  float A_All=(*asymAllZ)[iZ*binningZ[iPa1].size()+iZ2].first;
	  float A_Cut=(*asymCharmZ)[iZ*binningZ[iPa1].size()+iZ2].first;
	  float A_RevCut=(*asymUdsZ)[iZ*binningZ[iPa1].size()+iZ2].first;
	  yAll[iZ2]=A_All;
	  yCut[iZ2]=A_Cut;
	  yRevCut[iZ2]=A_RevCut;
	  double allInTestSample=passAsCharmZ[iZ][iZ2]+noPassAsCharmZ[iZ][iZ2]+passAsCharmZUds[iZ][iZ2]+noPassAsCharmZUds[iZ][iZ2];
	  double allPassedAsCharm=passAsCharmZ[iZ][iZ2]+passAsCharmZUds[iZ][iZ2];
	  double allNotPassedAsCharm=noPassAsCharmZ[iZ][iZ2]+noPassAsCharmZUds[iZ][iZ2];

	  //xvalsNum the right number?
	  double scaleFact=(double)1/allInTestSample; //amount as compared to all
	  double scaleFact2=(double)1/allPassedAsCharm; //amount as compared to all
	  double scaleFact3=(double)1/allNotPassedAsCharm;

	  double k2=passAsCharmZ[iZ][iZ2]*scaleFact2; //number of events expected in all sample
	  double k=(passAsCharmZ[iZ][iZ2]+noPassAsCharmZ[iZ][iZ2])*scaleFact; //num charm in all sample is num of events in charm mc (no matter if surv. the cut)
	  double kNoPass=noPassAsCharmZ[iZ][iZ2]*scaleFact3;

	  double l=(passAsCharmZUds[iZ][iZ2]+noPassAsCharmZUds[iZ][iZ2])*scaleFact;
	  double l2=passAsCharmZUds[iZ][iZ2]*scaleFact2;
	  double lNoPass=noPassAsCharmZUds[iZ][iZ2]*scaleFact3;

	  cout <<"k+l: " << k+l << ", k2+l2: " << k2+l2 <<endl;
	  ex[iZ2]=0;
	  xUds[iZ2]=xVals[zBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]/xValsNum[zBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
	  xCharm[iZ2]=xUds[iZ2]+0.01;
	  xAll[iZ2]=xUds[iZ2]+0.02;
	  xCut[iZ2]=xUds[iZ2]+0.03;
	  xRevCut[iZ2]=xUds[iZ2]+0.04;
	  //	  yUds[iZ2]=A_Cut-(k2*A_All /(k*l2-k2*l2));
	  ///	  yUds[iZ2]=(k*A_Cut-k2*A_All)/(k*l2-k2*l);
	  ///l->lNoPass, k->kNoPass, A_All->ARevCut
	  yUds[iZ2]=(kNoPass*A_Cut-k2*A_RevCut)/(kNoPass*l2-k2*lNoPass);
	  //	  yCharm[iZ2]=(l2*A_All-l*A_Cut)/(k*l2-k2*l);
	  ///	  yCharm[iZ2]=A_All/k -l*(A_Cut-(k2/k)*A_All)/(k*l2-k2*l);
	  yCharm[iZ2]=A_RevCut/kNoPass -lNoPass*(A_Cut-(k2/kNoPass)*A_RevCut)/(kNoPass*l2-k2*lNoPass);
	  cout <<"a_uds: " << yUds[iZ2] <<" charm: " << yCharm[iZ2] << endl;
	  //	  double el=sqrt(l);//not true since binomial  //variance for binomial is n*P*(1-P), P=l, n=1/scalefact
	  double el=sqrt(allInTestSample*(double)l*(1-l))/allInTestSample;
	  double el2=sqrt(allPassedAsCharm*(double)l2*(1-l2))/allPassedAsCharm;


	  double ek=sqrt(allInTestSample*(double)k*(1-k))/allInTestSample;
	  double ek2=sqrt(allPassedAsCharm*(double)k2*(1-k2))/allPassedAsCharm; //divided by N again...
	  double ekNoPass=sqrt(allNotPassedAsCharm*(double)kNoPass*(1-kNoPass))/allNotPassedAsCharm; //divided by N again...
	  double elNoPass=sqrt(allNotPassedAsCharm*(double)lNoPass*(1-lNoPass))/allNotPassedAsCharm;
	  cout <<"l: " << l <<" l2: " << l2<< " k: " << k <<" k2: " << k2 <<" el: " << el << " el2: " << el2 << " ek: "<< ek << " ek2: " << ek2 <<endl; 

	  double eA_All=(*asymAllZ)[iZ*binningZ[iPa1].size()+iZ2].second;
	  double eA_Cut=(*asymCharmZ)[iZ*binningZ[iPa1].size()+iZ2].second;
	  double eA_RevCut=(*asymUdsZ)[iZ*binningZ[iPa1].size()+iZ2].second;
	  cout <<"A_All: " << A_All <<"A_revCut: " << A_RevCut <<" A_cut: " << A_Cut <<" eAll: " << eA_All<<" eARevCut: " << eA_RevCut  << " eACut: " << eA_Cut <<endl;
	  eyAll[iZ2]=eA_All;
	  eyCut[iZ2]=eA_Cut;
	  eyRevCut[iZ2]=eA_RevCut;
	  /*	  double lowFact=k*l2-k2*l;
	    double dAuds_dAcut=k/lowFact;
	    double dAuds_dAall=k2/lowFact;
	    double dAuds_dk=(A_Cut*lowFact-(k*A_Cut-k2*A_All)*l2)/(lowFact*lowFact);
	    double dAuds_dk2=(-1)*(A_All*lowFact+(k*A_Cut-k2*A_All)*l)/(lowFact*lowFact);
	    double dAuds_dl=k2*(k*A_Cut-k2*A_All)/(lowFact*lowFact);
	    double dAuds_dl2=(-1)*k*(k*A_Cut-k2*A_All)/(lowFact*lowFact);

	    double dAcharm_dAcut=(-1)*l/lowFact;
	    double dAcharm_dAall=(double)1/k-l*k2/(k*lowFact);
	    double lowFact2=k*k*l2-k*k2*l;
	    double dAcharm_dk=(-1)*A_All/(k*k)-(l*A_Cut*(k*k*l2-k*k2*l)-(2*k*l2-k2*l)*(l*k*A_Cut-k*k2*A_All))/(lowFact2*lowFact2);
	    double dAcharm_dk2=(l*A_All*lowFact2-k*l*(l*k*A_Cut-l*k2*A_All))/(lowFact2*lowFact2);
	    double dAcharm_dl=(-1)*((k*A_Cut-k2*A_All)*lowFact+k*k2*(l*k*A_Cut-l*k2*A_All))/(lowFact2*lowFact2);
	    double dAcharm_dl2=k*k*(l*k*A_Cut-l*k2*A_All)/(lowFact2*lowFact2);*/

	  //lets do it again with k+l=1, should lead to smaller errors
	  double lowFact=(k-k2)*(k-k2);

	  double lowFact2=(k*k-k*k2)*(k*k-k*k2);
	  double lowFactNoPass=(kNoPass-k2)*(kNoPass-k2);
	  double lowFactNoPass2=(kNoPass*kNoPass-kNoPass*k2)*(kNoPass*kNoPass-kNoPass*k2);
	  ///double dAuds_dAcut=k/(double)(k-k2);
	  double dAuds_dAcut=kNoPass/(double)(kNoPass-k2);
///	  double dAuds_dAall=-k2/(double)(k-k2)
	  double dAuds_dARevCut=-k2/(double)(kNoPass-k2);


	  ///	  double dAuds_dk=A_Cut*(k-k2)-(k*A_Cut-k2*A_All)/lowFact;
	  double dAuds_dkNoPass=A_Cut*(kNoPass-k2)-(kNoPass*A_Cut-k2*A_RevCut)/lowFactNoPass;
	  ///	  double dAuds_dk2=((-1)*A_All*(k-k2)+(k*A_Cut-k2*A_All))/lowFact;
	  double dAuds_dk2=((-1)*A_RevCut*(kNoPass-k2)+(kNoPass*A_Cut-k2*A_RevCut))/lowFactNoPass;

	  ///	  double dAcharm_dAcut=(-1)/(k-k2);
	  double dAcharm_dAcut=(-1)/(kNoPass-k2);
	  ///	  double dAcharm_dAall=k2/(k*(k-k2));
	  double dAcharm_dARevCut=k2/(kNoPass*(kNoPass-k2));
	  ///	  double dAcharm_dk=(-1)/(k*k)*A_All-(A_Cut*(k*k-k*k2)-(2*k-k2)*(k*A_Cut-k2*A_All))/(lowFact2);
	  double dAcharm_dkNoPass=(-1)/(kNoPass*kNoPass)*A_RevCut-(A_Cut*(kNoPass*kNoPass-kNoPass*k2)-(2*kNoPass-k2)*(kNoPass*A_Cut-k2*A_RevCut))/(lowFactNoPass2);
	  ///	  double dAcharm_dk2=(-A_All*(k*k-k*k2)+k*k*A_Cut-k*k2*A_All)/lowFact2;
	  double dAcharm_dk2=(-A_RevCut*(kNoPass*kNoPass-kNoPass*k2)+kNoPass*kNoPass*A_Cut-kNoPass*k2*A_RevCut)/lowFactNoPass2;

	  dAuds_dAcut*=(dAuds_dAcut*eA_Cut*eA_Cut);
	  ///	  dAuds_dAall*=(dAuds_dAall*eA_All*eA_All);
	  dAuds_dARevCut*=(dAuds_dARevCut*eA_RevCut*eA_RevCut);
	  ///	  dAuds_dk*=(dAuds_dk*ek*ek);
	  dAuds_dkNoPass*=(dAuds_dkNoPass*ekNoPass*ekNoPass);
	  dAuds_dk2*=(dAuds_dk2*ek2*ek2);
	  //	  dAuds_dl*=(dAuds_dl*el*el);
	  //	  dAuds_dl2*=(dAuds_dl2*el2*el2);
	  dAcharm_dAcut*=(dAcharm_dAcut*eA_Cut*eA_Cut);
	  ///	  dAcharm_dAall*=(dAcharm_dAall*eA_All*eA_All);
	  dAcharm_dARevCut*=(dAcharm_dARevCut*eA_RevCut*eA_RevCut);
	  //	  dAcharm_dk*=(dAcharm_dk*ek*ek);
	  dAcharm_dkNoPass*=(dAcharm_dkNoPass*ekNoPass*ekNoPass);
	  dAcharm_dk2*=(dAcharm_dk2*ek2*ek2);
	  //	  dAcharm_dl*=(dAcharm_dl*el*el);
	  //	  dAcharm_dl2*=(dAcharm_dl2*el2*el2);
	  //	  cout <<"lowfact: " << lowFact << "lowFact2: " << lowFact2 << " d1: " << dAuds_dAcut <<" d2: " << dAuds_dAall <<" d3: " << dAuds_dk;
	  cout <<" d4: " << dAuds_dk2 << endl;
	  //	  cout <<"charm, cut: " << dAcharm_dAcut << " daall: " << dAcharm_dAall << " dk: " << dAcharm_dk <<" kd2: " << dAcharm_dk2 <<endl;
	  //	  eyUds[iZ2]=sqrt(dAuds_dAcut+dAuds_dAall+dAuds_dk+dAuds_dk2+dAuds_dl+dAuds_dl2);
	  //	  eyCharm[iZ2]=sqrt(dAcharm_dAcut+dAcharm_dAall+dAcharm_dk+dAcharm_dk2+dAcharm_dl+dAcharm_dl2);
	  ///eyUds[iZ2]=sqrt(dAuds_dAcut+dAuds_dAall+dAuds_dk+dAuds_dk2);
	  eyUds[iZ2]=sqrt(dAuds_dAcut+dAuds_dARevCut+dAuds_dkNoPass+dAuds_dk2);
	  ///eyCharm[iZ2]=sqrt(dAcharm_dAcut+dAcharm_dAall+dAcharm_dk+dAcharm_dk2);
	  eyCharm[iZ2]=sqrt(dAcharm_dAcut+dAcharm_dARevCut+dAcharm_dkNoPass+dAcharm_dk2);
	  oFile <<iZ2<<" uds: " << eyUds[iZ2]  <<" charm: " << eyCharm[iZ2]<<endl;
	  errorSumUdsZ+=eyUds[iZ2];
	  errorSumCharmZ+=eyCharm[iZ2];

	  wUds=1/(eyUds[iZ2]*eyUds[iZ2]);
	  wCharm=1/(eyCharm[iZ2]*eyCharm[iZ2]);
	  wAll=1/(eyAll[iZ2]*eyAll[iZ2]);
	  wRevCut=1/(eyRevCut[iZ2]*eyRevCut[iZ2]);
	  wCut=1/(eyCut[iZ2]*eyCut[iZ2]);

	  wUdsSum[iZ2]+=wUds;
	  wCharmSum[iZ2]+=wCharm;
	  wAllSum[iZ2]+=wAll;
	  wRevCutSum[iZ2]+=wRevCut;
	  wCutSum[iZ2]+=wCut;

	  yUdsS[iZ2]+=wUds*yUds[iZ2];
	  yCharmS[iZ2]+=wCharm*yCharm[iZ2];
	  yAllS[iZ2]+=wAll*yAll[iZ2];
	  yRevCutS[iZ2]+=wRevCut*yRevCut[iZ2];
	  yCutS[iZ2]+=wCut*yCut[iZ2];
	  xAllS[iZ2]+=wAll*xAll[iZ2];
	  xRevCutS[iZ2]+=wRevCut*xRevCut[iZ2];
	  xCutS[iZ2]+=wCut*xCut[iZ2];
	  xUdsS[iZ2]+=wUds*xUds[iZ2];
	  cout <<"adding at " << iZ2 <<" " << wUds << " * " << xUds[iZ2] <<endl;
	  xCharmS[iZ2]+=wCharm*xCharm[iZ2];
	  //	  eyCharm[iZ2]=sqrt(eFirst*eFirst+eSecond*eSecond+eThird*eThird);
	  //	  eyUds[iZ2]=sqSum(eUdsUp/(k*A_Cut-k2*A_All),eLowFact/(k*l2-k2*l2))*yUds[iZ2];

	  double eFirst=sqSum(eA_All/A_All,ek/k)*(A_All/k);
	  double eLowFact1=sqSum(ek/k,el2/l2)*k*l2;
	  double eLowFact2=sqSum(ek2/k2,el/l)*k2*l;
	  double eLowFact=sqSum(eLowFact1,eLowFact2);
	  double eSecond=sqSum(el/l,eLowFact/(k*l2-k2*l),eA_Cut/A_Cut);; //error on hte factor before A_cut
	  eSecond*=A_Cut*l/(k*l2-k2*l); //no div, because we omitted the mult earlier
	  double eThird=sqSum(ek2/k2,ek/k,el/l,eA_All/A_All,eLowFact/(k*l2-k2*l));
	  eThird*=(l*k2/k)*A_All/(k*l2-k2*l);
	  //for uds:
	  double eUdsUp= sqSum(sqSum(ek/k,eA_Cut/A_Cut)*k*A_Cut,sqSum(ek2/k2,eA_All,A_All)*k2*A_All);
	  cout <<"Z binning: eall: " << eA_All <<" eaCut: " << eA_Cut <<" first: " << eFirst << ": elowfact1: " << eLowFact1 << " lowfact2: " << eLowFact2 << " eLowFact " << eLowFact <<" esec: " << eSecond  << " third: " << eThird << " eUdsUp : " << eUdsUp <<endl;
	  //	  eyCharm[iZ2]=sqrt(eFirst*eFirst+eSecond*eSecond+eThird*eThird);
	  //	  eyUds[iZ2]=sqSum(eUdsUp/(k*A_Cut-k2*A_All),eLowFact/(k*l2-k2*l2))*yUds[iZ2];
	  //	  eyUds[iZ2]=(*asymCharmZ)[k].second;
	  //	  eyCharm[iZ2]=(*asymCharmZ)[k].second; //not correct, have to take all the factors into account...

	  //      outFile  <<"Bin " << k << " Asymmetry: " << (*asymAllZ)[k].first << " +- " <<  (*asymAllZ)[k].second <<endl;

	}
      stringstream sstU;
      stringstream sstC;
      sstU <<"Both5_asym_iZ"<<iZ;
      TGraphErrors tgU(binningZ[iPa1].size(),xUds,yUds,ex,eyUds);
      TGraphErrors tgC(binningZ[iPa1].size(),xCharm,yCharm,ex,eyCharm);
      TGraphErrors tgAll(binningZ[iPa1].size(),xAll,yAll,ex,eyAll);
      TGraphErrors tgCut(binningZ[iPa1].size(),xCut,yCut,ex,eyCut);
      TGraphErrors tgRevCut(binningZ[iPa1].size(),xRevCut,yRevCut,ex,eyRevCut);
      tgC.SetMarkerSize(1.2);
      tgC.SetMarkerColor(kRed);
      tgU.SetMarkerSize(1.2);
      tgU.SetMarkerColor(kBlue);
      tgAll.SetMarkerSize(1.2);
      tgAll.SetMarkerColor(kBlack);
      tgCut.SetMarkerSize(1.2);
      tgCut.SetMarkerColor(kGreen);
      tgRevCut.SetMarkerSize(1.2);
      tgRevCut.SetMarkerColor(kPink);
      tgRevCut.SetMarkerStyle(kFullTriangleDown);

      TCanvas cU(sstU.str().c_str(),sstU.str().c_str(),10,50,500,800);
      float gMin=0;
      float gMax=0;
      /*	  if(tgU.GetMaximum()> tgC.GetMaximum())
	gMax=tgU.GetMaximum()+0.01;
	else
	gMax=tgC.GetMaximum()+0.01;
	if(tgU.GetMinimum()< tgC.GetMinimum())
	gMin=tgU.GetMinimum()-0.01;
	else
	gMin=tgC.GetMinimum()+0.01;*/
      gMin=-0.3;
      gMax=0.3; // to have space for the legend
      tgRevCut.GetYaxis()->SetRangeUser(gMin,gMax);
      tgRevCut.GetXaxis()->SetRangeUser(0.0,1.0);
      tgRevCut.GetXaxis()->SetLimits(0.0,1.0);
      tgRevCut.Draw("AP*");
      tgU.Draw("P* SAME");
      tgC.Draw("P* SAME");
      tgAll.Draw("P* SAME");
      tgCut.Draw("P* SAME");
      TLegend leg(0.1,0.7,0.48,0.9);
      //	  leg.SetHeader("Legend");
      leg.AddEntry(&tgRevCut,"Reverse Cut","p");
      leg.AddEntry(&tgU,"UDS","p");
      leg.AddEntry(&tgC,"charm","p");
      leg.AddEntry(&tgAll,"All","p");
      leg.AddEntry(&tgCut,"Cut","p");
      leg.Draw();

#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
      cU.SaveAs(("AsUdsCharmOnlyDist/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharmOnlyDist/"+sstU.str()+c_tmp+".eps").c_str());
#else
      cU.SaveAs(("AsUdsCharmWOKaon/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharmWOKaon/"+sstU.str()+c_tmp+".eps").c_str());
#endif	
#else
      cU.SaveAs(("AsUdsCharm/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharm/"+sstU.str()+c_tmp+".eps").c_str());
#endif
#else
      cU.SaveAs(("AsUdsCharmOnlyMass/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharmOnlyMass/"+sstU.str()+c_tmp+".eps").c_str());
#endif	
    }
  for(int iZ=0;iZ<numBins;iZ++)
    {
      yUdsInt[iZ]=yUdsS[iZ]/wUdsSum[iZ];
      yCharmInt[iZ]=yCharmS[iZ]/wCharmSum[iZ];
      yAllInt[iZ]=yAllS[iZ]/wAllSum[iZ];
      yCutInt[iZ]=yCutS[iZ]/wCutSum[iZ];
      yRevCutInt[iZ]=yRevCutS[iZ]/wRevCutSum[iZ];

      eyCharmInt[iZ]=sqrt(1/wCharmSum[iZ]);
      eyUdsInt[iZ]=sqrt(1/wUdsSum[iZ]);
      eyAllInt[iZ]=sqrt(1/wAllSum[iZ]);
      eyCutInt[iZ]=sqrt(1/wCutSum[iZ]);
      eyRevCutInt[iZ]=sqrt(1/wRevCutSum[iZ]);

      xUdsInt[iZ]=xUdsS[iZ]/wUdsSum[iZ];
      xCharmInt[iZ]=xCharmS[iZ]/wCharmSum[iZ];
      xAllInt[iZ]=xAllS[iZ]/wAllSum[iZ];
      xCutInt[iZ]=xCutS[iZ]/wCutSum[iZ];
      xRevCutInt[iZ]=xRevCutS[iZ]/wRevCutSum[iZ];
    }

  TGraphErrors tgUAll(binningZ[iPa1].size(),xUdsInt,yUdsInt,ex,eyUdsInt);
  TGraphErrors tgCAll(binningZ[iPa1].size(),xCharmInt,yCharmInt,ex,eyCharmInt);
  TGraphErrors tgAllAllZ(binningZ[iPa1].size(),xAllInt,yAllInt,ex,eyAllInt);
  TGraphErrors tgCutAllZ(binningZ[iPa1].size(),xCutInt,yCutInt,ex,eyCutInt);
  TGraphErrors tgRevCutAllZ(binningZ[iPa1].size(),xRevCutInt,yRevCutInt,ex,eyRevCutInt);


  tgCAll.SetMarkerSize(1.2);
  tgCAll.SetMarkerColor(kRed);
  tgUAll.SetMarkerSize(1.2);
  tgUAll.SetMarkerColor(kBlue);
  tgAllAllZ.SetMarkerSize(1.2);
  tgAllAllZ.SetMarkerColor(kBlack);
  tgCutAllZ.SetMarkerSize(1.2);
  tgCutAllZ.SetMarkerColor(kGreen);
  tgRevCutAllZ.SetMarkerSize(1.2);
  tgRevCutAllZ.SetMarkerColor(kPink);
  tgRevCutAllZ.SetMarkerStyle(kFullTriangleDown);
  stringstream sstA;
  sstA <<"All_asym_Z";
  TCanvas cAll(sstA.str().c_str(),sstA.str().c_str(),10,50,500,800);
  float gMin=0;
  float gMax=0;
  gMin=-0.1;
  gMax=0.1; // to have space for the legend
  tgRevCutAllZ.GetYaxis()->SetRangeUser(gMin,gMax);
  tgRevCutAllZ.GetXaxis()->SetRangeUser(0.0,1.0);
  tgRevCutAllZ.GetXaxis()->SetLimits(0.0,1.0);
  tgRevCutAllZ.Draw("AP*");
  tgUAll.Draw("P* SAME");
  tgCAll.Draw("P* SAME");
  tgAllAllZ.Draw("P* SAME");
  tgCutAllZ.Draw("P* SAME");
  TLegend leg(0.1,0.7,0.48,0.9);
  leg.AddEntry(&tgRevCutAllZ,"Reverse Cut","p");
  leg.AddEntry(&tgUAll,"UDS","p");
  leg.AddEntry(&tgCAll,"charm","p");
  leg.AddEntry(&tgAllAllZ,"All","p");
  leg.AddEntry(&tgCutAllZ,"Cut","p");
  leg.Draw();

#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
  cAll.SaveAs(("AsUdsCharmOnlyDist/"+sstA.str()+c_tmp+".png").c_str());
  cAll.SaveAs(("AsUdsCharmOnlyDist/"+sstA.str()+c_tmp+".eps").c_str());
  cAll.SaveAs(("AsUdsCharmOnlyDist/"+sstA.str()+c_tmp+".C").c_str());
#else
  cAll.SaveAs(("AsUdsCharmWOKaon/"+sstA.str()+c_tmp+".png").c_str());
  cAll.SaveAs(("AsUdsCharmWOKaon/"+sstA.str()+c_tmp+".eps").c_str());
  cAll.SaveAs(("AsUdsCharmWOKaon/"+sstA.str()+c_tmp+".C").c_str());
#endif
#else
  cAll.SaveAs(("AsUdsCharm/"+sstA.str()+c_tmp+".png").c_str());
  cAll.SaveAs(("AsUdsCharm/"+sstA.str()+c_tmp+".eps").c_str());
  cAll.SaveAs(("AsUdsCharm/"+sstA.str()+c_tmp+".C").c_str());
#endif
#else
  cAll.SaveAs(("AsUdsCharmOnlyMass/"+sstA.str()+c_tmp+".png").c_str());
  cAll.SaveAs(("AsUdsCharmOnlyMass/"+sstA.str()+c_tmp+".eps").c_str());
  cAll.SaveAs(("AsUdsCharmOnlyMass/"+sstA.str()+c_tmp+".C").c_str());
#endif
  drawStats(passAsCharmZ,noPassAsCharmZ, passAsCharmZUds, noPassAsCharmZUds,xVals[zBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],xValsNum[zBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningZ[iPa1].size(), (string("ZBinning")+c_tmp).c_str());
  drawStats(passAsCharmM,noPassAsCharmM, passAsCharmMUds, noPassAsCharmMUds, xVals[mBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],xValsNum[mBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningM[iPa1].size(),(string("MBinning")+c_tmp).c_str());
  drawStats(passAsCharmMTheta,noPassAsCharmMTheta, passAsCharmMThetaUds, noPassAsCharmMThetaUds, xValsTheta[mBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],xValsNumTheta[mBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningTheta[iPa1].size(),(string("MThetaBinning")+c_tmp).c_str());
  drawStats(passAsCharmZTheta,noPassAsCharmZTheta, passAsCharmZThetaUds, noPassAsCharmZThetaUds, xValsTheta[zBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],xValsNumTheta[zBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningTheta[iPa1].size(),(string("ZThetaBinning")+c_tmp).c_str());

  numBins=binningM[iPa1].size();
  drawMeanKins(meanKin,meanKinPassCuts,meanKinFailCuts,maxKin,c_tmp);
  oFile <<"M binning"  <<endl;
  setZero(xUdsInt,xCharmInt,yUdsInt,exInt,eyUdsInt,eyCharmInt,numBins);
  setZero(wUdsSum,wCharmSum,xUdsS,yUdsS,eyUdsS2,xCharmS,numBins);
  setZero(wAllSum,wCutSum,xAllS,yAllS,xCutS,yCutS,numBins);
  setZero(yCharmS,eyCharmS2,ex,numBins);
  setZero(xRevCutInt,yRevCutInt,yRevCutS,wRevCutSum,xRevCutS,numBins);

  for(int iM=0;iM<numBins;iM++)
    {
      float wUds=0;
      float wCharm=0;
      float wAll=0;
      float wCut=0;
      float wRevCut=0;

      float xUds[numBins];
      float xCharm[numBins];
      float yUds[numBins];
      float yCharm[numBins];
      float ex[numBins];
      float eyUds[numBins];
      float eyCharm[numBins];
      oFile <<"bin " << iM <<endl;

      float yAll[numBins];
      float yCut[numBins];
      float yRevCut[numBins];
      float xAll[numBins];
      float xCut[numBins];
      float xRevCut[numBins];
      float eyAll[numBins];
      float eyCut[numBins];
      float eyRevCut[numBins];

      for(int iM2=0;iM2<numBins;iM2++)
	{
	  float A_All=(*asymAllM)[iM*binningM[iPa1].size()+iM2].first;
	  float A_Cut=(*asymCharmM)[iM*binningM[iPa1].size()+iM2].first;
	  float A_RevCut=(*asymUdsM)[iM*binningM[iPa1].size()+iM2].first;
	  yAll[iM2]=A_All;
	  yCut[iM2]=A_Cut;
	  yRevCut[iM2]=A_RevCut;
	  double allInTestSample=passAsCharmM[iM][iM2]+noPassAsCharmM[iM][iM2]+passAsCharmMUds[iM][iM2]+noPassAsCharmMUds[iM][iM2];
	  double allPassedAsCharm=passAsCharmM[iM][iM2]+passAsCharmMUds[iM][iM2];
	  double allNotPassedAsCharm=noPassAsCharmM[iM][iM2]+noPassAsCharmMUds[iM][iM2];
	  //xvalsNum the right number?
	  cout <<"xValsNum: " << xValsNum[mBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2] << " all: " << allInTestSample <<endl;
	  cout <<"yValsNum: " << yValsNum[mBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2] <<endl;
	  cout <<"xValsNumCharm: " << xValsNumCharm[mBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2] << " allCharm: " << allInTestSample <<" allpassed as charm: " << allPassedAsCharm << "  charm" << passAsCharmM[iM][iM2] << " in uds: " << passAsCharmMUds[iM][iM2] <<endl;
	  cout <<" noPass in charm: " << noPassAsCharmM[iM][iM2] << " in uds: " << noPassAsCharmMUds[iM][iM2] <<endl;
	  double scaleFact=(double)1/allInTestSample; //amount as compared to all
	  double scaleFact2=(double)1/allPassedAsCharm; //amount as compared to all
	  double scaleFact3=(double)1/allNotPassedAsCharm;

	  double k=(passAsCharmM[iM][iM2]+noPassAsCharmM[iM][iM2])*scaleFact; //num charm in all sample is num of events in charm mc (no matter if surv. the cut)

	  double kNoPass=noPassAsCharmZ[iM][iM2]*scaleFact3;

	  double l=(passAsCharmMUds[iM][iM2]+noPassAsCharmMUds[iM][iM2])*scaleFact;
	  double k2=passAsCharmM[iM][iM2]*scaleFact2; //number of events expected in all sample
	  double l2=passAsCharmMUds[iM][iM2]*scaleFact2;
	  double lNoPass=noPassAsCharmMUds[iM][iM2]*scaleFact3;

	  cout <<"k+l: " << k+l << ", k2+l2: " << k2+l2 <<endl;
	  cout <<"scaleFact: " << scaleFact <<" sc2: " << scaleFact2 << " l: " << l <<" k: " << k << " l2: "<< l2 <<" k2: " << k2 <<endl;
	  ex[iM2]=0;
	  xUds[iM2]=xVals[mBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]/xValsNum[mBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
	  xCharm[iM2]=xUds[iM2]+0.01;
	  xAll[iM2]=xUds[iM2]+0.02;
	  xCut[iM2]=xUds[iM2]+0.03;
	  xRevCut[iM2]=xUds[iM2]+0.04;
	  cout <<"yUds[iM2]: " << yUds[iM2] << endl;
	  ///	  yUds[iM2]=(k*A_Cut-k2*A_All)/(k*l2-k2*l);
	  yUds[iM2]=(kNoPass*A_Cut-k2*A_RevCut)/(kNoPass*l2-k2*lNoPass);
	  cout <<"yUds after set: " << yUds[iM2] << "lower fact: " << k*l2-k2*l <<" upper fact: " << k*A_Cut-k2*A_All <<endl;
	  //	  yCharm[iM2]=(l2*A_All-l*A_Cut)/(k*l2-k2*l);
	  ///	  yCharm[iM2]=A_All/k -l*(k*A_Cut-k2*A_All)/(k*(k*l2-k2*l));
	  yCharm[iM2]=A_RevCut/kNoPass -lNoPass*(kNoPass*A_Cut-k2*A_RevCut)/(kNoPass*(kNoPass*l2-k2*lNoPass));

	  cout <<"im: " << iM <<" im2: " << iM2 << " yuds: " << yUds[iM2] << " ycharm: " << yCharm[iM2] << " a_all: " << A_All <<" A_Cut: " << A_Cut <<" error all: " << (*asymAllM)[iM*binningM[iPa1].size()+iM2].second << " error charm: " <<(*asymCharmM)[iM*binningM[iPa1].size()+iM2].second <<endl;
	  //	  double el=sqrt(l);
	  //	  double el2=sqrt(l2);
	  //	  double ek=sqrt(k);
	  //	  double ek2=sqrt(k2);
	  double el=sqrt(allInTestSample*(double)l*(1-l))/allInTestSample;
	  double el2=sqrt(allPassedAsCharm*(double)l2*(1-l2))/allPassedAsCharm;
	  double ek=sqrt(allInTestSample*(double)k*(1-k))/allInTestSample;
	  double ek2=sqrt(allPassedAsCharm*(double)k2*(1-k2))/allPassedAsCharm;
	  double ekNoPass=sqrt(allNotPassedAsCharm*(double)kNoPass*(1-kNoPass))/allNotPassedAsCharm; //divided by N again...
	  double elNoPass=sqrt(allNotPassedAsCharm*(double)lNoPass*(1-lNoPass))/allNotPassedAsCharm;
	  cout <<"l: " << l <<" l2: " << l2<< " k: " << k <<" k2: " << k2 <<" el: " << el << " el2: " << el2 << " ek: "<< ek << " ek2: " << ek2 <<endl; 
	  double eA_All=(*asymAllM)[iM*binningM[iPa1].size()+iM2].second;
	  double eA_Cut=(*asymCharmM)[iM*binningM[iPa1].size()+iM2].second;
	  double eA_RevCut=(*asymUdsM)[iM*binningM[iPa1].size()+iM2].second;
	  eyAll[iM2]=eA_All;
	  eyCut[iM2]=eA_Cut;
	  eyRevCut[iM2]=eA_RevCut;

	  double eFirst=sqSum(eA_All/A_All,ek/k)*(A_All/k);
	  double eLowFact1=sqSum(ek/k,el2/l2)*k*l2;
	  double eLowFact2=sqSum(ek2/k2,el/l)*k2*l;
	  double eLowFact=sqSum(eLowFact1,eLowFact2);
	  double eSecond=sqSum(el/l,eLowFact/(k*l2-k2*l),eA_Cut/A_Cut);; //error on hte factor before A_cut
	  eSecond*=A_Cut*l/(k*l2-k2*l); //no div, because we omitted the mult earlier
	  double eThird=sqSum(ek2/k2,ek/k,el/l,eA_All/A_All,eLowFact/(k*l2-k2*l));
	  eThird*=(l*k2/k)*A_All/(k*l2-k2*l);
	  //for uds:
	  double eUdsUp= sqSum(sqSum(ek/k,eA_Cut/A_Cut)*k*A_Cut,sqSum(ek2/k2,eA_All,A_All)*k2*A_All);
	  cout <<"eall: " << eA_All <<" eaCut: " << eA_Cut <<" first: " << eFirst << ": elowfact1: " << eLowFact1 << " lowfact2: " << eLowFact2 << " eLowFact " << eLowFact <<" esec: " << eSecond  << " third: " << eThird << " eUdsUp : " << eUdsUp <<endl;


	  /*	  double lowFact=k*l2-k2*l;
	    double dAuds_dAcut=k/lowFact;
	    double dAuds_dAall=k2/lowFact;
	    double dAuds_dk=(A_Cut*lowFact-(k*A_Cut-k2*A_All)*l2)/(lowFact*lowFact);
	    double dAuds_dk2=(-1)*(A_All*lowFact+(k*A_Cut-k2*A_All)*l)/(lowFact*lowFact);
	    double dAuds_dl=k2*(k*A_Cut-k2*A_All)/(lowFact*lowFact);
	    double dAuds_dl2=(-1)*k*(k*A_Cut-k2*A_All)/(lowFact*lowFact);*/


	  /*	  dAuds_dAcut*=(dAuds_dAcut*eA_Cut*eA_Cut);
	    dAuds_dAall*=(dAuds_dAall*eA_All*eA_All);
	    dAuds_dk*=(dAuds_dk*ek*ek);
	    dAuds_dk2*=(dAuds_dk2*ek2*ek2);
	    dAuds_dl*=(dAuds_dl*el*el);
	    dAuds_dl2*=(dAuds_dl2*el2*el2);

	    dAcharm_dAcut*=(dAcharm_dAcut*eA_Cut*eA_Cut);
	    dAcharm_dAall*=(dAcharm_dAall*eA_All*eA_All);
	    dAcharm_dk*=(dAcharm_dk*ek*ek);
	    dAcharm_dk2*=(dAcharm_dk2*ek2*ek2);
	    dAcharm_dl*=(dAcharm_dl*el*el);
	    dAcharm_dl2*=(dAcharm_dl2*el2*el2);*/

	  /*	  double lowFact=(k-k2)*(k-k2);
	    double lowFact2=(k*k-k*k2)*(k*k-k*k2);
	    double dAcharm_dAcut=(-1)*l/lowFact;
	    double dAcharm_dAall=(double)1/k-l*k2/(k*lowFact);
	    double lowFact2=k*k*l2-k*k2*l;
	    double dAcharm_dk=(-1)*A_All/(k*k)-(l*A_Cut*(k*k*l2-k*k2*l)-(2*k*l2-k2*l)*(l*k*A_Cut-k*k2*A_All))/(lowFact2*lowFact2);
	    double dAcharm_dk2=(l*A_All*lowFact2-k*l*(l*k*A_Cut-l*k2*A_All))/(lowFact2*lowFact2);
	    double dAcharm_dl=(-1)*((k*A_Cut-k2*A_All)*lowFact+k*k2*(l*k*A_Cut-l*k2*A_All))/(lowFact2*lowFact2);
	    double dAcharm_dl2=k*k*(l*k*A_Cut-l*k2*A_All)/(lowFact2*lowFact2);*/


	  /*	  dAuds_dAcut*=(dAuds_dAcut*eA_Cut*eA_Cut);
	    dAuds_dAall*=(dAuds_dAall*eA_All*eA_All);
	    dAuds_dk*=(dAuds_dk*ek*ek);
	    dAuds_dk2*=(dAuds_dk2*ek2*ek2);*/

	  double lowFact=(k-k2)*(k-k2);
	  double lowFact2=(k*k-k*k2)*(k*k-k*k2);
	  double lowFactNoPass=(kNoPass-k2)*(kNoPass-k2);
	  double lowFactNoPass2=(kNoPass*kNoPass-kNoPass*k2)*(kNoPass*kNoPass-kNoPass*k2);

	  ///	  double dAuds_dAcut=k/(double)(k-k2);
	  double dAuds_dAcut=kNoPass/(double)(kNoPass-k2);
	  ///	  double dAuds_dAall=-k2/(double)(k-k2);
	  double dAuds_dARevCut=-k2/(double)(kNoPass-k2);

	  ///	  double dAuds_dk=A_Cut*(k-k2)-(k*A_Cut-k2*A_All)/lowFact;
	  double dAuds_dkNoPass=A_Cut*(kNoPass-k2)-(kNoPass*A_Cut-k2*A_RevCut)/lowFactNoPass;
	  ///	  double dAuds_dk2=((-1)*A_All*(k-k2)+(k*A_Cut-k2*A_All))/lowFact;
	  double dAuds_dk2=((-1)*A_RevCut*(kNoPass-k2)+(kNoPass*A_Cut-k2*A_RevCut))/lowFactNoPass;

	  ///	  double dAcharm_dAcut=(-1)/(k-k2);
	  double dAcharm_dAcut=(-1)/(kNoPass-k2);
	  ///	  double dAcharm_dAall=k2/(k*(k-k2));
	  double dAcharm_dARevCut=k2/(kNoPass*(kNoPass-k2));
	  ///	  double dAcharm_dk=(-1)/(k*k)*A_All-(A_Cut*(k*k-k*k2)-(2*k-k2)*(k*A_Cut-k2*A_All))/(lowFact2);
	  double dAcharm_dkNoPass=(-1)/(kNoPass*kNoPass)*A_RevCut-(A_Cut*(kNoPass*kNoPass-kNoPass*k2)-(2*kNoPass-k2)*(kNoPass*A_Cut-k2*A_All))/(lowFactNoPass2);
	  ///	  double dAcharm_dk2=(-A_All*(k*k-k*k2)+k*k*A_Cut-k*k2*A_All)/lowFact2;
	  double dAcharm_dk2=(-A_RevCut*(kNoPass*kNoPass-kNoPass*k2)+kNoPass*kNoPass*A_Cut-kNoPass*k2*A_RevCut)/lowFactNoPass2;

	  //	  dAuds_dl*=(dAuds_dl*el*el);
	  //	  dAuds_dl2*=(dAuds_dl2*el2*el2);
	  dAuds_dAcut*=(dAuds_dAcut*eA_Cut*eA_Cut);

	  dAuds_dARevCut*=(dAuds_dARevCut*eA_RevCut*eA_RevCut);
	  //	  dAuds_dAall*=(dAuds_dAall*eA_All*eA_All);
	  ///	  dAuds_dk*=(dAuds_dk*ek*ek);
	  dAuds_dkNoPass*=(dAuds_dkNoPass*ekNoPass*ekNoPass);
	  dAuds_dk2*=(dAuds_dk2*ek2*ek2);
	  dAcharm_dAcut*=(dAcharm_dAcut*eA_Cut*eA_Cut);
	  ///	  dAcharm_dAall*=(dAcharm_dAall*eA_All*eA_All);
	  dAcharm_dARevCut*=(dAcharm_dARevCut*eA_RevCut*eA_RevCut);
	  ///	  dAcharm_dk*=(dAcharm_dk*ek*ek);
	  dAcharm_dkNoPass*=(dAcharm_dkNoPass*ekNoPass*ekNoPass);
	  dAcharm_dk2*=(dAcharm_dk2*ek2*ek2);

	  //	  cout <<"M lowfact: " << lowFact << "lowFact2: " << lowFact2 << " d1: " << dAuds_dAcut <<" d2: " << dAuds_dAall <<" d3: " << dAuds_dk;
	  cout <<" d4: " << dAuds_dk2 <<endl;
	  //	  cout <<"charm, cut: " << dAcharm_dAcut << " daall: " << dAcharm_dAall << " dk: " << dAcharm_dk <<" kd2: " << dAcharm_dk2 <<endl;
	  //	  eyUds[iM2]=sqrt(dAuds_dAcut+dAuds_dAall+dAuds_dk+dAuds_dk2+dAuds_dl+dAuds_dl2);
	  //	  eyCharm[iM2]=sqrt(dAcharm_dAcut+dAcharm_dAall+dAcharm_dk+dAcharm_dk2+dAcharm_dl+dAcharm_dl2);
	  ///	  eyUds[iM2]=sqrt(dAuds_dAcut+dAuds_dAall+dAuds_dk+dAuds_dk2);
	  eyUds[iM2]=sqrt(dAuds_dAcut+dAuds_dARevCut+dAuds_dkNoPass+dAuds_dk2);
	  ///	  eyCharm[iM2]=sqrt(dAcharm_dAcut+dAcharm_dAall+dAcharm_dk+dAcharm_dk2);
	  eyCharm[iM2]=sqrt(dAcharm_dAcut+dAcharm_dARevCut+dAcharm_dkNoPass+dAcharm_dk2);
	  oFile <<iM2<<" uds: " << eyUds[iM2]  <<" charm: " << eyCharm[iM2]<<endl;
	  errorSumUdsM+=eyUds[iM2];
	  errorSumCharmM+=eyCharm[iM2];
	  //	  eyCharm[iM2]=sqrt(eFirst*eFirst+eSecond*eSecond+eThird*eThird);
	  //	  eyUds[iM2]=sqSum(eUdsUp/(k*A_Cut-k2*A_All),eLowFact/(k*l2-k2*l2))*yUds[iM2];
	  //	  eyUds[iM2]=(*asymCharmM)[k].second;
	  //	  eyCharm[iM2]=(*asymCharmM)[k].second; //not correct, have to take all the factors into account...
	  //      outFile  <<"Bin " << k << " Asymmetry: " << (*asymAllZ)[k].first << " +- " <<  (*asymAllZ)[k].second <<endl;
	  wUds=1/(eyUds[iM2]*eyUds[iM2]);
	  wCharm=1/(eyCharm[iM2]*eyCharm[iM2]);
	  wAll=1/(eyAll[iM2]*eyAll[iM2]);
	  wCut=1/(eyCut[iM2]*eyCut[iM2]);
	  wRevCut=1/(eyRevCut[iM2]*eyRevCut[iM2]);
	  cout <<"iM2: " << iM2 <<endl;
	  cout <<"wcut: " << wCut <<" wAll: " << wAll <<endl;

	  wUdsSum[iM2]+=wUds;
	  wCharmSum[iM2]+=wCharm;
	  wAllSum[iM2]+=wAll;
	  wCutSum[iM2]+=wCut;
	  wRevCutSum[iM2]+=wRevCut;

	  yUdsS[iM2]+=wUds*yUds[iM2];
	  yCharmS[iM2]+=wCharm*yCharm[iM2];
	  cout <<"yAllS before: " << yCutS[iM2] <<endl;
	  yAllS[iM2]+=wAll*yAll[iM2];
	  cout <<"adding " << wAll<<" * " << yAll[iM2] << " = " << yAllS[iM2]<<endl;
	  cout <<"ycuts before: " << yCutS[iM2] <<endl;
	  yCutS[iM2]+=wCut*yCut[iM2];
	  yRevCutS[iM2]+=wRevCut*yRevCut[iM2];
	  cout <<"adding " << wCut<<" * " << yCut[iM2] << " = " << yCutS[iM2]<<endl;
	  xUdsS[iM2]+=wUds*xUds[iM2];
	  cout <<"adding at M" << iM2 <<" " << wUds << " * " << xUds[iM2] <<endl;
	  xCharmS[iM2]+=wCharm*xCharm[iM2];
	  xAllS[iM2]+=wAll*xAll[iM2];
	  xCutS[iM2]+=wCut*xCut[iM2];
	  xRevCutS[iM2]+=wRevCut*xRevCut[iM2];

	}
      stringstream sstU;
      stringstream sstC;
      sstU <<"Both5_asym_iM"<<iM;;
      TGraphErrors tgU(binningM[iPa1].size(),xUds,yUds,ex,eyUds);
      TGraphErrors tgC(binningM[iPa1].size(),xCharm,yCharm,ex,eyCharm);
      TGraphErrors tgAll(binningM[iPa1].size(),xAll,yAll,ex,eyAll);
      TGraphErrors tgCut(binningM[iPa1].size(),xCut,yCut,ex,eyCut);
      TGraphErrors tgRevCut(binningM[iPa1].size(),xRevCut,yRevCut,ex,eyRevCut);


      tgRevCut.SetMarkerSize(1.2);
      tgRevCut.SetMarkerColor(kPink);
      tgRevCut.SetMarkerStyle(kFullTriangleDown);
      tgC.SetMarkerSize(1.2);
      tgC.SetMarkerColor(kRed);
      tgU.SetMarkerSize(1.2);
      tgU.SetMarkerColor(kBlue);
      tgAll.SetMarkerSize(1.2);
      tgAll.SetMarkerColor(kBlack);
      tgCut.SetMarkerSize(1.2);
      tgCut.SetMarkerColor(kGreen);
      TCanvas cU(sstU.str().c_str(),sstU.str().c_str(),10,50,500,800);
      float gMin=0;
      float gMax=0;
      /*	  if(tgU.GetMaximum()> tgC.GetMaximum())
	gMax=tgU.GetMaximum()+0.01;
	else
	gMax=tgC.GetMaximum()+0.01;
	if(tgU.GetMinimum()< tgC.GetMinimum())
	gMin=tgU.GetMinimum()-0.01;
	else
	gMin=tgC.GetMinimum()-0.01;*/
      gMin=-0.3;
      gMax=0.3;
      ///      tgU.GetYaxis()->SetRangeUser(gMin,gMax);
      ///      tgU.GetXaxis()->SetLimits(0.0,2.5);
      ///      tgU.GetXaxis()->SetRangeUser(0.0,2.5);
      tgRevCut.GetYaxis()->SetRangeUser(gMin,gMax);
      tgRevCut.GetXaxis()->SetRangeUser(0.0,2.0);
      tgRevCut.GetXaxis()->SetLimits(0.0,2.0);
      tgRevCut.Draw("AP*");
      tgU.Draw("P* SAME");
      tgC.Draw("P*SAME");
      tgAll.Draw("P*SAME");
      tgCut.Draw("P*SAME");
      TLegend leg(0.1,0.7,0.48,0.9);
      //	  leg.SetHeader("Legend");
      leg.AddEntry(&tgU,"UDS","p");
      leg.AddEntry(&tgC,"charm","p");
      leg.AddEntry(&tgCut,"Cut","p");
      leg.AddEntry(&tgAll,"All","p");
      leg.AddEntry(&tgRevCut,"Reverse Cut","p");
      leg.Draw();


#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
      cU.SaveAs(("AsUdsCharmOnlyDist/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharmOnlyDist/"+sstU.str()+c_tmp+".eps").c_str());
      cU.SaveAs(("AsUdsCharmOnlyDist/"+sstU.str()+c_tmp+".C").c_str());
#else
      cU.SaveAs(("AsUdsCharmWOKaon/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharmWOKaon/"+sstU.str()+c_tmp+".eps").c_str());
      cU.SaveAs(("AsUdsCharmWOKaon/"+sstU.str()+c_tmp+".C").c_str());
#endif
#else
      cU.SaveAs(("AsUdsCharm/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharm/"+sstU.str()+c_tmp+".eps").c_str());
      cU.SaveAs(("AsUdsCharm/"+sstU.str()+c_tmp+".C").c_str());
#endif
#else
      cU.SaveAs(("AsUdsCharmOnlyMass/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharmOnlyMass/"+sstU.str()+c_tmp+".eps").c_str());
      cU.SaveAs(("AsUdsCharmOnlyMass/"+sstU.str()+c_tmp+".C").c_str());
#endif
    }

  for(int iM=0;iM<numBins;iM++)
    {
      yUdsInt[iM]=yUdsS[iM]/wUdsSum[iM];
      yCharmInt[iM]=yCharmS[iM]/wCharmSum[iM];
      yAllInt[iM]=yAllS[iM]/wAllSum[iM];
      yCutInt[iM]=yCutS[iM]/wCutSum[iM];
      yRevCutInt[iM]=yRevCutS[iM]/wRevCutSum[iM];
      cout <<"setting yCutint: " << yCutInt[iM] <<" = " << yCutS[iM] <<" / " << wCutSum[iM] <<endl;
      cout <<"setting yCutint: " << yAllInt[iM] <<" = " << yAllS[iM] <<" / " << wAllSum[iM] <<endl;
      eyCharmInt[iM]=sqrt(1/wCharmSum[iM]);
      eyUdsInt[iM]=sqrt(1/wUdsSum[iM]);
      eyAllInt[iM]=sqrt(1/wAllSum[iM]);
      eyCutInt[iM]=sqrt(1/wCutSum[iM]);
      eyRevCutInt[iM]=sqrt(1/wRevCutSum[iM]);

      xUdsInt[iM]=xUdsS[iM]/wUdsSum[iM];
      cout <<"wudssum: " << wUdsSum[iM] <<" leading to: " << xUdsInt[iM] <<endl;
      xCharmInt[iM]=xCharmS[iM]/wCharmSum[iM];
      xAllInt[iM]=xAllS[iM]/wAllSum[iM];
      xCutInt[iM]=xCutS[iM]/wCutSum[iM];
      xRevCutInt[iM]=xRevCutS[iM]/wRevCutSum[iM];
    }
  TGraphErrors tgUAllM(binningM[iPa1].size(),xUdsInt,yUdsInt,ex,eyUdsInt);
  TGraphErrors tgCAllM(binningM[iPa1].size(),xCharmInt,yCharmInt,ex,eyCharmInt);
  TGraphErrors tgAllAllM(binningM[iPa1].size(),xAllInt,yAllInt,ex,eyAllInt);
  TGraphErrors tgCutAllM(binningM[iPa1].size(),xCutInt,yCutInt,ex,eyCutInt);
  TGraphErrors tgRevCutAllM(binningM[iPa1].size(),xRevCutInt,yRevCutInt,ex,eyRevCutInt);

  tgRevCutAllM.SetMarkerSize(1.2);
  tgRevCutAllM.SetMarkerColor(kPink);
  tgRevCutAllM.SetMarkerStyle(kFullTriangleDown);
  tgCAllM.SetMarkerSize(1.2);
  tgCAllM.SetMarkerColor(kRed);
  tgUAllM.SetMarkerSize(1.2);
  tgUAllM.SetMarkerColor(kBlue);
  tgAllAllM.SetMarkerSize(1.2);
  tgAllAllM.SetMarkerColor(kBlack);
  tgCutAllM.SetMarkerSize(1.2);
  tgCutAllM.SetMarkerColor(kGreen);
  stringstream sstAM;
  sstAM <<"All_asym_M";
  TCanvas cAllM(sstAM.str().c_str(),sstAM.str().c_str(),10,50,500,800);
  gMin=-0.2;
  gMax=0.3; // to have space for the legend
  ///  tgUAllM.GetYaxis()->SetRangeUser(gMin,gMax);
  tgRevCutAllM.GetYaxis()->SetRangeUser(gMin,gMax);
  tgRevCutAllM.GetXaxis()->SetRangeUser(0.0,1.0);
  tgRevCutAllM.GetXaxis()->SetLimits(0.0,2.5);
  tgRevCutAllM.Draw("AP*");
  tgUAllM.Draw("P* SAME");
  tgCAllM.Draw("P* SAME");
  tgAllAllM.Draw("P* SAME");
  tgCutAllM.Draw("P* SAME");
  TLegend legM(0.1,0.7,0.48,0.9);
  legM.AddEntry(&tgUAllM,"UDS","p");
  legM.AddEntry(&tgCAllM,"charm","p");
  legM.AddEntry(&tgAllAllM,"All","p");
  legM.AddEntry(&tgCutAllM,"Cut","p");
  legM.AddEntry(&tgRevCutAllM,"ReverseCut","p");
  legM.Draw();

#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
  cAllM.SaveAs(("AsUdsCharmOnlyDist/"+sstAM.str()+c_tmp+".png").c_str());
  cAllM.SaveAs(("AsUdsCharmOnlyDist/"+sstAM.str()+c_tmp+".eps").c_str());
#else
  cAllM.SaveAs(("AsUdsCharmWOKaon/"+sstAM.str()+c_tmp+".png").c_str());
  cAllM.SaveAs(("AsUdsCharmWOKaon/"+sstAM.str()+c_tmp+".eps").c_str());
#endif
#else
  cAllM.SaveAs(("AsUdsCharm/"+sstAM.str()+c_tmp+".png").c_str());
  cAllM.SaveAs(("AsUdsCharm/"+sstAM.str()+c_tmp+".eps").c_str());
#endif
#else
  cAllM.SaveAs(("AsUdsCharmOnlyMass/"+sstAM.str()+c_tmp+".png").c_str());
  cAllM.SaveAs(("AsUdsCharmOnlyMass/"+sstAM.str()+c_tmp+".eps").c_str());
#endif
  ///-----------------------------------------

  setZero(xUdsInt,xCharmInt,yUdsInt,exInt,eyUdsInt,eyCharmInt,numBins);
  setZero(wUdsSum,wCharmSum,xUdsS,yUdsS,eyUdsS2,xCharmS,numBins);
  setZero(wAllSum,wCutSum,xAllS,yAllS,xCutS,yCutS,numBins);
  setZero(yCharmS,eyCharmS2,ex,numBins);

  for(int iM=0;iM<numBins;iM++)
    {
      float wUds=0;
      float wCharm=0;
      float wAll=0;
      float wCut=0;

      float xUds[numBins];
      float xCharm[numBins];
      float yUds[numBins];
      float yCharm[numBins];
      float ex[numBins];
      float eyUds[numBins];
      float eyCharm[numBins];
      oFile <<"bin " << iM <<endl;

      float yAll[numBins];
      float yCut[numBins];
      float xAll[numBins];
      float xCut[numBins];
      float eyAll[numBins];
      float eyCut[numBins];

      for(int iT=0;iT<numBins;iT++)
	{
	  float A_All=(*asymAllMTheta)[iM*binningM[iPa1].size()+iT].first;
	  float A_Cut=(*asymCharmMTheta)[iM*binningM[iPa1].size()+iT].first;
	  yAll[iT]=A_All;
	  yCut[iT]=A_Cut;
	  double allInTestSample=passAsCharmMTheta[iM][iT]+noPassAsCharmMTheta[iM][iT]+passAsCharmMThetaUds[iM][iT]+noPassAsCharmMThetaUds[iM][iT];
	  double allPassedAsCharm=passAsCharmMTheta[iM][iT]+passAsCharmMThetaUds[iM][iT];
	  //xvalsNum the right number?
	  double scaleFact=(double)1/allInTestSample; //amount as compared to all
	  double scaleFact2=(double)1/allPassedAsCharm; //amount as compared to all
	  double k=(passAsCharmMTheta[iM][iT]+noPassAsCharmMTheta[iM][iT])*scaleFact; //num charm in all sample is num of events in charm mc (no matter if surv. the cut)
	  double l=(passAsCharmMThetaUds[iM][iT]+noPassAsCharmMThetaUds[iM][iT])*scaleFact;
	  double k2=passAsCharmMTheta[iM][iT]*scaleFact2; //number of events expected in all sample
	  double l2=passAsCharmMThetaUds[iM][iT]*scaleFact2;
	  ex[iT]=0;
	  //confusing
	  xUds[iT]=xValsTheta[mBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iT]/xValsNumTheta[mBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iT];
	  xCharm[iT]=xUds[iT]+0.01;
	  xAll[iT]=xUds[iT]+0.02;
	  xCut[iT]=xUds[iT]+0.03;
	  yUds[iT]=(k*A_Cut-k2*A_All)/(k*l2-k2*l);
	  yCharm[iT]=A_All/k -l*(k*A_Cut-k2*A_All)/(k*(k*l2-k2*l));
	  cout <<"Theta: im: " << iM <<" im2: " << iT << " yuds: " << yUds[iT] << " ycharm: " << yCharm[iT] << " a_all: " << A_All <<" A_Cut: " << A_Cut <<" error all: " << (*asymAllMTheta)[iM*binningM[iPa1].size()+iT].second << " error charm: " <<(*asymCharmMTheta)[iM*binningM[iPa1].size()+iT].second <<endl;
	  double el=sqrt(allInTestSample*(double)l*(1-l))/allInTestSample;
	  double el2=sqrt(allPassedAsCharm*(double)l2*(1-l2))/allPassedAsCharm;
	  double ek=sqrt(allInTestSample*(double)k*(1-k))/allInTestSample;
	  double ek2=sqrt(allPassedAsCharm*(double)k2*(1-k2))/allPassedAsCharm;
	  cout <<"l: " << l <<" l2: " << l2<< " k: " << k <<" k2: " << k2 <<" el: " << el << " el2: " << el2 << " ek: "<< ek << " ek2: " << ek2 <<endl; 
	  double eA_All=(*asymAllMTheta)[iM*binningM[iPa1].size()+iT].second;
	  double eA_Cut=(*asymCharmMTheta)[iM*binningM[iPa1].size()+iT].second;
	  eyAll[iT]=eA_All;
	  eyCut[iT]=eA_Cut;

	  double eFirst=sqSum(eA_All/A_All,ek/k)*(A_All/k);
	  double eLowFact1=sqSum(ek/k,el2/l2)*k*l2;
	  double eLowFact2=sqSum(ek2/k2,el/l)*k2*l;
	  double eLowFact=sqSum(eLowFact1,eLowFact2);
	  double eSecond=sqSum(el/l,eLowFact/(k*l2-k2*l),eA_Cut/A_Cut);; //error on hte factor before A_cut
	  eSecond*=A_Cut*l/(k*l2-k2*l); //no div, because we omitted the mult earlier
	  double eThird=sqSum(ek2/k2,ek/k,el/l,eA_All/A_All,eLowFact/(k*l2-k2*l));
	  eThird*=(l*k2/k)*A_All/(k*l2-k2*l);
	  //for uds:
	  double eUdsUp= sqSum(sqSum(ek/k,eA_Cut/A_Cut)*k*A_Cut,sqSum(ek2/k2,eA_All,A_All)*k2*A_All);
	  cout <<"eall: " << eA_All <<" eaCut: " << eA_Cut <<" first: " << eFirst << ": elowfact1: " << eLowFact1 << " lowfact2: " << eLowFact2 << " eLowFact " << eLowFact <<" esec: " << eSecond  << " third: " << eThird << " eUdsUp : " << eUdsUp <<endl;

	  double lowFact=(k-k2)*(k-k2);
	  double lowFact2=(k*k-k*k2)*(k*k-k*k2);
	  double dAuds_dAcut=k/(double)(k-k2);
	  double dAuds_dAall=-k2/(double)(k-k2);

	  double dAuds_dk=A_Cut*(k-k2)-(k*A_Cut-k2*A_All)/lowFact;
	  double dAuds_dk2=((-1)*A_All*(k-k2)+(k*A_Cut-k2*A_All))/lowFact;

	  double dAcharm_dAcut=(-1)/(k-k2);
	  double dAcharm_dAall=k2/(k*(k-k2));
	  double dAcharm_dk=(-1)/(k*k)*A_All-(A_Cut*(k*k-k*k2)-(2*k-k2)*(k*A_Cut-k2*A_All))/(lowFact2);
	  double dAcharm_dk2=(-A_All*(k*k-k*k2)+k*k*A_Cut-k*k2*A_All)/lowFact2;

	  dAuds_dAcut*=(dAuds_dAcut*eA_Cut*eA_Cut);
	  dAuds_dAall*=(dAuds_dAall*eA_All*eA_All);
	  dAuds_dk*=(dAuds_dk*ek*ek);
	  dAuds_dk2*=(dAuds_dk2*ek2*ek2);
	  dAcharm_dAcut*=(dAcharm_dAcut*eA_Cut*eA_Cut);
	  dAcharm_dAall*=(dAcharm_dAall*eA_All*eA_All);
	  dAcharm_dk*=(dAcharm_dk*ek*ek);
	  dAcharm_dk2*=(dAcharm_dk2*ek2*ek2);
	  cout <<"MTheta lowfact: " << lowFact << "lowFact2: " << lowFact2 << " d1: " << dAuds_dAcut <<" d2: " << dAuds_dAall <<" d3: " << dAuds_dk;
	  cout <<" d4: " << dAuds_dk2 <<endl;
	  cout <<"charm, cut: " << dAcharm_dAcut << " daall: " << dAcharm_dAall << " dk: " << dAcharm_dk <<" kd2: " << dAcharm_dk2 <<endl;
	  eyUds[iT]=sqrt(dAuds_dAcut+dAuds_dAall+dAuds_dk+dAuds_dk2);
	  eyCharm[iT]=sqrt(dAcharm_dAcut+dAcharm_dAall+dAcharm_dk+dAcharm_dk2);
	  oFile <<iT<<" uds: " << eyUds[iT]  <<" charm: " << eyCharm[iT]<<endl;
	  errorSumUdsM+=eyUds[iT];
	  errorSumCharmM+=eyCharm[iT];

	  wUds=1/(eyUds[iT]*eyUds[iT]);
	  wCharm=1/(eyCharm[iT]*eyCharm[iT]);
	  wAll=1/(eyAll[iT]*eyAll[iT]);
	  wCut=1/(eyCut[iT]*eyCut[iT]);

	  wUdsSum[iT]+=wUds;
	  wCharmSum[iT]+=wCharm;
	  wAllSum[iT]+=wAll;
	  wCutSum[iT]+=wCut;

	  yUdsS[iT]+=wUds*yUds[iT];
	  yCharmS[iT]+=wCharm*yCharm[iT];
	  yAllS[iT]+=wAll*yAll[iT];
	  yCutS[iT]+=wCut*yCut[iT];

	  xUdsS[iT]+=wUds*xUds[iT];
	  cout <<"adding at M" << iT <<" " << wUds << " * " << xUds[iT] <<endl;
	  xCharmS[iT]+=wCharm*xCharm[iT];
	  xAllS[iT]+=wAll*xAll[iT];
	  xCutS[iT]+=wCut*xCut[iT];
	}
      stringstream sstU;
      stringstream sstC;
      sstU <<"Both5_asym_iMTheta"<<iM;;
      TGraphErrors tgU(binningTheta[iPa1].size(),xUds,yUds,ex,eyUds);
      TGraphErrors tgC(binningTheta[iPa1].size(),xCharm,yCharm,ex,eyCharm);
      TGraphErrors tgAll(binningTheta[iPa1].size(),xAll,yAll,ex,eyAll);
      TGraphErrors tgCut(binningTheta[iPa1].size(),xCut,yCut,ex,eyCut);

      tgC.SetMarkerSize(1.2);
      tgC.SetMarkerColor(kRed);
      tgU.SetMarkerSize(1.2);
      tgU.SetMarkerColor(kBlue);
      tgAll.SetMarkerSize(1.2);
      tgAll.SetMarkerColor(kBlack);
      tgCut.SetMarkerSize(1.2);
      tgCut.SetMarkerColor(kGreen);
      TCanvas cU(sstU.str().c_str(),sstU.str().c_str(),10,50,800,1200);
      float gMin=0;
      float gMax=0;

      gMin=-0.7;
      gMax=0.7;
      tgU.GetYaxis()->SetRangeUser(gMin,gMax);
      tgU.GetXaxis()->SetLimits(0.0,1.2);
      tgU.GetXaxis()->SetRangeUser(0.0,1.2);
      tgU.Draw("AP*");
      tgC.Draw("P*SAME");
      tgAll.Draw("P*SAME");
      tgCut.Draw("P*SAME");
      TLegend leg(0.1,0.7,0.48,0.9);
      //	  leg.SetHeader("Legend");
      leg.AddEntry(&tgU,"UDS","p");
      leg.AddEntry(&tgC,"charm","p");
      leg.AddEntry(&tgCut,"Cut","p");
      leg.AddEntry(&tgAll,"All","p");
      leg.Draw();


#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
      cU.SaveAs(("AsUdsCharmOnlyDist/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharmOnlyDist/"+sstU.str()+c_tmp+".eps").c_str());
      cU.SaveAs(("AsUdsCharmOnlyDist/"+sstU.str()+c_tmp+".C").c_str());
#else
      cU.SaveAs(("AsUdsCharmWOKaon/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharmWOKaon/"+sstU.str()+c_tmp+".eps").c_str());
      cU.SaveAs(("AsUdsCharmWOKaon/"+sstU.str()+c_tmp+".C").c_str());
#endif
#else
      cU.SaveAs(("AsUdsCharm/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharm/"+sstU.str()+c_tmp+".eps").c_str());
      cU.SaveAs(("AsUdsCharm/"+sstU.str()+c_tmp+".C").c_str());
#endif
#else
      cU.SaveAs(("AsUdsCharmOnlyMass/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharmOnlyMass/"+sstU.str()+c_tmp+".eps").c_str());
      cU.SaveAs(("AsUdsCharmOnlyMass/"+sstU.str()+c_tmp+".C").c_str());
#endif
    }

  for(int iT=0;iT<numBins;iT++)
    {
      yUdsInt[iT]=yUdsS[iT]/wUdsSum[iT];
      yCharmInt[iT]=yCharmS[iT]/wCharmSum[iT];
      yAllInt[iT]=yAllS[iT]/wAllSum[iT];
      yCutInt[iT]=yCutS[iT]/wCutSum[iT];

      eyCharmInt[iT]=sqrt(1/wCharmSum[iT]);
      eyUdsInt[iT]=sqrt(1/wUdsSum[iT]);
      eyAllInt[iT]=sqrt(1/wAllSum[iT]);
      eyCutInt[iT]=sqrt(1/wCutSum[iT]);

      xUdsInt[iT]=xUdsS[iT]/wUdsSum[iT];
      cout <<"wudssum: " << wUdsSum[iT] <<" leading to: " << xUdsInt[iT] <<endl;
      xCharmInt[iT]=xCharmS[iT]/wCharmSum[iT];
      xAllInt[iT]=xAllS[iT]/wAllSum[iT];
      xCutInt[iT]=xCutS[iT]/wCutSum[iT];
      
      isnan(yUdsInt[iT]) ? yUdsInt[iT]=0 : true;
      isnan(xUdsInt[iT]) ? xUdsInt[iT]=0 : true;
      isnan(yCharmInt[iT]) ? yCharmInt[iT]=0 : true;
      isnan(xCharmInt[iT]) ? xCharmInt[iT]=0 : true;
      isnan(yAllInt[iT]) ? yAllInt[iT]=0 : true;
      isnan(xAllInt[iT]) ? xAllInt[iT]=0 : true;
      isnan(yCutInt[iT]) ? yCutInt[iT]=0 : true;
      isnan(xCutInt[iT]) ? xCutInt[iT]=0 : true;
      isnan(eyCutInt[iT]) ? eyCutInt[iT]=0 : true;
      isnan(eyAllInt[iT]) ? eyAllInt[iT]=0 : true;
      isnan(eyCharmInt[iT]) ? eyCharmInt[iT]=0 : true;
      isnan(eyUdsInt[iT]) ? eyUdsInt[iT]=0 : true;
    }
  TGraphErrors tgUAllMTheta(binningTheta[iPa1].size(),xUdsInt,yUdsInt,ex,eyUdsInt);
  TGraphErrors tgCAllMTheta(binningTheta[iPa1].size(),xCharmInt,yCharmInt,ex,eyCharmInt);
  TGraphErrors tgAllAllMTheta(binningTheta[iPa1].size(),xAllInt,yAllInt,ex,eyAllInt);
  TGraphErrors tgCutAllMTheta(binningTheta[iPa1].size(),xCutInt,yCutInt,ex,eyCutInt);
  tgCAllMTheta.SetMarkerSize(1.2);
  tgCAllMTheta.SetMarkerColor(kRed);
  tgUAllMTheta.SetMarkerSize(1.2);
  tgUAllMTheta.SetMarkerColor(kBlue);
  tgAllAllMTheta.SetMarkerSize(1.2);
  tgAllAllMTheta.SetMarkerColor(kBlack);
  tgCutAllMTheta.SetMarkerSize(1.2);
  tgCutAllMTheta.SetMarkerColor(kGreen);
  stringstream sstAMTheta;
  sstAMTheta <<"All_asym_MTheta";
  TCanvas cAllMTheta(sstAMTheta.str().c_str(),sstAMTheta.str().c_str(),10,50,800,1200);
  gMin=-0.7;
  gMax=0.7; // to have space for the legend
  tgUAllMTheta.GetYaxis()->SetRangeUser(gMin,gMax);
  tgUAllMTheta.Draw("AP*");
  tgCAllMTheta.Draw("P* SAME");
  tgAllAllMTheta.Draw("P* SAME");
  tgCutAllMTheta.Draw("P* SAME");
  TLegend legMTheta(0.1,0.7,0.48,0.9);
  legMTheta.AddEntry(&tgUAllMTheta,"UDS","p");
  legMTheta.AddEntry(&tgCAllMTheta,"charm","p");
  legMTheta.AddEntry(&tgAllAllMTheta,"All","p");
  legMTheta.AddEntry(&tgCutAllMTheta,"Cut","p");
  legMTheta.Draw();

#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
  cAllMTheta.SaveAs(("AsUdsCharmOnlyDist/"+sstAMTheta.str()+c_tmp+".png").c_str());
  cAllMTheta.SaveAs(("AsUdsCharmOnlyDist/"+sstAMTheta.str()+c_tmp+".eps").c_str());
#else
  cAllMTheta.SaveAs(("AsUdsCharmWOKaon/"+sstAMTheta.str()+c_tmp+".png").c_str());
  cAllMTheta.SaveAs(("AsUdsCharmWOKaon/"+sstAMTheta.str()+c_tmp+".eps").c_str());
#endif
#else
  cAllMTheta.SaveAs(("AsUdsCharm/"+sstAMTheta.str()+c_tmp+".png").c_str());
  cAllMTheta.SaveAs(("AsUdsCharm/"+sstAMTheta.str()+c_tmp+".eps").c_str());
#endif
#else
  cAllMTheta.SaveAs(("AsUdsCharmOnlyMass/"+sstAMTheta.str()+c_tmp+".png").c_str());
  cAllMTheta.SaveAs(("AsUdsCharmOnlyMass/"+sstAMTheta.str()+c_tmp+".eps").c_str());
#endif
  ///==========================================================================
  setZero(xUdsInt,xCharmInt,yUdsInt,exInt,eyUdsInt,eyCharmInt,numBins);
  setZero(wUdsSum,wCharmSum,xUdsS,yUdsS,eyUdsS2,xCharmS,numBins);
  setZero(wAllSum,wCutSum,xAllS,yAllS,xCutS,yCutS,numBins);
  setZero(yCharmS,eyCharmS2,ex,numBins);

  for(int iZ=0;iZ<numBins;iZ++)
    {
      float wUds=0;
      float wCharm=0;
      float wAll=0;
      float wCut=0;

      float xUds[numBins];
      float xCharm[numBins];
      float yUds[numBins];
      float yCharm[numBins];
      float ex[numBins];
      float eyUds[numBins];
      float eyCharm[numBins];
      oFile <<"bin " << iZ <<endl;

      float yAll[numBins];
      float yCut[numBins];
      float xAll[numBins];
      float xCut[numBins];
      float eyAll[numBins];
      float eyCut[numBins];

      for(int iT=0;iT<numBins;iT++)
	{
	  float A_All=(*asymAllZTheta)[iZ*binningZ[iPa1].size()+iT].first;
	  float A_Cut=(*asymCharmZTheta)[iZ*binningZ[iPa1].size()+iT].first;
	  yAll[iT]=A_All;
	  yCut[iT]=A_Cut;
	  double allInTestSample=passAsCharmZTheta[iZ][iT]+noPassAsCharmZTheta[iZ][iT]+passAsCharmZThetaUds[iZ][iT]+noPassAsCharmZThetaUds[iZ][iT];
	  double allPassedAsCharm=passAsCharmZTheta[iZ][iT]+passAsCharmZThetaUds[iZ][iT];
	  //xvalsNum the right number?
	  double scaleFact=(double)1/allInTestSample; //amount as compared to all
	  double scaleFact2=(double)1/allPassedAsCharm; //amount as compared to all
	  double k=(passAsCharmZTheta[iZ][iT]+noPassAsCharmZTheta[iZ][iT])*scaleFact; //num charm in all sample is num of events in charm mc (no matter if surv. the cut)
	  double l=(passAsCharmZThetaUds[iZ][iT]+noPassAsCharmZThetaUds[iZ][iT])*scaleFact;
	  double k2=passAsCharmZTheta[iZ][iT]*scaleFact2; //number of events expected in all sample
	  double l2=passAsCharmZThetaUds[iZ][iT]*scaleFact2;
	  ex[iT]=0;
	  //confusing
	  xUds[iT]=xValsTheta[zBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iT]/xValsNumTheta[zBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iT];
	  xCharm[iT]=xUds[iT]+0.01;
	  xAll[iT]=xUds[iT]+0.02;
	  xCut[iT]=xUds[iT]+0.03;
	  yUds[iT]=(k*A_Cut-k2*A_All)/(k*l2-k2*l);

	  yCharm[iT]=A_All/k -l*(k*A_Cut-k2*A_All)/(k*(k*l2-k2*l));
	  cout << " for it: " << iT <<" k: " << k << " l: " << l <<" k2: " << k2 << " l2: " << l2 <<endl;
	  cout <<"im: " << iZ <<" im2: " << iT << " yuds: " << yUds[iT] << " ycharm: " << yCharm[iT] << " a_all: " << A_All <<" A_Cut: " << A_Cut <<" error all: " << (*asymAllZTheta)[iZ*binningZ[iPa1].size()+iT].second << " error charm: " <<(*asymCharmZTheta)[iZ*binningZ[iPa1].size()+iT].second <<endl;
	  double el=sqrt(allInTestSample*(double)l*(1-l))/allInTestSample;
	  double el2=sqrt(allPassedAsCharm*(double)l2*(1-l2))/allPassedAsCharm;
	  double ek=sqrt(allInTestSample*(double)k*(1-k))/allInTestSample;
	  double ek2=sqrt(allPassedAsCharm*(double)k2*(1-k2))/allPassedAsCharm;
	  cout <<"l: " << l <<" l2: " << l2<< " k: " << k <<" k2: " << k2 <<" el: " << el << " el2: " << el2 << " ek: "<< ek << " ek2: " << ek2 <<endl; 
	  double eA_All=(*asymAllZTheta)[iZ*binningZ[iPa1].size()+iT].second;
	  double eA_Cut=(*asymCharmZTheta)[iZ*binningZ[iPa1].size()+iT].second;
	  eyAll[iT]=eA_All;
	  eyCut[iT]=eA_Cut;

	  double eFirst=sqSum(eA_All/A_All,ek/k)*(A_All/k);
	  double eLowFact1=sqSum(ek/k,el2/l2)*k*l2;
	  double eLowFact2=sqSum(ek2/k2,el/l)*k2*l;
	  double eLowFact=sqSum(eLowFact1,eLowFact2);
	  double eSecond=sqSum(el/l,eLowFact/(k*l2-k2*l),eA_Cut/A_Cut);; //error on hte factor before A_cut
	  eSecond*=A_Cut*l/(k*l2-k2*l); //no div, because we omitted the mult earlier
	  double eThird=sqSum(ek2/k2,ek/k,el/l,eA_All/A_All,eLowFact/(k*l2-k2*l));
	  eThird*=(l*k2/k)*A_All/(k*l2-k2*l);
	  //for uds:
	  double eUdsUp= sqSum(sqSum(ek/k,eA_Cut/A_Cut)*k*A_Cut,sqSum(ek2/k2,eA_All,A_All)*k2*A_All);
	  cout <<"eall: " << eA_All <<" eaCut: " << eA_Cut <<" first: " << eFirst << ": elowfact1: " << eLowFact1 << " lowfact2: " << eLowFact2 << " eLowFact " << eLowFact <<" esec: " << eSecond  << " third: " << eThird << " eUdsUp : " << eUdsUp <<endl;

	  double lowFact=(k-k2)*(k-k2);
	  double lowFact2=(k*k-k*k2)*(k*k-k*k2);
	  double dAuds_dAcut=k/(double)(k-k2);
	  double dAuds_dAall=-k2/(double)(k-k2);

	  double dAuds_dk=A_Cut*(k-k2)-(k*A_Cut-k2*A_All)/lowFact;
	  double dAuds_dk2=((-1)*A_All*(k-k2)+(k*A_Cut-k2*A_All))/lowFact;

	  double dAcharm_dAcut=(-1)/(k-k2);
	  double dAcharm_dAall=k2/(k*(k-k2));
	  double dAcharm_dk=(-1)/(k*k)*A_All-(A_Cut*(k*k-k*k2)-(2*k-k2)*(k*A_Cut-k2*A_All))/(lowFact2);
	  double dAcharm_dk2=(-A_All*(k*k-k*k2)+k*k*A_Cut-k*k2*A_All)/lowFact2;


	  dAuds_dAcut*=(dAuds_dAcut*eA_Cut*eA_Cut);
	  dAuds_dAall*=(dAuds_dAall*eA_All*eA_All);
	  dAuds_dk*=(dAuds_dk*ek*ek);
	  dAuds_dk2*=(dAuds_dk2*ek2*ek2);
	  dAcharm_dAcut*=(dAcharm_dAcut*eA_Cut*eA_Cut);
	  dAcharm_dAall*=(dAcharm_dAall*eA_All*eA_All);
	  dAcharm_dk*=(dAcharm_dk*ek*ek);
	  dAcharm_dk2*=(dAcharm_dk2*ek2*ek2);
	  cout <<"ZTheta lowfact: " << lowFact << "lowFact2: " << lowFact2 << " d1: " << dAuds_dAcut <<" d2: " << dAuds_dAall <<" d3: " << dAuds_dk;
	  cout <<" d4: " << dAuds_dk2 <<endl;
	  cout <<"charm, cut: " << dAcharm_dAcut << " daall: " << dAcharm_dAall << " dk: " << dAcharm_dk <<" kd2: " << dAcharm_dk2 <<endl;
	  eyUds[iT]=sqrt(dAuds_dAcut+dAuds_dAall+dAuds_dk+dAuds_dk2);
	  eyCharm[iT]=sqrt(dAcharm_dAcut+dAcharm_dAall+dAcharm_dk+dAcharm_dk2);
	  oFile <<iT<<" uds: " << eyUds[iT]  <<" charm: " << eyCharm[iT]<<endl;
	  errorSumUdsZ+=eyUds[iT];
	  errorSumCharmZ+=eyCharm[iT];

	  wUds=1/(eyUds[iT]*eyUds[iT]);
	  wCharm=1/(eyCharm[iT]*eyCharm[iT]);
	  wAll=1/(eyAll[iT]*eyAll[iT]);
	  wCut=1/(eyCut[iT]*eyCut[iT]);

	  wUdsSum[iT]+=wUds;
	  wCharmSum[iT]+=wCharm;
	  wAllSum[iT]+=wAll;
	  wCutSum[iT]+=wCut;

	  yUdsS[iT]+=wUds*yUds[iT];
	  yCharmS[iT]+=wCharm*yCharm[iT];
	  yAllS[iT]+=wAll*yAll[iT];
	  yCutS[iT]+=wCut*yCut[iT];

	  xUdsS[iT]+=wUds*xUds[iT];
	  cout <<"adding at Z" << iT <<" " << wUds << " * " << xUds[iT] <<endl;
	  xCharmS[iT]+=wCharm*xCharm[iT];
	  xAllS[iT]+=wAll*xAll[iT];
	  xCutS[iT]+=wCut*xCut[iT];
	}
      stringstream sstU;
      stringstream sstC;
      sstU <<"Both5_asym_iZTheta"<<iZ;;
      TGraphErrors tgU(binningTheta[iPa1].size(),xUds,yUds,ex,eyUds);
      TGraphErrors tgC(binningTheta[iPa1].size(),xCharm,yCharm,ex,eyCharm);
      TGraphErrors tgAll(binningTheta[iPa1].size(),xAll,yAll,ex,eyAll);
      TGraphErrors tgCut(binningTheta[iPa1].size(),xCut,yCut,ex,eyCut);

      tgC.SetMarkerSize(1.2);
      tgC.SetMarkerColor(kRed);
      tgU.SetMarkerSize(1.2);
      tgU.SetMarkerColor(kBlue);
      tgAll.SetMarkerSize(1.2);
      tgAll.SetMarkerColor(kBlack);
      tgCut.SetMarkerSize(1.2);
      tgCut.SetMarkerColor(kGreen);
      TCanvas cU(sstU.str().c_str(),sstU.str().c_str(),10,50,800,1200);
      float gMin=0;
      float gMax=0;

      gMin=-0.7;
      gMax=0.7;
      tgU.GetYaxis()->SetRangeUser(gMin,gMax);
      tgU.GetXaxis()->SetLimits(0.0,1.2);
      tgU.GetXaxis()->SetRangeUser(0.0,1.2);
      tgU.Draw("AP*");
      tgC.Draw("P*SAME");
      tgAll.Draw("P*SAME");
      tgCut.Draw("P*SAME");
      TLegend leg(0.1,0.7,0.48,0.9);
      //	  leg.SetHeader("Legend");
      leg.AddEntry(&tgU,"UDS","p");
      leg.AddEntry(&tgC,"charm","p");
      leg.AddEntry(&tgCut,"Cut","p");
      leg.AddEntry(&tgAll,"All","p");
      leg.Draw();


#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
      cU.SaveAs(("AsUdsCharmOnlyDist/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharmOnlyDist/"+sstU.str()+c_tmp+".eps").c_str());
      cU.SaveAs(("AsUdsCharmOnlyDist/"+sstU.str()+c_tmp+".C").c_str());
#else
      cU.SaveAs(("AsUdsCharmWOKaon/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharmWOKaon/"+sstU.str()+c_tmp+".eps").c_str());
      cU.SaveAs(("AsUdsCharmWOKaon/"+sstU.str()+c_tmp+".C").c_str());
#endif
#else
      cU.SaveAs(("AsUdsCharm/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharm/"+sstU.str()+c_tmp+".eps").c_str());
      cU.SaveAs(("AsUdsCharm/"+sstU.str()+c_tmp+".C").c_str());
#endif
#else
      cU.SaveAs(("AsUdsCharmOnlyMass/"+sstU.str()+c_tmp+".png").c_str());
      cU.SaveAs(("AsUdsCharmOnlyMass/"+sstU.str()+c_tmp+".eps").c_str());
      cU.SaveAs(("AsUdsCharmOnlyMass/"+sstU.str()+c_tmp+".C").c_str());
#endif
    }

  for(int iT=0;iT<numBins;iT++)
    {
      yUdsInt[iT]=yUdsS[iT]/wUdsSum[iT];
      yCharmInt[iT]=yCharmS[iT]/wCharmSum[iT];
      yAllInt[iT]=yAllS[iT]/wAllSum[iT];
      yCutInt[iT]=yCutS[iT]/wCutSum[iT];

      eyCharmInt[iT]=sqrt(1/wCharmSum[iT]);
      eyUdsInt[iT]=sqrt(1/wUdsSum[iT]);
      eyAllInt[iT]=sqrt(1/wAllSum[iT]);
      eyCutInt[iT]=sqrt(1/wCutSum[iT]);

      xUdsInt[iT]=xUdsS[iT]/wUdsSum[iT];
      cout <<"wudssum: " << wUdsSum[iT] <<" leading to: " << xUdsInt[iT] <<endl;
      xCharmInt[iT]=xCharmS[iT]/wCharmSum[iT];
      xAllInt[iT]=xAllS[iT]/wAllSum[iT];
      xCutInt[iT]=xCutS[iT]/wCutSum[iT];

      isnan(yUdsInt[iT]) ? yUdsInt[iT]=0 : true;
      isnan(xUdsInt[iT]) ? xUdsInt[iT]=0 : true;
      isnan(yCharmInt[iT]) ? yCharmInt[iT]=0 : true;
      isnan(xCharmInt[iT]) ? xCharmInt[iT]=0 : true;
      isnan(yAllInt[iT]) ? yAllInt[iT]=0 : true;
      isnan(xAllInt[iT]) ? xAllInt[iT]=0 : true;
      isnan(yCutInt[iT]) ? yCutInt[iT]=0 : true;
      isnan(xCutInt[iT]) ? xCutInt[iT]=0 : true;
      isnan(eyCutInt[iT]) ? eyCutInt[iT]=0 : true;
      isnan(eyAllInt[iT]) ? eyAllInt[iT]=0 : true;
      isnan(eyCharmInt[iT]) ? eyCharmInt[iT]=0 : true;
      isnan(eyUdsInt[iT]) ? eyUdsInt[iT]=0 : true;
    }
  TGraphErrors tgUAllZTheta(binningTheta[iPa1].size(),xUdsInt,yUdsInt,ex,eyUdsInt);
  TGraphErrors tgCAllZTheta(binningTheta[iPa1].size(),xCharmInt,yCharmInt,ex,eyCharmInt);
  TGraphErrors tgAllAllZTheta(binningTheta[iPa1].size(),xAllInt,yAllInt,ex,eyAllInt);
  TGraphErrors tgCutAllZTheta(binningTheta[iPa1].size(),xCutInt,yCutInt,ex,eyCutInt);
  tgCAllZTheta.SetMarkerSize(1.2);
  tgCAllZTheta.SetMarkerColor(kRed);
  tgUAllZTheta.SetMarkerSize(1.2);
  tgUAllZTheta.SetMarkerColor(kBlue);
  tgAllAllZTheta.SetMarkerSize(1.2);
  tgAllAllZTheta.SetMarkerColor(kBlack);
  tgCutAllZTheta.SetMarkerSize(1.2);
  tgCutAllZTheta.SetMarkerColor(kGreen);
  stringstream sstAZTheta;
  sstAZTheta <<"All_asym_ZTheta";
  TCanvas cAllZTheta(sstAZTheta.str().c_str(),sstAZTheta.str().c_str(),10,50,800,1200);
  gMin=-0.7;
  gMax=0.7; // to have space for the legend
  tgUAllZTheta.GetYaxis()->SetRangeUser(gMin,gMax);
  tgUAllZTheta.Draw("AP*");
  tgCAllZTheta.Draw("P* SAME");
  tgAllAllZTheta.Draw("P* SAME");
  tgCutAllZTheta.Draw("P* SAME");
  TLegend legZTheta(0.1,0.7,0.48,0.9);
  legZTheta.AddEntry(&tgUAllZTheta,"UDS","p");
  legZTheta.AddEntry(&tgCAllZTheta,"charm","p");
  legZTheta.AddEntry(&tgAllAllZTheta,"All","p");
  legZTheta.AddEntry(&tgCutAllZTheta,"Cut","p");
  legZTheta.Draw();

#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
  cAllZTheta.SaveAs(("AsUdsCharmOnlyDist/"+sstAZTheta.str()+c_tmp+".png").c_str());
  cAllZTheta.SaveAs(("AsUdsCharmOnlyDist/"+sstAZTheta.str()+c_tmp+".eps").c_str());
#else
  cAllZTheta.SaveAs(("AsUdsCharmWOKaon/"+sstAZTheta.str()+c_tmp+".png").c_str());
  cAllZTheta.SaveAs(("AsUdsCharmWOKaon/"+sstAZTheta.str()+c_tmp+".eps").c_str());
#endif
#else
  cAllZTheta.SaveAs(("AsUdsCharm/"+sstAZTheta.str()+c_tmp+".png").c_str());
  cAllZTheta.SaveAs(("AsUdsCharm/"+sstAZTheta.str()+c_tmp+".eps").c_str());
#endif
#else
  cAllZTheta.SaveAs(("AsUdsCharmOnlyMass/"+sstAZTheta.str()+c_tmp+".png").c_str());
  cAllZTheta.SaveAs(("AsUdsCharmOnlyMass/"+sstAZTheta.str()+c_tmp+".eps").c_str());
#endif

  /////============================================


  long countAgain=0;
  for(int iM=0;iM<numBins;iM++)
    { 
      for(int iM2=0;iM2<numBins;iM2++)
	{ 
	  int countLoc=0;
	  for(int i=0;i <xValsNum[mBinning][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];i++)
	    {
	      countLoc++;
	      countAgain++;
	    }
	  cout <<"for iM: " << iM <<" iM2: " << iM2 << " counts: " << countLoc <<endl;
	}
    }
  cout <<"count Hads: " << allHads <<" and again: " << countAgain <<"pipPN: " << allHadsPiPiPN << "  thrust  " << allHadsThrust << " opening; " << allHadsOpening << " highCuts; " << allHadsHighCuts << " before fill: " << allHadsBeforeFill <<" fillM " << allHadsBeforeFillM<<endl;

  cout <<"numCharm: " << numCharm <<" numUds:  " << numUds <<endl;
  cout <<"validmva: " << mvaCounter <<" nonvalid: " <<nonValidMvaCounter<<endl;
  cout <<"cutted ev: " << cuttedEv <<endl;
  cout <<"accEv: " << accEv <<endl;
  oFile <<" errorUdsSumZ: " << errorSumUdsZ <<" errorSumUdsM"  << errorSumUdsM <<" errorSumCharmZ: " << errorSumCharmZ <<" errorSumCharmM: " << errorSumCharmM << endl;
  oFile << " errorAllZ: " << errorSumUdsZ+errorSumCharmZ <<" errorAllM" << errorSumUdsM+errorSumCharmM << " errorAllUds: " << errorSumUdsZ+errorSumUdsM <<" errorAllCharm: " << errorSumCharmZ+errorSumCharmM <<endl;
  oFile <<flush;
  oFile.close();

}
