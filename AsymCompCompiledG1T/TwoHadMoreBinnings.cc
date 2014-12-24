//#define MC
//#define PAW_STYLE
//#define WITH_PURITIES
#define WITH_KIN_CORR


#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
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
using namespace std;
ofstream* errorVglFile;
vector<float>* binningM;
vector<float>* binningZ;
vector<float>* binningZSpect;
vector<float>* decThetaBinning;
vector<float>* dilBinning;

vector<float>* binningMReduced;
vector<float>* binningZReduced;


//to compare errors from the fits with the naive statistical errors sqrt(N)/N
#include "TwoHadAsymsCommons.h"
//#define MAX_EVENTS 10000

#define AMP_PHIR 0.1
#define OPENING_CUT 0.8

#include "Instances.h"

//dil factor: A(y) (1/2 -y -y^2) = (CMS) 1/4 (1+cos^2(\theta))
//B(y) = y(1-y) =(CMS) 1/4(sin^2\theta)
//dilution fact: b/a=sin^2\theta/(1+cos^2\theta)

//modulation of the counts with the mc angle
float getInc(float phiR, double mass1, double mass2,  double z1, double z2,int pt1, int pt2, int charge1,int charge2)
{
  //  cout <<"getinc: " << 1+AMP_PHIR*cos(phiR)<<endl;
  return 1+AMP_PHIR*cos(phiR);
}

double******* kinSpectraPars;
double****** kinSpectraReducedPars;
double****** kinSpectra;
double***** kinSpectraReduced;
TH2D***** kinSpectraH;
TH1D***** kinSpectra1H;

float getInc2(int binning, int index1, int index2, int mBin1, int mBin2, float zVal1, float zVal2)
{
  //  return 1;//test
  if(mBin1<0 || mBin2 <0 || zVal1 < 0 || zVal2 <0)
    {
      cout <<"ret 0 "<< mBin1 <<" b2: " << mBin2 <<" zval1: " << zVal1 <<" val2: " << zVal2<<endl;
      return 0;
    }
  //  cout << "spect: " <<kinSpectra[index1][index2][zBin1][zBin2][mBin1][mBin2] << endl;
  //  cout <<"getinc: " << 1+AMP_PHIR*cos(phiR)<<endl;
  //  cout << "returning: " <<(float)1/kinSpectra[index1][index2][zBin1][zBin2][mBin1][mBin2] <<endl;
  double val=0;
  switch(binning)
    {
    case zBinning:
      break;

    case mBinning:
      break;
    }
  //  cout <<"ind1: " << index1 <<", ind2: " << index2 <<" mbin1: " << mBin1 <<" m2: " << mBin2 <<endl;
  //  cout <<"pointer is : " << kinSpectra1H[index1][index2][mBin1][mBin2] <<endl;
  //  printf("or: %p \n",kinSpectra1H[index1][index2][mBin1][mBin2]);
  float normFactor=(float)(kinSpectra1H[index1][index2][mBin1][mBin2]->GetEntries());

  //  cout <<"entries; " << normFactor <<endl;
  //normalize to first spectrum
  double nf2=kinSpectra1H[0][0][0][0]->FindBin(zVal1)*kinSpectra1H[0][0][0][0]->FindBin(zVal2);
  normFactor*=nf2;
  normFactor/=kinSpectra1H[index1][index2][mBin1][mBin2]->GetNbinsX();

  int bin1=kinSpectra1H[index1][index2][mBin1][mBin2]->FindBin(zVal1);
  int bin2=kinSpectra1H[index1][index2][mBin1][mBin2]->FindBin(zVal2);
  float binVal1=kinSpectra1H[index1][index2][mBin1][mBin2]->GetBinContent(bin1);
  float binVal2=kinSpectra1H[index1][index2][mBin1][mBin2]->GetBinContent(bin2);
  if(binVal1==0|| binVal2==0)
    return 2;
  //  cout <<"m1: " << mBin1 <<" m2: " <<mBin2 <<"ret norm: " << normFactor <<" fact1: " << binVal1 << ", fact2: " << binVal2 << " bin1: " << bin1 <<" bin2: " << bin2 << " for z1: " << zVal1 <<" 2: " << zVal2 << " fact: " <<normFactor*normFactor/(binVal1*binVal2)<<endl
  ;
  if(normFactor*normFactor/(binVal1*binVal2)>10)
    return 10;
  return normFactor*normFactor/(binVal1*binVal2);
  //  return (float)1/kinSpectra[index1][index2][zBin1][zBin2][mBin1][mBin2];
}

//variables for the purity computation
namespace purities
{
  int kinBin1, kinBin2, binning; //kinematic bin of the computed purity
  int chargeType1;    //charge & particle type for the purity
  int particleType1;
  int chargeType2;    //charge & particle type for the purity
  int particleType2;
  float purity;
  float purityError;
}

namespace asymmetries
{
  int binning;
  int kinBin1,kinBin2, kinBin3,kinBin4;
  int chargeType1;
  int particleType1;
  int chargeType2;
  int particleType2;
  float asymmetry;
  float asError;
  float meanKinVal1;
  float meanKinVal2;
  float meanKinVal3;
  float meanKinVal4;
  float meanTheta;
  float chi2Fit;
  int ndfFit;
  int expNr;
  int isOnRes;
}


int main(int argc, char** argv)
{
  long cuttedEv=0;
  long accEv=0;
  long nanEv=0;
  vector<TH1D*> massHistos;
  TH1D* massPiPi_PN=new TH1D("massPiPi_PN","Invariant Mass #pi^{+}/#pi^{-} pairs",1000,0,3);
  TH1D* massPiPi_PZ=new TH1D("massPiPi_PZ","Invariant Mass #pi^{+}/#pi^{0} pairs",1000,0,3);
  TH1D* massPiPi_ZN=new TH1D("massPiPi_ZN","Invariant Mass #pi^{0}/#pi^{-} pairs",1000,0,3);
  TH1D* massPiK_PN=new TH1D("massPiK_PN","Invariant Mass #pi^{+}/K^{-} pairs",1000,0,3);

  TH1D hPhiR("cosPhi12","cosPhi12",16,-1,1);
  TH1D hPhiZ("cosPhi20","cosPhi20",16,-1,1);

  TGraphErrors* tgPurPiPi_PN_Z;
  TGraphErrors* tgPurPiPi_PN_M;

  TGraphErrors* tgPurPiPi_PZ_Z;
  TGraphErrors* tgPurPiPi_PZ_M;

  TGraphErrors* tgPurPiPi_ZN_Z;
  TGraphErrors* tgPurPiPi_ZN_M;

  TGraphErrors* tgPurPiK_PN_Z;
  TGraphErrors* tgPurPiK_PN_M;

  gStyle->SetOptFit(1);
  gROOT->SetStyle("Plain");
  ofstream outputFile("outputFile");
  errorVglFile=new ofstream("errorVgl");
  //continuum luminosities
  int numLumi=24;
  float luminositiesCont[]={0.594,0.,1.211,1.203,1.402,0.853,3.562,0.,1.416,1.671,3.745,2.393,2.722,
			    1.944,6.078,6.315,5.657,6.524,2.315,3.438,2.586,4.825,0.,7.821};//luminosities for exp 7-55 continuum, 
  float luminositiesOnRes[]={5.928,4.44,8.132,10.739,12.682,11.181,24.953,4.375,6.266,25.741,25.427,17.827,17.619,16.733,61.658,43.639,59.989,56.989,13.048,
			     37.577,27.293,38.935,0.0,73.514};
  float ZCut_for_mPlot=0.2;

  float errFact=0;
  for(int i=0;i<numLumi-1;i++)
    {
      errFact+=luminositiesCont[i];
      //      errFact+=luminositiesOnRes[i];
    }
  outputFile <<"int Lumi: " << errFact <<endl;
  errFact=sqrt(errFact/(luminositiesCont[numLumi-1]+luminositiesOnRes[numLumi-1]));//this is the factor that compares error from exp 55 to the exp error for all experiments
  cout <<"errFact: " << errFact <<endl;

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
  if(folderName.find("on_resonance")!=string::npos)
    onResonance=true;

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
  ofstream outFile("Asymmetries");
#ifdef PAW_STYLE
  //for paw style....
  TChain* chI=new TChain("h7");//stupid paw...
  TChain* chF=new TChain("h8");
#else
  cout <<"no paw style" <<endl;
  TChain* chAll=new TChain("DataTree");
#endif
  //why is the z distribution so weird??
  TH1D hZ1("hZ1","hZ1",200,0,5);
  TH1D hZ2("hZ2","hZ2",200,0,5);

  TH1D hAsDiff("AsDiff","AsDiff",200,-10,10);

  TH1D hThetaDiff11("ThetaDiff11","ThetaDiff11",200,-4,4);
  TH1D hThetaDiff12("ThetaDiff12","ThetaDiff12",200,-4,4);
  TH1D hThetaDiff21("ThetaDiff21","ThetaDiff21",200,-4,4);
  TH1D hThetaDiff22("ThetaDiff22","ThetaDiff22",200,-4,4);

  TH1D hThrustThetaRes("Thrust_Rec-MC","Thrust_Rec-MC",200,-3,3);


  TH1D hR("R","R",200,-4,4);

  //  const float zCut=1.2;//is also the max for the z binning. Rought estimate, should follow from resolution of z
  const float zCut=1.0;//is also the max for the z binning. Rought estimate, should follow from resolution of z
  const float mCut=1000;
  const float thetaCut=1000;

  vector<string> fieldNamesF;
  vector<string> fieldNamesI;

  binningM=new vector<float>[NumParticle];
  binningZ=new vector<float>[NumParticle];
  binningMReduced=new vector<float>[NumParticle];
  binningZReduced=new vector<float>[NumParticle];
  decThetaBinning=new vector<float>[NumParticle];
  dilBinning=new vector<float>[NumParticle];

  decThetaBinning[PiPi].push_back(0.1);
  decThetaBinning[PiPi].push_back(10);
  dilBinning[PiPi].push_back(0.1);
  dilBinning[PiPi].push_back(10);

  loadBinning(binningM, binningZ, binningMReduced, binningZReduced);
  int NumSingleBins=2;
  int NumFourBins=1;
  vector<float>** singleBinning=new vector<float>*[NumSingleBins];
  singleBinning[0]=decThetaBinning;
  singleBinning[1]=dilBinning;

  vector<float>*** fourBinning=new vector<float>**[NumFourBins];
  for(int i=0;i<NumFourBins;i++)
    {
      fourBinning[i]=new vector<float>*[2];
      fourBinning[i][0]=binningMReduced;
      fourBinning[i][1]=binningZReduced;
    }


  binningZSpect=new vector<float>[NumParticle];

  float* mLowCuts=new float[NumParticle];
  float* mHighCuts=new float[NumParticle];

  float* zLowCuts=new float[NumParticle];
  float* zHighCuts=new float[NumParticle];

  vector<float> binningAng;

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


  binningZSpect[PiPi].push_back(0.25);
  binningZSpect[PiPi].push_back(0.3);
  binningZSpect[PiPi].push_back(0.35);
  binningZSpect[PiPi].push_back(0.4);
  binningZSpect[PiPi].push_back(0.45);
  binningZSpect[PiPi].push_back(0.55);
  binningZSpect[PiPi].push_back(0.65);
  binningZSpect[PiPi].push_back(0.75);
  binningZSpect[PiPi].push_back(1.0);

  cpyVect(binningZSpect[PiK],binningZ[PiK]);
  cpyVect(binningZSpect[KPi],binningZ[KPi]);
  cpyVect(binningZSpect[KK],binningZ[KK]);
  cpyVect(binningZSpect[UNKNOWN],binningZ[UNKNOWN]);

  int numAngBins=16;
  //see promo thesis, this is the factor for 16 bin, one-D fit. One has to divide by this factor  e.g. a_fit/a = corrfact
  float finiteBinCorrFact=0.9936;
  //  float finiteBinCorrFact=1;

  double pi=3.14159265;
  for(int i=0;i<numAngBins;i++)
    {
      binningAng.push_back((i+1)*2*pi/((float)numAngBins)-pi);
      cout <<"ang binning: " << binningAng[i]<<endl;
    }

  int maxKin=-1;
  int maxKinSpect=binningZSpect[PiPi].size();
  //pipi has probably the most
  if(binningM[PiPi].size()>binningZ[PiPi].size())
    {
      maxKin=binningM[PiPi].size()+3;
    }
  else
    {
      maxKin=binningZ[PiPi].size()+3;
    }
  //the +3 for safety
  vector<float*> memLocsF;
  vector<int*> memLocsI;

  int z1Counter;
  float z1[2000];
  int z2Counter;
  float z2[2000];
  int mass1Counter;
  float mass1[2000];
  int mass2Counter;

  int z1_mcCounter;
  float z1_mc[2000];
  int z2_mcCounter;
  float z2_mc[2000];
  int mass1_mcCounter;
  float mass1_mc[2000];
  int mass2_mcCounter;
  float mass2_mc[2000];

  float mass2[2000];
  int phiRCounter;
  float phiRSum[2000];

  int phiR1Counter;
  float phiR1[2000];

  int twoPhiZeroCounter;
  float twoPhiZero[2000];

  int phiZero1Counter;
  float phiZero1[2000];

  float theta11[2000];
  int theta11Counter;
  float theta12[2000];
  int theta12Counter;
  float theta21[2000];
  int theta21Counter;
  float theta22[2000];
  int theta22Counter;

  float thrustTheta_mc;
  float thrustTheta;

  float z1Ratio[2000];
  float z2Ratio[2000];
  int z1RatioCounter;
  int z2RatioCounter;

  float z1Ratio_mc[2000];
  float z2Ratio_mc[2000];
  int z1RatioCounter_mc;
  int z2RatioCounter_mc;

  float phiRSum_mc[2000];
  int phiRSum_mcCounter;

  int charge1Counter;
  int chargeType1[2000];
  int particle1Counter;
  int particleType1[2000];

  int charge2Counter;
  int chargeType2[2000];
  int particle2Counter;
  int particleType2[2000];

  int charge1Counter_mc;
  int chargeType1_mc[2000];
  int particle1Counter_mc;
  int particleType1_mc[2000];

  int charge2Counter_mc;
  int chargeType2_mc[2000];
  int particle2Counter_mc;
  int particleType2_mc[2000];

  float thetaEThrust=pi/(float)2; //for now until the new data is generated //correction factor is 1/4 sin^2(theta)
  int decayTheta1Counter;
  float decayTheta1[2000];
  int decayTheta2Counter;
  float decayTheta2[2000];

  float thrustProj11[2000];
  int thrustProj11Counter;
  float thrustProj12[2000];
  int thrustProj12Counter;
  float thrustProj21[2000];
  int thrustProj21Counter;
  float thrustProj22[2000];
  int thrustProj22Counter;

  int z1Loc=0;
  int z2Loc=1;
  int mass1Loc=2;
  int mass2Loc=3;
  int phiRSumLoc=4;


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

  fieldNamesI.push_back("phiR1Counter");
  fieldNamesF.push_back("phiR1");
  memLocsI.push_back(&phiR1Counter);
  memLocsF.push_back(phiR1);
  
  fieldNamesI.push_back("twoPhiZeroCounter");
  fieldNamesF.push_back("twoPhiZero");
  memLocsI.push_back(&twoPhiZeroCounter);
  memLocsF.push_back(twoPhiZero);

  fieldNamesI.push_back("phiZero1Counter");
  fieldNamesF.push_back("phiZero1");
  memLocsI.push_back(&phiZero1Counter);
  memLocsF.push_back(phiZero1);

  fieldNamesI.push_back("theta11Counter");
  fieldNamesF.push_back("theta11");
  memLocsI.push_back(&theta11Counter);
  memLocsF.push_back(theta11);

  fieldNamesI.push_back("theta12Counter");
  fieldNamesF.push_back("theta12");
  memLocsI.push_back(&theta12Counter);
  memLocsF.push_back(theta12);

  fieldNamesI.push_back("theta21Counter");
  fieldNamesF.push_back("theta21");
  memLocsI.push_back(&theta21Counter);
  memLocsF.push_back(theta21);

  fieldNamesI.push_back("theta22Counter");
  fieldNamesF.push_back("theta22");
  memLocsI.push_back(&theta22Counter);
  memLocsF.push_back(theta22);

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

  fieldNamesI.push_back("z1RatioCounter");
  fieldNamesF.push_back("z1Ratio");
  memLocsI.push_back(&z1RatioCounter);
  memLocsF.push_back(z1Ratio);

  fieldNamesI.push_back("z2RatioCounter");
  fieldNamesF.push_back("z2Ratio");
  memLocsI.push_back(&z2RatioCounter);
  memLocsF.push_back(z2Ratio);


  fieldNamesI.push_back("decayTheta1Counter");
  memLocsI.push_back(&decayTheta1Counter);
  fieldNamesF.push_back("decayTheta1");
  memLocsF.push_back(decayTheta1);

  fieldNamesI.push_back("decayTheta2Counter");
  memLocsI.push_back(&decayTheta2Counter);
  fieldNamesF.push_back("decayTheta2");
  memLocsF.push_back(decayTheta2);
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

#ifdef MC
  fieldNamesF.push_back("thrustTheta_mc");
  memLocsF.push_back(&thrustTheta_mc);

  fieldNamesF.push_back("thetaEThrust");
  memLocsF.push_back(&thetaEThrust);
  fieldNamesF.push_back("phiRSum_mc");
  memLocsF.push_back(phiRSum_mc);
  fieldNamesI.push_back("phiRSum_mcCounter");
  memLocsI.push_back(&phiRSum_mcCounter);

  fieldNamesI.push_back("chargeType1_mcCounter");
  fieldNamesI.push_back("chargeType1_mc");
  memLocsI.push_back(&charge1Counter_mc);
  memLocsI.push_back(chargeType1_mc);


  fieldNamesI.push_back("particleType1_mcCounter");
  fieldNamesI.push_back("particleType1_mc");
  memLocsI.push_back(&particle1Counter_mc);
  memLocsI.push_back(particleType1_mc);

  fieldNamesI.push_back("chargeType2_mcCounter");
  fieldNamesI.push_back("chargeType2_mc");
  memLocsI.push_back(&charge2Counter_mc);
  memLocsI.push_back(chargeType2_mc);


  fieldNamesI.push_back("particleType2_mcCounter");
  fieldNamesI.push_back("particleType2_mc");
  memLocsI.push_back(&particle2Counter_mc);
  memLocsI.push_back(particleType2_mc);

  fieldNamesI.push_back("z1Ratio_mcCounter");
  fieldNamesF.push_back("z1Ratio_mc");
  memLocsI.push_back(&z1RatioCounter_mc);
  memLocsF.push_back(z1Ratio_mc);

  fieldNamesI.push_back("z2Ratio_mcCounter");
  fieldNamesF.push_back("z2Ratio_mc");
  memLocsI.push_back(&z2RatioCounter_mc);
  memLocsF.push_back(z2Ratio_mc);

  fieldNamesI.push_back("mass1_mcCounter");
  fieldNamesF.push_back("mass1_mc");
  memLocsI.push_back(&mass1_mcCounter);
  memLocsF.push_back(mass1_mc);

  fieldNamesI.push_back("mass2_mcCounter");
  fieldNamesF.push_back("mass2_mc");
  memLocsI.push_back(&mass2_mcCounter);
  memLocsF.push_back(mass2_mc);

  fieldNamesI.push_back("z1_mcCounter");
  fieldNamesF.push_back("z1_mc");
  memLocsI.push_back(&z1_mcCounter);
  memLocsF.push_back(z1_mc);

  fieldNamesI.push_back("z2_mcCounter");
  fieldNamesF.push_back("z2_mc");
  memLocsI.push_back(&z2_mcCounter);
  memLocsF.push_back(z2_mc);


  TFile m_file("purities.root","recreate");
  TTree m_tree("PurityTree","PurityTree");

  m_tree.Branch("kinBin1",&purities::kinBin1,"kinBin1/I");
  m_tree.Branch("kinBin2",&purities::kinBin2,"kinBin2/I");
  m_tree.Branch("binning",&purities::binning,"binning/I");
  m_tree.Branch("chargeType1",&purities::chargeType1,"chargeType1/I");
  m_tree.Branch("particleType1",&purities::particleType1,"particleType1/I");
  m_tree.Branch("chargeType2",&purities::chargeType2,"chargeType2/I");
  m_tree.Branch("particleType2",&purities::particleType2,"particleType2/I");
  m_tree.Branch("purity",&purities::purity,"purity/F");
  //m_tree.Branch("purity",&purities::purity,"purity/F");
  m_tree.Branch("purityError",&purities::purityError,"purityError/F");
#endif
  long nevents=0;
#ifdef PAW_STYLE
  chI->Add((string(rootPath)+"/singlerootfiles/*.hr.root").c_str());
  chF->Add((string(rootPath)+"/singlerootfiles/*.hr.root").c_str());

  for(int i=0;i<fieldNamesF.size();i++)
    {
      cout <<"trying to branch" <<endl;
      float* memLoc=memLocsF[i];
      chF->SetBranchAddress(fieldNamesF[i].c_str(),memLoc);
    }
  for(int i=0;i<fieldNamesI.size();i++)
    {
      int* memLoc2=memLocsI[i];
      chI->SetBranchAddress(fieldNamesI[i].c_str(),memLoc2);
    }
  long neventsF=chF->GetEntries();
  long neventsI=chI->GetEntries();

  if(neventsF!=neventsI)
    {
      cout << "event mismatch" <<endl;
      exit(0);
    }
  nevents=neventsF;
#else
  chAll->Add((string(rootPath)+"/singlerootfiles/*.mr.root").c_str());
  for(int i=0;i<fieldNamesF.size();i++)
    {
      cout <<"trying to branch on " << fieldNamesF[i]<<endl;
      float* memLoc=memLocsF[i];
      chAll->SetBranchAddress(fieldNamesF[i].c_str(),memLoc);
    }
  for(int i=0;i<fieldNamesI.size();i++)
    {
      int* memLoc2=memLocsI[i];
      chAll->SetBranchAddress(fieldNamesI[i].c_str(),memLoc2);
    }
  nevents=chAll->GetEntries();
#endif
  pair<TH1D*,TH1D*> zSpectH[25];
  for(int i=0;i<5;i++)
    {
      for(int j=0;j<5;j++)
	{
	  stringstream ss;
	  ss <<"z1SpectForMBin"<<i<<"_"<<j;
	  stringstream ss2;
	  ss2 <<"z2SpectForMBin"<<i<<"_"<<j;
	   
	  zSpectH[i*5+j].first=new TH1D(ss.str().c_str(),ss.str().c_str(),1000,0,1);
	  zSpectH[i*5+j].second=new TH1D(ss2.str().c_str(),ss2.str().c_str(),1000,0,1);
	}
    }
  pair<TH1D*,TH1D*> zSpectWeighted[25];
  for(int i=0;i<5;i++)
    {
      for(int j=0;j<5;j++)
	{
	  stringstream ss;
	  stringstream ss2;
	  ss <<"z1SpectWeightedMBin"<<i<<"_"<<j;
	  ss2 <<"z2SpectWeightedMBin"<<i<<"_"<<j;
	  zSpectWeighted[i*5+j].first=new TH1D(ss.str().c_str(),ss.str().c_str(),1000,0,1);;
	  zSpectWeighted[i*5+j].second=new TH1D(ss2.str().c_str(),ss2.str().c_str(),1000,0,1);;
	}
    }
    
  //save counts differential in Type of binning, charge , had type, had1 kin, had2 kin, angular bin

  cout<<" numbin " << NumBin <<" numcharge: " << NumCharge <<", NumParticle: " << NumParticle <<", maxKin: " << maxKin <<", numAngBins: " << numAngBins<<endl;

  float* singleBinVal=new float[6];
  int****** countsFour=allocateArray<int>(NumBin,maxKin,maxKin,maxKin,maxKin,numAngBins);
  int****** counts=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  int***** countsSingle=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,numAngBins);
  //combine hads in the same hemi --- look at dihadana to see how they are combined... might be only pn or all..
  int****** countsSameHemi=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  int****** countsDiffEvt=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  float****** countsF=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  float****** countsF2=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  //relates z and m spectra
  kinSpectraPars=allocateArray<double>((NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKinSpect,maxKinSpect, maxKin,maxKin,NumPars);
  kinSpectraReducedPars=allocateArray<double>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKinSpect,maxKinSpect,NumPars);
  //more z bins to capture structure in the beginning...
  kinSpectra=allocateArray<double>((NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKinSpect,maxKinSpect, maxKin,maxKin);
  kinSpectraReduced=allocateArray<double>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKinSpect,maxKinSpect);
  kinSpectraH=allocHistos((NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKinSpect,maxKinSpect); //one th2d for all mbins
  gHCounter=0;
  kinSpectra1H=allocHistos1D((NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKinSpect,maxKinSpect); //one th1d for all mbins



  //avg over value2, conforming to notation in paper
  float***** xVals=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** xValsNum=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** yVals=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** yValsNum=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);

  float***** xValsFour=allocateArray<float>(NumBin,maxKin,maxKin,maxKin,maxKin);
  int***** xValsNumFour=allocateArray<int>(NumBin,maxKin,maxKin,maxKin,maxKin);
  float***** yValsFour=allocateArray<float>(NumBin,maxKin,maxKin,maxKin,maxKin);
  int***** yValsNumFour=allocateArray<int>(NumBin,maxKin,maxKin,maxKin,maxKin);
  float***** aValsFour=allocateArray<float>(NumBin,maxKin,maxKin,maxKin,maxKin);
  int***** aValsNumFour=allocateArray<int>(NumBin,maxKin,maxKin,maxKin,maxKin);
  float***** bValsFour=allocateArray<float>(NumBin,maxKin,maxKin,maxKin,maxKin);
  int***** bValsNumFour=allocateArray<int>(NumBin,maxKin,maxKin,maxKin,maxKin);

  int**** xValsNumSingle=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin);
  float**** xValsSingle=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin);

  float***** purities=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** puritiesErrors=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  //the sin2(theta) of thrust and beam axis
  float***** kinCorrFact=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);

#ifdef MC
  unsigned long***** numCorrIdent=allocateArray<unsigned long>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  unsigned long***** numFalseIdent=allocateArray<unsigned long>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
#else
  cout <<"loading purities" <<endl;
#ifdef WITH_PURITIES
  loadPurities(purities,puritiesErrors,"purities.root",LAST_PARTICLE_COMB,LAST_CHARGE_COMB,binningZ[PiPi].size(),binningM[PiPi].size());
#endif
#ifdef WITH_KIN_CORR
  //  loadKinematics(kinSpectraReducedPars,kinSpectraPars,kinSpectraReduced,kinSpectra,"kinematics.root",LAST_PARTICLE_COMB, LAST_CHARGE_COMB,binningZSpect[PiPi].size(), binningM[PiPi].size());
  loadKinPars1D(kinSpectra1H);
  cout<<"test cont ..." <<endl;
  cout <<"name; " << kinSpectra1H[0][0][0][0]->GetName() <<endl;
  cout <<"numBins; " << kinSpectra1H[0][0][0][0]->GetNbinsX() <<endl;

#endif

  //  cout <<"done" <<endl;
#endif 

  cout <<"there are " <<nevents << " events"<<endl;
  for(long i=0;i<nevents;i++)
    {
      if(!(i%100000))
	cout <<"processing event nr " << i << " of " << nevents << "(" << 100*i/(float)nevents<< "% )"<<endl;
#ifdef PAW_STYLE
      chF->GetEntry(i);
      chI->GetEntry(i);
#else
      chAll->GetEntry(i);
#endif

#ifdef MAX_EVENTS
      if(i>MAX_EVENTS)
	break;
#endif


      for(int iHadsInEv=1;iHadsInEv<z1Counter;iHadsInEv++)
	{
	  //	  if((z1Ratio[iHadsInEv]==z1Ratio[iHadsInEv-1])
	    
	}


      //      cout <<"after ge" <<endl;
      //all counters should be the same, so it is a waste of space
      //      cout <<"z1Counter; " << z1Counter <<", z1Mem: " << memLocsF[z1Loc]<<endl;
      for(int iHadsInEv=0;iHadsInEv<z1Counter;iHadsInEv++)
	{
	  if(fabs(thrustProj11[iHadsInEv])<OPENING_CUT || fabs(thrustProj12[iHadsInEv])<OPENING_CUT ||fabs(thrustProj21[iHadsInEv])<OPENING_CUT ||fabs(thrustProj22[iHadsInEv])<OPENING_CUT)
	    {
	      continue;
	    }
	  int chargeIndx=ind(chargeType1[iHadsInEv],chargeType2[iHadsInEv],NumCharge);
	  int mIndex=ind(particleType1[iHadsInEv],particleType2[iHadsInEv],NumParticle);	  
	  hThetaDiff11.Fill(theta11[iHadsInEv]-thrustTheta);
	  hThetaDiff12.Fill(theta12[iHadsInEv]-thrustTheta);
	  hThetaDiff21.Fill(theta21[iHadsInEv]-thrustTheta);
	  hThetaDiff22.Fill(theta22[iHadsInEv]-thrustTheta);
	  hThrustThetaRes.Fill(thrustTheta-thrustTheta_mc);
	  if(((float*)(memLocsF[z1Loc]))[iHadsInEv]>ZCut_for_mPlot && ((float*)(memLocsF[z2Loc]))[iHadsInEv]>ZCut_for_mPlot)
	    {
	      if(particleType1[iHadsInEv]==PiPi&&chargeType1[iHadsInEv]==PN&&particleType2[iHadsInEv]==PiPi&&chargeType2[iHadsInEv]==PN )
		{
		  hPhiR.Fill(cos(((float*)(memLocsF[phiRSumLoc]))[iHadsInEv]));
		  hPhiZ.Fill(cos(twoPhiZero[iHadsInEv]));
		  massPiPi_PN->Fill(((float*)(memLocsF[mass1Loc]))[iHadsInEv]);
		  massPiPi_PN->Fill(((float*)(memLocsF[mass2Loc]))[iHadsInEv]);
		}
	      if(particleType1[iHadsInEv]==PiPi&&chargeType1[iHadsInEv]==PZ&&particleType2[iHadsInEv]==PiPi&&chargeType2[iHadsInEv]==PZ)
		{
		  massPiPi_PZ->Fill(((float*)(memLocsF[mass1Loc]))[iHadsInEv]);
		  massPiPi_PZ->Fill(((float*)(memLocsF[mass2Loc]))[iHadsInEv]);
		}
	      if(particleType1[iHadsInEv]==PiPi&&chargeType1[iHadsInEv]==ZN &&particleType2[iHadsInEv]==PiPi&&chargeType2[iHadsInEv]==ZN)
		{
		  massPiPi_ZN->Fill(((float*)(memLocsF[mass1Loc]))[iHadsInEv]);
		  massPiPi_ZN->Fill(((float*)(memLocsF[mass2Loc]))[iHadsInEv]);
		}
	      if(particleType1[iHadsInEv]==PiK&&chargeType1[iHadsInEv]==PN && particleType2[iHadsInEv]==PiK&&chargeType2[iHadsInEv]==PN)
		{
		  massPiK_PN->Fill(((float*)(memLocsF[mass1Loc]))[iHadsInEv]);
		  massPiK_PN->Fill(((float*)(memLocsF[mass2Loc]))[iHadsInEv]);
		}
	    }

	  float zValue1=((float*)(memLocsF[z1Loc]))[iHadsInEv];
	  float mValue1=((float*)(memLocsF[mass1Loc]))[iHadsInEv];
	  float zValue2=((float*)(memLocsF[z2Loc]))[iHadsInEv];
	  float mValue2=((float*)(memLocsF[mass2Loc]))[iHadsInEv];

	  if(mValue1 > 2.0 || mValue2 > 2.0 || zValue1 > 1.0 || zValue2 > 1.0)
	    continue;
	  //	  kinematics[iBin][ind(chargeType1[iHadsInEv],chargeType2[iHadsInEv],NumCharge)][ind(particleType1[iHadsInEv],particleType2[iHadsInEv],NumParticle)][getBin(binningZ[particleType1[iHadsInEv]],zValue)][getBin(binningM[particleType1[iHadsInEv]],mValue)];

	  //doesn't make any sense, since the one is in cms, the other in lab
	  //	  if(fabs(theta11[iHadsInEv]-thrustTheta) > thetaCut || fabs(theta12[iHadsInEv]-thrustTheta) >thetaCut || fabs(theta21[iHadsInEv]-thrustTheta)>thetaCut || fabs(theta22[iHadsInEv]-thrustTheta)>thetaCut)
	  //	    continue;
	  int zbin1=getBin(binningZ[particleType1[iHadsInEv]],zValue1);
	  int zbin2=getBin(binningZ[particleType1[iHadsInEv]],zValue2);
	  int zSpectbin1=getBin(binningZSpect[particleType1[iHadsInEv]],zValue1);
	  int zSpectbin2=getBin(binningZSpect[particleType1[iHadsInEv]],zValue2);
	  int mbin1=getBin(binningM[particleType1[iHadsInEv]],mValue1);
	  int mbin2=getBin(binningM[particleType1[iHadsInEv]],mValue2);

	  int m_chType1=chargeType1[iHadsInEv];
	  int m_pt1=particleType1[iHadsInEv];

	  int m_chType2=chargeType2[iHadsInEv];
	  int m_pt2=particleType2[iHadsInEv];
	  float value1=((float*)(memLocsF[z1Loc]))[iHadsInEv];
	  if(isnan(value1))
	    {
	      cout <<"value 1 nan" <<endl;
	      nanEv++;
	      continue;
	    }
	  //		    cout <<"1+1 " <<endl;
	  float value2=((float*)(memLocsF[z2Loc]))[iHadsInEv];
	  if(isnan(value2))
	    {
	      cout <<"value 2 nan" <<endl;
	      nanEv++;
	      continue;
	    }
	  float value3=((float*)(memLocsF[phiRSumLoc]))[iHadsInEv];
	  if(isnan(value3))
	    {
	      cout <<"value 3 nan" <<endl;
	      nanEv++;
	      continue;
	    }
	  //single binnings, for current release, so only pipi (probably going to regret this in the future)
	  if(PiPi==particleType1[iHadsInEv]&& PiPi==particleType2[iHadsInEv])
	    {
	      singleBinVal[thetaBinning]=0.25*sin(thetaEThrust)*sin(thetaEThrust);
	      singleBinVal[dilFactBinning]=1;
	      for(int iBin=thetaBinning;iBin<=dilFactBinning;iBin++)
		{
		    fillBorders(singleBinning[iBin][particleType1[iHadsInEv]],binningAng,singleBinVal[iBin],value3,countsSingle[iBin][chargeIndx][mIndex],xValsSingle[iBin][chargeIndx][mIndex],xValsNumSingle[iBin][chargeIndx][mIndex]);

		}
	      for(int iBin=fourBinningOne;iBin<=fourBinningOne;iBin++)
		{
		  fillBorders4(fourBinning[iBin][0][particleType1[iHadsInEv]],fourBinning[iBin][1][particleType1[iHadsInEv]],binningAng,mValue1,mValue2,zValue1,zValue2,value3,countsFour[iBin],xValsFour[iBin],xValsNumFour[iBin],yValsFour[iBin],yValsNumFour[iBin],aValsFour[iBin],aValsNumFour[iBin],bValsFour[iBin],bValsNumFour[iBin]);
		}
	    }
	  for(int iBin=zBinning;iBin<=mBinning;iBin++)
	    {
	      switch(iBin)
		{
		case zBinning:
		  {
		    hZ1.Fill(value1);
		    hZ2.Fill(value2);

		    if(value2>zHighCuts[particleType1[iHadsInEv]] || value1 > zHighCuts[particleType1[iHadsInEv]])
		      {
			cuttedEv++;
			continue;
		      }
		    accEv++;
#ifndef WITH_KIN_CORR
		    if(zSpectbin1>=0 && zSpectbin2>=0 && mbin1>=0 && mbin2>=0)
		      {
			kinSpectra[chargeIndx][mIndex][zSpectbin1][zSpectbin2][mbin1][mbin2]++;
			cout <<"Filling: " << chargeIndx<<", " <<mIndex<< " " << mbin1 <<" " << mbin2 <<endl;
			kinSpectraH[chargeIndx][mIndex][mbin1][mbin2]->Fill(value1,value2);
			kinSpectra1H[chargeIndx][mIndex][mbin1][mbin2]->Fill(value1);
			kinSpectra1H[chargeIndx][mIndex][mbin1][mbin2]->Fill(value2);
		      }
#endif

#ifdef MC
		    if(particleType1[iHadsInEv]==particleType1_mc[iHadsInEv] && chargeType1[iHadsInEv]==chargeType1_mc[iHadsInEv]&& particleType2[iHadsInEv]==particleType2_mc[iHadsInEv] && chargeType2[iHadsInEv]==chargeType2_mc[iHadsInEv])
		      incBorders(binningZ[particleType1[iHadsInEv]],value1,value2,numCorrIdent[iBin][chargeIndx][mIndex]);
		    else
		      incBorders(binningZ[particleType1[iHadsInEv]],value1,value2,numFalseIdent[iBin][chargeIndx][mIndex]);
#endif
		    incBorders(binningZ[particleType1[iHadsInEv]],value1,value2,kinCorrFact[iBin][chargeIndx][mIndex],0.25*sin(thetaEThrust)*sin(thetaEThrust));
		    fillBorders(binningZ[particleType1[iHadsInEv]],binningAng,value1,value2,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);


#ifdef WITH_KIN_CORR

		    if((chargeType1[iHadsInEv]==PN || chargeType1[iHadsInEv]==PZ || chargeType1[iHadsInEv]==ZN)&&(chargeType2[iHadsInEv]==PN || chargeType2[iHadsInEv]==PZ || chargeType2[iHadsInEv]==ZN))
		      {
			//	cout << "ch1: " << m_chType1 <<" ch2: " << m_chType2 << " pt1: " << m_pt1 << " pt2:  " << m_pt2 << " zbin1: " << zbin1 << " zbin2; " << zbin2 << " mbin1: " << mbin1 <<" m2: " << mbin2 <<endl;
			//	fillBorders(binningZ[particleType1[iHadsInEv]],binningAng,value1,value2,value3,countsF2[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex],getInc2(iBin,chargeIndx,mIndex,mbin1,mbin2,zbin1, zbin2));
		      }
#endif
#ifdef MC
		    if(phiRSum_mc[iHadsInEv]!=-1)
		      {
			hR.Fill(phiRSum_mc[iHadsInEv]-((float*)memLocsF[phiRSumLoc])[iHadsInEv]);
			fillBorders(binningZ[particleType1[iHadsInEv]],binningAng,value1,value2,value3,countsF[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex],getInc(phiRSum_mc[iHadsInEv],mass1_mc[iHadsInEv],mass2_mc[iHadsInEv],z1_mc[iHadsInEv], z2_mc[iHadsInEv],particleType1_mc[iHadsInEv],particleType2_mc[iHadsInEv],chargeType1_mc[iHadsInEv],chargeType2_mc[iHadsInEv]));
		      }
#endif
		    break;
		  }
		case mBinning:
		  {
		    if(PiPi==particleType1[iHadsInEv]&& PiPi==particleType2[iHadsInEv])
		      {
			if(mbin1>=0 && mbin2 >=0)
			  {
			    zSpectH[mbin1*5+mbin2].first->Fill(value1);
			    zSpectH[mbin1*5+mbin2].second->Fill(value2);
			  }
			//			zSpectWeighted[mbin1][mbin2].first->Fill(value1);
			//			zSpectWeighted[mbin1][mbin2].second->Fill(value2);
		      }


		    float value1=((float*)(memLocsF[mass1Loc]))[iHadsInEv];
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
		      continue;

		    accEv++;
#ifdef MC
		    if(particleType1[iHadsInEv]==particleType1_mc[iHadsInEv] && chargeType1[iHadsInEv]==chargeType1_mc[iHadsInEv]&&particleType2[iHadsInEv]==particleType2_mc[iHadsInEv] && chargeType2[iHadsInEv]==chargeType2_mc[iHadsInEv])
		      incBorders(binningM[particleType1[iHadsInEv]],value1,value2,numCorrIdent[iBin][chargeIndx][mIndex]);
		    else
		      incBorders(binningM[particleType1[iHadsInEv]],value1,value2,numFalseIdent[iBin][chargeIndx][mIndex]);
#endif
		    incBorders(binningM[particleType1[iHadsInEv]],value1,value2,kinCorrFact[iBin][chargeIndx][mIndex],0.25*sin(thetaEThrust)*sin(thetaEThrust));
		    fillBorders(binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);
#ifdef WITH_KIN_CORR
		    if((chargeType1[iHadsInEv]==PN || chargeType1[iHadsInEv]==PZ || chargeType1[iHadsInEv]==ZN)&&(chargeType2[iHadsInEv]==PN || chargeType2[iHadsInEv]==PZ || chargeType2[iHadsInEv]==ZN))
		      {
			//	cout << "ch1: " << m_chType1 <<" ch2: " << m_chType2 << " pt1: " << m_pt1 << " pt2:  " << m_pt2 << " zbin1: " << zbin1 << " zbin2; " << zbin2 << " mbin1: " << mbin1 <<" m2: " << mbin2 <<endl;
			fillBorders(binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,countsF2[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex],getInc2(iBin,chargeIndx,mIndex,mbin1,mbin2,zValue1, zValue2));
		      }
#endif
#ifdef MC		    
		    if(phiRSum_mc[iHadsInEv]!=-1)
		      fillBorders(binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,countsF[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex],getInc(phiRSum_mc[iHadsInEv],mass1_mc[iHadsInEv],mass2_mc[iHadsInEv],z1_mc[iHadsInEv], z2_mc[iHadsInEv],particleType1_mc[iHadsInEv],particleType2_mc[iHadsInEv],chargeType1_mc[iHadsInEv],chargeType2_mc[iHadsInEv]));
#endif
		    break;
		  }
		default:
		  cout << "binning not recognized" << endl;
		  exit(0);
		}
	    }
	}
    }
  //and now the fit...
  cout<<"done with filling arryas, fitting..."<<endl;
  vector<vector<pair<float,float> >* > zAsyms;
  vector<vector<pair<float,float> >* >mAsyms;
  vector<vector<pair<float,float> >*** > fourAsyms;
  vector<vector<pair<float,float> >* > singleAsyms;
  vector<vector<pair<float,float> >* > zAsymsMod;
  vector<vector<pair<float,float> >* >mAsymsMod;
  vector<vector<pair<float,float> >* > zAsymsMod2;
  vector<vector<pair<float,float> >* >mAsymsMod2;
  int numBins;
  string strCharge;
  string strParticle;
  vector<pair<float ,float> > *** tmpAsymFour;
  vector<pair<float ,float> > * tmpAsymSingle;
  vector<pair<float ,float> > * tmpAsym;
  vector<pair<float ,float> > * tmpAsymMod;
  vector<pair<float ,float> > * tmpAsymMod2;
  vector<vector<float>* > v_chi2Z;
  vector<vector<int>* > v_ndfsZ;
  vector<vector<float>* > v_chi2M;
  vector<vector<int>* > v_ndfsM;
  //
  //     Fit with the four kinematic bins
  //
  for(int iBin=fourBinningOne;iBin<=fourBinningOne;iBin++)
    {
      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	{
	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	      if(!(iCh1==PN))
		continue;
	      if(!(iCh2==PN))
		continue;
	      strCharge=iToStrCharge(iCh1)+"_"+iToStrCharge(iCh2);
	      cout <<"iCh: " << iCh1<<"ich2: " << iCh2 <<endl;
	      for(int iPa1=PiPi;iPa1<=LAST_PARTICLE_COMB;iPa1++)
		//	  for(int iPa=PiPi;iPa<=UNKNOWN;iPa++)
		{
		  for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
		    {
		      if(iPa1!=PiPi || iPa2!=PiPi)
			continue;
		      strParticle=iToStrPart(iPa1)+"_"+iToStrPart(iPa2);
		      cout << "ipa: " << iPa1 << " ipa2: " <<iPa2<< endl;
		      vector<float>* t_chi2=new vector<float>;
		      vector<int>* t_ndfs=new vector<int>;
		      outFile << strCharge <<" " << strParticle << " zBinning: " <<endl;
		      numBins=binningZ[iPa1].size();
		      int redMbinningNum=fourBinning[iBin][0][iPa1].size();
		      int redZbinningNum=fourBinning[iBin][0][iPa1].size();
		      cout <<"for fitting: " << redMbinningNum <<" and z: " << redZbinningNum <<endl;
		      tmpAsymFour=fitTheFour(countsFour[iBin],binningAng,redMbinningNum, redZbinningNum, errorVglFile);
		      fourAsyms.push_back(tmpAsymFour);
			}
		    }
		}
	    }
    }

//
//  Fit Single Kinematic bins
//
for(int iBin=thetaBinning;iBin<=dilFactBinning;iBin++)
  {
    for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
      {
	for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	  {
	    if(!(iCh1==PN))
	      continue;
	    if(!(iCh2==PN))
	      continue;
	    strCharge=iToStrCharge(iCh1)+"_"+iToStrCharge(iCh2);
	    cout <<"iCh: " << iCh1<<"ich2: " << iCh2 <<endl;
	    for(int iPa1=PiPi;iPa1<=LAST_PARTICLE_COMB;iPa1++)
	      {
		for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
		  {
		    if(iPa1!=PiPi || iPa2!=PiPi)
		      continue;
		    strParticle=iToStrPart(iPa1)+"_"+iToStrPart(iPa2);
		    cout << "ipa: " << iPa1 << " ipa2: " <<iPa2<< endl;
		    vector<float>* t_chi2=new vector<float>;
		    vector<int>* t_ndfs=new vector<int>;
		    outFile << strCharge <<" " << strParticle << " zBinning: " <<endl;
		    numBins=binningZ[iPa1].size();
		    int numSingleBin=singleBinning[iBin][iPa1].size();
		    tmpAsymSingle=fitTheSingle(countsSingle[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numSingleBin, errorVglFile);
		    singleAsyms.push_back(tmpAsymSingle);
		  }
	      }
	  }
      }
  }

  for(int iBin=zBinning;iBin<=mBinning;iBin++)
    {
      //      for(int iCh=PN;iCh<=ZZ;iCh++)
      cout <<"iBin: " << iBin <<endl;
      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	{
	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	      /*      if(!(iCh1==PN || iCh1==PZ || iCh1==ZN || iCh1==PNNP || iCh1==ZNNZ))
		      continue;
		      if(!(iCh2==PN || iCh2==PZ || iCh2==ZN || iCh2==PNNP || iCh2==ZNNZ))
		      continue;*/
	      if(!(iCh1==PN || iCh1==PZ || iCh1==ZN))
		continue;
	      if(!(iCh2==PN || iCh2==PZ || iCh2==ZN))
		continue;
	      strCharge=iToStrCharge(iCh1)+"_"+iToStrCharge(iCh2);
	      cout <<"iCh: " << iCh1<<"ich2: " << iCh2 <<endl;
	      for(int iPa1=PiPi;iPa1<=LAST_PARTICLE_COMB;iPa1++)
		//	  for(int iPa=PiPi;iPa<=UNKNOWN;iPa++)
		{
		  for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
		    {
		      strParticle=iToStrPart(iPa1)+"_"+iToStrPart(iPa2);
		      cout << "ipa: " << iPa1 << " ipa2: " <<iPa2<< endl;
		      vector<float>* t_chi2=new vector<float>;
		      vector<int>* t_ndfs=new vector<int>;
		      switch(iBin)
			{
			case zBinning:
			  {
			    outFile << strCharge <<" " << strParticle << " zBinning: " <<endl;
			    numBins=binningZ[iPa1].size();


			    tmpAsym=fitTheSh__(counts[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2,t_ndfs);
			    for(int iChi=0;iChi<t_chi2->size();iChi++)
			      {
				cout <<"chi2: " << (*t_chi2)[iChi] << ", ndf: " << (*t_ndfs)[iChi]<<", chi/ndf:" << (*t_chi2)[iChi]/(*t_ndfs)[iChi] <<endl;
			      }

			    v_chi2Z.push_back(t_chi2);
			    v_ndfsZ.push_back(t_ndfs);
			    zAsyms.push_back(tmpAsym);
			    tmpAsymMod2=fitTheSh__(countsF2[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins,errorVglFile);
			    zAsymsMod2.push_back(tmpAsymMod2);
#ifdef MC
			    tmpAsymMod=fitTheSh__(countsF[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile);
			    zAsymsMod.push_back(tmpAsymMod);
#endif
			    cout <<"size of asy vec: " << tmpAsym->size()<<endl;
			    for(int k=0;k<tmpAsym->size();k++)
			      {
				outFile  <<"Bin " << k << " Asymmetry: " << (*tmpAsym)[k].first << " +- " <<  (*tmpAsym)[k].second <<endl;
#ifdef MC
				outFile  <<"Bin " << k << " Asymmetry: " << (*tmpAsymMod)[k].first << " +- " <<  (*tmpAsymMod)[k].second <<endl;
#endif
			      }
			    cout <<"wrote to file " <<endl;
			    break;
			  }
			case mBinning:
			  {
			    outFile << strCharge <<" " << strParticle << " M_inv Binning: " <<endl;
			    numBins=binningM[iPa1].size();
			    tmpAsym=fitTheSh__(counts[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2,t_ndfs);
			    v_chi2M.push_back(t_chi2);
			    v_ndfsM.push_back(t_ndfs);
			    mAsyms.push_back(tmpAsym);
			    tmpAsymMod2=fitTheSh__(countsF2[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins,errorVglFile);
			    mAsymsMod2.push_back(tmpAsymMod2);
#ifdef MC
			    tmpAsymMod=fitTheSh__(countsF[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins,errorVglFile);
			    mAsymsMod.push_back(tmpAsym);
#endif
			    for(int k=0;k<tmpAsym->size();k++)
			      {
				outFile  <<"Bin " << k << " Asymmetry: " << (*tmpAsym)[k].first << " +- " <<  (*tmpAsym)[k].second <<endl;
#ifdef MC
				outFile  <<"Bin " << k << " Asymmetry: " << (*tmpAsymMod)[k].first << " +- " <<  (*tmpAsymMod)[k].second <<endl;
#endif
			      }
			    cout <<"after fit m" <<endl;
			    break;
			  }
			default:
			  cout << "binning not recognized" << endl;
			  exit(0);
			}

		    }
		}
	    }
	}
    }

  vector<int> vecBinning;
  //draw the asymmetries...
  vector<TGraphErrors*> vecGraphs;
  vector<TGraphErrors*> vecGraphsKinCorr;
  vector<TGraphErrors*> vecGraphsProjections;
  vector<TGraphErrors*> vecPurities;
  vector<string> graphTitles;
  vector<string> graphTitlesKinCorr;
  vector<string> graphTitlesPurities;
  vector<string> graphTitlesProj;
  int asymIndexZ=0;
  int asymIndexM=0;
  int asymIndexSingle=0;
  int asymIndexFour=0;

  cout << "maxKin : " << maxKin <<endl;
  float* mX=new float[maxKin];
  float* mY=new float[maxKin];
  float* mYKinCorr=new float[maxKin];
  float* mY_Pr=new float[maxKin];
  float* mXErr=new float[maxKin];
  float* mYErr=new float[maxKin];
  float* mYErrKinCorr=new float[maxKin];
  float* mYErr_Pr=new float[maxKin];

  //as xCheck, sum all asyms and check if they are zero:
  int aCounter=0;
  float asymSum=0;
  float asymSumErr=0;
  for(int i=0;i<zAsyms.size();i++)
    {
      for(int j=0;j<(*zAsyms[i]).size();j++)
	{
	  if(!zAsyms[i]||(*zAsyms[i])[j].second > 3) //no counts
	    continue;
	  asymSum+=(*zAsyms[i])[j].first;
	  asymSumErr+=(*zAsyms[i])[j].second*(*zAsyms[i])[j].second;

	  outputFile << "adding Asymmetry (z): " << (*zAsyms[i])[j].first <<" +- " << (*zAsyms[i])[j].second <<endl;
	  aCounter++;
	}
    }
  cout <<"zAsymAvg: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"zAsymAvg: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"num Asyms: " << aCounter <<endl;
  aCounter=0;
  asymSum=0;

  for(int i=0;i<mAsyms.size();i++)
    {
      for(int j=0;j<(*mAsyms[i]).size();j++)
	{
	  if((*mAsyms[i])[j].second > 3) //no counts
	    continue;
	  asymSum+=(*mAsyms[i])[j].first;
	  asymSumErr+=(*mAsyms[i])[j].second*(*mAsyms[i])[j].second;
	  outputFile << "adding Asymmetry (m): " << (*mAsyms[i])[j].first <<" +- " << (*mAsyms[i])[j].second <<endl;
	  aCounter++;
	}
    }
  cout <<"mAsymAvg: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"mAsymAvg: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"num Asyms: " << aCounter <<endl;

#ifdef MC
  asymSum=0;
  aCounter=0;
  for(int i=0;i<zAsymsMod.size();i++)
    {
      for(int j=0;j<(*zAsymsMod[i]).size();j++)
	{
	  if((*zAsymsMod[i])[j].second > 3) //no counts
	    continue;
	  asymSum+=(*zAsymsMod[i])[j].first;
	  asymSumErr+=(*zAsymsMod[i])[j].second*(*zAsymsMod[i])[j].second;
	  hAsDiff.Fill(((*zAsymsMod[i])[j].first-AMP_PHIR)/(*zAsymsMod[i])[j].second);
	  outputFile << "adding AsymmetryMod: " << (*zAsymsMod[i])[j].first <<" +- " << (*zAsymsMod[i])[j].second <<endl;
	  aCounter++;
	}
    }
  cout <<"zAsymAvgMod: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"zAsymAvgMod: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"num ModAsyms: " << aCounter <<endl;
  aCounter=0;
  asymSum=0;

  for(int i=0;i<mAsymsMod.size();i++)
    {
      for(int j=0;j<(*mAsymsMod[i]).size();j++)
	{
	  if((*mAsymsMod[i])[j].second > 3) //no counts
	    continue;
	  asymSum+=(*mAsymsMod[i])[j].first;
	  asymSumErr+=(*mAsymsMod[i])[j].second*(*mAsymsMod[i])[j].second;
	  hAsDiff.Fill(((*mAsymsMod[i])[j].first-AMP_PHIR)/(*mAsymsMod[i])[j].second);
	  outputFile << "adding AsymmetryMod: " << (*mAsymsMod[i])[j].first <<" +- " << (*mAsymsMod[i])[j].second <<endl;
	  aCounter++;
	}
    }
  cout <<"mAsymAvgMod: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"mAsymAvgMod: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"num ModAsyms: " << aCounter <<endl;
#endif

  asymSum=0;
  aCounter=0;
  for(int i=0;i<zAsymsMod2.size();i++)
    {
      for(int j=0;j<(*zAsymsMod2[i]).size();j++)
	{
	  if((*zAsymsMod2[i])[j].second > 3) //no counts
	    continue;
	  asymSum+=(*zAsymsMod2[i])[j].first;
	  asymSumErr+=(*zAsymsMod2[i])[j].second*(*zAsymsMod2[i])[j].second;
	  hAsDiff.Fill(((*zAsymsMod2[i])[j].first-AMP_PHIR)/(*zAsymsMod2[i])[j].second);
	  outputFile << "adding AsymmetryMod: " << (*zAsymsMod2[i])[j].first <<" +- " << (*zAsymsMod2[i])[j].second <<endl;
	  aCounter++;
	}
    }
  cout <<"zAsymAvgMod2: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"zAsymAvgMod2: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"num ModAsyms2: " << aCounter <<endl;
  aCounter=0;
  asymSum=0;

  for(int i=0;i<mAsymsMod.size();i++)
    {
      for(int j=0;j<(*mAsymsMod[i]).size();j++)
	{
	  if((*mAsymsMod[i])[j].second > 3) //no counts
	    continue;
	  asymSum+=(*mAsymsMod[i])[j].first;
	  asymSumErr+=(*mAsymsMod[i])[j].second*(*mAsymsMod[i])[j].second;
	  hAsDiff.Fill(((*mAsymsMod[i])[j].first-AMP_PHIR)/(*mAsymsMod[i])[j].second);
	  outputFile << "adding AsymmetryMod: " << (*mAsymsMod[i])[j].first <<" +- " << (*mAsymsMod[i])[j].second <<endl;
	  aCounter++;
	}
    }
  cout <<"mAsymAvgMod2: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"mAsymAvgMod2: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"num ModAsyms2: " << aCounter <<endl;
  saveKinPars(kinSpectraH, binningM);
  saveKinPars(kinSpectra1H, binningM);
  ///-----end mod2
  TFile asFile(("AsymmetriesRFiles/"+folderName+".root").c_str(),"recreate");
  TTree as_tree("AsymmetriesTree","AsymmetriesTree");
  TTree as_treeFour("AsymmetriesTreeFour","AsymmetriesTreeFour");
  TTree as_treeSingle("AsymmetriesTreeSingle","AsymmetriesTreeSingle");

  TTree* pTrees[]={&as_tree,&as_treeFour,&as_treeSingle};
  for(int i=0;i<3;i++)
    {
      pTrees[i]->Branch("kinBin1",&asymmetries::kinBin1,"kinBin1/I");
      pTrees[i]->Branch("kinBin2",&asymmetries::kinBin2,"kinBin2/I");
      pTrees[i]->Branch("kinBin3",&asymmetries::kinBin3,"kinBin3/I");
      pTrees[i]->Branch("kinBin4",&asymmetries::kinBin4,"kinBin4/I");
      pTrees[i]->Branch("binning",&asymmetries::binning,"binning/I");
      pTrees[i]->Branch("chargeType1",&asymmetries::chargeType1,"chargeType1/I");
      pTrees[i]->Branch("particleType1",&asymmetries::particleType1,"particleType1/I");
      pTrees[i]->Branch("chargeType2",&asymmetries::chargeType2,"chargeType2/I");
      pTrees[i]->Branch("particleType2",&asymmetries::particleType2,"particleType2/I");
      pTrees[i]->Branch("asymmetry",&asymmetries::asymmetry,"asymmetry/F");
      pTrees[i]->Branch("asError",&asymmetries::asError,"asError/F");
      pTrees[i]->Branch("meanKinVal1",&asymmetries::meanKinVal1,"meanKinVal1/F");
      pTrees[i]->Branch("meanKinVal2",&asymmetries::meanKinVal2,"meanKinVal2/F");
      pTrees[i]->Branch("meanKinVal3",&asymmetries::meanKinVal3,"meanKinVal3/F");
      pTrees[i]->Branch("meanKinVal4",&asymmetries::meanKinVal4,"meanKinVal4/F");
      pTrees[i]->Branch("meanTheta",&asymmetries::meanTheta,"meanTheta/F");
      pTrees[i]->Branch("chi2Fit",&asymmetries::chi2Fit,"chi2Fit/F");
      pTrees[i]->Branch("ndfFit",&asymmetries::ndfFit,"ndfFit/I");
      //not necessary to save for every value but... what the heck
      pTrees[i]->Branch("expNr",&asymmetries::expNr,"expNr/I");
      pTrees[i]->Branch("isOnRes",&asymmetries::isOnRes,"isOnRes/I");
    }


  //singles
  for(int iBin=thetaBinning;iBin<=dilFactBinning;iBin++)
    {
      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	//      for(int iCh=PN;iCh<=ZZ;iCh++)
	{
	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	      if(!(iCh1==PN))
		continue;
	      if(!(iCh2==PN))
		continue;
	      strCharge=iToStrCharge(iCh1)+"_"+iToStrCharge(iCh2);
	      for(int iPa1=PiPi;iPa1<=LAST_PARTICLE_COMB;iPa1++)
		{
		  for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
		    {
		      if(iPa1!=PiPi)
			continue;
		      if(iPa2!=PiPi)
			continue;
		      int numSBins=singleBinning[iBin][iCh1].size();
		      for(int iSBin=0;iSBin<numSBins;iSBin++)
			{
			  mX[iSBin]=xValsSingle[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iSBin]/xValsNumSingle[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iSBin];
			  mXErr[iSBin]=0;
			  mY[iSBin]=(*singleAsyms[asymIndexSingle])[iSBin].first;
			  mYErr[iSBin]=(*singleAsyms[asymIndexSingle])[iSBin].second;
			  mXErr[iSBin]=0;
			  asymmetries::binning=iBin;
			  asymmetries::kinBin1=iSBin;
			  asymmetries::kinBin2=0;
			  asymmetries::kinBin3=0;
			  asymmetries::kinBin4=0;
			  asymmetries::chargeType1=PN;
			  asymmetries::chargeType2=PN;
			  asymmetries::particleType1=PiPi;
			  asymmetries::particleType2=PiPi;
			  asymmetries::asymmetry=mY[iSBin];
			  asymmetries::asError=mYErr[iSBin];
			  asymmetries::meanKinVal1=mX[iSBin];
			  asymmetries::meanKinVal2=0;
			  asymmetries::meanKinVal3=0;
			  asymmetries::meanKinVal4=0;
			  asymmetries::chi2Fit=0;
			  asymmetries::ndfFit=0;
			  asymmetries::expNr=expNumber;
			  asymmetries::isOnRes=onResonance;
			  as_treeSingle.Fill();
			}
		      asymIndexSingle++;
		    }
		}
	    }
	}
    }
  //four  binning
  for(int iBin=fourBinningOne;iBin<=fourBinningOne;iBin++)
    {
      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	//      for(int iCh=PN;iCh<=ZZ;iCh++)
	{
	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	      //important for the index counter, so that the correct asymmetry is associated with the bins etc...
	      if(!(iCh1==PN))
		continue;
	      if(!(iCh2==PN))
		continue;
	      strCharge=iToStrCharge(iCh1)+"_"+iToStrCharge(iCh2);
	      //we only have pipi, pn anyways for this binning...
	      int redMbinningNum=fourBinning[iBin][0][PiPi].size();
	      int redZbinningNum=fourBinning[iBin][0][PiPi].size();
	      cout <<" redMbinningNum: " << redMbinningNum <<", zbinningNum: " << redZbinningNum <<endl;
	      for(int iM1=0;iM1<redMbinningNum;iM1++)
		{
		  for(int iM2=0;iM2<redMbinningNum;iM2++)
		    {
		      for(int iZ1=0;iZ1<redZbinningNum;iZ1++)
			{
			  for(int iZ2=0;iZ2<redZbinningNum;iZ2++)
			    {
			      mX[iZ2]=xValsFour[iBin][iM1][iM2][iZ1][iZ2]/xValsNum[iBin][iM1][iM2][iZ1][iZ2];
			      mXErr[iZ2]=0;
			      cout <<"asymIndexFour: " << asymIndexFour << "im: " <<iM1 << " , iM2: " << iM2 << ", iZ1 " << iZ1 << " iz2: " << iZ2 <<endl;
			      cout <<"num in fourasyms.size(): " << fourAsyms.size() <<endl;
			      cout <<" size2: " << (*fourAsyms[asymIndexFour])[iM1][iM2].size()<<endl;
			      mY[iZ2]=(*(fourAsyms[asymIndexFour])[iM1][iM2])[iZ1*redZbinningNum+iZ2].first;
			      cout <<"my: " << mY[iZ2] <<endl;
			      mYErr[iZ2]=(*(fourAsyms[asymIndexFour])[iM1][iM2])[iZ1*redZbinningNum+iZ2].second;
			      asymmetries::binning=iBin;
			      asymmetries::kinBin1=iM1;
			      asymmetries::kinBin2=iM2;
			      asymmetries::kinBin3=iZ1;
			      asymmetries::kinBin4=iZ2;
			      asymmetries::chargeType1=PN;
			      asymmetries::chargeType2=PN;
			      asymmetries::particleType1=PiPi;
			      asymmetries::particleType2=PiPi;
			      asymmetries::asymmetry=mY[iZ2];
			      asymmetries::asError=mYErr[iZ2];
			      asymmetries::meanKinVal1=yValsFour[iBin][iM1][iM2][iZ1][iZ2]/yValsNumFour[iBin][iM1][iM2][iZ1][iZ2];
			      asymmetries::meanKinVal2=mX[iZ2];
			      asymmetries::meanKinVal3=aValsFour[iBin][iM1][iM2][iZ1][iZ2]/aValsNumFour[iBin][iM1][iM2][iZ1][iZ2];
			      asymmetries::meanKinVal4=bValsFour[iBin][iM1][iM2][iZ1][iZ2]/bValsNumFour[iBin][iM1][iM2][iZ1][iZ2];
			      asymmetries::chi2Fit=0;
			      asymmetries::ndfFit=0;
			      asymmetries::expNr=expNumber;
			      asymmetries::isOnRes=onResonance;
			      as_treeFour.Fill();
	
			    }
			}
		    }
		}
	      asymIndexFour++;
	    }
	}
    }
  


  for(int iBin=zBinning;iBin<=mBinning;iBin++)
    {
      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	//      for(int iCh=PN;iCh<=ZZ;iCh++)
	{
	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	      /*	  if(!(iCh1==PN || iCh1==PZ || iCh1==ZN || iCh1==PNNP || iCh1==ZNNZ))
		continue;
		if(!(iCh2==PN || iCh2==PZ || iCh2==ZN || iCh2==PNNP || iCh2==ZNNZ))
		continue;*/
	      //we want to have all possible combinations
	      /*	 not quite sure 24.1.09
		if(!(iCh1==PN || iCh1==PZ || iCh1==ZN || iCh1==NP || iCh1==ZP || iCh1==NZ))
		continue;
		if(!(iCh2==PN || iCh2==PZ || iCh2==ZN || iCh2==NP || iCh2==ZP||  iCh2==NZ))
		continue;*/
	      if(!(iCh1==PN || iCh1==PZ || iCh1==ZN))
		continue;
	      if(!(iCh2==PN || iCh2==PZ || iCh2==ZN))
		continue;
	      strCharge=iToStrCharge(iCh1)+"_"+iToStrCharge(iCh2);
	      for(int iPa1=PiPi;iPa1<=LAST_PARTICLE_COMB;iPa1++)
		//for(int iPa=PiPi;iPa<=UNKNOWN;iPa++)
		{
		  for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
		    {

		      strParticle=iToStrPart(iPa1)+"_"+iToStrPart(iPa2);
		      switch(iBin)
			{
			case zBinning:
			  {
			    vecBinning.push_back(binningZ[iPa1].size());
			    for(int iZ=0;iZ<binningZ[iPa1].size();iZ++)
			      {
				stringstream str;
				str <<"_z1_from_";
				if(iZ==0)
				  str << 0;
				else
				  str << (binningZ[iPa1])[iZ-1];
				str << "_to_" << binningZ[iPa1][iZ];


				for(int iZ2=0;iZ2<binningZ[iPa1].size();iZ2++)
				  {
				    mX[iZ2]=xVals[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]/xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    mXErr[iZ2]=0;
				    mY[iZ2]=(*zAsyms[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].first;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    //			    cout <<"my: " << mY[iZ2] <<endl;
				    mYKinCorr[iZ2]=(*zAsymsMod2[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].first;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    //			    cout <<"myCorr: " << mYKinCorr[iZ2] <<endl;
				    mY_Pr[iZ2]=0; //for projections
				    mYErr[iZ2]=(*zAsyms[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].second;///kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];

				    asymmetries::binning=iBin;
				    asymmetries::kinBin1=iZ;
				    asymmetries::kinBin2=iZ2;
				    asymmetries::chargeType1=iCh1;
				    asymmetries::chargeType2=iCh2;
				    asymmetries::particleType1=iPa1;
				    asymmetries::particleType2=iPa2;
				    asymmetries::asymmetry=mY[iZ2];
				    asymmetries::asError=mYErr[iZ2];
				    asymmetries::meanKinVal1=yVals[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]/yValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    asymmetries::meanKinVal2=mX[iZ2];
				    asymmetries::chi2Fit=(*v_chi2Z[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2];
				    asymmetries::ndfFit=(*v_ndfsZ[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2];
				    asymmetries::expNr=expNumber;
				    asymmetries::isOnRes=onResonance;
				    as_tree.Fill();
				    mYErrKinCorr[iZ2]=(*zAsymsMod2[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].second;///kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    mYErr_Pr[iZ2]=(*zAsyms[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].second/errFact;//(kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]t)*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];

				    outputFile<<"projected error (z): " <<mYErr_Pr[iZ2] <<endl;


#ifdef MC
				    //there is a special case for pnnp. There the particle type does not have to match....
				    if(numFalseIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]+numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]!=0)
				      {
					purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]=numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]/(float)(numFalseIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]+numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]);
					puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]=getError(numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2],numFalseIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]);
				      }
				    cout <<"purity for " << strCharge+"_" + strParticle + "_" +str.str() << " z bin  : " << iZ2 << ": " <<purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]<<endl;
#else
#ifdef WITH_PURITIES
				    if(purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]>0)
				      {
					mY[iZ2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					mYErr[iZ2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					mYKinCorr[iZ2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					mYErrKinCorr[iZ2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					mYErr_Pr[iZ2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					outputFile <<"using purity of " << purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2] << " for "  << strCharge+"_" + strParticle + "_" +str.str() << " z bin  : " << iZ2;
					outputFile <<" purity Error: " << puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					outputFile <<" and correction fact of "  << kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]/(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]<<endl;
					outputFile <<"resulting Asym: " <<  mY[iZ2];
					outputFile << " error: " << mYErr[iZ2] <<endl;
				      }
				    else
				      {
					outputFile << " purity for " << strCharge+"_" + strParticle + "_" +str.str() << " m bin  : " << iZ2 << " is " <<purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]<<endl;
				      }
#endif
#endif
				  }
				graphTitles.push_back(strCharge+"_" + strParticle + "_" +str.str());
				graphTitlesKinCorr.push_back(strCharge+"_" + strParticle + "_" +str.str()+"KinCorr");
				graphTitlesProj.push_back(strCharge+"_" + strParticle + "_" +str.str()+"Proj");
				graphTitlesPurities.push_back(strCharge+"_"+strParticle+"_"+str.str()+"Purity");
				TGraphErrors* pTG=new TGraphErrors(binningZ[iPa1].size(),mX,mY,mXErr,mYErr);
				//TGraphErrors* pTG=new TGraphErrors(binningZ[iPa1].size(),mX,mYKinCorr,mXErr,mYErrKinCorr);
				TGraphErrors* pTGKinCorr=new TGraphErrors(binningZ[iPa1].size(),mX,mYKinCorr,mXErr,mYErrKinCorr);
				TGraphErrors* pTG_Proj=new TGraphErrors(binningZ[iPa1].size(),mX,mY_Pr,mXErr,mYErr_Pr);
				TGraphErrors* pGP=new TGraphErrors(binningZ[iPa1].size(),mX,purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ],mXErr,puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ]);

				pGP->GetYaxis()->SetRangeUser(0,1);
				pGP->GetYaxis()->SetNdivisions(6,true);
				pGP->GetXaxis()->SetNdivisions(6,true);
				pGP->SetMarkerStyle(8); 
				pGP->SetMarkerSize(1.2);
				pGP->SetMarkerColor(kRed);
				/*			pGP->GetYaxis()->SetLabelSize(0.07);
				  pGP->GetYaxis()->SetTitleSize(0.07);
				  pGP->GetXaxis()->SetTitleSize(0.07);
				  pGP->GetXaxis()->SetLabelSize(0.07);*/
				pGP->GetYaxis()->SetTitle("Purity");
				pGP->GetXaxis()->SetTitle("z");
				vecPurities.push_back(pGP);
				if(iCh1==PN && iPa1==PiPi&&iCh2==PN && iPa2==PiPi)
				  tgPurPiPi_PN_Z=pGP;
				if(iCh1==PZ && iPa1==PiPi&&iCh2==PZ && iPa2==PiPi)
				  tgPurPiPi_PZ_Z=pGP;
				if(iCh1==ZN && iPa1==PiPi&&iCh2==ZN && iPa2==PiPi)
				  tgPurPiPi_ZN_Z=pGP;
				if(iCh1==PN && iPa1==PiK && iCh2==PN && iPa2==PiK)
				  tgPurPiK_PN_Z=pGP;

				vecGraphs.push_back(pTG);
				vecGraphsKinCorr.push_back(pTGKinCorr);
				pTG_Proj->GetYaxis()->SetLabelSize(0.07);
				pTG_Proj->GetYaxis()->SetTitleSize(0.07);
				pTG_Proj->GetXaxis()->SetTitleSize(0.07);
				pTG_Proj->GetXaxis()->SetLabelSize(0.07);
				pTG_Proj->GetYaxis()->SetNdivisions(6,true);
				pTG_Proj->GetXaxis()->SetNdivisions(6,true);
				pTG_Proj->SetMarkerStyle(8); 
				pTG_Proj->SetMarkerSize(1.2);
				pTG_Proj->SetMarkerColor(kRed);
				pTG_Proj->GetYaxis()->SetTitle("A_{(#Phi_{1}+#Phi_{2})}");
				pTG_Proj->GetXaxis()->SetTitle("z");
				vecGraphsProjections.push_back(pTG_Proj);
				//			pTG_Proj->SaveAs((strCharge+" " + strParticle + " " +str.str())+".C").c_str());
			      }
			    asymIndexZ++;
			    break;
			  }

			case mBinning:
			  {
			    vecBinning.push_back(binningM[iPa1].size());
			    for(int iM=0;iM<binningM[iPa1].size();iM++)
			      {
				stringstream str;
				str <<"m1_from_";
				if(iM==0)
				  str << 0;
				else
				  str << binningM[iPa1][iM-1];
				str << "_to_" << binningM[iPa1][iM];
				for(int iM2=0;iM2<binningM[iPa1].size();iM2++)
				  {
				    mX[iM2]=xVals[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]/xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    mXErr[iM2]=0;
				    mY[iM2]=(*mAsyms[asymIndexM])[iM*binningM[iPa1].size()+iM2].first;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];

				    mYKinCorr[iM2]=(*mAsymsMod2[asymIndexM])[iM*binningM[iPa1].size()+iM2].first;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    mY_Pr[iM2]=0;
				    mYErr[iM2]=(*mAsyms[asymIndexM])[iM*binningM[iPa1].size()+iM2].second;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];

				    asymmetries::binning=iBin;
				    asymmetries::kinBin1=iM;
				    asymmetries::kinBin2=iM2;
				    asymmetries::chargeType1=iCh1;
				    asymmetries::chargeType2=iCh2;
				    asymmetries::particleType1=iPa1;
				    asymmetries::particleType2=iPa2;
				    asymmetries::asymmetry=mY[iM2];
				    asymmetries::asError=mYErr[iM2];
				    asymmetries::meanKinVal1=yVals[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]/yValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    asymmetries::meanKinVal2=mX[iM2];
				    asymmetries::chi2Fit=(*v_chi2M[asymIndexM])[iM*binningM[iPa1].size()+iM2];
				    asymmetries::ndfFit=(*v_ndfsM[asymIndexM])[iM*binningM[iPa1].size()+iM2];
				    asymmetries::expNr=expNumber;
				    asymmetries::isOnRes=onResonance;
				    as_tree.Fill();
				    mYErrKinCorr[iM2]=(*mAsymsMod2[asymIndexM])[iM*binningM[iPa1].size()+iM2].second;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    mYErr_Pr[iM2]=(*mAsyms[asymIndexM])[iM*binningM[iPa1].size()+iM2].second/errFact;///(kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*errFact)*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];;
				    outputFile<<"projected error (m): " <<mYErr_Pr[iM2] <<endl;
				    cout <<"myErr: " << mYErr[iM2] <<" myErr+Pr: " << mYErr_Pr[iM2] << " errFact: " << errFact <<endl;
#ifdef MC
				    if(numFalseIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]+numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]!=0)
				      {
					purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]=numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]/(float)(numFalseIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]+numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]);
					puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]=getError(numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2],numFalseIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]);
				      }
				    cout <<"purity for " << strCharge+"_" + strParticle + "_" +str.str() << " m bin  : " << iM2 << ": " <<purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]<<endl;
#else
#ifdef WITH_PURITIES
				    if(purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]>0)
				      {
					mY[iM2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
					mYErr[iM2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
					mYKinCorr[iM2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
					mYErrKinCorr[iM2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
					mYErr_Pr[iM2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
					outputFile <<"using purity of " << purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2] << " for "  << strCharge+"_" + strParticle + "_" +str.str() << " m bin  : " << iM2;
					outputFile <<" purity Error: " << puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
					outputFile <<" and correction fact of "  << kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]/(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]<<endl;
					outputFile <<"resulting Asym: " <<  mY[iM2];
					outputFile << " error: " << mYErr[iM2] <<endl;
				      }
				    else
				      {
					outputFile << " purity for " << strCharge+"_" + strParticle + "_" +str.str() << " m bin  : " << iM2 << " is " <<purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]<<endl;
				      }
#endif

#endif
				  }
				graphTitles.push_back(strCharge+"_" + strParticle + "_" +str.str());
				graphTitlesKinCorr.push_back(strCharge+"_" + strParticle + "_" +str.str()+"KinCorr");
				graphTitlesProj.push_back(strCharge+"_" + strParticle + "_" +str.str()+ "Proj");
				graphTitlesPurities.push_back(strCharge+"_"+strParticle+"_"+str.str()+"Purity");
				TGraphErrors* pTG=new TGraphErrors(binningM[iPa1].size(),mX,mY,mXErr,mYErr);
				//			TGraphErrors* pTG=new TGraphErrors(binningM[iPa1].size(),mX,mYKinCorr,mXErr,mYErr);
				TGraphErrors* pTGKinCorr=new TGraphErrors(binningM[iPa1].size(),mX,mYKinCorr,mXErr,mYErrKinCorr);
				TGraphErrors* pGP=new TGraphErrors(binningM[iPa1].size(),mX,purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM],mXErr,puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM]);


				pGP->GetYaxis()->SetRangeUser(0,1);
				pGP->GetYaxis()->SetLabelSize(0.07);
				pGP->GetYaxis()->SetTitleSize(0.07);
				pGP->GetXaxis()->SetTitleSize(0.07);
				pGP->GetXaxis()->SetLabelSize(0.07);
				pGP->GetYaxis()->SetNdivisions(6,true);
				pGP->GetXaxis()->SetNdivisions(6,true);
				pGP->SetMarkerStyle(8); 
				pGP->SetMarkerSize(1.2);
				pGP->SetMarkerColor(kRed);
				pGP->GetYaxis()->SetTitle("Purity");
				pGP->GetXaxis()->SetTitle("M_{Inv} [GeV]");
				/*			pTG->GetYaxis()->SetRangeUser(0,1);
				  pTG->GetYaxis()->SetLabelSize(0.07);
				  pTG->GetYaxis()->SetTitleSize(0.07);
				  pTG->GetXaxis()->SetTitleSize(0.07);
				  pTG->GetXaxis()->SetLabelSize(0.07);
				  pTG->GetYaxis()->SetNdivisions(6,true);
				  pTG->GetXaxis()->SetNdivisions(6,true);
				  pTG->SetMarkerStyle(8); 
				  pTG->SetMarkerSize(1.2);
				  pTG->SetMarkerColor(kRed);
				  pTG->GetYaxis()->SetTitle("A");
				  pTG->GetXaxis()->SetTitle("M_{Inv} [GeV]");*/

				vecPurities.push_back(pGP);
				if(iCh1==PN && iPa1==PiPi&& iCh2==PN && iPa2==PiPi)
				  tgPurPiPi_PN_M=pGP;
				if(iCh1==PZ && iPa1==PiPi&& iCh2==PZ && iPa2==PiPi)
				  tgPurPiPi_PZ_M=pGP;
				if(iCh1==ZN && iPa1==PiPi&& iCh2==ZN && iPa2==PiPi)
				  tgPurPiPi_ZN_M=pGP;
				if(iCh1==PN && iPa1==PiK && iCh2==PN && iPa2==PiK)
				  tgPurPiK_PN_M=pGP;
				vecGraphs.push_back(pTG);
				vecGraphsKinCorr.push_back(pTGKinCorr);
				for(int li=0;li<binningM[iPa1].size();li++)
				  {
				    cout << "myErr_Pr["<<li<<"]= " << mYErr_Pr[li] <<endl;
				  }
				TGraphErrors* pTG_Proj=new TGraphErrors(binningM[iPa1].size(),mX,mY_Pr,mXErr,mYErr_Pr);
				/*			pTG_Proj->GetYaxis()->SetLabelSize(0.07);
				  pTG_Proj->GetYaxis()->SetTitleSize(0.07);
				  pTG_Proj->GetXaxis()->SetTitleSize(0.07);
				  pTG_Proj->GetXaxis()->SetLabelSize(0.07);*/
				pTG_Proj->GetYaxis()->SetNdivisions(6,true);
				pTG_Proj->GetXaxis()->SetNdivisions(6,true);
				pTG_Proj->SetMarkerStyle(8); 
				pTG_Proj->SetMarkerSize(1.2);
				pTG_Proj->SetMarkerColor(kRed);
				pTG_Proj->GetYaxis()->SetTitle("A_{(#Phi_{1}+#Phi_{2})}");
				pTG_Proj->GetXaxis()->SetTitle("M_{Inv} [GeV]");
				vecGraphsProjections.push_back(pTG_Proj);
				//			pTG_Proj->SaveAs((strCharge+" " + strParticle + " " +str.str())+".C").c_str());
			      }
			    asymIndexM++;
			    break;
			  }
			default:
			  cout << "binning not recognized" << endl;
			  exit(0);
			}
		    }
		}
	    }
	}
    }
  asFile.Write();
  asFile.Close();


  TCanvas c1("c1","c1",10,20,500,200);
  c1.Divide(1,2);
  c1.cd(1);
  hZ1.Draw();
  c1.cd(2);
  hZ2.Draw();
  c1.SaveAs("ZCanvas.ps");
  TCanvas c2("c2","c2",10,20,500,200);

  c2.Divide(1,2);
  c2.cd(1);
  hR.Draw();
  c2.cd(2);
  hAsDiff.Draw();
  c2.SaveAs("AsDiffs.ps");

  TH2D twoDZ("twoDZ","twoDZ",10,0,10,5,0,5);

  for(int i=0;i<5;i++)
    {
      for(int j=0;j<5;j++)
	{
	  stringstream ss;
	  stringstream ss2;
	  ss <<"Z1_m"<<i <<"_m"<<j;
	  ss2<<"zSpectra/";
	  ss2 <<"Z1_m"<<i <<"_m"<<j<<".pdf";
	  TCanvas c1(ss.str().c_str(),ss.str().c_str(),10,20,500,200);
	  c1.Divide(1,2);
	  c1.cd(1);
	  zSpectH[i*5+j].first->Draw();
	  c2.cd(2);
	  zSpectH[i*5+j].second->Draw();
	  c1.SaveAs(ss2.str().c_str());
	  twoDZ.SetBinContent(2*i,j,zSpectH[i*5+j].first->GetMean());
	  twoDZ.SetBinContent(2*i+1,j,zSpectH[i*5+j].second->GetMean());
	}
    }

  TCanvas cTDZ("ctwoDZ","ctwoDZ",10,20,500,200);
  twoDZ.Draw("col");
  cTDZ.SaveAs("zSpectra/twoDZ.pdf");


  TCanvas cPhiR;
  hPhiR.Draw();
  cPhiR.SaveAs("PhiR.C");

  TCanvas cPhiZ;
  hPhiZ.Draw();
  cPhiZ.SaveAs("PhiZ.C");

  vector<TH1D*> v_thetaDiff;
  v_thetaDiff.push_back(&hThetaDiff11);
  v_thetaDiff.push_back(&hThetaDiff12);
  v_thetaDiff.push_back(&hThetaDiff21);
  v_thetaDiff.push_back(&hThetaDiff22);
  saveAs("thetaDiff.ps",v_thetaDiff);
  vector<TH1D*> v_thrustRes;
  v_thrustRes.push_back(&hThrustThetaRes);
  saveAs("thrustRes.ps",v_thrustRes);
  saveAs("AsKinCorr.ps",vecGraphsKinCorr,graphTitlesKinCorr,outputFile, vecBinning);
  saveAs("Asymmetries.ps",vecGraphs,graphTitles,outputFile,vecBinning);
  outputFile << "saving " << vecGraphsProjections.size() <<" projections having " << graphTitlesProj.size() <<" titles"<<endl;
  saveAs("Pr.C",vecGraphsProjections,graphTitlesProj,outputFile,vecBinning);
  saveAs("Pr.ps",vecGraphsProjections,graphTitlesProj,outputFile,vecBinning);
  outputFile << "done "<<endl;
#ifdef WITH_PURITIES
  saveAs(".C",vecPurities,graphTitlesPurities, outputFile,vecBinning);
  saveAs(".ps",vecPurities,graphTitlesPurities, outputFile, vecBinning);
#endif
#ifdef MC
  for(int iBin=zBinning;iBin<=mBinning;iBin++)
    {
      purities::binning=iBin;
      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	//      for(int iCh=PN;iCh<=ZZ;iCh++)
	{
	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	      /*	  if(!(iCh1==PN || iCh1==PZ || iCh1==ZN || iCh1==PNNP || iCh1==ZNNZ))
		continue;
		if(!(iCh2==PN || iCh2==PZ || iCh2==ZN || iCh2==PNNP || iCh2==ZNNZ))
		continue;*/

	      if(!(iCh1==PN || iCh1==PZ || iCh1==ZN))
		continue;
	      if(!(iCh2==PN || iCh2==PZ || iCh2==ZN))
		continue;
	      strCharge=iToStrCharge(iCh1)+"_"+iToStrCharge(iCh2);
	      purities::chargeType1=iCh1;
	      purities::chargeType2=iCh2;
	      for(int iPa1=PiPi;iPa1<=LAST_PARTICLE_COMB;iPa1++)
		//for(int iPa=PiPi;iPa<=UNKNOWN;iPa++)
		{
		  for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
		    {
		      strParticle=iToStrPart(iPa1)+"_"+iToStrPart(iPa2);
		      purities::particleType1=iPa1;
		      purities::particleType2=iPa2;
		      switch(iBin)
			{
			case zBinning:
			  {

			    for(int iZ=0;iZ<binningZ[iPa1].size();iZ++)
			      {
				for(int iZ2=0;iZ2<binningZ[iPa1].size();iZ2++)
				  {
				    purities::kinBin1=iZ;
				    purities::kinBin2=iZ2;
				    purities::purity=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    purities::purityError=puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
				    m_tree.Fill();
				  }
			      }
			    break;
			  }
			case mBinning:
			  {
			    for(int iM=0;iM<binningM[iPa1].size();iM++)
			      {
				for(int iM2=0;iM2<binningM[iPa1].size();iM2++)
				  {
				    purities::kinBin1=iM;
				    purities::kinBin2=iM2;
				    purities::purity=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    purities::purityError=puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    m_tree.Fill();
				  }
			      }
			    break;
			  }
			default:
			  cout << "last binning not recognized" << iBin<<endl;
			  exit(0);
			}
		    }
		}
	    }
	}
    }
  m_file.Write();
  m_file.Close();
#endif
#ifndef WITH_KIN_CORR
  cout <<"writing kins" <<endl;
  TFile fKinematics("kinematics.root","recreate");
  TTree m_treeKin("KinematicsTree","KinematicsTree");

  m_treeKin.Branch("zBin1",&kinematics::zBin1,"zBin1/I");
  m_treeKin.Branch("zBin2",&kinematics::zBin2,"zBin2/I");
  m_treeKin.Branch("mBin1",&kinematics::mBin1,"mBin1/I");
  m_treeKin.Branch("mBin2",&kinematics::mBin2,"mBin2/I");
  m_treeKin.Branch("chargeType1",&kinematics::chargeType1,"chargeType1/I");
  m_treeKin.Branch("particleType1",&kinematics::particleType1,"particleType1/I");
  m_treeKin.Branch("chargeType2",&kinematics::chargeType2,"chargeType2/I");
  m_treeKin.Branch("particleType2",&kinematics::particleType2,"particleType2/I");
  m_treeKin.Branch("counts",&kinematics::counts,"counts/I");
  for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
    {
      kinematics::chargeType1=iCh1;
      for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	{
	  if(!(iCh1==PN || iCh1==PZ || iCh1==ZN))
	    continue;
	  if(!(iCh2==PN || iCh2==PZ || iCh2==ZN))
	    continue;
	  kinematics::chargeType2=iCh2;
	  for(int iPa1=PiPi;iPa1<=LAST_PARTICLE_COMB;iPa1++)
	    {
	      kinematics::particleType1=iPa1;
	      for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
		{
		  kinematics::particleType2=iPa2;
		  for(int iZ=0;iZ<binningZSpect[iPa1].size();iZ++)
		    {
		      kinematics::zBin1=iZ;
		      for(int iZ2=0;iZ2<binningZSpect[iPa1].size();iZ2++)
			{
			  kinematics::zBin2=iZ2;
			  for(int iM=0;iM<binningM[iPa1].size();iM++)
			    {
			      kinematics::mBin1=iM;
			      for(int iM2=0;iM2<binningM[iPa1].size();iM2++)
				{
				  //			      cout <<"about to write for ich1: " << iCh1 <<" ich2: " << iCh2 <<" ipa1: " << iPa1 <<" ipa2: " << iPa2 << " iz1: " << iZ << " iZ2:  " << iZ2 <<" im1: " << iM <<" im2: " << iM2 <<endl;
				  //		      cout <<"writing for ich1: " << iCh1 <<" ich2: " << iCh2 <<" ipa1: " << iPa1 <<" ipa2: " << iPa2 << " iz1: " << iZ << " iZ2:  " << iZ2 <<" im1: " << iM <<" im2: " << iM2 <<" : " << kinSpectra[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2][iM][iM2]<<endl;
				  kinematics::mBin2=iM2;
				  kinematics::counts=kinSpectra[ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2][iM][iM2];
				  m_treeKin.Fill();
				}
			    }
			}
		    }
		}
	    }
	}
    }
  fKinematics.Write();
  fKinematics.Close();




#endif
  TCanvas cMPiPi_PN("cMPiPi_PN","M_PiPi_PN",10,20,200,500);
  massPiPi_PN->Draw();
  cMPiPi_PN.SaveAs("cMPiPi_PN.C");
  TCanvas cMPiPi_PZ("cMPiPi_PZ","M_PiPi_PZ",10,20,200,500);
  massPiPi_PZ->Draw();
  cMPiPi_PZ.SaveAs("cMPiPi_PZ.C");
  TCanvas cMPiPi_ZN("cMPiPi_ZN","M_PiPi_ZN",10,20,200,500);
  massPiPi_ZN->Draw();
  cMPiPi_ZN.SaveAs("cMPiPi_ZN.C");
  TCanvas cMPiK_PN("cMPiK_PN","M_PiK_PN",10,20,200,500);
  massPiK_PN->Draw();
  cMPiK_PN.SaveAs("cMPiK_PN.C");

  TCanvas cPPiPi_PN_Z("cPPiPi_PN_Z","P_PiPi_PN_Z",10,20,200,500);
  tgPurPiPi_PN_Z->Draw();
  tgPurPiPi_PN_Z->GetYaxis()->SetRangeUser(0,1);
  cPPiPi_PN_Z.SaveAs("purityPiPi_PN_Z.C");
  cPPiPi_PN_Z.SaveAs("purityPiPi_PN_Z.ps");
  TCanvas cPPiPi_PN_M("cPPiPi_PN_M","P_PiPi_PN_M",10,20,200,500);
  tgPurPiPi_PN_M->Draw();
  tgPurPiPi_PN_M->GetYaxis()->SetRangeUser(0,1);
  cPPiPi_PN_M.SaveAs("purityPiPi_PN_M.C");
  cPPiPi_PN_M.SaveAs("purityPiPi_PN_M.ps");

  TCanvas cPPiPi_PZ_Z("cPPiPi_PZ_Z","P_PiPi_PZ_Z",10,20,200,500);
  tgPurPiPi_PZ_Z->Draw();
  tgPurPiPi_PZ_Z->GetYaxis()->SetRangeUser(0,1);
  cPPiPi_PZ_Z.SaveAs("purityPiPi_PZ_Z.C");
  cPPiPi_PZ_Z.SaveAs("purityPiPi_PZ_Z.ps");
  TCanvas cPPiPi_PZ_M("cPPiPi_PZ_M","P_PiPi_PZ_M",10,20,200,500);
  tgPurPiPi_PZ_M->Draw();
  tgPurPiPi_PZ_M->GetYaxis()->SetRangeUser(0,1);
  cPPiPi_PZ_M.SaveAs("purityPiPi_PZ_M.C");
  cPPiPi_PZ_M.SaveAs("purityPiPi_PZ_M.ps");

  TCanvas cPPiPi_ZN_Z("cPPiPi_ZN_Z","P_PiPi_ZN_Z",10,20,200,500);
  tgPurPiPi_ZN_Z->Draw();
  tgPurPiPi_ZN_Z->GetYaxis()->SetRangeUser(0,1);
  cPPiPi_ZN_Z.SaveAs("purityPiPi_ZN_Z.C");
  cPPiPi_ZN_Z.SaveAs("purityPiPi_ZN_Z.ps");
  TCanvas cPPiPi_ZN_M("cPPiPi_ZN_M","P_PiPi_ZN_M",10,20,200,500);
  tgPurPiPi_ZN_M->Draw();
  tgPurPiPi_ZN_M->GetYaxis()->SetRangeUser(0,1);
  cPPiPi_ZN_M.SaveAs("purityPiPi_ZN_M.C");
  cPPiPi_ZN_M.SaveAs("purityPiPi_ZN_M.ps");

  TCanvas cPPiK_PN_Z("cPPiK_PN_Z","P_PiK_PN_Z",10,20,200,500);
  tgPurPiK_PN_Z->Draw();
  tgPurPiK_PN_Z->GetYaxis()->SetRangeUser(0,1);
  cPPiK_PN_Z.SaveAs("purityPiK_PN_Z.C");
  cPPiK_PN_Z.SaveAs("purityPiK_PN_Z.ps");
  TCanvas cPPiK_PN_M("cPPiK_PN_M","P_PiK_PN_M",10,20,200,500);
  tgPurPiK_PN_M->Draw();
  tgPurPiK_PN_M->GetYaxis()->SetRangeUser(0,1);
  cPPiK_PN_M.SaveAs("purityPiK_PN_M.C");
  cPPiK_PN_M.SaveAs("purityPiK_PN_M.ps");

  cout <<"cuttedEv: " << cuttedEv <<" nanEv: " << nanEv << " accEv: " << accEv <<endl;
}

