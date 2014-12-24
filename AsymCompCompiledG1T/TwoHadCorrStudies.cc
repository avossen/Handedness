 #define MC
//#define PAW_STYLE
//#define WITH_PURITIES
//#define WITH_HIGHER_HARMONICS


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


//to compare errors from the fits with the naive statistical errors sqrt(N)/N
#include "TwoHadAsymsCommons.h"
//#define MAX_EVENTS 10000

#define AMP_PHIR 0.05

#include "Instances.h"

//modulation of the counts with the mc angle
float getInc(float phiR, double mass1, double mass2,  double z1, double z2,int pt1, int pt2, int charge1,int charge2)
{
  //  cout <<"getinc: " << 1+AMP_PHIR*cos(phiR)<<", phir: " << phiR<<endl;
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
  int kinBin1,kinBin2;
  int chargeType1;
  int particleType1;
  int chargeType2;
  int particleType2;
  int asymmetryCounter;
  float asymmetry[4];
  float asError[4];
  float meanKinVal1;
  float meanKinVal2;
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

  TChain* chAll=new TChain("DataTree");
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
  loadBinning(binningM, binningZ);
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
      maxKin=binningM[PiPi].size();
    }
  else
    {
      maxKin=binningZ[PiPi].size();
    }
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
  int****** counts=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
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
  float***** purities=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** puritiesErrors=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  //the sin2(theta) of thrust and beam axis
  float***** kinCorrFact=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);

#ifdef MC
  unsigned long***** numCorrIdent=allocateArray<unsigned long>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  unsigned long***** numFalseIdent=allocateArray<unsigned long>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
#else
#endif 
  int runCount=0;
  int numEvtsPerRun=500000; //split in runs with a specific number of events to check for correlations
  // int  numEvtsPerRun=nevents-1;
  while((runCount+1)*numEvtsPerRun< nevents)
    {
      cout <<"runCount: " << runCount <<endl;
      for(long i=0;i<numEvtsPerRun;i++)
	{
	  if(!(i%100000))
	    cout <<"processing event nr " << i << " of " << nevents << "(" << 100*i/(float)nevents<< "% )"<<endl;

	  chAll->GetEntry(runCount*numEvtsPerRun+i);

#ifdef MAX_EVENTS
	  if(i>MAX_EVENTS)
	    break;
#endif

	  for(int iHadsInEv=0;iHadsInEv<z1Counter;iHadsInEv++)
	    {
	      int chargeIndx=ind(chargeType1[iHadsInEv],chargeType2[iHadsInEv],NumCharge);
	      int mIndex=ind(particleType1[iHadsInEv],particleType2[iHadsInEv],NumParticle);	  

	      float zValue1=((float*)(memLocsF[z1Loc]))[iHadsInEv];
	      float mValue1=((float*)(memLocsF[mass1Loc]))[iHadsInEv];
	      float zValue2=((float*)(memLocsF[z2Loc]))[iHadsInEv];
	      float mValue2=((float*)(memLocsF[mass2Loc]))[iHadsInEv];
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
	      //	      float value3=((float*)(memLocsF[phiRSumLoc]))[iHadsInEv];//coult be that the loc ist not correct anymore
	      float value3=phiRSum[iHadsInEv];
	      if(isnan(value3))
		{
		  cout <<"value 3 nan" <<endl;
		  nanEv++;
		  continue;
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

#ifdef MC
			if(particleType1[iHadsInEv]==particleType1_mc[iHadsInEv] && chargeType1[iHadsInEv]==chargeType1_mc[iHadsInEv]&& particleType2[iHadsInEv]==particleType2_mc[iHadsInEv] && chargeType2[iHadsInEv]==chargeType2_mc[iHadsInEv])
			  incBorders(binningZ[particleType1[iHadsInEv]],value1,value2,numCorrIdent[iBin][chargeIndx][mIndex]);
			else
			  incBorders(binningZ[particleType1[iHadsInEv]],value1,value2,numFalseIdent[iBin][chargeIndx][mIndex]);
#endif
			//			incBorders(binningZ[particleType1[iHadsInEv]],value1,value2,kinCorrFact[iBin][chargeIndx][mIndex],0.25*sin(thetaEThrust)*sin(thetaEThrust));
			fillBorders(binningZ[particleType1[iHadsInEv]],binningAng,value1,value2,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);


#ifdef MC
			if(phiRSum_mc[iHadsInEv]!=-1 && phiRSum_mc[iHadsInEv]!=0 && phiRSum[iHadsInEv]!=-1)
			  {
			    hR.Fill(phiRSum_mc[iHadsInEv]-phiRSum[iHadsInEv]);
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
				//no 5x5 binning anymore..
				//				zSpectH[mbin1*5+mbin2].first->Fill(value1);
				//				zSpectH[mbin1*5+mbin2].second->Fill(value2);

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
			//			float value3=((float*)(memLocsF[phiRSumLoc]))[iHadsInEv];
			value3=phiRSum[iHadsInEv];
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
			//			incBorders(binningM[particleType1[iHadsInEv]],value1,value2,kinCorrFact[iBin][chargeIndx][mIndex],0.25*sin(thetaEThrust)*sin(thetaEThrust));
			fillBorders(binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);



#ifdef MC		    
			if(phiRSum_mc[iHadsInEv]!=-1&& phiRSum_mc[iHadsInEv]!=0&& phiRSum[iHadsInEv]!=-1)
			  {
			    fillBorders(binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,countsF[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex],getInc(phiRSum_mc[iHadsInEv],mass1_mc[iHadsInEv],mass2_mc[iHadsInEv],z1_mc[iHadsInEv], z2_mc[iHadsInEv],particleType1_mc[iHadsInEv],particleType2_mc[iHadsInEv],chargeType1_mc[iHadsInEv],chargeType2_mc[iHadsInEv]));
			  }
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
      vector<vector<pair<float,float> >* > zAsymsMod;
      vector<vector<pair<float,float> >* >mAsymsMod;
      vector<vector<pair<float,float> >* > zAsymsMod2;
      vector<vector<pair<float,float> >* >mAsymsMod2;


      vector<vector<pair<float,float> >* > zAsymsPar2;
      vector<vector<pair<float,float> >* >mAsymsPar2;
      vector<vector<pair<float,float> >* > zAsymsModPar2;
      vector<vector<pair<float,float> >* >mAsymsModPar2;
      vector<vector<pair<float,float> >* > zAsymsMod2Par2;
      vector<vector<pair<float,float> >* >mAsymsMod2Par2;

      vector<vector<pair<float,float> >* > zAsymsPar3;
      vector<vector<pair<float,float> >* >mAsymsPar3;
      vector<vector<pair<float,float> >* > zAsymsModPar3;
      vector<vector<pair<float,float> >* >mAsymsModPar3;
      vector<vector<pair<float,float> >* > zAsymsMod2Par3;
      vector<vector<pair<float,float> >* >mAsymsMod2Par3;

      vector<vector<pair<float,float> >* > zAsymsPar4;
      vector<vector<pair<float,float> >* >mAsymsPar4;
      vector<vector<pair<float,float> >* > zAsymsModPar4;
      vector<vector<pair<float,float> >* >mAsymsModPar4;
      vector<vector<pair<float,float> >* > zAsymsMod2Par4;
      vector<vector<pair<float,float> >* >mAsymsMod2Par4;

      vector<vector<pair<float,float> >* >* allZAsyms[]={&zAsyms, &zAsymsPar2, &zAsymsPar3, &zAsymsPar4};
      vector<vector<pair<float,float> >* >* allMAsyms[]={&mAsyms, &mAsymsPar2, &mAsymsPar3, &mAsymsPar4};
      vector<vector<pair<float,float> >* >* allZAsymsMod[]={&zAsymsMod, &zAsymsModPar2, &zAsymsModPar3, &zAsymsModPar4};
      vector<vector<pair<float,float> >* >* allMAsymsMod[]={&mAsymsMod, &mAsymsModPar2, &mAsymsModPar3, &mAsymsModPar4};
      vector<vector<pair<float,float> >* >* allZAsymsMod2[]={&zAsymsMod2, &zAsymsMod2Par2, &zAsymsMod2Par3, &zAsymsMod2Par4};
      vector<vector<pair<float,float> >* >* allMAsymsMod2[]={&mAsymsMod2, &mAsymsMod2Par2, &mAsymsMod2Par3, &mAsymsMod2Par4};


      int numBins;
      string strCharge;
      string strParticle;
      vector<pair<float ,float> > * tmpAsym;
      vector<pair<float ,float> > * tmpAsymMod;
      vector<pair<float ,float> > * tmpAsymMod2;
      vector<vector<float>* > v_chi2Z;
      vector<vector<int>* > v_ndfsZ;
      vector<vector<float>* > v_chi2M;
      vector<vector<int>* > v_ndfsM;


      int lastPar=1;
#ifdef WITH_HIGHER_HARMONICS
      lastPar=4;
#endif

      for(int iPar=1;iPar<=lastPar;iPar++)
	{
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
			  if(iPa1!=PiPi||iPa2!=PiPi)
			    continue;
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

				
				tmpAsym=fitTheSh__(counts[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2,t_ndfs,iPar);

				v_chi2Z.push_back(t_chi2);
				v_ndfsZ.push_back(t_ndfs);
				//				zAsyms.push_back(tmpAsym);
				allZAsyms[iPar-1]->push_back(tmpAsym);
				tmpAsymMod2=fitTheSh__(countsF2[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins,errorVglFile,iPar);
				//				zAsymsMod2.push_back(tmpAsymMod2);
				allZAsymsMod2[iPar-1]->push_back(tmpAsymMod2);
#ifdef MC
				tmpAsymMod=fitTheSh__(countsF[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile,iPar);
				//				zAsymsMod.push_back(tmpAsymMod);
				allZAsymsMod[iPar-1]->push_back(tmpAsymMod);
#endif
				cout <<"size of asy vec: " << tmpAsym->size()<<endl;
				for(int k=0;k<tmpAsym->size();k++)
				  {
				    outFile  <<"Bin " << k << " Asymmetry: " << (*tmpAsym)[k].first << " +- " <<  (*tmpAsym)[k].second <<endl;
#ifdef MC
				    outFile  <<"Bin " << k << " AsymmetryMod: " << (*tmpAsymMod)[k].first << " +- " <<  (*tmpAsymMod)[k].second <<endl;
#endif
				  }
				cout <<"wrote to file " <<endl;
				break;
			      }
			    case mBinning:
			      {
				outFile << strCharge <<" " << strParticle << " M_inv Binning: " <<endl;
				numBins=binningM[iPa1].size();
				tmpAsym=fitTheSh__(counts[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2,t_ndfs,iPar);
				allMAsyms[iPar-1]->push_back(tmpAsym);
				//				mAsyms.push_back(tmpAsym);
				v_chi2M.push_back(t_chi2);
				v_ndfsM.push_back(t_ndfs);

				tmpAsymMod2=fitTheSh__(countsF2[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins,errorVglFile,iPar);
				allMAsymsMod2[iPar-1]->push_back(tmpAsymMod2);
				//				mAsymsMod2.push_back(tmpAsymMod2);
#ifdef MC
				tmpAsymMod=fitTheSh__(countsF[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins,errorVglFile,iPar);
				//				mAsymsMod.push_back(tmpAsymMod);
				allMAsymsMod[iPar-1]->push_back(tmpAsymMod);
#endif
				for(int k=0;k<tmpAsym->size();k++)
				  {
				    outFile  <<"Bin " << k << " Asymmetry: " << (*tmpAsym)[k].first << " +- " <<  (*tmpAsym)[k].second <<endl;
#ifdef MC
				    outFile  <<"Bin " << k << " AsymmetryMod: " << (*tmpAsymMod)[k].first << " +- " <<  (*tmpAsymMod)[k].second <<endl;
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
      for(int i=0;i<allZAsyms[0]->size();i++)
	{
	  for(int j=0;j<(*(*allZAsyms[0])[i]).size();j++)
	    {
	      if(!(*allZAsyms[0])[i]||(*(*allZAsyms[0])[i])[j].second > 3) //no counts
		continue;
	      asymSum+=(*(*allZAsyms[0])[i])[j].first;
	      asymSumErr+=(*(*allZAsyms[0])[i])[j].second*(*(*allZAsyms[0])[i])[j].second;

	      outputFile << "adding Asymmetry (z): " << (*(*allZAsyms[0])[i])[j].first <<" +- " << (*(*allZAsyms[0])[i])[j].second <<endl;
	      aCounter++;
	    }
	}
      cout <<"zAsymAvg: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
      outputFile <<"zAsymAvg: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
      outputFile <<"num Asyms: " << aCounter <<endl;
      aCounter=0;
      asymSum=0;

      for(int i=0;i<(*allMAsyms[0]).size();i++)
	{
	  for(int j=0;j<(*(*allMAsyms[0])[i]).size();j++)
	    {
	      if((*(*allMAsyms[0])[i])[j].second > 3) //no counts
		continue;
	      asymSum+=(*(*allMAsyms[0])[i])[j].first;
	      asymSumErr+=(*(*allMAsyms[0])[i])[j].second*(*(*allMAsyms[0])[i])[j].second;
	      outputFile << "adding Asymmetry (m): " << (*(*allMAsyms[0])[i])[j].first <<" +- " << (*(*allMAsyms[0])[i])[j].second <<endl;
	      aCounter++;
	    }
	}
      cout <<"mAsymAvg: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
      outputFile <<"mAsymAvg: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
      outputFile <<"num Asyms: " << aCounter <<endl;

#ifdef MC
      asymSum=0;
      aCounter=0;
      for(int i=0;i<(*allZAsymsMod[0]).size();i++)
	{
	  for(int j=0;j<(*(*allZAsymsMod[0])[i]).size();j++)
	    {
	      if((*(*allZAsymsMod[0])[i])[j].second > 3) //no counts
		continue;
	      asymSum+=(*(*allZAsymsMod[0])[i])[j].first;
	      asymSumErr+=(*(*allZAsymsMod[0])[i])[j].second*(*(*allZAsymsMod[0])[i])[j].second;
	      hAsDiff.Fill(((*(*allZAsymsMod[0])[i])[j].first-AMP_PHIR)/(*(*allZAsymsMod[0])[i])[j].second);
	      outputFile << "adding AsymmetryMod: " << (*(*allZAsymsMod[0])[i])[j].first <<" +- " << (*(*allZAsymsMod[0])[i])[j].second <<endl;
	      aCounter++;
	    }
	}
      cout <<"zAsymAvgMod: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
      outputFile <<"zAsymAvgMod: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
      outputFile <<"num ModAsyms: " << aCounter <<endl;
      aCounter=0;
      asymSum=0;

      for(int i=0;i<(*allMAsymsMod[0]).size();i++)
	{
	  for(int j=0;j<(*(*allMAsymsMod[0])[i]).size();j++)
	    {
	      if((*(*allMAsymsMod[0])[i])[j].second > 3) //no counts
		continue;
	      asymSum+=(*(*allMAsymsMod[0])[i])[j].first;
	      asymSumErr+=(*(*allMAsymsMod[0])[i])[j].second*(*(*allMAsymsMod[0])[i])[j].second;
	      hAsDiff.Fill(((*(*allMAsymsMod[0])[i])[j].first-AMP_PHIR)/(*(*allMAsymsMod[0])[i])[j].second);
	      outputFile << "adding AsymmetryMod: " << (*(*allMAsymsMod[0])[i])[j].first <<" +- " << (*(*allMAsymsMod[0])[i])[j].second <<endl;
	      aCounter++;
	    }
	}
      cout <<"mAsymAvgMod: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
      outputFile <<"mAsymAvgMod: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
      outputFile <<"num ModAsyms: " << aCounter <<endl;
#endif

      asymSum=0;
      aCounter=0;
      for(int i=0;i<(*allZAsymsMod2[0]).size();i++)
	{
	  for(int j=0;j<(*(*allZAsymsMod2[0])[i]).size();j++)
	    {
	      if((*(*allZAsymsMod2[0])[i])[j].second > 3) //no counts
		continue;
	      asymSum+=(*(*allZAsymsMod2[0])[i])[j].first;
	      asymSumErr+=(*(*allZAsymsMod2[0])[i])[j].second*(*(*allZAsymsMod2[0])[i])[j].second;
	      hAsDiff.Fill(((*(*allZAsymsMod2[0])[i])[j].first-AMP_PHIR)/(*(*allZAsymsMod2[0])[i])[j].second);
	      outputFile << "adding AsymmetryMod: " << (*(*allZAsymsMod2[0])[i])[j].first <<" +- " << (*(*allZAsymsMod2[0])[i])[j].second <<endl;
	      aCounter++;
	    }
	}
      cout <<"zAsymAvgMod2: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
      outputFile <<"zAsymAvgMod2: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
      outputFile <<"num ModAsyms2: " << aCounter <<endl;
      aCounter=0;
      asymSum=0;

      for(int i=0;i<(*allMAsymsMod[0]).size();i++)
	{
	  for(int j=0;j<(*(*allMAsymsMod[0])[i]).size();j++)
	    {
	      if((*(*allMAsymsMod[0])[i])[j].second > 3) //no counts
		continue;
	      asymSum+=(*(*allMAsymsMod[0])[i])[j].first;
	      asymSumErr+=(*(*allMAsymsMod[0])[i])[j].second*(*(*allMAsymsMod[0])[i])[j].second;
	      hAsDiff.Fill(((*(*allMAsymsMod[0])[i])[j].first-AMP_PHIR)/(*(*allMAsymsMod[0])[i])[j].second);
	      outputFile << "adding AsymmetryMod: " << (*(*allMAsymsMod[0])[i])[j].first <<" +- " << (*(*allMAsymsMod[0])[i])[j].second <<endl;
	      aCounter++;
	    }
	}
      cout <<"mAsymAvgMod2: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
      outputFile <<"mAsymAvgMod2: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
      outputFile <<"num ModAsyms2: " << aCounter <<endl;
      saveKinPars(kinSpectraH, binningM);
      saveKinPars(kinSpectra1H, binningM);
      ///-----end mod2
      float inpAsym=AMP_PHIR;
      stringstream sstr;
      sstr <<"AsymmetriesCorrStudies/"<<folderName+"_"<<runCount<<".root";
      TFile asFile(sstr.str().c_str(),"recreate");
      TTree as_tree("AsymmetriesTree","AsymmetriesTree");
      as_tree.Branch("kinBin1",&asymmetries::kinBin1,"kinBin1/I");
      as_tree.Branch("kinBin2",&asymmetries::kinBin2,"kinBin2/I");
      as_tree.Branch("binning",&asymmetries::binning,"binning/I");
      as_tree.Branch("chargeType1",&asymmetries::chargeType1,"chargeType1/I");
      as_tree.Branch("particleType1",&asymmetries::particleType1,"particleType1/I");
      as_tree.Branch("chargeType2",&asymmetries::chargeType2,"chargeType2/I");
      as_tree.Branch("particleType2",&asymmetries::particleType2,"particleType2/I");
      as_tree.Branch("asymmetryCounter",&asymmetries::asymmetryCounter,"asymmetryCounter/I");
      as_tree.Branch("asymmetry",asymmetries::asymmetry,"asymmetry[asymmetryCounter]/F");
      as_tree.Branch("asError",asymmetries::asError,"asError[asymmetryCounter]/F");
      as_tree.Branch("meanKinVal1",&asymmetries::meanKinVal1,"meanKinVal1/F");
      as_tree.Branch("meanKinVal2",&asymmetries::meanKinVal2,"meanKinVal2/F");
      as_tree.Branch("meanTheta",&asymmetries::meanTheta,"meanTheta/F");
      as_tree.Branch("chi2Fit",&asymmetries::chi2Fit,"chi2Fit/F");
      as_tree.Branch("ndfFit",&asymmetries::ndfFit,"ndfFit/I");
      //not necessary to save for every value but... what the heck
      as_tree.Branch("expNr",&asymmetries::expNr,"expNr/I");
      as_tree.Branch("isOnRes",&asymmetries::isOnRes,"isOnRes/I");
      as_tree.Branch("RunNum",&runCount,"RunNum/I");
      as_tree.Branch("InputAsymmetry",&inpAsym,"InputAsymmetry/F");

      for(int iBin=zBinning;iBin<=mBinning;iBin++)
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
		    //for(int iPa=PiPi;iPa<=UNKNOWN;iPa++)
		    {
		      for(int iPa2=PiPi;iPa2<=iPa1;iPa2++)
			{
			  if(iPa1!=PiPi || iPa2!=PiPi)
			    continue;
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
					mY[iZ2]=(*(*allZAsyms[0])[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].first;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					//			    cout <<"my: " << mY[iZ2] <<endl;
					mYKinCorr[iZ2]=(*(*allZAsymsMod2[0])[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].first;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					//			    cout <<"myCorr: " << mYKinCorr[iZ2] <<endl;
					mY_Pr[iZ2]=0; //for projections
					mYErr[iZ2]=(*(*allZAsyms[0])[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].second;///kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					//changed this to the mod2, to have the weighted in there
					for(int iPar=1;iPar<=lastPar;iPar++)
					  {
					    asymmetries::asymmetry[iPar-1]=(*(*allZAsymsMod[iPar-1])[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].first;
					    asymmetries::asError[iPar-1]=(*(*allZAsymsMod[iPar-1])[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].second;
					  }
					asymmetries::asymmetryCounter=lastPar;
					asymmetries::binning=iBin;
					asymmetries::kinBin1=iZ;
					asymmetries::kinBin2=iZ2;
					asymmetries::chargeType1=iCh1;
					asymmetries::chargeType2=iCh2;
					asymmetries::particleType1=iPa1;
					asymmetries::particleType2=iPa2;
					//					asymmetries::asymmetry=mY[iZ2];
					//					asymmetries::asError=mYErr[iZ2];
					asymmetries::meanKinVal1=yVals[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]/yValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					asymmetries::meanKinVal2=mX[iZ2];
					asymmetries::chi2Fit=(*v_chi2Z[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2];
					asymmetries::ndfFit=(*v_ndfsZ[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2];
					asymmetries::expNr=expNumber;
					asymmetries::isOnRes=onResonance;
					as_tree.Fill();
					mYErrKinCorr[iZ2]=(*(*allZAsymsMod2[0])[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].second;///kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					mYErr_Pr[iZ2]=(*(*allZAsyms[0])[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2].second/errFact;//(kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]t)*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];

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
					mY[iM2]=(*(*allMAsyms[0])[asymIndexM])[iM*binningM[iPa1].size()+iM2].first;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];

					mYKinCorr[iM2]=(*(*allMAsymsMod2[0])[asymIndexM])[iM*binningM[iPa1].size()+iM2].first;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
					mY_Pr[iM2]=0;
					mYErr[iM2]=(*(*allMAsyms[0])[asymIndexM])[iM*binningM[iPa1].size()+iM2].second;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
					for(int iPar=1;iPar<=lastPar;iPar++)
					  {
					    asymmetries::asymmetry[iPar-1]=(*(*allMAsymsMod[iPar-1])[asymIndexM])[iM*binningM[iPa1].size()+iM2].first;
					    asymmetries::asError[iPar-1]=(*(*allMAsymsMod[iPar-1])[asymIndexM])[iM*binningM[iPa1].size()+iM2].second;
					  }
					asymmetries::binning=iBin;
					asymmetries::kinBin1=iM;
					asymmetries::kinBin2=iM2;
					asymmetries::chargeType1=iCh1;
					asymmetries::chargeType2=iCh2;
					asymmetries::particleType1=iPa1;
					asymmetries::particleType2=iPa2;
					//					asymmetries::asymmetry=mY[iM2];
					//					asymmetries::asError=mYErr[iM2];
					asymmetries::meanKinVal1=yVals[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]/yValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
					asymmetries::meanKinVal2=mX[iM2];
					asymmetries::chi2Fit=(*v_chi2M[asymIndexM])[iM*binningM[iPa1].size()+iM2];
					asymmetries::ndfFit=(*v_ndfsM[asymIndexM])[iM*binningM[iPa1].size()+iM2];
					asymmetries::expNr=expNumber;
					asymmetries::isOnRes=onResonance;
					as_tree.Fill();
					mYErrKinCorr[iM2]=(*(*allMAsymsMod2[0])[asymIndexM])[iM*binningM[iPa1].size()+iM2].second;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
					mYErr_Pr[iM2]=(*(*allMAsyms[0])[asymIndexM])[iM*binningM[iPa1].size()+iM2].second/errFact;///(kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*errFact)*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];;
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
      runCount++;
    }
}

