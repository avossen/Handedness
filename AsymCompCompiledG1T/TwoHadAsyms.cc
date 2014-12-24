//#define MC

//#define WITH_PURITIES  //have to generate purities.root new, if binning is changed
//#define PURITIES_ONLY  //might not generate the plots

//#define WITH_HIGHER_HARMONICS

//had before

#define D0Mass 1.865
//essentially resolution, just a guess
#define D0Width 0.6
#define D0Lower  D0Mass-D0Width/2


//going after the resolution of the az angle...
//#define MAX_THRUST_THETA_LAB 2.4
//#define MIN_THRUST_THETA_LAB 0.7

#define MAX_THRUST_THETA_LAB 2.0
#define MIN_THRUST_THETA_LAB 1.0

//#define PHIR1_MINUS_PHIR2//seems to have signal
#define TWO_TIMES_PHIR1_MINUS_PHIR2
#define MIN_LONG_PROJ 0.3

#define D0Upper  D0Mass+D0Width/2


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

#include <set>

using namespace std;
ofstream* errorVglFile;
vector<float>* binningM;
vector<float>* binningZ;
vector<float>* binningZSpect;
vector<float> binningThetaPiPi;


//to compare errors from the fits with the naive statistical errors sqrt(N)/N
#include "TwoHadAsymsCommons.h"
//#define MAX_EVENTS 100000
#define AMP_PHIR 0.1
//#define AMP_PHIR 0.0

//only starting at 0.9 ok
#define OPENING_CUT 0.8
#define JET_CUT_ETA 0.7
#define JET_CUT_PHI 0.7
#define JET_CUT_R 0.5

#define MAX_OPENING_CUT 2.0

#include "Instances.h"


//modulation of the counts with the mc angle
float getInc(float phiR, double mass1, double mass2,  double z1, double z2,int pt1, int pt2, int charge1,int charge2)
{
  //  cout <<"getinc: " << 1+AMP_PHIR*cos(phiR)<<endl;
  return 1+AMP_PHIR*cos(phiR);
}


float getIncT(float phiR, float z1, float z2, float mass1, float mass2,float kinVal)
{
  if(z1> 1 || z2 > 1 || mass1 > 3 || mass2 > 3)
    return 0.0;
  //  cout <<"kinVal: " << kinVal << " phiR: " << phiR <<endl;
  //  cout << "returning: " << 1+kinVal*AMP_PHIR*cos(phiR) <<endl;
  //linear in  z1, z2, mass1, mass2 with 1*z1 and 0.5 mass
  //  cout <<"z1: " << z1 << " z2: " << z2 << " m1: " << mass1 << " m2: " << mass2 << endl;
  return 1+z1*z2*mass1*mass2*4*kinVal*AMP_PHIR*cos(phiR);
}


float getInc2(int binning, int index1, int index2, int mBin1, int mBin2, float zVal1, float zVal2)
{
  //  return 1;//test
  if(mBin1<0 || mBin2 <0 || zVal1 < 0 || zVal2 <0)
    {
      cout <<"ret 0 "<< mBin1 <<" b2: " << mBin2 <<" zval1: " << zVal1 <<" val2: " << zVal2<<endl;
      return 0;
    }

  double val=0;
  switch(binning)
    {
    case zBinning:
      break;

    case mBinning:
      break;
    }

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
  float asymmetry;
  float asError;
  float meanKinVal1;
  float meanKinVal2;
  float meanTheta;
  float meanThetaSin;
  float meanThetaCos;
  float chi2Fit;
  int ndfFit;
  int expNr;
  int isOnRes;
}


int main(int argc, char** argv)
{
  long overallCorrIdentPN=0; //for the paper, average purity for PiPi PN
  long overallFalseIdentPN=0;

  long numPiPiPN_before_cut=0;
  long numPiPiPN_after_cut=0;


  long angTot=0;
  long gr06=0;
  long gr07=0;
  long gr08=0;

  long ihads=0;
  long ihadsC=0;

  gROOT->SetStyle("Plain");
  long cuttedEv=0;
  long accEv=0;
  long nanEv=0;
  //see if the cuts for charm enrichted or depleted samples work
  long  numCharmEv=0;
  long allEvts=0;
  long keptEvts=0;
  long nonCharmEv=0;

  long allCombs=0;
  long keptCombs=0;

  vector<TH1D*> massHistos;
  TH1D* massPiPi_PN=new TH1D("massPiPi_PN","Invariant Mass #pi^{+}/#pi^{-} pairs",1000,0,3);
  TH1D* massPiPi_PZ=new TH1D("massPiPi_PZ","Invariant Mass #pi^{+}/#pi^{0} pairs",1000,0,3);
  TH1D* massPiPi_ZN=new TH1D("massPiPi_ZN","Invariant Mass #pi^{0}/#pi^{-} pairs",1000,0,3);
  TH1D* massPiK_PN=new TH1D("massPiK_PN","Invariant Mass #pi^{+}/K^{-} pairs",1000,0,3);

  TH1D* thrustProj_PN=new TH1D("thrustProjPN","thrustProjPN",1000,-2,2);




  //  TH1D* hPhiComp=new TH1D("PhiR_Diff","PhiR_Diff",500,-1,1);

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
  //gStyle->SetOptTitle(kFALSE);
  //  gStyle->SetOptStat(kFALSE);
  gStyle->SetPalette(1,0);

  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptTitle(0);

  gROOT->ForceStyle();
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadBottomMargin(0.14);

  gStyle->SetLabelOffset(0.00,"x");
  gStyle->SetTitleOffset(1.,"x");

  gStyle->SetMarkerSize(1.3);
  gStyle->SetLineColor(1);

  gStyle->SetLabelSize(0.06);
  gStyle->SetTitleXSize(0.07);

  gStyle->SetTitleYSize(0.07);

  gStyle->SetTextSize(0.07);

 
  ofstream outputFile("outputFile");
  ofstream xcheckFile("xcheck");
  errorVglFile=new ofstream("errorVgl");
  //continuum luminosities
  int numLumi=24;
  float luminositiesCont[]={0.594,0.,1.211,1.203,1.402,0.853,3.562,0.,1.416,1.671,3.745,2.393,2.722,   //sum: 67.771
			    1.944,6.078,6.315,5.657,6.524,2.315,3.438,2.586,4.825,0.,7.821};//luminosities for exp 7-55 continuum, 

  float luminositiesOnRes[]={5.928,4.44,8.132,10.739,12.682,11.181,24.953,4.375,6.266,25.741,25.427,17.827,17.619,16.733,61.658,43.639,59.989,56.989,13.048,
			     37.577,27.293,38.935,0.0,73.514}; //604,685
  float ZCut_for_mPlot=0.2;

  float errFact=0;
  for(int i=0;i<numLumi-1;i++)
    {
      errFact+=luminositiesCont[i];
      errFact+=luminositiesOnRes[i];
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
  TH1D hZ1("hZ1","hZ1",200,0,1);
  TH1D hZ2("hZ2","hZ2",200,0,1);
  TH1D hZCharm("hZCharm","hZCharm",200,0,5);

  TH1D hAsDiff("AsDiff","AsDiff",200,-10,10);

  TH1D hThetaDiff11("ThetaDiff11","ThetaDiff11",200,-4,4);
  TH1D hThetaDiff12("ThetaDiff12","ThetaDiff12",200,-4,4);
  TH1D hThetaDiff21("ThetaDiff21","ThetaDiff21",200,-4,4);
  TH1D hThetaDiff22("ThetaDiff22","ThetaDiff22",200,-4,4);

  TH1D hThrustThetaRes("Thrust_Rec-MC","Thrust_Rec-MC",200,-3,3);
  TH1I hNumD0FirstHemi("numD0FirstHemi","numD0FirstHemi",10,0,9);
  TH1I hNumD0SecondHemi("numD0SecondHemi","numD0SecondHemi",10,0,9);

  TH1D hR("R_diff","R_diff",500,-0.5,0.5);

  //  const float zCut=1.2;//is also the max for the z binning. Rought estimate, should follow from resolution of z
  const float zCut=1.0;//is also the max for the z binning. Rought estimate, should follow from resolution of z
  const float mCut=1000;
  const float thetaCut=1000;

  vector<string> fieldNamesF;
  vector<string> fieldNamesI;

  binningM=new vector<float>[NumParticle];
  binningZ=new vector<float>[NumParticle];
  //  binningMix1=new vector<float>[NumParticle];
  //  binningMix2=new vector<float>[NumParticle];
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
  binningThetaPiPi.push_back(0.2);
  binningThetaPiPi.push_back(0.4);
  binningThetaPiPi.push_back(0.5);
  binningThetaPiPi.push_back(0.6);
  binningThetaPiPi.push_back(0.7);
  binningThetaPiPi.push_back(0.8);
  binningThetaPiPi.push_back(0.9);
  binningThetaPiPi.push_back(1.0);
  binningThetaPiPi.push_back(1.2);

  float** thetaMeanRatio=allocateArray<float>(binningM[PiPi].size(),binningM[PiPi].size());//mean of ratios (see aurores mail from June 29)
  float** thetaMeanSin=allocateArray<float>(binningM[PiPi].size(),binningM[PiPi].size());//ratio of means
  float** thetaMeanCos=allocateArray<float>(binningM[PiPi].size(),binningM[PiPi].size());//ratio of means
  float** thetaCount=allocateArray<float>(binningM[PiPi].size(),binningM[PiPi].size());//ratio of means

  float** thetaMeanRatioZ=allocateArray<float>(binningZ[PiPi].size(),binningZ[PiPi].size());//mean of ratios (see aurores mail from June 29)
  float** thetaMeanSinZ=allocateArray<float>(binningZ[PiPi].size(),binningZ[PiPi].size());//ratio of means
  float** thetaMeanCosZ=allocateArray<float>(binningZ[PiPi].size(),binningZ[PiPi].size());//ratio of means
  float** thetaCountZ=allocateArray<float>(binningZ[PiPi].size(),binningZ[PiPi].size());//ratio of means

  float** thetaMeanRatioMix=allocateArray<float>(binningZ[PiPi].size(),binningM[PiPi].size());//mean of ratios (see aurores mail from June 29)
  float** thetaMeanSinMix=allocateArray<float>(binningZ[PiPi].size(),binningM[PiPi].size());//ratio of means
  float** thetaMeanCosMix=allocateArray<float>(binningZ[PiPi].size(),binningM[PiPi].size());//ratio of means
  float** thetaCountMix=allocateArray<float>(binningZ[PiPi].size(),binningM[PiPi].size());//ratio of means

  float** meanZforThetaMix=allocateArray<float>(binningZ[PiPi].size(),binningM[PiPi].size());
  float** meanMforThetaMix=allocateArray<float>(binningZ[PiPi].size(),binningM[PiPi].size());

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

  float thrustProj11[2000];
  int thrustProj11Counter;
  float thrustProj12[2000];
  int thrustProj12Counter;
  float thrustProj21[2000];
  int thrustProj21Counter;
  float thrustProj22[2000];
  int thrustProj22Counter;

  float theta11[2000];
  int theta11Counter;
  float theta12[2000];
  int theta12Counter;
  float theta21[2000];
  int theta21Counter;
  float theta22[2000];
  int theta22Counter;



  float eta11[2000];
  int eta11Counter;
  float eta12[2000];
  int eta12Counter;
  float eta21[2000];
  int eta21Counter;
  float eta22[2000];
  int eta22Counter;

  float phi11[2000];
  int phi11Counter;
  float phi12[2000];
  int phi12Counter;
  float phi21[2000];
  int phi21Counter;
  float phi22[2000];
  int phi22Counter;

  float thrustTheta_mc;
  float thrustTheta;
  float thrustThetaLab;

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

  int chargeCounter_mc;
  int chargeType_mc[2000];

  int charge1Counter_mc;
  int chargeType1_mc[2000];
  int particle1Counter_mc;
  int particleType1_mc[2000];

  int charge2Counter_mc;
  int chargeType2_mc[2000];
  int particle2Counter_mc;
  int particleType2_mc[2000];

  float thetaEThrust=pi/(float)2; //for now until the new data is generated //correction factor is 1/4 sin^2(theta)
  float thetaEThrust_mc=pi/(float)2; //for now until the new data is generated //correction factor is 1/4 sin^2(theta)
  float phiEThrust;
  float phiThrust;

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


  fieldNamesI.push_back("eta11Counter");
  fieldNamesF.push_back("eta11");
  memLocsI.push_back(&eta11Counter);
  memLocsF.push_back(eta11);

  fieldNamesI.push_back("eta12Counter");
  fieldNamesF.push_back("eta12");
  memLocsI.push_back(&eta12Counter);
  memLocsF.push_back(eta12);

  fieldNamesI.push_back("eta21Counter");
  fieldNamesF.push_back("eta21");
  memLocsI.push_back(&eta21Counter);
  memLocsF.push_back(eta21);

  fieldNamesI.push_back("eta22Counter");
  fieldNamesF.push_back("eta22");
  memLocsI.push_back(&eta22Counter);
  memLocsF.push_back(eta22);

  fieldNamesI.push_back("phi11Counter");
  fieldNamesF.push_back("phi11");
  memLocsI.push_back(&phi11Counter);
  memLocsF.push_back(phi11);

  fieldNamesI.push_back("phi12Counter");
  fieldNamesF.push_back("phi12");
  memLocsI.push_back(&phi12Counter);
  memLocsF.push_back(phi12);

  fieldNamesI.push_back("phi21Counter");
  fieldNamesF.push_back("phi21");
  memLocsI.push_back(&phi21Counter);
  memLocsF.push_back(phi21);

  fieldNamesI.push_back("phi22Counter");
  fieldNamesF.push_back("phi22");
  memLocsI.push_back(&phi22Counter);
  memLocsF.push_back(phi22);



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
  fieldNamesF.push_back("thrustThetaLab");
  memLocsF.push_back(&thrustThetaLab);

  fieldNamesF.push_back("thrustPhi");
  memLocsF.push_back(&phiThrust);

  fieldNamesF.push_back("phiEThrust");
  memLocsF.push_back(&phiEThrust);

  fieldNamesI.push_back("z1RatioCounter");
  fieldNamesF.push_back("z1Ratio");
  memLocsI.push_back(&z1RatioCounter);
  memLocsF.push_back(z1Ratio);

  fieldNamesI.push_back("z2RatioCounter");
  fieldNamesF.push_back("z2Ratio");
  memLocsI.push_back(&z2RatioCounter);
  memLocsF.push_back(z2Ratio);

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
  fieldNamesF.push_back("thetaEThrust_mc");
  memLocsF.push_back(&thetaEThrust_mc);
  fieldNamesF.push_back("phiRSum_mc");
  memLocsF.push_back(phiRSum_mc);
  fieldNamesI.push_back("phiRSum_mcCounter");
  memLocsI.push_back(&phiRSum_mcCounter);

  fieldNamesI.push_back("chargeType1_mcCounter");
  fieldNamesI.push_back("chargeType1_mc");
  memLocsI.push_back(&charge1Counter_mc);
  memLocsI.push_back(chargeType1_mc);

  fieldNamesI.push_back("chargeType_mcCounter");
  fieldNamesI.push_back("chargeType_mc");
  memLocsI.push_back(&chargeCounter_mc);
  memLocsI.push_back(chargeType_mc);

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
  cout <<"about to branch" <<endl;
  m_tree.Branch("kinBin1" ,&purities::kinBin1,"kinBin1/I");
  m_tree.Branch("kinBin2",&purities::kinBin2,"kinBin2/I");
  m_tree.Branch("binning",&purities::binning,"binning/I");
  m_tree.Branch("chargeType1",&purities::chargeType1,"chargeType1/I");
  m_tree.Branch("particleType1",&purities::particleType1,"particleType1/I");
  m_tree.Branch("chargeType2",&purities::chargeType2,"chargeType2/I");
  m_tree.Branch("particleType2",&purities::particleType2,"particleType2/I");
  m_tree.Branch("purity",&purities::purity,"purity/F");
  //m_tree.Branch("purity",&purities::purity,"purity/F");
  m_tree.Branch("purityError",&purities::purityError,"purityError/F");
  cout <<"done" <<endl;
#endif //mc
  long nevents=0;
  chAll->Add((string(rootPath)+"/*.root").c_str());
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
  cout <<"done branching "<<endl;
  nevents=chAll->GetEntries();


  int****** counts=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  int** countsTheta=allocateArray<int>(binningThetaPiPi.size(),numAngBins);
  int** countsThetaOpeningCut=allocateArray<int>(binningThetaPiPi.size(),numAngBins);
  float** countsThetaF=allocateArray<float>(binningThetaPiPi.size(),numAngBins);
  float** countsThetaOpeningCutF=allocateArray<float>(binningThetaPiPi.size(),numAngBins);
  int** countsThetaJetCut=allocateArray<int>(binningThetaPiPi.size(),numAngBins);

  //combine hads in the same hemi --- look at dihadana to see how they are combined... might be only pn or all..
  int****** countsSameHemi=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  int****** countsDiffEvt=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  float****** countsF=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  float****** countsF2=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  int** countsZBin=allocateArray<int>(maxKin,maxKin);
  int** countsMBin=allocateArray<int>(maxKin,maxKin);


  gHCounter=0;

  //avg over value2, conforming to notation in paper
  float***** xVals=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** xValsNum=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** yVals=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** yValsNum=allocateArray<int>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** purities=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** puritiesErrors=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);

  float* thetaVal=new float[binningThetaPiPi.size()];
  int* thetaValNum=new int[binningThetaPiPi.size()];
  float* thetaValF=new float[binningThetaPiPi.size()];
  int* thetaValNumF=new int[binningThetaPiPi.size()];

  float* thetaValOpeningCut=new float[binningThetaPiPi.size()];
  int* thetaValNumOpeningCut=new int[binningThetaPiPi.size()];
  float* thetaValOpeningCutF=new float[binningThetaPiPi.size()];
  int* thetaValNumOpeningCutF=new int[binningThetaPiPi.size()];

  float* thetaValJetCut=new float[binningThetaPiPi.size()];
  int* thetaValNumJetCut=new int[binningThetaPiPi.size()];

  //to gain statistics integrate over one side
  float**** puritiesInt=allocateArray<float>(NumBin,NumCharge,NumParticle,maxKin);
  float**** puritiesErrorsInt=allocateArray<float>(NumBin,NumCharge,NumParticle,maxKin);

  //the sin2(theta) of thrust and beam axis
  float***** kinCorrFact=allocateArray<float>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);

#ifdef MC
  unsigned long***** numCorrIdent=allocateArray<unsigned long>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  unsigned long**** numCorrIdentInt=allocateArray<unsigned long>(NumBin,NumCharge,NumParticle,maxKin);

  unsigned long***** numFalseIdent=allocateArray<unsigned long>(NumBin,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  unsigned long**** numFalseIdentInt=allocateArray<unsigned long>(NumBin,NumCharge,NumParticle,maxKin);

#else

  cout <<"lastcharge: " << LAST_CHARGE_COMB<<endl;
  //  cout <<"done" <<endl;
#endif  //no mc

  cout <<"there are " <<nevents << " events"<<endl;
  for(long i=0;i<nevents;i++)
    {
      if(!(i%100000))
	cout <<"processing event nr " << i << " of " << nevents << "(" << 100*i/(float)nevents<< "% )"<<endl;
      chAll->GetEntry(i);

      float mySin=sin(thetaEThrust);
      float myCos=cos(thetaEThrust);

      float kinFact=mySin*mySin/(1+myCos*myCos);// transverse spin projection
      //and convert to longitudinal spin projection      
      kinFact=sin(acos(kinFact));
      if(kinFact<MIN_LONG_PROJ)
	continue;
      if(thrustThetaLab>MAX_THRUST_THETA_LAB || thrustThetaLab < MIN_THRUST_THETA_LAB)
	continue;

#ifdef MAX_EVENTS
      if(i>MAX_EVENTS)
	break;
#endif
      allEvts++;

      int locIHads=0;

      for(int iHadsInEv=0;iHadsInEv<z1Counter;iHadsInEv++)
	{
	  double eta11T=-log(tan((theta11[iHadsInEv])/2));
	  double eta12T=-log(tan((theta12[iHadsInEv])/2));
	  double eta21T=-log(tan((theta21[iHadsInEv])/2));
	  double eta22T=-log(tan((theta22[iHadsInEv])/2));
	  double etaT=-log(tan((thrustTheta)/2));

	  double etaL11=-log(tan((theta11[iHadsInEv]+thrustTheta)/2));
	  double etaL12=-log(tan((theta12[iHadsInEv]+thrustTheta)/2));
	  double etaL21=-log(tan((theta21[iHadsInEv]+thrustTheta)/2));
	  double etaL22=-log(tan((theta22[iHadsInEv]+thrustTheta)/2));

	  //	  cout <<"eta11: " << eta11 <<" eta12 " << eta12<<" eta21: " << eta21 << " eta22: " << eta22 <<endl;

	  keptCombs++;
	  if(particleType1[iHadsInEv]==PiPi&&chargeType1[iHadsInEv]==PN&&particleType2[iHadsInEv]==PiPi&&chargeType2[iHadsInEv]==PN )
	    {
	      locIHads++;



	      //	      cout <<"kinFact: " << kinFact <<endl;
	      fillBorders(binningThetaPiPi,binningAng,kinFact,((float*)(memLocsF[phiRSumLoc]))[iHadsInEv],countsTheta,thetaVal, thetaValNum);
#ifdef MC 
	      float kinFact_mc=sin(thetaEThrust_mc)*sin(thetaEThrust_mc)/(1+cos(thetaEThrust_mc)*cos(thetaEThrust_mc));
	      //convert in log:
	      kinFact_mc=sin(acos(kinFact_mc));

	      //	      cout <<" thrust theta is: " << thetaEThrust_mc << " reconstructed: " << thetaEThrust << " kin fact: " << kinFact_mc<<" otherone: "<< kinFact <<endl;
	      fillBorders(binningThetaPiPi,binningAng,kinFact,((float*)(memLocsF[phiRSumLoc]))[iHadsInEv],countsThetaF, thetaValF,thetaValNumF,getIncT(phiRSum_mc[iHadsInEv],z1_mc[iHadsInEv], z2_mc[iHadsInEv], mass1_mc[iHadsInEv], mass2_mc[iHadsInEv], kinFact_mc));
#endif 
	    }
	  
	  if(fabs(thrustProj11[iHadsInEv])>MAX_OPENING_CUT || fabs(thrustProj12[iHadsInEv])>MAX_OPENING_CUT ||fabs(thrustProj21[iHadsInEv])>MAX_OPENING_CUT ||fabs(thrustProj22[iHadsInEv])>MAX_OPENING_CUT)
	    {
	      cout <<"opening cut" <<endl;
	      continue;
	    }

	  float r11=sqrt((etaT-eta11[iHadsInEv])*(etaT-eta11[iHadsInEv])+(phiThrust-phi11[iHadsInEv])*(phiThrust-phi11[iHadsInEv]));
	  float r12=sqrt((etaT-eta12[iHadsInEv])*(etaT-eta12[iHadsInEv])+(phiThrust-phi12[iHadsInEv])*(phiThrust-phi12[iHadsInEv]));
	  float r21=sqrt(((-1)*etaT-eta21[iHadsInEv])*((-1)*etaT-eta21[iHadsInEv])+(phiThrust-phi21[iHadsInEv]+pi)*(phiThrust-phi21[iHadsInEv]+pi));
	  float r22=sqrt(((-1)*etaT-eta22[iHadsInEv])*((-1)*etaT-eta22[iHadsInEv])+(phiThrust-phi22[iHadsInEv]+pi)*(phiThrust-phi22[iHadsInEv]+pi));
	  // cout <<"r11: " << r11 << " r12: " << r12 << " r21 : " << r21 << " r22: " << r22 <<endl;
	  if((r11 < JET_CUT_R)  && (r12 < JET_CUT_R) && (r21 < JET_CUT_R) &&  (r22 < JET_CUT_R))
	    {
	      //    cout <<"jet cut ok..." <<endl;
	  /*	  if(fabs(etaT-eta11[iHadsInEv])<JET_CUT_ETA || fabs(etaT-eta12[iHadsInEv])<JET_CUT_ETA ||fabs(etaT-eta21[iHadsInEv])<JET_CUT_ETA ||fabs(etaT-eta22[iHadsInEv])<JET_CUT_ETA)
	    {
	      if(fabs(phiThrust-phi11[iHadsInEv])<JET_CUT_PHI || fabs(phiThrust-phi12[iHadsInEv])<JET_CUT_PHI ||fabs(phiThrust-phi21[iHadsInEv])<JET_CUT_PHI ||fabs(phiThrust-phi22[iHadsInEv])<JET_CUT_PHI)
		{
	  */


	      //	      cout <<"angTot: " << angTot << " gr06: " << gr06 << " gr 07: " << gr07 <<" gr 08: " << gr08 <<endl;
			  //	      cout <<"kinFact: " << kinFact <<endl;
		  fillBorders(binningThetaPiPi,binningAng,kinFact,((float*)(memLocsF[phiRSumLoc]))[iHadsInEv],countsThetaJetCut,thetaValJetCut, thetaValNumJetCut);
		  //		}
	    }
	  
	  if(fabs(thrustProj11[iHadsInEv])<OPENING_CUT || fabs(thrustProj12[iHadsInEv])<OPENING_CUT ||fabs(thrustProj21[iHadsInEv])<OPENING_CUT ||fabs(thrustProj22[iHadsInEv])<OPENING_CUT)
	    {
	      //	      cout <<"thrust" <<endl;
	      continue;
	    }
	  if(particleType1[iHadsInEv]==PiPi&&chargeType1[iHadsInEv]==PN&&particleType2[iHadsInEv]==PiPi&&chargeType2[iHadsInEv]==PN )
	    {
	      //	      cout <<"kinFact: " << kinFact <<endl;
	      fillBorders(binningThetaPiPi,binningAng,kinFact,((float*)(memLocsF[phiRSumLoc]))[iHadsInEv],countsThetaOpeningCut,thetaValOpeningCut, thetaValNumOpeningCut);
#ifdef MC
	      float kinFact_mc=sin(thetaEThrust_mc)*sin(thetaEThrust_mc)/(1+cos(thetaEThrust_mc)*cos(thetaEThrust_mc));
	      //convert into long:
	      kinFact_mc=sin(acos(kinFact_mc));
	      fillBorders(binningThetaPiPi,binningAng,kinFact,((float*)(memLocsF[phiRSumLoc]))[iHadsInEv],countsThetaOpeningCutF,thetaValOpeningCutF,thetaValNumOpeningCutF,getIncT(phiRSum_mc[iHadsInEv],z1_mc[iHadsInEv],z2_mc[iHadsInEv], mass1_mc[iHadsInEv], mass2_mc[iHadsInEv],kinFact_mc));
#endif 

	    }
	  //since we use the 'all' collection in diHadAna, we have to set the charge type for the quadruple manually
	  if(chargeType1[iHadsInEv]==PN && chargeType2[iHadsInEv]==PN && ((particleType1[iHadsInEv]==PiK&&particleType2[iHadsInEv]==KPi)|| (particleType1[iHadsInEv]==KPi&&particleType2[iHadsInEv]==PiK)))
	    {
	      //cout <<"found PNNP.." <<endl;
	      	      chargeType[iHadsInEv]=PNNP;
	    }
	  if(chargeType[iHadsInEv]==PNNP || chargeType[iHadsInEv]==ZNNZ)
	    {
	      chargeType1[iHadsInEv]=chargeType[iHadsInEv];
	      chargeType2[iHadsInEv]=chargeType[iHadsInEv];
	      particleType1[iHadsInEv]=PiK;
	      particleType2[iHadsInEv]=PiK;
	    }
#ifdef MC
	  if(chargeType_mc[iHadsInEv]==PNNP || chargeType_mc[iHadsInEv]==ZNNZ)
	    {
	      chargeType1_mc[iHadsInEv]=chargeType_mc[iHadsInEv];
	      chargeType2_mc[iHadsInEv]=chargeType_mc[iHadsInEv];
	      particleType1_mc[iHadsInEv]=PiK;
	      particleType2_mc[iHadsInEv]=PiK;
	    }
#endif
	  int chargeIndx=ind(chargeType1[iHadsInEv],chargeType2[iHadsInEv],NumCharge);
	  int mIndex=ind(particleType1[iHadsInEv],particleType2[iHadsInEv],NumParticle);	  
	  if(((float*)(memLocsF[z1Loc]))[iHadsInEv]>ZCut_for_mPlot && ((float*)(memLocsF[z2Loc]))[iHadsInEv]>ZCut_for_mPlot)
	    {
	      if(particleType1[iHadsInEv]==PiPi&&chargeType1[iHadsInEv]==PN&&particleType2[iHadsInEv]==PiPi&&chargeType2[iHadsInEv]==PN )
		{
		  hPhiR.Fill(cos(((float*)(memLocsF[phiRSumLoc]))[iHadsInEv]));
		  hPhiZ.Fill(cos(twoPhiZero[iHadsInEv]));
		  massPiPi_PN->Fill(((float*)(memLocsF[mass1Loc]))[iHadsInEv]);
		  massPiPi_PN->Fill(((float*)(memLocsF[mass2Loc]))[iHadsInEv]);

		  thrustProj_PN->Fill(thrustProj11[iHadsInEv]);
		  thrustProj_PN->Fill(thrustProj21[iHadsInEv]);
		  thrustProj_PN->Fill(thrustProj22[iHadsInEv]);
		  thrustProj_PN->Fill(thrustProj12[iHadsInEv]);
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
	  //let's look only at the last mBins, for arguments sake in first hemi
	  //  if(mValue1< 1.2)
	  //    continue;


	  if(mValue1<mLowCuts[particleType1[iHadsInEv]] || mValue2<mLowCuts[particleType1[iHadsInEv]])
	    {
	      cuttedEv++;
	      continue;
	    }
	  if(mValue1>mHighCuts[particleType1[iHadsInEv]] || mValue2>mHighCuts[particleType1[iHadsInEv]])
	    continue;
	  
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
#ifdef PHIR1_MINUS_PHIR2
	   value3=value3-2*phiR1[iHadsInEv];
	   if(value3 > pi)
	     value3-=2*pi;
	   if(value3 < -pi)
	     value3+=2*pi;
#endif
#ifdef TWO_TIMES_PHIR1_MINUS_PHIR2
	   //	   cout <<"value3 : " << value3 <<endl;
	   value3=2*(value3-2*phiR1[iHadsInEv]); //max 6pi -epsionl
	   //	   cout <<"value3: " <<value3 << " " << phiR1[iHadsInEv] <<endl;
	   if(value3 > pi)
	     value3-=2*pi;
	   if(value3 < -pi)
	     value3+=2*pi;
	   if(value3 > pi)
	     value3-=2*pi;
	   if(value3 < -pi)
	     value3+=2*pi;
	   if(value3 > pi)
	     value3-=2*pi;
	   if(value3 < -pi)
	     value3+=2*pi;
#endif
	  if(isnan(value3))
	    {
	      cout <<"value 3 nan" <<endl;
	      nanEv++;
	      continue;
	    }

	  for(int iBin=zBinning;iBin<=mixedBinning;iBin++)
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
		    if(PiPi==particleType1[iHadsInEv]&& PiPi==particleType2[iHadsInEv] && m_chType1==PN &&  m_chType2==PN)
		      {
			countsZBin[zbin1][zbin2]++;
		      }

		    accEv++;


#ifdef MC
		    if(particleType1[iHadsInEv]==particleType1_mc[iHadsInEv] && chargeType1[iHadsInEv]==chargeType1_mc[iHadsInEv]&& particleType2[iHadsInEv]==particleType2_mc[iHadsInEv] && chargeType2[iHadsInEv]==chargeType2_mc[iHadsInEv])
		      {
			incBorders(binningZ[particleType1[iHadsInEv]],value1,value2,numCorrIdent[iBin][chargeIndx][mIndex]);
		      }
		    else
		      {
			incBorders(binningZ[particleType1[iHadsInEv]],value1,value2,numFalseIdent[iBin][chargeIndx][mIndex]);
		      }
		    if(particleType1[iHadsInEv]==PiPi&& chargeType1[iHadsInEv]==PN)
		      {
			/*			if(value1 > 0.825)
			  cout <<"first pari high z "<<endl;
			if(value2 > 0.825)
			cout <<"second pari high z "<<endl;*/
		      }

		    if(particleType1[iHadsInEv]==particleType1_mc[iHadsInEv] && chargeType1[iHadsInEv]==chargeType1_mc[iHadsInEv])
		      {
			numCorrIdentInt[iBin][chargeType1[iHadsInEv]][particleType1[iHadsInEv]][getBin(binningZ[particleType1[iHadsInEv]],value1)]++;
 			if(particleType1[iHadsInEv]==PiPi && chargeType1[iHadsInEv]==PN)
			  {
			    // cout <<"inc overall " << overallCrrIdentPN <<endl;
			    overallCorrIdentPN++;
			  }
		      }
		    else
		      {
			numFalseIdentInt[iBin][chargeType1[iHadsInEv]][particleType1[iHadsInEv]][getBin(binningZ[particleType1[iHadsInEv]],value1)]++;
			if(particleType1[iHadsInEv]==PiPi && chargeType1[iHadsInEv]==PN)
			  {
			    //			    cout <<"inc overall false " << overallFalseIdentPN <<endl;
			    overallFalseIdentPN++;
			  }
		      }

		    if(particleType2[iHadsInEv]==particleType2_mc[iHadsInEv] && chargeType2[iHadsInEv]==chargeType2_mc[iHadsInEv])
		      {
		      numCorrIdentInt[iBin][chargeType2[iHadsInEv]][particleType2[iHadsInEv]][getBin(binningZ[particleType2[iHadsInEv]],value2)]++;
 			if(particleType1[iHadsInEv]==PiPi && chargeType1[iHadsInEv]==PN)
			  overallCorrIdentPN++;
		      }
		    else
		      {
			numFalseIdentInt[iBin][chargeType2[iHadsInEv]][particleType2[iHadsInEv]][getBin(binningZ[particleType2[iHadsInEv]],value2)]++;
			if(particleType1[iHadsInEv]==PiPi && chargeType1[iHadsInEv]==PN)
			  overallFalseIdentPN++;
		      }

#endif
		    incBorders(binningZ[particleType1[iHadsInEv]],value1,value2,kinCorrFact[iBin][chargeIndx][mIndex],0.25*sin(thetaEThrust)*sin(thetaEThrust));
		    fillBorders(binningZ[particleType1[iHadsInEv]],binningAng,value1,value2,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);
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

#ifdef PHIR1_MINUS_PHIR2
	   value3=value3-2*phiR1[iHadsInEv];
	   if(value3 > pi)
	     value3-=2*pi;
	   if(value3 < -pi)
	     value3+=2*pi;
#endif
#ifdef TWO_TIMES_PHIR1_MINUS_PHIR2
	   value3=2*(value3-2*phiR1[iHadsInEv]);
	   if(value3 > pi)
	     value3-=2*pi;
	   if(value3 < -pi)
	     value3+=2*pi;
	   if(value3 > pi)
	     value3-=2*pi;
	   if(value3 < -pi)
	     value3+=2*pi;
   if(value3 > pi)
	     value3-=2*pi;
	   if(value3 < -pi)
	     value3+=2*pi;
#endif
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

		    if(PiPi==particleType1[iHadsInEv]&& PiPi==particleType2[iHadsInEv])
		      {
			countsMBin[mbin1][mbin2]++;
		      }
#ifdef MC
		    if(particleType1[iHadsInEv]==particleType1_mc[iHadsInEv] && chargeType1[iHadsInEv]==chargeType1_mc[iHadsInEv]&&particleType2[iHadsInEv]==particleType2_mc[iHadsInEv] && chargeType2[iHadsInEv]==chargeType2_mc[iHadsInEv])
		      incBorders(binningM[particleType1[iHadsInEv]],value1,value2,numCorrIdent[iBin][chargeIndx][mIndex]);
		    else
		      incBorders(binningM[particleType1[iHadsInEv]],value1,value2,numFalseIdent[iBin][chargeIndx][mIndex]);
		    if(particleType1[iHadsInEv]==particleType1_mc[iHadsInEv] && chargeType1[iHadsInEv]==chargeType1_mc[iHadsInEv])
		      numCorrIdentInt[iBin][chargeType1[iHadsInEv]][particleType1[iHadsInEv]][getBin(binningM[particleType1[iHadsInEv]],value1)]++;
		    else
		      numFalseIdentInt[iBin][chargeType1[iHadsInEv]][particleType1[iHadsInEv]][getBin(binningM[particleType1[iHadsInEv]],value1)]++;

		    if(particleType2[iHadsInEv]==particleType2_mc[iHadsInEv] && chargeType2[iHadsInEv]==chargeType2_mc[iHadsInEv])
		      numCorrIdentInt[iBin][chargeType2[iHadsInEv]][particleType2[iHadsInEv]][getBin(binningM[particleType2[iHadsInEv]],value2)]++;
		    else
		      numFalseIdentInt[iBin][chargeType2[iHadsInEv]][particleType2[iHadsInEv]][getBin(binningM[particleType2[iHadsInEv]],value2)]++;
#endif
		    incBorders(binningM[particleType1[iHadsInEv]],value1,value2,kinCorrFact[iBin][chargeIndx][mIndex],0.25*sin(thetaEThrust)*sin(thetaEThrust));
		    fillBorders(binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);

#ifdef MC		    
		    if(phiRSum_mc[iHadsInEv]!=-1)
		      fillBorders(binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,countsF[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex],getInc(phiRSum_mc[iHadsInEv],mass1_mc[iHadsInEv],mass2_mc[iHadsInEv],z1_mc[iHadsInEv], z2_mc[iHadsInEv],particleType1_mc[iHadsInEv],particleType2_mc[iHadsInEv],chargeType1_mc[iHadsInEv],chargeType2_mc[iHadsInEv]));
#endif
		    break;
		  }
		case multBinning: 
		  {
		    break;
		  }
		case mixedBinning:
		  {
		    float value1=((float*)(memLocsF[z1Loc]))[iHadsInEv];
		    //switchted...
		    float value1_sw=((float*)(memLocsF[z2Loc]))[iHadsInEv];
		    if(isnan(value1))
		      {
			cout <<"value 1 z mix nan" <<endl;
			nanEv++;
			continue;
		      }
		    float value2=((float*)(memLocsF[mass2Loc]))[iHadsInEv];
		    float value2_sw=((float*)(memLocsF[mass1Loc]))[iHadsInEv];
		    if(isnan(value2))
		      {
			cout <<"value 2 m nan" <<endl;
			nanEv++;
			continue;
		      }
		    if(value2>mCut || value1 > zCut || value1_sw> zCut || value2_sw > mCut)
		      {
			//			cout <<"z1 value: " << value1 << " m2 value: " << value2 <<endl;
			cuttedEv++;
			continue;
		      }
		    float value3=((float*)(memLocsF[phiRSumLoc]))[iHadsInEv];
#ifdef PHIR1_MINUS_PHIR2
	   value3=value3-2*phiR1[iHadsInEv];
	   if(value3 > pi)
	     value3-=2*pi;
	   if(value3 < -pi)
	     value3+=2*pi;
#endif
#ifdef TWO_TIMES_PHIR1_MINUS_PHIR2
	   value3=2*(value3-2*phiR1[iHadsInEv]);
	   if(value3 > pi)
	     value3-=2*pi;
	   if(value3 < -pi)
	     value3+=2*pi;
	   if(value3 > pi)
	     value3-=2*pi;
	   if(value3 < -pi)
	     value3+=2*pi;
   if(value3 > pi)
	     value3-=2*pi;
	   if(value3 < -pi)
	     value3+=2*pi;
#endif
		    if(isnan(value3))
		      {
			nanEv++;
			continue;
		      }
		    if(value1<zLowCuts[particleType1[iHadsInEv]] || value2<mLowCuts[particleType1[iHadsInEv]] || value1_sw < zLowCuts[particleType1[iHadsInEv]] || value2_sw < mLowCuts[particleType1[iHadsInEv]])
		      {
			cuttedEv++;
			continue;
		      }
		    if(value1>zHighCuts[particleType1[iHadsInEv]] || value2>mHighCuts[particleType1[iHadsInEv]] || value1_sw > zHighCuts[particleType1[iHadsInEv]]|| value2_sw > mHighCuts[particleType1[iHadsInEv]])
		      continue;
		    //mix
		    accEv++;

#ifdef MC
		    if(particleType1[iHadsInEv]==particleType1_mc[iHadsInEv] && chargeType1[iHadsInEv]==chargeType1_mc[iHadsInEv]&&particleType2[iHadsInEv]==particleType2_mc[iHadsInEv] && chargeType2[iHadsInEv]==chargeType2_mc[iHadsInEv])
		      {
			//			cout <<"inc mixed "<<endl;
			incBorders(binningZ[particleType1[iHadsInEv]], binningM[particleType1[iHadsInEv]],value1,value2,numCorrIdent[iBin][chargeIndx][mIndex]);
			incBorders(binningZ[particleType1[iHadsInEv]], binningM[particleType1[iHadsInEv]],value1_sw,value2_sw,numCorrIdent[iBin][chargeIndx][mIndex]);
		      }
		    else
		      {
			incBorders(binningZ[particleType1[iHadsInEv]],binningM[particleType1[iHadsInEv]],value1,value2,numFalseIdent[iBin][chargeIndx][mIndex]);
			incBorders(binningZ[particleType1[iHadsInEv]],binningM[particleType1[iHadsInEv]],value1_sw,value2_sw,numFalseIdent[iBin][chargeIndx][mIndex]);
		      }
#endif
		    //filling mixed counts
		    fillBordersMix(binningZ[particleType1[iHadsInEv]],binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);
		    fillBordersMix(binningZ[particleType1[iHadsInEv]],binningM[particleType1[iHadsInEv]],binningAng,value1_sw, value2_sw,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);

#ifdef MC		    
		    if(phiRSum_mc[iHadsInEv]!=-1)
		      {
		      fillBordersMix(binningZ[particleType1[iHadsInEv]],binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,countsF[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex],getInc(phiRSum_mc[iHadsInEv],mass1_mc[iHadsInEv],mass2_mc[iHadsInEv],z1_mc[iHadsInEv], z2_mc[iHadsInEv],particleType1_mc[iHadsInEv],particleType2_mc[iHadsInEv],chargeType1_mc[iHadsInEv],chargeType2_mc[iHadsInEv]));
		      fillBordersMix(binningZ[particleType1[iHadsInEv]],binningM[particleType1[iHadsInEv]],binningAng,value1_sw, value2_sw,value3,countsF[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex],getInc(phiRSum_mc[iHadsInEv],mass1_mc[iHadsInEv],mass2_mc[iHadsInEv],z1_mc[iHadsInEv], z2_mc[iHadsInEv],particleType1_mc[iHadsInEv],particleType2_mc[iHadsInEv],chargeType1_mc[iHadsInEv],chargeType2_mc[iHadsInEv]));
		      }
#endif
		    break;
		  }
		default:
		  cout << "binning not recognized1" << endl;
		  exit(0);
		}
	    }
	}
      ihads+=locIHads;
      ihadsC++;
    }

  cout <<"num pipi pn before opening cut: " <<  numPiPiPN_before_cut << " after cut: " << numPiPiPN_after_cut <<endl;
  float overallPurity=overallCorrIdentPN/(double)(overallFalseIdentPN+overallCorrIdentPN);
  cout <<" overallCorr: " << overallCorrIdentPN <<" false: " << overallFalseIdentPN <<endl;
  cout <<"overallPurity is " <<overallPurity <<endl;

  cout <<"mean ihads: "<< ihads << " / " << ihadsC <<", "  << (double)ihads/(double)ihadsC <<endl;


  cout <<"evt ratio: " << allEvts << " / " << keptEvts << " = " << allEvts/(double)keptEvts <<endl;
  cout <<"keptCombs: " << keptCombs <<endl;
  //and now the fit...
  cout<<"done with filling arryas, fitting..."<<endl;
  vector<vector<pair<float,float> >* > zAsyms;
  vector<vector<pair<float,float> >* >mAsyms;
  vector<vector<pair<float,float> >* > zAsymsMod;
  vector<vector<pair<float,float> >* >mAsymsMod;
  vector<vector<pair<float,float> >* > zAsymsMod2;
  vector<vector<pair<float,float> >* >mAsymsMod2;
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
  vector<vector<float>* > v_chi2Mixed;
  vector<vector<int>* > v_ndfsMixed;



  for(int iBin=zBinning;iBin<=mixedBinning;iBin++)
    {
      //      for(int iCh=PN;iCh<=ZZ;iCh++)
      cout <<"iBin: " << iBin <<endl;
      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	{
	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	      if(!(iCh1==PN || iCh1==PZ || iCh1==ZN || iCh1==PNNP || iCh1==ZNNZ))
		continue;
	      if(!(iCh2==PN || iCh2==PZ || iCh2==ZN || iCh2==PNNP || iCh2==ZNNZ))
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
			    mAsymsMod.push_back(tmpAsymMod);
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
			case multBinning:
			  {
			    break;
			  }
			case mixedBinning:
			  {
			    outFile << strCharge <<" " << strParticle << " M_inv Binning: " <<endl;
			    numBins=binningM[iPa1].size();
#ifndef PURITIES_ONLY
			    tmpAsym=fitTheSh__(counts[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins, errorVglFile, t_chi2,t_ndfs);
			    v_chi2Mixed.push_back(t_chi2);
			    v_ndfsMixed.push_back(t_ndfs);
			    mAsyms.push_back(tmpAsym);
			    tmpAsymMod2=fitTheSh__(countsF2[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins,errorVglFile);
			    mAsymsMod2.push_back(tmpAsymMod2);
#ifdef MC
			    tmpAsymMod=fitTheSh__(countsF[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins,errorVglFile);
			    mAsymsMod.push_back(tmpAsymMod);
#endif
			    for(int k=0;k<tmpAsym->size();k++)
			      {
				outFile  <<"Bin " << k << " Asymmetry: " << (*tmpAsym)[k].first << " +- " <<  (*tmpAsym)[k].second <<endl;
#ifdef MC
				outFile  <<"Bin " << k << " Asymmetry: " << (*tmpAsymMod)[k].first << " +- " <<  (*tmpAsymMod)[k].second <<endl;
#endif
			      }
			    cout <<"after fit mixed" <<endl;
#endif//only purities
			    break;
			  }
			default:
			  cout << "binning not recognized2" << endl;
			  exit(0);
			}

		    }
		}
	    }
	}
    }


  vector<float>* t_chi2T=new vector<float>;
  vector<int>* t_ndfsT=new vector<int>;

  vector<pair<float ,float> > * tmpThetaAsOpeningCut;
  vector<pair<float ,float> > * tmpThetaAs;
  vector<pair<float ,float> > * tmpThetaAsOpeningCutF;
  vector<pair<float ,float> > * tmpThetaAsF;
  vector<pair<float ,float> > * tmpThetaAsJetCut;

  tmpThetaAs=fitTheSingle(countsTheta,binningAng,binningThetaPiPi.size(), errorVglFile, t_chi2T,t_ndfsT);
  
  float yT[30];
  float xT[30];
  float yET[30];
  float xET[30];
  for(int i=0;i<binningThetaPiPi.size();i++)
    {
      xT[i]=0.0;
      xET[i]=0.0;
      yET[i]=0.0;
      yT[i]=0.0;
      if(thetaValNum[i]>0)
	{
	  xT[i]=thetaVal[i]/thetaValNum[i];
	  yT[i]=(*tmpThetaAs)[i].first;
	  yET[i]=(*tmpThetaAs)[i].second;
	}
    }
  TGraphErrors tgTheta(binningThetaPiPi.size(),xT,yT,xET,yET);
  TCanvas cThetaAs;
  tgTheta.Draw("AP*");
  cThetaAs.SaveAs("kinFactAsym.png");
  cThetaAs.Write();

  tmpThetaAsF=fitTheSingle(countsThetaF,binningAng,binningThetaPiPi.size(), errorVglFile, t_chi2T,t_ndfsT);
  
  for(int i=0;i<binningThetaPiPi.size();i++)
    {
      xT[i]=0.0;
      xET[i]=0.0;
      yET[i]=0.0;
      yT[i]=0.0;
      if(thetaValNumF[i]>0)
	{
	  xT[i]=thetaValF[i]/thetaValNumF[i];
	  yT[i]=(*tmpThetaAsF)[i].first;
	  yET[i]=(*tmpThetaAsF)[i].second;
	}
    }
  TGraphErrors tgThetaF(binningThetaPiPi.size(),xT,yT,xET,yET);
  TCanvas cThetaAsF;
  tgThetaF.Draw("AP*");
  tgThetaF.Write();

  tmpThetaAsOpeningCut=fitTheSingle(countsThetaOpeningCut,binningAng,binningThetaPiPi.size(), errorVglFile, t_chi2T,t_ndfsT);
  for(int i=0;i<binningThetaPiPi.size();i++)
    {
      xT[i]=0.0;
      xET[i]=0.0;
      yET[i]=0.0;
      yT[i]=0.0;
      if(thetaValNumOpeningCut[i]>0)
	{
	  xT[i]=thetaValOpeningCut[i]/thetaValNumOpeningCut[i];
	  yT[i]=(*tmpThetaAsOpeningCut)[i].first;
	  yET[i]=(*tmpThetaAsOpeningCut)[i].second;
	}
    }

  TGraphErrors tgThetaOC(binningThetaPiPi.size(),xT,yT,xET,yET);
  TCanvas cThetaAsOC;
  tgThetaOC.Draw("AP*");
  tgThetaOC.Write();

  tmpThetaAsOpeningCutF=fitTheSingle(countsThetaOpeningCutF,binningAng,binningThetaPiPi.size(), errorVglFile, t_chi2T,t_ndfsT);
  for(int i=0;i<binningThetaPiPi.size();i++)
    {
      xT[i]=0.0;
      xET[i]=0.0;
      yET[i]=0.0;
      yT[i]=0.0;
      if(thetaValNumOpeningCutF[i]>0)
	{
	  xT[i]=thetaValOpeningCutF[i]/thetaValNumOpeningCutF[i];
	  yT[i]=(*tmpThetaAsOpeningCutF)[i].first;
	  yET[i]=(*tmpThetaAsOpeningCutF)[i].second;
	}
    }

  TGraphErrors tgThetaOCF(binningThetaPiPi.size(),xT,yT,xET,yET);
  TCanvas cThetaAsOCF;
  tgThetaOCF.Draw("AP*");
  tgThetaOCF.Write();




tmpThetaAsJetCut=fitTheSingle(countsThetaJetCut,binningAng,binningThetaPiPi.size(), errorVglFile, t_chi2T,t_ndfsT);
  for(int i=0;i<binningThetaPiPi.size();i++)
    {
      xT[i]=0.0;
      xET[i]=0.0;
      yET[i]=0.0;
      yT[i]=0.0;
      if(thetaValNumJetCut[i]>0)
	{
	  xT[i]=thetaValJetCut[i]/thetaValNumJetCut[i];
	  yT[i]=(*tmpThetaAsJetCut)[i].first;
	  yET[i]=(*tmpThetaAsJetCut)[i].second;
	}
    }

  TGraphErrors tgThetaJC(binningThetaPiPi.size(),xT,yT,xET,yET);
  TCanvas cThetaAsJC;
  tgThetaJC.Draw("AP*");
  cThetaAsJC.SaveAs("thetaAsJetCut.png");




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
  int asymIndexMixed=0;

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

#ifndef PURITIES_ONLY
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


#endif //purities_only

#ifdef WITH_HIGHER_HARMONICS
  TFile asFile(("AsRFilesHighHarm/"+folderName+".root").c_str(),"recreate");
#else
  TFile asFile(("AsymmetriesRFilesForA/"+folderName+".root").c_str(),"recreate");
#endif

#ifndef PURITIES_ONLY

  TTree as_tree("AsymmetriesTree","AsymmetriesTree");
  as_tree.Branch("kinBin1",&asymmetries::kinBin1,"kinBin1/I");
  as_tree.Branch("kinBin2",&asymmetries::kinBin2,"kinBin2/I");
  as_tree.Branch("binning",&asymmetries::binning,"binning/I");
  as_tree.Branch("chargeType1",&asymmetries::chargeType1,"chargeType1/I");
  as_tree.Branch("particleType1",&asymmetries::particleType1,"particleType1/I");
  as_tree.Branch("chargeType2",&asymmetries::chargeType2,"chargeType2/I");
  as_tree.Branch("particleType2",&asymmetries::particleType2,"particleType2/I");
  as_tree.Branch("asymmetry",&asymmetries::asymmetry,"asymmetry/F");
  as_tree.Branch("asError",&asymmetries::asError,"asError/F");
  as_tree.Branch("meanKinVal1",&asymmetries::meanKinVal1,"meanKinVal1/F");
  as_tree.Branch("meanKinVal2",&asymmetries::meanKinVal2,"meanKinVal2/F");
  as_tree.Branch("meanTheta",&asymmetries::meanTheta,"meanTheta/F");
  //hack for for aurore
  as_tree.Branch("meanThetaSin",&asymmetries::meanThetaSin,"meanThetaSin/F");
  as_tree.Branch("meanThetaCos",&asymmetries::meanThetaCos,"meanThetaCos/F");


  as_tree.Branch("chi2Fit",&asymmetries::chi2Fit,"chi2Fit/F");
  as_tree.Branch("ndfFit",&asymmetries::ndfFit,"ndfFit/I");
  //not necessary to save for every value but... what the heck
  as_tree.Branch("expNr",&asymmetries::expNr,"expNr/I");
  as_tree.Branch("isOnRes",&asymmetries::isOnRes,"isOnRes/I");


  for(int iBin=zBinning;iBin<=mixedBinning;iBin++)
    {
      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	//      for(int iCh=PN;iCh<=ZZ;iCh++)
	{
	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	        if(!(iCh1==PN || iCh1==PZ || iCh1==ZN || iCh1==PNNP || iCh1==ZNNZ))
		continue;
		if(!(iCh2==PN || iCh2==PZ || iCh2==ZN || iCh2==PNNP || iCh2==ZNNZ))
		continue;
	      //we want to have all possible combinations
	      /*	 not quite sure 24.1.09
		if(!(iCh1==PN || iCh1==PZ || iCh1==ZN || iCh1==NP || iCh1==ZP || iCh1==NZ))
		continue;
		if(!(iCh2==PN || iCh2==PZ || iCh2==ZN || iCh2==NP || iCh2==ZP||  iCh2==NZ))
		continue;*/
		/*	      if(!(iCh1==PN || iCh1==PZ || iCh1==ZN))
		continue;
	      if(!(iCh2==PN || iCh2==PZ || iCh2==ZN))
	      continue;*/
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
				    outputFile <<" my: iz " << iZ2 << " " << mY[iZ2]<<endl;
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
				    asymmetries::ndfFit= (*v_ndfsZ[asymIndexZ])[iZ*binningZ[iPa1].size()+iZ2];
				    cout <<"w ndf: " << asymmetries::ndfFit<<endl;
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
				    //there is a special case for pnnp. There the particle type does not have to match....
#else
#ifdef WITH_PURITIES
				    //				    if(iCh1==PNNP || iCh2 ==PNNP)
				    //				      cout <<"purity: " << purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2] <<endl;

				    if(purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2]>0)
				      {
					mY[iZ2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					mYErr[iZ2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					mYKinCorr[iZ2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					mYErrKinCorr[iZ2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					outputFile << "myerr_pr: " << mYErr_Pr[iZ2];
					mYErr_Pr[iZ2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iZ2];
					outputFile <<" and now: " << mYErr_Pr[iZ2] <<endl;
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
#endif//mc
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
				pGP->SetMarkerStyle(23); 
				pGP->SetMarkerSize(1.2);
				pGP->SetMarkerColor(kRed);
				/*			pGP->GetYaxis()->SetLabelSize(0.07);
				  pGP->GetYaxis()->SetTitleSize(0.07);
				  pGP->GetXaxis()->SetTitleSize(0.07);
				  pGP->GetXaxis()->SetLabelSize(0.07);*/
				if(iCh1==PN && iPa1==PiPi&&iCh2==PN && iPa2==PiPi)
				  pGP->GetYaxis()->SetRangeUser(0.8,1.0);
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
				outputFile <<" pushing back ptg" <<endl;
				vecGraphs.push_back(pTG);
				vecGraphsKinCorr.push_back(pTGKinCorr);
				//				pTG_Proj->GetYaxis()->SetLabelSize(0.07);
				//				pTG_Proj->GetYaxis()->SetTitleSize(0.07);
				//				pTG_Proj->GetXaxis()->SetTitleSize(0.07);
				//				pTG_Proj->GetXaxis()->SetLabelSize(0.07);
				pTG_Proj->GetYaxis()->SetNdivisions(6,true);
				pTG_Proj->GetXaxis()->SetNdivisions(6,true);
				//				pTG_Proj->SetMarkerStyle(23); 
				//				pTG_Proj->SetMarkerSize(1.2);
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
				    cout <<"ndf: " << asymmetries::ndfFit<<endl;
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
					outputFile << "myerr_pr: " << mYErr_Pr[iM2];
					mYErr_Pr[iM2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
					outputFile <<" and now: " << mYErr_Pr[iM2] <<endl;
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
				pGP->SetMarkerStyle(23); 
				pGP->SetMarkerSize(1.2);
				pGP->SetMarkerColor(kRed);
				pGP->GetYaxis()->SetTitle("Purity");
				if(iCh1==PN && iPa1==PiPi&&iCh2==PN && iPa2==PiPi)
				  pGP->GetYaxis()->SetRangeUser(0.8,1.0);
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
				//				pTG_Proj->SetMarkerStyle(23); 
				//				pTG_Proj->SetMarkerSize(1.2);
				pTG_Proj->SetMarkerColor(kRed);
				pTG_Proj->GetYaxis()->SetTitle("A_{(#Phi_{1}+#Phi_{2})}");
				pTG_Proj->GetXaxis()->SetTitle("M_{Inv} [GeV]");
				vecGraphsProjections.push_back(pTG_Proj);
				//			pTG_Proj->SaveAs((strCharge+" " + strParticle + " " +str.str())+".C").c_str());
			      }
			    asymIndexM++;
			    break;
			  }
			case multBinning:
			  {
			    break;
			  }
			case mixedBinning:
			  {	
			    cout <<"-1-" <<endl;
			    vecBinning.push_back(binningZ[iPa1].size());
			    for(int iZ=0;iZ<binningZ[iPa1].size();iZ++)
			      {
				cout <<"iZ: " << iZ <<endl;
				stringstream str;
				str <<"z_from_";
				if(iZ==0)
				  str << 0;
				else
				  str << binningZ[iPa1][iZ-1];
				str << "_to_" << binningZ[iPa1][iZ];
				for(int iM2=0;iM2<binningM[iPa1].size();iM2++)
				  {
				    cout <<"im2: " << iM2 <<endl;
				    mX[iM2]=xVals[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]/xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2];
				    mXErr[iM2]=0;
				    mY[iM2]=(*mAsyms[asymIndexMixed])[iZ*binningZ[iPa1].size()+iM2].first;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];

				    mYKinCorr[iM2]=(*mAsymsMod2[asymIndexMixed])[iZ*binningZ[iPa1].size()+iM2].first;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    mY_Pr[iM2]=0;
				    mYErr[iM2]=(*mAsyms[asymIndexMixed])[iZ*binningZ[iPa1].size()+iM2].second;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    cout <<"2" <<endl;
				    cout <<"asbin for mixed: " << iZ <<" " << iM2 << endl;
				    asymmetries::binning=iBin;
				    asymmetries::kinBin1=iZ;
				    asymmetries::kinBin2=iM2;
				    asymmetries::chargeType1=iCh1;
				    asymmetries::chargeType2=iCh2;
				    asymmetries::particleType1=iPa1;
				    asymmetries::particleType2=iPa2;
				    asymmetries::asymmetry=mY[iM2];
				    asymmetries::asError=mYErr[iM2];
				    cout <<"3" <<endl;
				    asymmetries::meanKinVal1=yVals[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]/yValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2];
				    cout <<"4" <<endl;
				    asymmetries::meanKinVal2=mX[iM2];
				    cout <<"5" <<endl;
				    asymmetries::chi2Fit=(*v_chi2Mixed[asymIndexMixed])[iZ*binningZ[iPa1].size()+iM2];
				    cout <<"6" <<endl;
				    asymmetries::ndfFit=(*v_ndfsMixed[asymIndexMixed])[iZ*binningZ[iPa1].size()+iM2];
				    cout <<"ndf: " << asymmetries::ndfFit<<endl;
				    asymmetries::expNr=expNumber;
				    asymmetries::isOnRes=onResonance;
				    asymmetries::meanTheta= thetaMeanRatioMix[iZ][iM2]/thetaCountMix[iZ][iM2];
				    asymmetries::meanThetaSin= thetaMeanSinMix[iZ][iM2]/thetaCountMix[iZ][iM2];
				    asymmetries::meanThetaCos=thetaMeanCosMix[iZ][iM2]/thetaCountMix[iZ][iM2];
				    as_tree.Fill();
				    mYErrKinCorr[iM2]=(*mAsymsMod2[asymIndexMixed])[iZ*binningZ[iPa1].size()+iM2].second;//kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];
				    mYErr_Pr[iM2]=(*mAsyms[asymIndexMixed])[iZ*binningZ[iPa1].size()+iM2].second/errFact;///(kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2]*errFact)*(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iM][iM2];;
				    outputFile<<"projected error (m): " <<mYErr_Pr[iM2] <<endl;
				    cout <<"myErr: " << mYErr[iM2] <<" myErr+Pr: " << mYErr_Pr[iM2] << " errFact: " << errFact <<endl;
#ifdef MC
				    cout <<"for mixed, numFalse: " << numFalseIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2] <<" numcorr: " <<numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2] <<" for iZ: " << iZ << " iM: " << iM2<<endl;
				    if(numFalseIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]+numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]!=0)
				      {
					purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]=numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]/(float)(numFalseIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]+numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]);
					cout <<"for mixed purities: bin: " << iZ <<", iM2 " <<iM2<<" purity: " <<purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2] <<endl;
					puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]=getError(numCorrIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2],numFalseIdent[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]);
				      }

				    cout <<"mixed purity for " << strCharge+"_" + strParticle + "_" +str.str() << " m bin  : " << iM2 << ": " <<purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]<<endl;
#else
#ifdef WITH_PURITIES
				    if(purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]>0)
				      {
					mY[iM2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2];
					mYErr[iM2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2];
					mYKinCorr[iM2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2];
					mYErrKinCorr[iM2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2];
					outputFile << "myerr_pr: " << mYErr_Pr[iM2];
					cout << "myerr_pr: " << mYErr_Pr[iM2];
					mYErr_Pr[iM2]/=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2];
					outputFile <<" and now: " << mYErr_Pr[iM2] <<endl;
					cout <<" and now: " << mYErr_Pr[iM2] <<endl;
					outputFile <<"using purity of " << purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2] << " for "  << strCharge+"_" + strParticle + "_" +str.str() << " m bin  : " << iM2;
					cout <<"using purity of " << purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2] << " for "  << strCharge+"_" + strParticle + "_" +str.str() << " m bin  : " << iM2;
					outputFile <<" purity Error: " << puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2];
					cout <<" purity Error: " << puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2];
					outputFile <<" and correction fact of "  << kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]/(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]<<endl;
					cout <<" and correction fact of "  << kinCorrFact[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]/(float)xValsNum[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]<<endl;
					outputFile <<"resulting Asym: " <<  mY[iM2];
					outputFile << " error: " << mYErr[iM2] <<endl;
					cout <<"resulting Asym: " <<  mY[iM2];
					cout << " error: " << mYErr[iM2] <<endl;
				      }
				    else
				      {
					cout << " mixed (fail) purity for " << strCharge+"_" + strParticle + "_" +str.str() << " m bin  : " << iM2 << " is " <<purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]<<endl;
					outputFile << " purity for " << strCharge+"_" + strParticle + "_" +str.str() << " m bin  : " << iM2 << " is " <<purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]<<endl;
				      }
#endif
#endif//mc
				  }
				graphTitles.push_back(strCharge+"_" + strParticle + "_" +str.str());
				graphTitlesKinCorr.push_back(strCharge+"_" + strParticle + "_" +str.str()+"KinCorr");
				graphTitlesProj.push_back(strCharge+"_" + strParticle + "_" +str.str()+ "Proj");
				graphTitlesPurities.push_back(strCharge+"_"+strParticle+"_"+str.str()+"MixedPurity");
				//mX goes over all inv mass
				TGraphErrors* pTG=new TGraphErrors(binningM[iPa1].size(),mX,mY,mXErr,mYErr);
				//			TGraphErrors* pTG=new TGraphErrors(binningM[iPa1].size(),mX,mYKinCorr,mXErr,mYErr);
				TGraphErrors* pTGKinCorr=new TGraphErrors(binningM[iPa1].size(),mX,mYKinCorr,mXErr,mYErrKinCorr);
				TGraphErrors* pGP=new TGraphErrors(binningM[iPa1].size(),mX,purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ],mXErr,puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ]);


				pGP->GetYaxis()->SetRangeUser(0,1);
				pGP->GetYaxis()->SetLabelSize(0.07);
				pGP->GetYaxis()->SetTitleSize(0.07);
				pGP->GetXaxis()->SetTitleSize(0.07);
				pGP->GetXaxis()->SetLabelSize(0.07);
				pGP->GetYaxis()->SetNdivisions(6,true);
				pGP->GetXaxis()->SetNdivisions(6,true);
				pGP->SetMarkerStyle(23); 
				pGP->SetMarkerSize(1.2);
				pGP->SetMarkerColor(kRed);
				pGP->GetYaxis()->SetTitle("Purity");
				if(iCh1==PN && iPa1==PiPi&&iCh2==PN && iPa2==PiPi)
				  pGP->GetYaxis()->SetRangeUser(0.8,1.0);
				pGP->GetXaxis()->SetTitle("Mixed M_{Inv} [GeV]");
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
				//			n	pTG_Proj->SetMarkerStyle(23); 
				//				pTG_Proj->SetMarkerSize(1.2);
				pTG_Proj->SetMarkerColor(kRed);
				pTG_Proj->GetYaxis()->SetTitle("A_{(#Phi_{1}+#Phi_{2})}");
				pTG_Proj->GetXaxis()->SetTitle("M_{Inv} [GeV]");
				vecGraphsProjections.push_back(pTG_Proj);
				//			pTG_Proj->SaveAs((strCharge+" " + strParticle + " " +str.str())+".C").c_str());
			      }
			    asymIndexMixed++;
			    break;
			  }
			default:
			  cout << "binning not recognized3" << endl;
			  exit(0);
			}
		    }
		}
	    }
	}
    }
  asFile.Write();
  asFile.Close();

#endif //purities_only
  for(int i=0;i<binningZ[PiPi].size();i++)
    {
      for(int j=0;j<binningZ[PiPi].size();j++)
	{
	  xcheckFile << countsZBin[i][j]<<endl;
	}
    }
  xcheckFile <<endl<<endl<<"inv mass counts " << endl <<endl;
  for(int i=0;i<binningM[PiPi].size();i++)
    {
      for(int j=0;j<binningM[PiPi].size();j++)
	{
	  xcheckFile << countsMBin[i][j]<<endl;
	}
    }


  TCanvas c1("c1","c1",10,20,500,200);
  c1.Divide(1,2);
  c1.cd(1);
  hZ1.Draw();
  c1.cd(2);
  hZ2.Draw();
  c1.SaveAs("ZCanvas.ps");
  c1.SaveAs("ZCanvas.png");
  TCanvas c2("c2","c2",10,20,500,200);

  c2.Divide(1,2);
  c2.cd(1);
  hR.Draw();
  c2.cd(2);
  hAsDiff.Draw();
  c2.SaveAs("AsDiffs.ps");
  c2.SaveAs("AsDiffs.pdf");

  TCanvas cD0("cD0","cD0",10,20,500,200);
  cD0.Divide(1,2);
  cD0.cd(1);
  hNumD0FirstHemi.Draw();
  cD0.cd(2);
  hNumD0SecondHemi.Draw();
  cD0.SaveAs("cD0.png");
  cD0.SaveAs("cD0.pdf");

  TH2D twoDZ("twoDZ","twoDZ",10,0,10,5,0,5);

  /*  for(int i=0;i<5;i++) //only works for 5x5 bins
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
  */ 
  TCanvas cTDZ("ctwoDZ","ctwoDZ",10,20,500,200);
  twoDZ.Draw("col");
  cTDZ.SaveAs("zSpectra/twoDZ.pdf");


  TCanvas cPhiR;
  hPhiR.Draw();
  cPhiR.SaveAs("PhiR.C");

  TCanvas cPhiZ;
  hPhiZ.Draw();
  cPhiZ.SaveAs("PhiZ.C");
#ifndef PURITIES_ONLY
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
  saveAs("Asymmetries.png",vecGraphs,graphTitles,outputFile,vecBinning);
  //    outputFile << "saving " << vecGraphsProjections.size() <<" projections having " << graphTitlesProj.size() <<" titles"<<endl;
  //  saveAs("Pr.C",vecGraphsProjections,graphTitlesProj,outputFile,vecBinning);
  //  saveAs("Pr.ps",vecGraphsProjections,graphTitlesProj,outputFile,vecBinning);
  //  saveAs("Pr.pdf",vecGraphsProjections,graphTitlesProj,outputFile,vecBinning);

  outputFile << "done "<<endl;
#ifdef WITH_PURITIES
  cout <<"my puri" <<endl;
  cout <<"saving " << vecPurities.size() <<"purity vectors to save " <<endl;
  for(int i=0;i<graphTitlesPurities.size();i++)
    {
      cout <<"saving " << graphTitlesPurities[i]<<endl;
    }
  saveAs(".pdf",vecPurities,graphTitlesPurities, outputFile, vecBinning);
  saveAs(".C",vecPurities,graphTitlesPurities, outputFile,vecBinning);
  saveAs(".eps",vecPurities,graphTitlesPurities, outputFile, vecBinning);
#endif
#endif
#ifdef MC
  for(int iBin=zBinning;iBin<=mixedBinning;iBin++)
    {
      purities::binning=iBin;
      for(int iCh1=PN;iCh1<=LAST_CHARGE_COMB;iCh1++)
	//      for(int iCh=PN;iCh<=ZZ;iCh++)
	{
	  for(int iCh2=PN;iCh2<=iCh1;iCh2++)
	    {
	      if(!(iCh1==PN || iCh1==PZ || iCh1==ZN || iCh1==PNNP || iCh1==ZNNZ))
		continue;
		if(!(iCh2==PN || iCh2==PZ || iCh2==ZN || iCh2==PNNP || iCh2==ZNNZ))
		continue;

		/*	      if(!(iCh1==PN || iCh1==PZ || iCh1==ZN))
		continue;
	      if(!(iCh2==PN || iCh2==PZ || iCh2==ZN))
	      continue;*/
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
			case multBinning:
			  {
			    break;
			  }
			case mixedBinning:
			  {
			    for(int iZ=0;iZ<binningZ[iPa1].size();iZ++)
			      {
				for(int iM2=0;iM2<binningM[iPa1].size();iM2++)
				  {
				    purities::kinBin1=iZ;
				    purities::kinBin2=iM2;
				    cout <<"writing purity for mixed: " << purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2]<<endl;
				    purities::purity=purities[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2];
				    purities::purityError=puritiesErrors[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)][iZ][iM2];
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
  cout <<"writing file " <<endl;
  m_file.Write();
  m_file.Close();
#endif



  TCanvas cThrustProj("thrustProj","thrustProj",10,20,200,500);
  thrustProj_PN->Draw();
  cThrustProj.SaveAs("cThurstProj.png");
  cThrustProj.SaveAs("cThurstProj.pdf");

  TCanvas cMPiPi_PN("cMPiPi_PN","M_PiPi_PN",10,20,200,500);
  massPiPi_PN->Draw();
  cMPiPi_PN.SaveAs("cMPiPi_PN.C");
  cMPiPi_PN.SaveAs("cMPiPi_PN.pdf");
  cMPiPi_PN.SaveAs("cMPiPi_PN.ps");
  cMPiPi_PN.SaveAs("cMPiPi_PN.png");
  TCanvas cMPiPi_PZ("cMPiPi_PZ","M_PiPi_PZ",10,20,200,500);
  massPiPi_PZ->Draw();
  cMPiPi_PZ.SaveAs("cMPiPi_PZ.C");
  cMPiPi_PZ.SaveAs("cMPiPi_PZ.pdf");
  cMPiPi_PZ.SaveAs("cMPiPi_PZ.ps");
  cMPiPi_PZ.SaveAs("cMPiPi_PZ.png");
  TCanvas cMPiPi_ZN("cMPiPi_ZN","M_PiPi_ZN",10,20,200,500);
  massPiPi_ZN->Draw();
  cMPiPi_ZN.SaveAs("cMPiPi_ZN.C");
  TCanvas cMPiK_PN("cMPiK_PN","M_PiK_PN",10,20,200,500);
  massPiK_PN->Draw();
  cMPiK_PN.SaveAs("cMPiK_PN.C");
  cMPiK_PN.SaveAs("cMPiK_PN.pdf");
  cMPiK_PN.SaveAs("cMPiK_PN.png");
  cMPiK_PN.SaveAs("cMPiK_PN.ps");

#ifdef MC
  for(int iBin=zBinning;iBin<=mixedBinning;iBin++)
    {
      for(int iCh=PN;iCh<=LAST_CHARGE_COMB;iCh++)
	{
	  if(!(iCh==PN || iCh==PZ || iCh==ZN || iCh==PNNP || iCh==ZNNZ))
	    continue;
	  strCharge=iToStrCharge(iCh);
	  for(int iPa=PiPi;iPa<=LAST_PARTICLE_COMB;iPa++)
	    {
	      if(iPa!=PiPi)
		continue;
	      strParticle=iToStrPart(iPa);
	      switch(iBin)
		{
		case zBinning:
		  {
		    //		    cout <<" numBins: " << binningZ[iPa].size()<<endl;
		    for(int iZ=0;iZ<binningZ[iPa].size();iZ++)
		      {
			//			cout <<"z: " << iZ <<" numFalse: " <<numFalseIdentInt[iBin][iCh][iPa][iZ] <<" " << numCorrIdentInt[iBin][iCh][iPa][iZ] <<endl;
			if(numFalseIdentInt[iBin][iCh][iPa][iZ]+numCorrIdentInt[iBin][iCh][iPa][iZ]!=0)
			  {
			    puritiesInt[iBin][iCh][iPa][iZ]=numCorrIdentInt[iBin][iCh][iPa][iZ]/(float)(numFalseIdentInt[iBin][iCh][iPa][iZ]+numCorrIdentInt[iBin][iCh][iPa][iZ]);
			    puritiesErrorsInt[iBin][iCh][iPa][iZ]=getError(numCorrIdentInt[iBin][iCh][iPa][iZ],numFalseIdentInt[iBin][iCh][iPa][iZ]);
			    cout <<"purity for " << strCharge+"_" + strParticle + "_"  << " z bin  : " << iZ << ": " <<puritiesInt[iBin][iCh][iPa][iZ]<<" " << puritiesErrorsInt[iBin][iCh][iPa][iZ]<<endl;
			  }
		      }
		    break;
		  }
		case mBinning:
		  {
		    for(int iM=0;iM<binningM[iPa].size();iM++)
		      {
			if(numFalseIdentInt[iBin][iCh][iPa][iM]+numCorrIdentInt[iBin][iCh][iPa][iM]!=0)
			  {
			    puritiesInt[iBin][iCh][iPa][iM]=numCorrIdentInt[iBin][iCh][iPa][iM]/(float)(numFalseIdentInt[iBin][iCh][iPa][iM]+numCorrIdentInt[iBin][iCh][iPa][iM]);
			      puritiesErrorsInt[iBin][iCh][iPa][iM]=getError(numCorrIdentInt[iBin][iCh][iPa][iM],numFalseIdentInt[iBin][iCh][iPa][iM]);
			      cout <<"purity for " << strCharge+"_" + strParticle + "_" << " m bin  : " << iM << ": " <<puritiesInt[iBin][iCh][iPa][iM]<<" " << puritiesErrorsInt[iBin][iCh][iPa][iM]<<endl;
			  }
		      }
		    break;
		  }
		}
	    }
	}
      
    }		

#endif    
  cout <<"num Charm: " << numCharmEv << " all: " << allEvts <<"  ratio: " << numCharmEv/(float)(allEvts) <<endl;
  cout <<"cuttedEv: " << cuttedEv <<" nanEv: " << nanEv << " accEv: " << accEv <<endl;
}
