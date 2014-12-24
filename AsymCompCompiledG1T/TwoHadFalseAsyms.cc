//#define MC
//#define PAW_STYLE
//#define WoA   //without acceptance
//#define WITH_PURITIES
//#define WITH_KIN_CORR


//#define PHIR1_MINUS_PHIR2
//#define TWO_TIMES_PHIR1_MINUS_PHIR2


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
#include <math.h>
#include <time.h>
using namespace std;
ofstream* errorVglFile;
vector<float>* binningM;
vector<float>* binningZ;
vector<float>* binningMult;
vector<float>* binningZSpect;
//to compare errors from the fits with the naive statistical errors sqrt(N)/N
#include "TwoHadAsymsCommons.h"
#include "modview.h"
//#define MAX_EVENTS 100

#define AMP_PHIR 0.1
#define OPENING_CUT 0.8

//modulation of the counts with the mc angle
float getInc(float phiR, double mass1, double mass2,  double z1, double z2,int pt1, int pt2, int charge1,int charge2);

double******* kinSpectraPars;
double****** kinSpectraReducedPars;
double****** kinSpectra;
double***** kinSpectraReduced;
TH2D***** kinSpectraH;
TH1D***** kinSpectra1H;

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
  float chi2Fit;
  float ndfFit;
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
  ofstream outFile("Asymmetries");
#ifdef PAW_STYLE
  //for paw style....
  TChain* chI=new TChain("h7");//stupid paw...
  TChain* chF=new TChain("h8");
#else
  cout <<"no paw style" <<endl;
#ifdef WoA
  TChain* chAll=new TChain("GenTree");
  cout <<"read GenTree" <<endl;
#else
  TChain* chAll=new TChain("DataTree");
#endif
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
  binningZSpect=new vector<float>[NumParticle];
  binningMult=new vector<float>[NumParticle];

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

  /*  binningM[PiPi].push_back(0.4);
  binningM[PiPi].push_back(0.55);
  binningM[PiPi].push_back(0.77);
  binningM[PiPi].push_back(1.2);
  binningM[PiPi].push_back(2.0);

  binningZ[PiPi].push_back(0.3);
  binningZ[PiPi].push_back(0.4);
  binningZ[PiPi].push_back(0.55);
  binningZ[PiPi].push_back(0.75);
  binningZ[PiPi].push_back(1.0);
  */

  binningMult[PiPi].push_back(3);
  binningMult[PiPi].push_back(6);
  binningMult[PiPi].push_back(9);
  binningMult[PiPi].push_back(14);
  binningMult[PiPi].push_back(1000);

  binningMult[PiK].push_back(3);
  binningMult[PiK].push_back(6);
  binningMult[PiK].push_back(9);
  binningMult[PiK].push_back(14);
  binningMult[PiK].push_back(1000);

  binningMult[KPi].push_back(3);
  binningMult[KPi].push_back(6);
  binningMult[KPi].push_back(9);
  binningMult[KPi].push_back(14);
  binningMult[KPi].push_back(1000);


  binningMult[KK].push_back(3);
  binningMult[KK].push_back(6);
  binningMult[KK].push_back(9);
  binningMult[KK].push_back(14);
  binningMult[KK].push_back(1000);

binningMult[UNKNOWN].push_back(3);
  binningMult[UNKNOWN].push_back(6);
  binningMult[UNKNOWN].push_back(9);
  binningMult[UNKNOWN].push_back(14);
  binningMult[UNKNOWN].push_back(1000);


  binningZSpect[PiPi].push_back(0.25);
  binningZSpect[PiPi].push_back(0.3);
  binningZSpect[PiPi].push_back(0.35);
  binningZSpect[PiPi].push_back(0.4);
  binningZSpect[PiPi].push_back(0.45);
  binningZSpect[PiPi].push_back(0.55);
  binningZSpect[PiPi].push_back(0.65);
  binningZSpect[PiPi].push_back(0.75);
  binningZSpect[PiPi].push_back(1.0);


  /*  binningM[PiK].push_back(0.7);
  binningM[PiK].push_back(0.9);
  binningM[PiK].push_back(1.2);
  binningM[PiK].push_back(1.6);
  binningM[PiK].push_back(2.0);

  binningZ[PiK].push_back(0.3);
  binningZ[PiK].push_back(0.5);
  binningZ[PiK].push_back(0.7);
  binningZ[PiK].push_back(1.2);



  binningM[KPi].push_back(0.7);
  binningM[KPi].push_back(0.85);
  binningM[KPi].push_back(1.0);
  binningM[KPi].push_back(1000.0);

  binningZ[KPi].push_back(0.3);
  binningZ[KPi].push_back(0.5);
  binningZ[KPi].push_back(0.7);
  binningZ[KPi].push_back(1.2);

  binningM[KK].push_back(1.2);
  binningM[KK].push_back(1.5);
  binningM[KK].push_back(1.7);
  binningM[KK].push_back(1000.0);

  binningZ[KK].push_back(0.3);
  binningZ[KK].push_back(0.5);
  binningZ[KK].push_back(0.7);
  binningZ[KK].push_back(1.2);

  binningM[UNKNOWN].push_back(0.4);
  binningM[UNKNOWN].push_back(0.55);
  binningM[UNKNOWN].push_back(0.77);
  binningM[UNKNOWN].push_back(1000.0);

  binningZ[UNKNOWN].push_back(0.3);
  binningZ[UNKNOWN].push_back(0.5);
  binningZ[UNKNOWN].push_back(0.7);
  binningZ[UNKNOWN].push_back(1.2);*/

  loadBinning(binningM, binningZ); // if the binning is done any other way, it is incompatible with modview

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
  float phiRDiff[2000];

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

  float thrustProj11[2000];
  int thrustProj11Counter;
  float thrustProj12[2000];
  int thrustProj12Counter;
  float thrustProj21[2000];
  int thrustProj21Counter;
  float thrustProj22[2000];
  int thrustProj22Counter;

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
  float phiRDiff_mc[2000];
  int phiRSum_mcCounter;

  int charge1Counter;
  int chargeType1[2000];
  int particle1Counter;
  int particleType1[2000];
  int multiplicity[2000];

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

#ifdef WoA
  for(int i=0;i<fieldNamesI.size();i++)
    {
      if(fieldNamesI[i].find("Counter")!=string::npos)
	{
	  int mpos=fieldNamesI[i].find("Counter");
	  string mst=fieldNamesI[i];
	  cout <<"first mst1: " << mst <<endl;
	  mst.replace(mpos,mpos+7,"_mcWoACounter");
	  cout <<"new mst: " << mst <<endl;
	  fieldNamesI[i]=mst;
	}
      else
	fieldNamesI[i]=((fieldNamesI[i]+string("_mcWoA")));
    }
  for(int i=0;i<fieldNamesF.size();i++)
    {
      fieldNamesF[i]=((fieldNamesF[i]+string("_mcWoA")));
    }
#endif

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
      cout <<"trying to branch" <<fieldNamesF[i]<<endl;
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
#ifdef WoA
  for(int i=0;i<fieldNamesI.size();i++)
    {
      if(fieldNamesI[i].find("Counter")!=string::npos)
	{
	  int mpos=fieldNamesI[i].find("Counter");
	  string mst=fieldNamesI[i];
	  cout <<"first mst: " << mst <<endl;

	  mst.replace(mpos,mpos+7,"_mcWoACounter");
	  cout <<"new mst: " << mst <<endl;
	  fieldNamesI[i]=mst;
	}
      else
	fieldNamesI[i]=(fieldNamesI[i]+string("_mcWoA"));
    }
  for(int i=0;i<fieldNamesF.size();i++)
    {
      fieldNamesF[i]=(fieldNamesF[i]+string("_mcWoA"));
    }

#endif

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
  cout <<nevents << " events " <<endl;
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
    int locNumBins=NumBin+1;//due to new mult binning
  cout<<" numbin " << locNumBins <<" numcharge: " << NumCharge <<", NumParticle: " << NumParticle <<", maxKin: " << maxKin <<", numAngBins: " << numAngBins<<endl;
  int****** counts=allocateArray<int>(locNumBins,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  //combine hads in the same hemi --- look at dihadana to see how they are combined... might be only pn or all..
  int****** countsSameHemi=allocateArray<int>(locNumBins,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  int****** countsDiffEvt=allocateArray<int>(locNumBins,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin,numAngBins);
  //relates z and m spectra
  gHCounter=0;

  //avg over value2, conforming to notation in paper
  float***** xVals=allocateArray<float>(locNumBins,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** meanAngDiff=allocateArray<float>(locNumBins,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** meanAngSum=allocateArray<float>(locNumBins,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** xValsNum=allocateArray<int>(locNumBins,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  float***** yVals=allocateArray<float>(locNumBins,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  int***** yValsNum=allocateArray<int>(locNumBins,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);

  cout << "yvals am anfang " <<yVals[0][0][0][0][0] <<" nm: " << yValsNum[0][0][0][0][0]<<endl;


#ifdef MC
  unsigned long***** numCorrIdent=allocateArray<unsigned long>(locNumBins,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
  unsigned long***** numFalseIdent=allocateArray<unsigned long>(locNumBins,(NumCharge*NumCharge+NumCharge)/2,(NumParticle*NumParticle+NumParticle)/2,maxKin,maxKin);
#else

#ifdef WITH_PURITIES
#endif
#ifdef WITH_KIN_CORR

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

      int numFalseCombs1=0;
      int numFalseCombs2=0;
      int numFalseCombs=0;
      float falseAngles1st[2000];
      float falseAngles2nd[2000];
      float falseAngles[2000];
      float falseMass1[2000];
      float falseMass2[2000];
      float falseMass[2000];
      float falseZ1[2000];
      float falseZ2[2000];
      float falseZRatio1[2000];
      float falseZRatio2[2000];
      float falseMult[2000];
      int falsePT1[2000];
      int falsePT2[2000];
      int falseCT1[2000];
      int falseCT2[2000];

      bool doneSecHalf=false;

      for(int iHadsInEv=1;iHadsInEv<z1Counter;iHadsInEv++)
	{
	  if(fabs(thrustProj11[iHadsInEv])<OPENING_CUT || fabs(thrustProj12[iHadsInEv])<OPENING_CUT ||fabs(thrustProj21[iHadsInEv])<OPENING_CUT ||fabs(thrustProj22[iHadsInEv])<OPENING_CUT)
	    {
	      //	      continue;  /this cut leads to weird efects, because then some combinations are not filled
	    }
	  //since the loop has the first pair as the outer loop
	  if(!(z1[iHadsInEv]==z1[iHadsInEv-1]))
	    {
	      //	      if(chargeType1[iHadsInEv-1]==PN && particleType1[iHadsInEv-1]==PiPi)
		{
		  falseAngles1st[numFalseCombs1]=phiR1[iHadsInEv-1];
		  falseZ1[numFalseCombs1]=z1[iHadsInEv-1];
		  falseZRatio1[numFalseCombs1]=z1Ratio[iHadsInEv-1];
		  //ratio: z of pair / z of one hadron.... wrong order here?
		  falseMass1[numFalseCombs1]=mass1[iHadsInEv-1];
		  falsePT1[numFalseCombs1]=particleType1[iHadsInEv-1];
		  falseCT1[numFalseCombs1]=chargeType1[iHadsInEv-1];
		  //  cout <<"falsecomb: " << numFalseCombs1 <<endl;
		  numFalseCombs1++;
		}
	      doneSecHalf=true;
	    }
	  else //same first hadron pair
	    {
	      if(!doneSecHalf)
		{
		  if(iHadsInEv==1)
		    {
			{
			  falseAngles2nd[numFalseCombs2]=phiRSum[iHadsInEv-1]-phiR1[iHadsInEv-1];
			  if(falseAngles2nd[numFalseCombs2]<-pi)
			    falseAngles2nd[numFalseCombs2]+=2*pi;
			  if(falseAngles2nd[numFalseCombs2]>pi)
			    falseAngles2nd[numFalseCombs2]-=2*pi;
			  
			  falseMass2[numFalseCombs2]=mass2[iHadsInEv-1];
			  falsePT2[numFalseCombs2]=particleType2[iHadsInEv-1];
			  falseCT2[numFalseCombs2]=chargeType2[iHadsInEv-1];
			  falseZ2[numFalseCombs2]=z2[iHadsInEv-1];
			  falseZRatio2[numFalseCombs2]=z2Ratio[iHadsInEv-1];
			  numFalseCombs2++;

			}
		    }
		    {
		      falseAngles2nd[numFalseCombs2]=phiRSum[iHadsInEv]-phiR1[iHadsInEv];
		      if(falseAngles2nd[numFalseCombs2]<-pi)
			falseAngles2nd[numFalseCombs2]+=2*pi;
		      if(falseAngles2nd[numFalseCombs2]>pi)
			falseAngles2nd[numFalseCombs2]-=2*pi;

		      falseZ2[numFalseCombs2]=z2[iHadsInEv];
		      falseZRatio2[numFalseCombs2]=z2Ratio[iHadsInEv];
		      falseMass2[numFalseCombs2]=mass2[iHadsInEv];
		      falsePT2[numFalseCombs2]=particleType2[iHadsInEv];
		      falseCT2[numFalseCombs2]=chargeType2[iHadsInEv];
		      numFalseCombs2++;
		    }
		}
	    }
	  if(iHadsInEv==(z1Counter-1))
	    {
		{
		    {
		      falseAngles1st[numFalseCombs1]=phiR1[iHadsInEv];
		      falseZ1[numFalseCombs1]=z1[iHadsInEv];
		      falseZRatio1[numFalseCombs1]=z1Ratio[iHadsInEv];
		      falseMass1[numFalseCombs1]=mass1[iHadsInEv-1];
		      falsePT1[numFalseCombs1]=particleType1[iHadsInEv-1];
		      falseCT1[numFalseCombs1]=chargeType1[iHadsInEv-1];
		      numFalseCombs1++;
		    }
		}
	    }
	}
      
      for(int i=0;i<numFalseCombs1-1;i++)
	{
	  for(int j=i+1;j<numFalseCombs1;j++)
	    {
	      //ratio is sum/z1 in one pair-> z1=sum/ratio z2 = sum -z1
	      float locZ1_1=falseZ1[i]/falseZRatio1[i];
	      float locZ1_2=falseZ1[i]-locZ1_1;
	      float locZ2_1=falseZ1[j]/falseZRatio1[j];
	      float locZ2_2=falseZ1[j]-locZ2_1;
	      //make sure that we don't have the same hadrons in each combination
	      if(locZ1_1!=locZ2_1 && locZ1_2!=locZ2_2)
		{
		  //	      cout <<"nachher: z1_1: " << locZ1_1 << " z1_2: " << locZ1_2 << " z2_1: " << locZ2_1 <<" z2_2: " << locZ2_2 <<endl;

		  //try only one angle...

		  //this is the angle that we look at, so if we take the difference here for Boer's stuff, it will be a misnomer... but so be it
 phiRSum[numFalseCombs]=falseAngles1st[i]+falseAngles1st[j];
#ifdef PHIR1_MINUS_PHIR2
		  phiRSum[numFalseCombs]=falseAngles1st[i]-falseAngles1st[j];
	   if(phiRSum[numFalseCombs] > pi)
	     phiRSum[numFalseCombs]-=2*pi;
	   if(phiRSum[numFalseCombs] < -pi)
	     phiRSum[numFalseCombs]+=2*pi;
#endif
#ifdef TWO_TIMES_PHIR1_MINUS_PHIR2
		  phiRSum[numFalseCombs]=2*(falseAngles1st[i]-falseAngles1st[j]);
	   if(phiRSum[numFalseCombs] > pi)
	     phiRSum[numFalseCombs]-=2*pi;
	   if(phiRSum[numFalseCombs] < -pi)
	     phiRSum[numFalseCombs]+=2*pi;
	   if(phiRSum[numFalseCombs] > pi)
	     phiRSum[numFalseCombs]-=2*pi;
	   if(phiRSum[numFalseCombs] < -pi)
	     phiRSum[numFalseCombs]+=2*pi;
#endif
		 
		  if(falseAngles1st[i]>falseAngles1st[j])
		    phiRDiff[numFalseCombs]=falseAngles1st[i]-falseAngles1st[j];
		  else
		    phiRDiff[numFalseCombs]=falseAngles1st[j]-falseAngles1st[i];

		  if(phiRSum[numFalseCombs]<-pi)
		    phiRSum[numFalseCombs]+=2*pi;
		  if(phiRSum[numFalseCombs]>pi)
		    phiRSum[numFalseCombs]-=2*pi;
		  //	      cout <<"phirsum: " << phiRSum[numFalseCombs]<<" for " << numFalseCombs<<" before " << falseAngles1st[i] <<" + " << falseAngles1st[j] << " = " <<falseAngles1st[i]+falseAngles1st[j]<<endl;
		  z2[numFalseCombs]=falseZ1[i];
		  		  z2[numFalseCombs]=falseZ1[j];
				  //z1[numFalseCombs]=0.25;
		  //		  printf("combining in first pair z1: %f z2: %f ang1: %f ang2: %f\n ",falseZ1[i],falseZ1[j],falseAngles1st[i],falseAngles1st[j]);
		  chargeType1[numFalseCombs]=falseCT1[i];
		  chargeType2[numFalseCombs]=falseCT1[j];
		  particleType1[numFalseCombs]=falsePT1[i];
		  if(falsePT2[i] > 20)
		    cout <<"falsept1: " << falsePT1[i] << " i: "  << i <<endl;
		  if(particleType1[numFalseCombs]>NumParticle-1)
		    cout <<"false pt: " << particleType1[numFalseCombs]<<endl;
		  particleType2[numFalseCombs]=falsePT1[j];
		  if(particleType2[numFalseCombs]>NumParticle-1)
		    cout <<"false pt: 2" << particleType2[numFalseCombs]<<endl;
		  multiplicity[numFalseCombs]=numFalseCombs1;
		  //		  cout << "false combs1: " << numFalseCombs1<<endl;
		  numFalseCombs++;
		}
	    }
	}
      for(int i=0;i<numFalseCombs2-1;i++)
	{
	  for(int j=i+1;j<numFalseCombs2;j++)
	    {
	      //ratio is sum/z1 in one pair-> z1=sum/ratio z2 = sum -z1
	      float locZ1_1=falseZ2[i]/falseZRatio2[i];
	      float locZ1_2=falseZ2[i]-locZ1_1;
	      float locZ2_1=falseZ2[j]/falseZRatio2[j];
	      float locZ2_2=falseZ2[j]-locZ2_1;
	      //make sure that we don't have the same hadrons in each combination
	      if(locZ1_1!=locZ2_1 && locZ1_2!=locZ2_2)
		{
		  phiRSum[numFalseCombs]=falseAngles2nd[i]+falseAngles2nd[j];
		  if(phiRSum[numFalseCombs]<-pi)
		    phiRSum[numFalseCombs]+=2*pi;
		  if(phiRSum[numFalseCombs]>pi)
		    phiRSum[numFalseCombs]-=2*pi;

		  if(falseAngles1st[i]>falseAngles1st[j])
		    phiRDiff[numFalseCombs]=falseAngles1st[i]-falseAngles1st[j];
		  else
		    phiRDiff[numFalseCombs]=falseAngles1st[j]-falseAngles1st[i];


		  //	      cout <<"phirsum2: " << phiRSum[numFalseCombs]<<" for " << numFalseCombs<<" before " << falseAngles2nd[i] <<" + " << falseAngles2nd[j] << " = " <<falseAngles2nd[i]+falseAngles2nd[j]<<endl;
		  z2[numFalseCombs]=falseZ2[i];
		  z2[numFalseCombs]=falseZ2[j];
				  //z1[numFalseCombs]=0.35;
		  //		  	      printf("combining in second pair z1: %f z2: %f ang1: %f ang2: %f\n",falseZ2[i],falseZ2[j],falseAngles2nd[i],falseAngles2nd[j]);
		  chargeType1[numFalseCombs]=falseCT2[i];
		  chargeType2[numFalseCombs]=falseCT2[j];
		  particleType1[numFalseCombs]=falsePT2[i];
		  if(falsePT2[i] > 20)
		    cout <<"falsept2: " << falsePT2[j] << " j: "  << j <<endl;
		  particleType2[numFalseCombs]=falsePT2[j];
		  multiplicity[numFalseCombs]=numFalseCombs2;
		  //		  cout << "false combs2: " << numFalseCombs2<<endl;
		  numFalseCombs++;
		}
	    }
	}
      if(numFalseCombs1==2 && numFalseCombs2==2)
	{
	  //	  cout <<"z1 counter: " << z1Counter <<endl;
	  //	  cout <<"num comb: "  << numFalseCombs << endl;
	  for(int i=0;i<numFalseCombs;i++)
	    {
	      //	      cout <<" sum angle <<" <<phiRSum[i]<<endl;     
	    }
	}
      //      cout <<"after ge" <<endl;
      //all counters should be the same, so it is a waste of space
      //      cout <<"z1Counter; " << z1Counter <<", z1Mem: " << memLocsF[z1Loc]<<endl;
      //      cout <<"false combs overall: " << numFalseCombs <<endl;
      for(int iHadsInEv=0;iHadsInEv<numFalseCombs;iHadsInEv++)
      //      for(int iHadsInEv=0;iHadsInEv<z1Counter;iHadsInEv++)
	{
	  //for some reason particle ids 90 or 99 are showing up...
	  int chargeIndx=ind(chargeType1[iHadsInEv],chargeType2[iHadsInEv],NumCharge);
	  int mIndex=ind(particleType1[iHadsInEv],particleType2[iHadsInEv],NumParticle);	  
	  if(chargeIndx > 100 || mIndex > 100)
	    {
	      cout <<  "ihads in : " << iHadsInEv << endl;
	      cout << "index screwd: " << chargeType1[iHadsInEv] <<" charge2: " << chargeType2[iHadsInEv] << " " <<  particleType1[iHadsInEv] <<" " << particleType2[iHadsInEv]<<endl;
	      continue;
	    }
	  if(particleType1[iHadsInEv] >=NumParticle || particleType2[iHadsInEv] >=NumParticle)
	    continue;
	  /*	  hThetaDiff11.Fill(theta11[iHadsInEv]-thrustTheta);
	  hThetaDiff12.Fill(theta12[iHadsInEv]-thrustTheta);
	  hThetaDiff21.Fill(theta21[iHadsInEv]-thrustTheta);
	  hThetaDiff22.Fill(theta22[iHadsInEv]-thrustTheta);
	  hThrustThetaRes.Fill(thrustTheta-thrustTheta_mc);*/
	  /*	  if(((float*)(memLocsF[z1Loc]))[iHadsInEv]>ZCut_for_mPlot && ((float*)(memLocsF[z2Loc]))[iHadsInEv]>ZCut_for_mPlot)
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
	  */
	  float zValue1=((float*)(memLocsF[z1Loc]))[iHadsInEv];
	  float mValue1=((float*)(memLocsF[mass1Loc]))[iHadsInEv];
	  float zValue2=((float*)(memLocsF[z2Loc]))[iHadsInEv];
	  float mValue2=((float*)(memLocsF[mass2Loc]))[iHadsInEv];
	  //	  kinematics[iBin][ind(chargeType1[iHadsInEv],chargeType2[iHadsInEv],NumCharge)][ind(particleType1[iHadsInEv],particleType2[iHadsInEv],NumParticle)][getBin(binningZ[particleType1[iHadsInEv]],zValue)][getBin(binningM[particleType1[iHadsInEv]],mValue)];

	  //doesn't make any sense, since the one is in cms, the other in lab
	  //	  if(fabs(theta11[iHadsInEv]-thrustTheta) > thetaCut || fabs(theta12[iHadsInEv]-thrustTheta) >thetaCut || fabs(theta21[iHadsInEv]-thrustTheta)>thetaCut || fabs(theta22[iHadsInEv]-thrustTheta)>thetaCut)
	  //	    continue;
	  //	  cout <<"getting zbin" << zValue1 <<"part: " << particleType1[iHadsInEv]<<endl;
	  int zbin1=getBin(binningZ[particleType1[iHadsInEv]],zValue1);
	  //	  	  cout <<"yup :"<< zbin1 <<endl;
	  int zbin2=getBin(binningZ[particleType1[iHadsInEv]],zValue2);
	  int zSpectbin1=getBin(binningZSpect[particleType1[iHadsInEv]],zValue1);
	  int zSpectbin2=getBin(binningZSpect[particleType1[iHadsInEv]],zValue2);
	  int mbin1=getBin(binningM[particleType1[iHadsInEv]],mValue1);
	  int mbin2=getBin(binningM[particleType1[iHadsInEv]],mValue2);
	    //	  int multbin1=getBin(binningMult[particleType1[iHadsInEv]],z1Counter);
	    //	  int multbin2=getBin(binningMult[particleType1[iHadsInEv]],z2Counter);

	  //	  int multbin1=getBin(binningMult[particleType1[iHadsInEv]],numFalseCombs1);
	  //	  int multbin2=getBin(binningMult[particleType1[iHadsInEv]],numFalseCombs2);
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
	  //	  for(int iBin=zBinning;iBin<=mBinning;iBin++)
	  for(int iBin=zBinning;iBin<=multBinning;iBin++)
	    {
	      switch(iBin)
		{
		case zBinning:
		  {

		    //		    hZ1.Fill(value1);
		    //		    hZ2.Fill(value2);

		    if(value2>zHighCuts[particleType1[iHadsInEv]] || value1 > zHighCuts[particleType1[iHadsInEv]])
		      {
			//			cout <<"z1 value: " << value1 << " z2 value: " << value2 <<endl;
			cuttedEv++;
			continue;
		      }
		    accEv++;

		    //		    cout <<"iHads: " << iHadsInEv <<endl;
		    //		    cout <<"val1: " << value1 <<" val2: " << value2 <<" value3 "<<value3 <<endl;
		    //		    cout <<"f1" <<endl;
		    //		    cout <<"fill with: " << particleType1[iHadsInEv] << " val1: " << value1 <<" val2: " << value2 << " val3: " << value3 <<" chargeIndex: " << chargeIndx <<" mIndex: " << mIndex <<endl;
		    fillBorders(binningZ[particleType1[iHadsInEv]],binningAng,value1,value2,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);
		    //		    cout <<"f2" <<endl;
#ifdef WITH_KIN_CORR

      if((chargeType1[iHadsInEv]==PN || chargeType1[iHadsInEv]==PZ || chargeType1[iHadsInEv]==ZN)&&(chargeType2[iHadsInEv]==PN || chargeType2[iHadsInEv]==PZ || chargeType2[iHadsInEv]==ZN))
      {
	//	cout << "ch1: " << m_chType1 <<" ch2: " << m_chType2 << " pt1: " << m_pt1 << " pt2:  " << m_pt2 << " zbin1: " << zbin1 << " zbin2; " << zbin2 << " mbin1: " << mbin1 <<" m2: " << mbin2 <<endl;

      }
#endif
#ifdef MC
		    if(phiRSum_mc[iHadsInEv]!=-1)
		      {
			hR.Fill(phiRSum_mc[iHadsInEv]-((float*)memLocsF[phiRSumLoc])[iHadsInEv]);
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
			    //more than 5 bin now...
			    //			    zSpectH[mbin1*5+mbin2].first->Fill(value1);
			    //			    zSpectH[mbin1*5+mbin2].second->Fill(value2);
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
		    fillBorders(binningM[particleType1[iHadsInEv]],binningAng,value1, value2,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);
#ifdef WITH_KIN_CORR
      if((chargeType1[iHadsInEv]==PN || chargeType1[iHadsInEv]==PZ || chargeType1[iHadsInEv]==ZN)&&(chargeType2[iHadsInEv]==PN || chargeType2[iHadsInEv]==PZ || chargeType2[iHadsInEv]==ZN))
      {
	//	cout << "ch1: " << m_chType1 <<" ch2: " << m_chType2 << " pt1: " << m_pt1 << " pt2:  " << m_pt2 << " zbin1: " << zbin1 << " zbin2; " << zbin2 << " mbin1: " << mbin1 <<" m2: " << mbin2 <<endl;

      }
#endif
#ifdef MC		    


#endif
		    break;
		  }
		case multBinning:
		  {
		    //		    float value1=z1Counter;
		    float value1=multiplicity[iHadsInEv];
		    if(isnan(value1))
		      {
			cout <<"value 1 m nan" <<endl;
			nanEv++;
		      continue;
		      }
		    float value2=multiplicity[iHadsInEv];
		    if(isnan(value2))
		      {
			cout <<"value 2 m nan" <<endl;
			nanEv++;
			continue;
		      }
		    float value3=((float*)(memLocsF[phiRSumLoc]))[iHadsInEv];
		    if(isnan(value3))
		      {
			cout <<"value 3 m nan" <<endl;
			nanEv++;
			continue;
		      }
		    if(particleType1[iHadsInEv]!=PiPi)
		      {
			//						cout << "had type not pip: " << particleType1[iHadsInEv]<<endl;
			//			continue;
		      }

		    //		    cout <<"1 value1: "<<value1 << " value2: " << value2 << " value3: " << value3<<endl;
		    //		    cout <<"filling xvals index: " << chargeIndx <<" " << mIndex <<" with "<< value1 <<" and " << value2 <<endl;
		    fillBorders(binningMult[particleType1[iHadsInEv]],binningAng,value1, value2,value3,counts[iBin][chargeIndx][mIndex],xVals[iBin][chargeIndx][mIndex],xValsNum[iBin][chargeIndx][mIndex],yVals[iBin][chargeIndx][mIndex],yValsNum[iBin][chargeIndx][mIndex]);
		    int multbin1=getBin(binningMult[particleType1[iHadsInEv]],value1);
		    int multbin2=getBin(binningMult[particleType1[iHadsInEv]],value2);
		    meanAngDiff[iBin][chargeIndx][mIndex][multbin1][multbin2]+=phiRDiff[iHadsInEv];
		    meanAngSum[iBin][chargeIndx][mIndex][multbin1][multbin2]+=phiRSum[iHadsInEv];
		    //		    cout <<"2" <<endl;
		    break;
		      }
		default:
		  cout << "binning not recognized" << endl;
		  //		  exit(0);
		}
	    }
	}
    }
  //and now the fit...
  cout<<"done with filling arryas, fitting..."<<endl;
  vector<vector<pair<float,float> >* > zAsyms;
  vector<vector<pair<float,float> >* >mAsyms;
  vector<vector<pair<float,float> >* >multAsyms;
  vector<vector<pair<float,float> >* > zAsymsMod;
  vector<vector<pair<float,float> >* >mAsymsMod;
  int numBins;
  string strCharge;
  string strParticle;
  vector<pair<float ,float> > * tmpAsym;
  vector<pair<float ,float> > * tmpAsymMod;
  vector<float> v_chi2;
  vector<int> v_ndfs;
  vector<vector<float>* > v_chi2Z;
  vector<vector<int>* > v_ndfsZ;
  vector<vector<float>* > v_chi2M;
  vector<vector<int>* > v_ndfsM;

  //for(int iBin=zBinning;iBin<=mBinning;iBin++)
  for(int iBin=zBinning;iBin<=multBinning;iBin++)
    {
      //      for(int iCh=PN;iCh<=ZZ;iCh++)
      //      cout <<"iBin: " << iBin <<endl;
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
	      //	      cout <<"last particlecomb: " << LAST_PARTICLE_COMB <<endl;
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
			    tmpAsym=fitTheSh__(counts[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins,errorVglFile,t_chi2,t_ndfs);
			    v_chi2Z.push_back(t_chi2);
			    v_ndfsZ.push_back(t_ndfs);
			    zAsyms.push_back(tmpAsym);

			    cout <<"size of asy vec: " << tmpAsym->size()<<endl;
			    for(int k=0;k<tmpAsym->size();k++)
			      {
				outFile  <<"Bin " << k << " Asymmetry: " << (*tmpAsym)[k].first << " +- " <<  (*tmpAsym)[k].second <<endl;

			      }
			    cout <<"wrote to file " <<endl;
			    break;
			  }
			case mBinning:
			  {
			    outFile << strCharge <<" " << strParticle << " M_inv Binning: " <<endl;
			    numBins=binningM[iPa1].size();

			    tmpAsym=fitTheSh__(counts[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins,errorVglFile, t_chi2,t_ndfs);
			    v_chi2M.push_back(t_chi2);
			    v_ndfsM.push_back(t_ndfs);
			    mAsyms.push_back(tmpAsym);

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
			    outFile << strCharge <<" " << strParticle << " mult Binning: " <<endl;
			    numBins=binningMult[iPa1].size();
			    tmpAsym=fitTheSh__(counts[iBin][ind(iCh1,iCh2,NumCharge)][ind(iPa1,iPa2,NumParticle)],binningAng,numBins,errorVglFile, t_chi2,t_ndfs);
			    v_chi2M.push_back(t_chi2);
			    v_ndfsM.push_back(t_ndfs);
			    multAsyms.push_back(tmpAsym);

			    for(int k=0;k<tmpAsym->size();k++)
			      {
				outFile  <<"Bin " << k << " Asymmetry: " << (*tmpAsym)[k].first << " +- " <<  (*tmpAsym)[k].second <<endl;
			      }
			    cout <<"after fit mult" <<endl;
			    break;
			  }
			default:
			  cout << "binning not recognized" << endl;
			  //			  exit(0);
			}
		    }
		}
	    }
	}
    }


  // maybe this part is not usable for the mult binning, where there is one binning more and only pip...
  modview myMod; 
  cout << "before " <<yVals[0][0][0][0][0] <<" nm: " << yValsNum[0][0][0][0][0]<<endl;
  myMod.getAsyms(zAsyms,mAsyms,xVals,yVals,xValsNum,yValsNum,v_chi2Z,v_chi2M,v_ndfsZ,v_ndfsM);
  myMod.doIt(zBinning,0,PN,PN,PiPi,PiPi,binningZ[PiPi].size(),0);
  TCanvas ct("ct","ct",10,20,500,200);
  myMod.pTG->Draw("AP*");
  ct.SaveAs(("FalseAsyms/"+myMod.name+".gif").c_str());
  TCanvas ct2("ct2","ct2",10,20,500,200);
  myMod.doIt(zBinning,1,PN,PN,PiPi,PiPi,binningZ[PiPi].size(),0);
  myMod.pTG->Draw("AP*");
  ct2.SaveAs(("FalseAsyms/"+myMod.name+".gif").c_str());
  TCanvas ct3("ct3","ct3",10,20,500,200);
  myMod.doIt(zBinning,2,PN,PN,PiPi,PiPi,binningZ[PiPi].size(),0);
  myMod.pTG->Draw("AP*");
  ct3.SaveAs(("FalseAsyms/"+myMod.name+".gif").c_str());
  TCanvas ct4("ct4","ct4",10,20,500,200);
  myMod.doIt(zBinning,3,PN,PN,PiPi,PiPi,binningZ[PiPi].size(),0);
  myMod.pTG->Draw("AP*");
  ct4.SaveAs(("FalseAsyms/"+myMod.name+".gif").c_str());
  TCanvas ct5("ct5","ct5",10,20,500,200);
  myMod.doIt(zBinning,4,PN,PN,PiPi,PiPi,binningZ[PiPi].size(),0);
  myMod.pTG->Draw("AP*");
  ct5.SaveAs(("FalseAsyms/"+myMod.name+".gif").c_str());
  
  float x[100];
  float y[100];
  float ex[100];
  float ey[100];

  float yPhiDiff[100];
  float eyPhiDiff[100];

  float yPhiSum[100];
  float eyPhiSum[100];

  TCanvas cMult("mult","mult",10,20,500,200);


      //  cMult.Divide(2,2);
  for(int i=0;i<binningMult[PiPi].size();i++)
    {
      cout <<"cd ing : " << i+1 <<endl;
      x[i]=xVals[multBinning][ind(PN,PN,NumCharge)][ind(PiPi,PiPi,NumParticle)][i][i]/xValsNum[multBinning][ind(PN,PN,NumCharge)][ind(PiPi,PiPi,NumParticle)][i][i];
      y[i]=(*multAsyms[0])[i*binningMult[PiPi].size()+i].first;
      yPhiDiff[i]=meanAngDiff[multBinning][ind(PN,PN,NumCharge)][ind(PiPi,PiPi,NumParticle)][i][i]/xValsNum[multBinning][ind(PN,PN,NumCharge)][ind(PiPi,PiPi,NumParticle)][i][i];
      cout <<"yphidiff: " << yPhiDiff[i]<<endl;
      yPhiSum[i]=meanAngSum[multBinning][ind(PN,PN,NumCharge)][ind(PiPi,PiPi,NumParticle)][i][i]/xValsNum[multBinning][ind(PN,PN,NumCharge)][ind(PiPi,PiPi,NumParticle)][i][i];
      eyPhiDiff[i]=0.0;
      eyPhiSum[i]=0.0;
      ex[i]=0;
      ey[i]=(*multAsyms[0])[i*binningMult[PiPi].size()+i].second;
      cout <<"x: " << x[i] <<", y[i] " << y[i] << ", ey: " << ey[i] <<endl;
    }
  TGraphErrors tg(binningMult[PiPi].size(),x,y,ex,ey);
  TGraphErrors tgSum(binningMult[PiPi].size(),x,yPhiSum,ex,eyPhiSum);
  TGraphErrors tgDiff(binningMult[PiPi].size(),x,yPhiDiff,ex,eyPhiDiff);

  //  tg.GetYaxis()->SetRangeUser(-0.1,0);
  //				tg.GetYaxis()->SetNdivisions(6,true);
  //				tg.GetXaxis()->SetNdivisions(6,true);
  tg.SetMarkerStyle(23); 
  tg.SetMarkerSize(1.2);
  tg.SetMarkerColor(kRed);

  tgSum.SetMarkerStyle(23); 
  tgSum.SetMarkerSize(1.2);
  tgSum.SetMarkerColor(kRed);
  
  tgDiff.SetMarkerStyle(23); 
  tgDiff.SetMarkerSize(1.2);
  tgDiff.SetMarkerColor(kRed);
  
      /*			tg.GetYaxis()->SetLabelSize(0.07);
				  tg.GetYaxis()->SetTitleSize(0.07);
				  tg.GetXaxis()->SetTitleSize(0.07);
				  tg.GetXaxis()->SetLabelSize(0.07);*/
  tg.GetYaxis()->SetTitle("A_{(#Phi_{1}+#Phi_{2})}");
  tg.GetXaxis()->SetTitle("Number of Combinations");
  tg.Draw("AP+");

  cMult.SaveAs("multAs.pdf");
  cMult.SaveAs("multAs.png");
  cMult.SaveAs("multAs.eps");
  TCanvas cMultPhiSum("multPhiSum","multPhiSum",10,20,500,200);
  tgSum.GetYaxis()->SetTitle("A_{(#Phi_{1}+#Phi_{2})}");
  tgSum.GetXaxis()->SetTitle("Number of Combinations");
  tgSum.Draw("AP+");

  cMultPhiSum.SaveAs("multAsSum.pdf");
  cMultPhiSum.SaveAs("multAsSum.png");
  cMultPhiSum.SaveAs("multAsSum.eps");
  TCanvas cMultPhiDiff("multPhiDiff","multPhiDiff",10,20,500,200);
  tgDiff.GetYaxis()->SetTitle("A_{(#Phi_{1}+#Phi_{2})}");
  tgDiff.GetXaxis()->SetTitle("Number of Combinations");
  tgDiff.Draw("AP+");

  cMultPhiDiff.SaveAs("multAsDiff.pdf");
  cMultPhiDiff.SaveAs("multAsDiff.png");
  cMultPhiDiff.SaveAs("multAsDiff.eps");
      
  vector<int> vecBinning;
  //draw the asymmetries...
  vector<TGraphErrors*> vecGraphs;
  vector<TGraphErrors*> vecGraphsKinCorr;
  vector<TGraphErrors*> vecGraphsProjections;
  vector<string> graphTitles;
  vector<string> graphTitlesKinCorr;
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

  TCanvas cZ("cz","cz",10,20,500,200);
  TH1D hZ("hZ","hZ",50,-5,5);

  for(int i=0;i<zAsyms.size();i++)
    {
      for(int j=0;j<(*zAsyms[i]).size();j++)
	{
	  if(!zAsyms[i]||(*zAsyms[i])[j].second > 3) //no counts
	    continue;
	  hZ.Fill((*zAsyms[i])[j].first/(*zAsyms[i])[j].second);
	  asymSum+=(*zAsyms[i])[j].first;
	  asymSumErr+=(*zAsyms[i])[j].second*(*zAsyms[i])[j].second;

	  outputFile << "adding Asymmetry (z): " << (*zAsyms[i])[j].first <<" +- " << (*zAsyms[i])[j].second <<endl;
	  aCounter++;
	}
    }
  hZ.Draw();
  cZ.SaveAs("FakeAsymZ.ps");
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
  outputFile <<"mAsymAvg: " << asymSum/(float)aCounter<< "+-" << sqrt(asymSumErr)/(float)aCounter<<endl;
  outputFile <<"num Asyms: " << aCounter <<endl;

 
  aCounter=0;
  asymSum=0;

  

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
  //  saveAs("Asymmetries.ps",vecGraphs,graphTitles,outputFile,vecBinning);
  outputFile << "saving " << vecGraphsProjections.size() <<" projections having " << graphTitlesProj.size() <<" titles"<<endl;

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


  cout <<"cuttedEv: " << cuttedEv <<" nanEv: " << nanEv << " accEv: " << accEv <<endl;
}



  
