//#define MAX_EVENTS 100000
#define KAON_LESS_CLASS //this has to be consistent with Ext2Charm... better as command line argument
#define ONLY_DIST
//#define ONLY_MASS
#define D0Mass 1.865
//essentially resolution, just a guess
#define D0Width 0.6
#define D0Lower  D0Mass-D0Width/2
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

// z bins in bacchetta: 0.2 < z<0.3 ; 0.3<z<0.4
//0.4<z<0.55; 0.55<z<0.75

using namespace std;

vector<float>* binningM;
vector<float>* binningZ;
vector<float>* binningTheta;
#define OPENING_CUT 0.8
#define MAX_OPENING_CUT 2.0


int main(int argc, char** argv)
{
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

  char mvaMethod[]="MLP method";
  //  float mvaCut=0.9;
  int maxKin=-1;
  vector<string> fieldNamesF;
  vector<string> fieldNamesI;

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
    }

  int** passAsCharmZ=new int*[maxKin];
  int** passAsCharmM=new int*[maxKin];
  int** noPassAsCharmZ=new int*[maxKin];
  int** noPassAsCharmM=new int*[maxKin];

  int** passAsCharmZTheta=new int*[maxKin];
  int** passAsCharmMTheta=new int*[maxKin];
  int** noPassAsCharmZTheta=new int*[maxKin];
  int** noPassAsCharmMTheta=new int*[maxKin];

  int** passAsCharmMZ=new int*[maxKin];
  int** noPassAsCharmMZ=new int*[maxKin];

  //mean kinematics of m/z in z/m bins 
  // numKinBins (m or z), (m orz)
  float*** meanKins=allocateArray<float>(2,2,maxKin);
  int** meanKinsNum=allocateArray<int>(2,maxKin);
  float*** meanKinsPassCuts=allocateArray<float>(2,2,maxKin);
  int** meanKinsNumPassCuts=allocateArray<int>(2,maxKin);

  float*** meanKinsFailCuts=allocateArray<float>(2,2,maxKin);
  int** meanKinsNumFailCuts=allocateArray<int>(2,maxKin);
  //  float** meanM2inZ=new float*[maxKin];
  //  float** meanZ1inM=new float*[maxKin]; //for each z bin, we 
  //  float** meanZ2inM=new float*[maxKin]; //for each z bin, we 

  //  int** meanMinZNum=new int*[maxKin];
  //  int** meanMinMNum=new int*[maxKin];

  for(int i=0;i<maxKin;i++)
    {
      passAsCharmM[i]=new int[maxKin];
      noPassAsCharmM[i]=new int[maxKin];
      passAsCharmZ[i]=new int[maxKin];
      noPassAsCharmZ[i]=new int[maxKin];
      passAsCharmMZ[i]=new int[maxKin];

      passAsCharmMTheta[i]=new int[maxKin];
      noPassAsCharmMTheta[i]=new int[maxKin];
      passAsCharmZTheta[i]=new int[maxKin];
      noPassAsCharmZTheta[i]=new int[maxKin];
      noPassAsCharmMZ[i]=new int[maxKin];

      /*  meanM1inZ[i]=new float[maxKin];
      meanM2inZ[i]=new float[maxKin];
      meanZ1inM[i]=new float[maxKin];
      meanZ2inM[i]=new float[maxKin];
      meanMinZNum[i]=new int[maxKin];
      meanZinMNum[i]=new int[maxKin];*/

      for(int j=0;j<maxKin;j++)
	{
	  passAsCharmM[i][j]=0;
	  passAsCharmZ[i][j]=0;
	  noPassAsCharmM[i][j]=0;
	  noPassAsCharmZ[i][j]=0;

	  passAsCharmMZ[i][j]=0;
	  noPassAsCharmMZ[i][j]=0;
	  passAsCharmMTheta[i][j]=0;
	  passAsCharmZTheta[i][j]=0;
	  noPassAsCharmMTheta[i][j]=0;
	  noPassAsCharmZTheta[i][j]=0;

	  /*	  meanM1inZ[i][j]=0;
	  meanM2inZ[i][j]=0;
	  meanZ1inM[i][j]=0;
	  meanZ2inM[i][j]=0;

	  meanMinZNum[i][j]=0;
	  meanZinMNum[i][j]=0;*/
	}
    }

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
  float mvaCut=0.4; //bdtd is pos/neg.. default
  if(argc==3)
    {
      char* c_mvaCut=argv2[2];
      mvaCut=atof(c_mvaCut);
    }
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
  if((folderName.find("charm")!=string::npos )||folderName.find("Charm")!=string::npos )
    isCharm=true;

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
  vector<string> varNamesF;
  vector<string> varNamesI;
  //respect order!!
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
  int noMvaCounter=0;
  for(long i=0;i<nevents;i++)
    {
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
	  noMvaCounter++;
#ifndef KAON_LESS_CLASS
	  continue;
#endif
	}

      if(!(i%100000))
	cout <<"processing event nr " << i << " of " << nevents << "(" << 100*i/(float)nevents<< "% )"<<endl;

      float shiftDr=log(fabs(dr));
      float LogmeanDist=log(meanDist);
      float LogmeanDistPN=log(meanDistPN);
      float LogmeanDistK=log(meanDistK);
      float LogmaxExtraVertDistK=log(maxExtraVertDistK);
      //      cout <<"maxIntraVertDistK: " << maxIntraVertDistK <<endl;
      float LogmaxIntraVertDistK=log(maxIntraVertDistK);
      /*      if(maxIntraVertDistK==0)
	{
	LogmaxIntraVertDistK=0;
	continue;//be safe
	}*/
      float LogmaxExtraVertDistPN=log(maxExtraVertDistPN);
      float LogmaxIntraVertDistPN=log(maxIntraVertDistPN);
      float LogminMassDiffToD=log(minMassDiffToD);
      float LogminMassDiffToD_PN=log(minMassDiffToD_PN);

      //test that all variables are within range

#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
      
      //      if(false)
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
	  if(!((particleType1[iHadsInEv]==PiPi) && (particleType2[iHadsInEv]==PiPi)))
	    continue;
	  if(!((chargeType1[iHadsInEv]==PN) && (chargeType2[iHadsInEv]==PN)))
	    continue;
	  if(fabs(thrustProj11[iHadsInEv])>MAX_OPENING_CUT || fabs(thrustProj12[iHadsInEv])>MAX_OPENING_CUT ||fabs(thrustProj21[iHadsInEv])>MAX_OPENING_CUT ||fabs(thrustProj22[iHadsInEv])>MAX_OPENING_CUT)
	    continue;
	  if(fabs(thrustProj11[iHadsInEv])<OPENING_CUT || fabs(thrustProj12[iHadsInEv])<OPENING_CUT ||fabs(thrustProj21[iHadsInEv])<OPENING_CUT ||fabs(thrustProj22[iHadsInEv])<OPENING_CUT)
	    {
	      continue;
	    }

	  int chargeIndx=ind(chargeType1[iHadsInEv],chargeType2[iHadsInEv],NumCharge);
	  int mIndex=ind(particleType1[iHadsInEv],particleType2[iHadsInEv],NumParticle);	  

	  float zValue1=((float*)(memLocsF[z1Loc]))[iHadsInEv];
	  float mValue1=((float*)(memLocsF[mass1Loc]))[iHadsInEv];
	  float zValue2=((float*)(memLocsF[z2Loc]))[iHadsInEv];
	  float mValue2=((float*)(memLocsF[mass2Loc]))[iHadsInEv];
	  float thetaValue=sin(thetaEThrust)*sin(thetaEThrust)/(1+cos(thetaEThrust)*cos(thetaEThrust));

	  if(mValue1>mHighCuts[particleType1[iHadsInEv]] || mValue2>mHighCuts[particleType1[iHadsInEv]])
	    continue;
	  if(zValue1 >1 || zValue1 <0 || zValue2>1 || zValue2 <0)
	    continue;
	  int zbin1=getBin(binningZ[particleType1[iHadsInEv]],zValue1);
	  int zbin2=getBin(binningZ[particleType1[iHadsInEv]],zValue2);
	  int mbin1=getBin(binningM[particleType1[iHadsInEv]],mValue1);
	  int mbin2=getBin(binningM[particleType1[iHadsInEv]],mValue2);
	  int thetaBin=getBin(binningTheta[particleType1[iHadsInEv]],thetaValue);

	  int m_chType1=chargeType1[iHadsInEv];
	  int m_pt1=particleType1[iHadsInEv];

	  int m_chType2=chargeType2[iHadsInEv];
	  int m_pt2=particleType2[iHadsInEv];
	  float value1=((float*)(memLocsF[z1Loc]))[iHadsInEv];
	  if(isnan(value1))
	    {
	      cout <<"value 1 nan" <<endl;
	      continue;
	    }
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
	  //	  valsF.push_back(LogmaxIntraVertDistK);
	  valsF.push_back(LogmaxExtraVertDistPN);
	  valsF.push_back(LogmaxIntraVertDistPN);
#ifndef KAON_LESS_CLASS
	  valsF.push_back(maxCosThetaK);
#endif
#ifndef ONLY_DIST
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
	  meanKins[zBinning][mBinning][zbin1]+=mValue1;
	  meanKins[zBinning][mBinning][zbin2]+=mValue2;
	  meanKins[zBinning][zBinning][zbin1]+=zValue1;
	  meanKins[zBinning][zBinning][zbin2]+=zValue2;
	  meanKins[mBinning][zBinning][mbin1]+=zValue1;
	  meanKins[mBinning][zBinning][mbin2]+=zValue2;
	  meanKins[mBinning][mBinning][mbin1]+=mValue1;
	  meanKins[mBinning][mBinning][mbin2]+=mValue2;
	  meanKinsNum[zBinning][zbin1]++;
	  meanKinsNum[zBinning][zbin2]++;
	  meanKinsNum[mBinning][mbin1]++;
	  meanKinsNum[mBinning][mbin2]++;
#ifdef ONLY_MASS
	  if(fabs(minMassDiffToD)< mvaCut*D0Width)
#else
	  if(myMva.evaluate(valsF,valsI,mvaMethod)>mvaCut)
#endif
	    {
	      passAsCharmMZ[mbin1][zbin1]++;
	      passAsCharmZ[zbin1][zbin2]++;
	      passAsCharmM[mbin1][mbin2]++;
	      passAsCharmZTheta[zbin1][thetaBin]++;
	      passAsCharmMTheta[mbin1][thetaBin]++;

	      meanKinsPassCuts[zBinning][mBinning][zbin1]+=mValue1;
	      meanKinsPassCuts[zBinning][mBinning][zbin2]+=mValue2;
	      meanKinsPassCuts[zBinning][zBinning][zbin1]+=zValue1;
	      meanKinsPassCuts[zBinning][zBinning][zbin2]+=zValue2;
	      meanKinsPassCuts[mBinning][zBinning][mbin1]+=zValue1;
	      meanKinsPassCuts[mBinning][zBinning][mbin2]+=zValue2;
	      meanKinsPassCuts[mBinning][mBinning][mbin1]+=mValue1;
	      meanKinsPassCuts[mBinning][mBinning][mbin2]+=mValue2;
	      meanKinsNumPassCuts[zBinning][zbin1]++;
	      meanKinsNumPassCuts[zBinning][zbin2]++;
	      meanKinsNumPassCuts[mBinning][mbin1]++;
	      meanKinsNumPassCuts[mBinning][mbin2]++;
	    }
	  else
	    {
	      noPassAsCharmMZ[mbin1][zbin1]++;
	      noPassAsCharmZ[zbin1][zbin2]++;
	      noPassAsCharmM[mbin1][mbin2]++;
	      noPassAsCharmZTheta[zbin1][thetaBin]++;
	      noPassAsCharmMTheta[mbin1][thetaBin]++;

	      meanKinsFailCuts[zBinning][mBinning][zbin1]+=mValue1;
	      meanKinsFailCuts[zBinning][mBinning][zbin2]+=mValue2;
	      meanKinsFailCuts[zBinning][zBinning][zbin1]+=zValue1;
	      meanKinsFailCuts[zBinning][zBinning][zbin2]+=zValue2;
	      meanKinsFailCuts[mBinning][zBinning][mbin1]+=zValue1;
	      meanKinsFailCuts[mBinning][zBinning][mbin2]+=zValue2;
	      meanKinsFailCuts[mBinning][mBinning][mbin1]+=mValue1;
	      meanKinsFailCuts[mBinning][mBinning][mbin2]+=mValue2;
	      meanKinsNumFailCuts[zBinning][zbin1]++;
	      meanKinsNumFailCuts[zBinning][zbin2]++;
	      meanKinsNumFailCuts[mBinning][mbin1]++;
	      meanKinsNumFailCuts[mBinning][mbin2]++;

	    }
	}
    }
  cout <<" z binning: " << endl;

  char c_tmp[100];
  sprintf(c_tmp,"%0.2f",mvaCut);
  replDot(c_tmp);
#ifndef ONLY_MASS
#ifdef KAON_LESS_CLASS
#ifdef ONLY_DIST
  TFile m_file((string("charmContribOnlyDist")+c_tmp+".root").c_str(),"UPDATE");
#else
  TFile m_file((string("charmContribKaonLess")+c_tmp+".root").c_str(),"UPDATE");
#endif
#else
  TFile m_file((string("charmContrib")+c_tmp+".root").c_str(),"UPDATE");
#endif
#else
  TFile m_file((string("charmContribOnlyMass")+c_tmp+".root").c_str(),"UPDATE");
#endif
  string s;
  if(isCharm)
    s="charmTree";
  else
    s="udsTree";
  TTree m_tree(s.c_str(),s.c_str());

  int kinBin1=0;
  int kinBin2=0;
  int binning=0;
  int numCorr=0;
  int numFalse=0;
  float cutValue=mvaCut;
  float meanM;
  float meanZ;
  float meanMFailCuts;
  float meanZFailCuts;
  float meanMPassCuts;
  float meanZPassCuts;

  m_tree.Branch("kinBin1",&kinBin1,"kinBin1/I");
  m_tree.Branch("kinBin2",&kinBin2,"kinBin2/I");
  m_tree.Branch("binning",&binning,"binning/I");
  m_tree.Branch("numCorr",&numCorr,"numCorr/I");
  m_tree.Branch("numFalse",&numFalse,"numFalse/I");
  m_tree.Branch("cutValue",&cutValue,"cutValue/F");
  m_tree.Branch("meanM",&meanM,"meanM/F");
  m_tree.Branch("meanZ",&meanZ,"meanZ/F");
  m_tree.Branch("meanMFailCuts",&meanMFailCuts,"meanMFailCuts/F");
  m_tree.Branch("meanZFailCuts",&meanZFailCuts,"meanZFailCuts/F");
  m_tree.Branch("meanMPassCuts",&meanMPassCuts,"meanMPassCuts/F");
  m_tree.Branch("meanZPassCuts",&meanZPassCuts,"meanZPassCuts/F");
  //ok, we 'mis'use the mult and mixed binning for m theta and z theta


  for(int i=0;i<maxKin;i++)
    {
      cout <<"mbin1: " << i <<endl;
      for(int j=0;j<maxKin;j++)
	{
	  cout <<"mbin2: " << j <<" numCounts: " << passAsCharmM[i][j]+noPassAsCharmM[i][j];
	}
    }

  int overall=0;
  for(int i=0;i<maxKin;i++)
    {
      cout <<"mbin: " << i <<endl;
      for(int j=0;j<maxKin;j++)
	{
	  cout <<"zbin: " << j <<" numCounts: " << passAsCharmMZ[i][j]+noPassAsCharmMZ[i][j];
	  overall+=passAsCharmMZ[i][j]+noPassAsCharmMZ[i][j];
	}
      cout <<endl;
    }

  cout <<"overall: " << overall<<endl;
  for(int iBin=zBinning;iBin<=mzBinning;iBin++)
    {
      binning=iBin;
      for(int i=0;i<maxKin;i++)
	{
	  for(int j=0;j<maxKin;j++)
	    {
	      kinBin1=i;
	      kinBin2=j;
	      cout <<"kinBin1: " << kinBin1 <<" kinBin2: " << kinBin2 <<endl;
	      switch(iBin)
		{
		case zBinning:
		  numCorr=passAsCharmZ[i][j];
		  numFalse=noPassAsCharmZ[i][j];
		  cout <<"writing to zbinning, numCorr: " << numCorr <<" numFalse: " << numFalse <<endl;
		  meanM=meanKins[zBinning][mBinning][kinBin1]/meanKinsNum[zBinning][kinBin1];
		  meanZ=meanKins[zBinning][zBinning][kinBin1]/meanKinsNum[zBinning][kinBin1];
		  meanMFailCuts=meanKinsFailCuts[zBinning][mBinning][kinBin1]/meanKinsNumFailCuts[zBinning][kinBin1];
		  meanZFailCuts=meanKinsFailCuts[zBinning][zBinning][kinBin1]/meanKinsNumFailCuts[zBinning][kinBin1];
		  meanMPassCuts=meanKinsPassCuts[zBinning][mBinning][kinBin1]/meanKinsNumPassCuts[zBinning][kinBin1];
		  meanZPassCuts=meanKinsPassCuts[zBinning][zBinning][kinBin1]/meanKinsNumPassCuts[zBinning][kinBin1];
		  break;
		case mBinning:
		  numCorr=passAsCharmM[i][j];
		  numFalse=noPassAsCharmM[i][j];
		  meanM=meanKins[mBinning][mBinning][kinBin1]/meanKinsNum[mBinning][kinBin1];
		  meanZ=meanKins[mBinning][zBinning][kinBin1]/meanKinsNum[mBinning][kinBin1];
		  meanMFailCuts=meanKinsFailCuts[mBinning][mBinning][kinBin1]/meanKinsNumFailCuts[mBinning][kinBin1];
		  meanZFailCuts=meanKinsFailCuts[mBinning][zBinning][kinBin1]/meanKinsNumFailCuts[mBinning][kinBin1];
		  meanMPassCuts=meanKinsPassCuts[mBinning][mBinning][kinBin1]/meanKinsNumPassCuts[mBinning][kinBin1];
		  meanZPassCuts=meanKinsPassCuts[mBinning][zBinning][kinBin1]/meanKinsNumPassCuts[mBinning][kinBin1];
		  cout <<"writing to mbinning, numCorr: " << numCorr <<" numFalse: " << numFalse <<endl;
		  break;
		case multBinning:
		  numCorr=passAsCharmMTheta[i][j];
		  numFalse=noPassAsCharmMTheta[i][j];
		case mixedBinning:
		  numCorr=passAsCharmZTheta[i][j];
		  numFalse=noPassAsCharmZTheta[i][j];
		  break;
		case mzBinning:
		  numCorr=passAsCharmMZ[i][j];
		  numFalse=noPassAsCharmMZ[i][j];
		  break;
		default:
		  cout <<"wrong binning!!!";
		}

	      m_tree.Fill();
	    }
	}
    }
  m_file.Write();
  m_file.Close();
  /*  for(int i=0;i<maxKin;i++)
    {
      cout << i << " " << passAsCharmZ[i] <<"  noPass: " << noPassAsCharmZ[i] <<" ratio: " << passAsCharmZ[i]/(float)noPassAsCharmZ[i]<<endl;
    }
  cout <<" m binning: " <<endl;
  for(int i=0;i<maxKin;i++)
    {
      cout << i << " " << passAsCharmM[i] <<"  noPass: " << noPassAsCharmM[i] <<" ratio: " << passAsCharmM[i]/(float)noPassAsCharmM[i]<<endl;
      }*/
  cout <<"mva: " << mvaCounter <<" nomva: " << noMvaCounter <<endl;
}
