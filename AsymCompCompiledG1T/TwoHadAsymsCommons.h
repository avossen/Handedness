#ifndef TWOHADASYMSCommons_H
#define TWOHADASYMSCommons_H

//#define LAST_CHARGE_COMB ZNNZ
#define LAST_CHARGE_COMB PNNP
//#define LAST_PARTICLE_COMB KK
#define LAST_PARTICLE_COMB PiK

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>
#include "TCanvas.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <TF2.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph2DErrors.h>
#include <TGraphErrors.h>


//#include <pair> 
using namespace std;
//enum TwoHadCharge{PN, NP,PP,NN,PZ,ZP,ZN,NZ,ZZ,PNNP,PZZP,ZNNZ,NA};//the unknown is for mc, when we
//have to order according to charge
enum TwoHadCharge{PN, NP,PP,NN,PZ,ZP,ZN,NZ,ZZ,PNNP,PZZP,ZNNZ,NA, hadChargeEnd};//the unknown is for mc, when we
//enum TwoHadCharge{PN, NP,PP,NN,PZ,ZP,ZN,NZ,ZZ,NA};//the unknown is for mc
enum kMCFlags{mcFlagNone,mcFlagMC,mcFlagWoA,mcAsData};
enum TwoHadPType{PiPi,PiK,KPi,KK,UNKNOWN, hadTypeEnd};
enum WeightingType{Simple,ThetaOpeningAngle,MZ};


enum ErrorCodes{
  smallRt
  };

//single binnings
enum binning1Types{thetaBinning, dilFactBinning};
//if mult binning is later, the false asymmetry code won't work anymore...
enum binningTypes{zBinning, mBinning, multBinning, mixedBinning,mzBinning};

enum binning4Types{fourBinningOne};

extern int NumCharge;
extern int NumPars; //number of parameters for kinematics fitting
extern int NumBin;
extern double pi;
extern int __numKinBins;
extern int NumParticle;
extern int gHCounter;

float getInc(float phiR, double mass1, double mass2,  double z1, double z2,int pt1, int pt2, int charge1,int charge2);
float getInc2(int binning, int index1, int index2, int mBin1, int mBin2, float zVal1, float zVal2);
void loadBinning(vector<float>* binningM, vector<float>* binningZ, vector<float>* binningMReduced=0, vector<float>* binningZReduced=0);

void cpyVect(vector<float>& v1, vector<float>& v2);

//gives index in charge, ptype  array. Since the combination does not depend on the order, we can save space by using an array of half the size
unsigned int ind(int i, int j, int max);

/*
to compute the error on the purities... hopefully the error is computed correctly..
 */


//void 

namespace kinematics
  {
    extern int zBin1, zBin2, mBin1, mBin2; //kinematic bin of the computed purity
    extern int chargeType1;    //charge & particle type for the purity
    extern int particleType1;
    extern int chargeType2;    //charge & particle type for the purity
    extern int particleType2;
    extern int counts;
  }

float getError(unsigned long numCorr, unsigned long  numFalse);

//loads purities, has of course to be synchronized with the root file
void loadPurities(float***** purities, float***** puritiesErrors,char* rFileName, int lastPT, int lastCharge, int numBinZ, int numBinM);

void loadKinPars1D(TH1D***** m_kinSpectraH);

void saveKinPars(TH1D***** kinSpectraH, vector<float>* binningM );

void saveKinPars(TH2D***** kinSpectraH, vector<float>* binningM);

//loads kinematics, has of course to be synchronized with the root file
void loadKinematics(double****** kinSpectraReducedPars,double******* kinSpectraPars,double***** kinSpectraReduced,double****** kinSpectra, char* rFileName, int lastPT, int lastCharge, int numBinZ, int numBinM, vector<float>* binningM);

string iToStrCharge(int i);

string iToStrPart(int i);

int getBin(vector<float>& b1, float value);
int getHBin(float val, int numBin, float maxVal);


//inc the kinematic border
template<class T> void incBorders(vector<float>& b1, float value1, float value2, T** counts);
void incBorders(vector<float>& b1, float value1, float value2, float** counts, float incr);
void incBorders(vector<float>& b1, vector<float>& b2, float value1, float value2, float** counts, float incr);
template<class T> void incBorders(vector<float>& b1, vector<float>& b2,float value1, float value2, T** counts);


void fillBorders(vector<float>& b1, float value1, int* counts);
void fillBorders(vector<float>& b1, vector<float>& b2, float value1, int* counts);

//b1 is the border for the kin binning, need two coordinates, since we have 2D binning, 
void fillBorders(vector<float>& b1, vector<float>& b2,  float value1, float value2, float value3,  int*** m_counts,float** xVals,int** xValsNum, float** yVals, int**yValsNum);
void fillBordersMix(vector<float>& b1,vector<float>& b2, vector<float>& b3,  float value1, float value2, float value3,  int*** m_counts,float** xVals,int** xValsNum, float** yVals, int** yValsNum);

void fillBorders(vector<float>& b1, float value1, int* counts);

void fillBorders(vector<float>& b1, float value1, float* counts, float inc);
//b1 is the border for the kin binning, need two coordinates, since we have 2D binning, 
void fillBorders(vector<float>& b1, vector<float>& b2,  float value1, float value2, float value3,  float*** m_counts,float** xVals,int** xValsNum, float** yVals, int** yValsNum, float inc);
void fillBordersMix(vector<float>& b1,vector<float>& b2, vector<float>& b3,  float value1, float value2, float value3,  float*** m_counts,float** xVals,int** xValsNum, float** yVals, int** yValsNum, float inc);

void fillBorders(vector<float>& b1, float value1, float* counts,float inc);
//for singles
void fillBorders(vector<float>& b1, vector<float>& b2,float value1, float value2, int** counts, float* xVals, int* xValsNum);
void fillBorders(vector<float>& b1, vector<float>& b2,float value1, float value2, float** counts, float* xVals, int* xValsNum, float inc);
void fillBorders(vector<float>& b1, vector<float>& b2,vector<float>& b3,float value1, float value2, int** counts, float* xVals, int* xValsNum);
 
void fillBorders4(vector<float>& b1, vector<float>& b2, vector<float>& b3,float value1, float value2, float value3, float value4,float value5,int***** counts, float**** xVals, int**** xValsNum, float**** yVals, int**** yValsNum,float**** aVals, int**** aValsNum,float**** bVals, int**** bValsNum );

vector<pair<float, float> >*** fitTheFour(int***** counts,vector<float>& binningAng,int numKinBin1, int numKinBins2, ofstream* errorVglFile, vector<float>* v_chi2=0,vector<int>* v_ndf=0);

vector<pair<float, float> >* fitTheSingle(int** counts,vector<float>& binningAng,int numKinBin1, ofstream* errorVglFile, vector<float>* v_chi2=0,vector<int>* v_ndf=0);
vector<pair<float, float> >* fitTheSingle(float** counts,vector<float>& binningAng,int numKinBin1, ofstream* errorVglFile, vector<float>* v_chi2=0,vector<int>* v_ndf=0);

//numPar is the parameter that one wants to retrieve (only useful with higher harmonics), 1 is the normal amplitude, 2, 3, 4 are the sin and  higher harmonics
vector<pair<float, float> >* fitTheSh__(int*** counts, vector<float>& binningAng,int numKinBins, ofstream* errorVglFile, vector<float>* v_chi2=0, vector<int>* v_ndf=0, int numPar=1);

vector<pair<float, float> >* fitTheSh__(float*** counts, vector<float>& binningAng,int numKinBins, ofstream* errorVglFile, int numPara=1);
   


TH1D*** allocHistos1D(int dim1, int dim2);

TH1D**** allocHistos1D(int dim1, int dim2, int dim3);

TH1D***** allocHistos1D(int dim1, int dim2, int dim3, int dim4);

TH1D****** allocHistos1D(int dim1, int dim2, int dim3, int dim4, int dim5);

//---

TH2D*** allocHistos(int dim1, int dim2);

TH2D**** allocHistos(int dim1, int dim2, int dim3);

TH2D***** allocHistos(int dim1, int dim2, int dim3, int dim4);

TH2D****** allocHistos(int dim1, int dim2, int dim3, int dim4, int dim5);

template<class T> T** allocateArray(int dim1, int dim2);
template<class T> T*** allocateArray(int dim1, int dim2, int dim3);
template<class T> T**** allocateArray(int dim1, int dim2, int dim3, int dim4);
template<class T> T***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5);
template<class T> T****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6);
template<class T> T******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7);



void saveAs(char* fileName,vector<TH1D*>& histos);

void replDot(char* title);

void replDash(char* title);

void saveAs(char* fileExt, vector<TGraphErrors*>& mGraphs, vector<string>& graphTitles, ofstream& oFile,vector<int>& vecKinBins);

void saveAs(char* fileExt, vector<TGraph*>& mGraphs, vector<string>& graphTitles,vector<int>& vecKinBins);

#endif
